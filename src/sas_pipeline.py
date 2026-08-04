#!/usr/bin/env python3
"""
SAS (Sufficient Alignment Support) Integrated Pipeline
Complete workflow for multi-technology sequencing data analysis

This script orchestrates the complete SAS pipeline integrating:
1. BE (Base Error) detection using short-read, HiFi, and ONT data
2. SE (Structural Error) detection using filtered ONT data

Features:
- Multi-technology data integration (Illumina, HiFi, ONT)
- Automatic reference index generation
- Optional BAM filtering skip
- Comprehensive error detection and classification
- Parallel processing support
- Detailed logging and reporting
- Chromosome-specific analysis and QV calculation
- Bed file masking support

License: MIT
"""

import os
import sys
import argparse
import subprocess
import logging
import shutil
import tempfile
import contextlib
from pathlib import Path
from datetime import datetime
import json
import sqlite3
from collections import OrderedDict, defaultdict
import re
from bisect import bisect_left

try:
    import pysam
except ModuleNotFoundError:
    pysam = None

try:
    from ._version import __version__
except ImportError:
    try:
        from src._version import __version__
    except ImportError:
        from _version import __version__


MAX_OPEN_CHROMOSOME_FILES = 64
MAX_COMMAND_LOG_BYTES = 64 * 1024


class CLIHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Keep examples and explicitly documented defaults readable."""


def positive_int(value):
    """Parse a strictly positive command-line integer."""
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be greater than 0")
    return parsed


def _same_existing_file(first_path, second_path):
    """Return whether two existing paths identify the same filesystem object.

    ``samefile`` itself is allowed to fail here: a missing/inaccessible path
    must continue through the normal rollback copy path so that any real
    restore failure is still reported to the caller.
    """
    try:
        return os.path.samefile(first_path, second_path)
    except OSError:
        return False


def _is_chromosome_header(columns):
    """Return whether columns begin with a chromosome/start/end header."""
    if isinstance(columns, str):
        columns = columns.split("\t")
    if len(columns) < 3:
        return False
    normalized = [column.strip().lower() for column in columns[:3]]
    return (
        normalized[0] in {"chromosome", "chrom", "chr", "seqname"}
        and normalized[1] in {"start", "window_start"}
        and normalized[2] in {"end", "stop", "window_end"}
    )


def _read_log_tail(log_file, limit=MAX_COMMAND_LOG_BYTES):
    """Read at most ``limit`` bytes from the end of a binary temporary log."""
    log_file.flush()
    size = log_file.tell()
    log_file.seek(max(0, size - limit))
    data = log_file.read(limit)
    return data.decode("utf-8", errors="replace"), size > limit


def _write_bed_regions_atomic(regions, output_bed):
    """Write BED triples and atomically publish the completed file."""
    output_bed = Path(output_bed)
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    temp_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=f".{output_bed.name}.",
            suffix=".tmp",
            dir=output_bed.parent,
            delete=False,
        ) as output:
            temp_name = output.name
            for chrom, start, end in regions:
                output.write(f"{chrom}\t{start}\t{end}\n")
        os.replace(temp_name, output_bed)
        temp_name = None
    finally:
        if temp_name is not None:
            try:
                os.unlink(temp_name)
            except FileNotFoundError:
                pass


def _iter_se_anomaly_records(anomalies_file):
    """Yield normalized SE anomaly records from the public TSV schema."""
    tag_index = None
    with open(anomalies_file, "r", encoding="utf-8") as input_file:
        for line_number, raw_line in enumerate(input_file, 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if _is_chromosome_header(fields):
                normalized = [field.strip().lower() for field in fields]
                try:
                    tag_index = normalized.index("anomalytags")
                except ValueError:
                    tag_index = None
                continue
            if len(fields) < 3:
                raise ValueError(
                    f"Invalid SE anomaly row {line_number} in {anomalies_file}: {line}"
                )
            try:
                one_based_start = int(fields[1])
                closed_end = int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"Invalid SE coordinates on row {line_number} in "
                    f"{anomalies_file}: {line}"
                ) from error
            if one_based_start < 1 or closed_end < one_based_start:
                raise ValueError(
                    f"Invalid SE interval on row {line_number} in "
                    f"{anomalies_file}: {line}"
                )
            tags_field = (
                fields[tag_index]
                if tag_index is not None and tag_index < len(fields)
                else fields[-1]
            )
            tags = tuple(
                tag.strip() for tag in tags_field.split(",") if tag.strip()
            )
            yield fields[0], one_based_start - 1, closed_end, tags


def _write_merged_se_bed_from_anomalies(anomalies_file, output_bed):
    """Create merged BED intervals from a (possibly masked) SE anomaly TSV."""
    regions_by_chromosome = OrderedDict()
    for chrom, start, end, _ in _iter_se_anomaly_records(anomalies_file):
        regions_by_chromosome.setdefault(chrom, []).append((start, end))

    merged_regions = []
    for chrom, intervals in regions_by_chromosome.items():
        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        merged_regions.extend((chrom, start, end) for start, end in merged)

    _write_bed_regions_atomic(merged_regions, output_bed)


def _write_se_anomaly_summary_from_tsv(
    anomalies_file, output_summary, chromosome_order=()
):
    """Create the standard SE tag-count summary from an anomaly TSV."""
    tag_types = ("low_depth", "high_depth", "indel", "clip", "non_spanning")
    counts = defaultdict(lambda: defaultdict(int))
    observed_chromosomes = []
    observed_set = set()

    for chrom, _, _, tags in _iter_se_anomaly_records(anomalies_file):
        if chrom not in observed_set:
            observed_set.add(chrom)
            observed_chromosomes.append(chrom)
        counts[chrom]["nr-SE"] += 1
        for tag in tags:
            if tag in tag_types:
                counts[chrom][tag] += 1

    chromosomes = sorted(dict.fromkeys(chromosome_order))
    known_chromosomes = set(chromosomes)
    for chrom in sorted(observed_chromosomes):
        if chrom not in known_chromosomes:
            chromosomes.append(chrom)
            known_chromosomes.add(chrom)

    output_summary = Path(output_summary)
    output_summary.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=f".{output_summary.name}.",
            suffix=".tmp",
            dir=output_summary.parent,
            delete=False,
        ) as output:
            temporary_name = output.name
            output.write(
                "Chromosome\tlow_depth\thigh_depth\tindel\tclip\t"
                "non-spanning\tnr-SE\n"
            )
            totals = defaultdict(int)
            for chrom in chromosomes:
                values = [counts[chrom][tag] for tag in tag_types]
                nr_se = counts[chrom]["nr-SE"]
                output.write(
                    "\t".join(
                        [chrom, *[str(value) for value in values], str(nr_se)]
                    )
                    + "\n"
                )
                for tag, value in zip(tag_types, values):
                    totals[tag] += value
                totals["nr-SE"] += nr_se
            output.write(
                "\t".join(
                    [
                        "Total",
                        *[str(totals[tag]) for tag in tag_types],
                        str(totals["nr-SE"]),
                    ]
                )
                + "\n"
            )
        os.replace(temporary_name, output_summary)
        temporary_name = None
    finally:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass


def _publish_file_set(replacements, backup_prefix=".publish-backup-"):
    """Publish staged files as one rollback-safe transaction.

    POSIX does not provide one atomic rename for multiple files.  Before any
    live pathname changes, retain each existing destination in a private
    same-filesystem directory.  If any replacement fails (including after the
    rename itself completed), restore every attempted destination; paths that
    did not exist before the transaction are removed again.

    All destinations must share a parent directory.  Sources are consumed on
    success and may be left missing after rollback; callers normally keep them
    in a temporary staging directory.
    """
    replacements = [
        (Path(source), Path(destination))
        for source, destination in replacements
    ]
    if not replacements:
        return

    output_dir = replacements[0][1].parent
    if any(
        destination.parent != output_dir
        for _, destination in replacements
    ):
        raise ValueError("Transactional destinations must share a directory")
    output_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(
        prefix=backup_prefix, dir=output_dir
    ) as backup_dir_name:
        backup_dir = Path(backup_dir_name)
        backups = {}

        # Complete all fallible backup copies before changing the live set.
        for _, destination in replacements:
            if not destination.exists():
                backups[destination] = None
                continue
            backup = backup_dir / destination.name
            try:
                # A hard link is constant-space and remains valid when the
                # live pathname is atomically replaced.
                os.link(destination, backup)
            except OSError:
                shutil.copy2(destination, backup)
            backups[destination] = backup

        attempted_destinations = []
        try:
            for source, destination in replacements:
                attempted_destinations.append(destination)
                os.replace(source, destination)
        except Exception as publish_error:
            rollback_errors = []
            for destination in reversed(attempted_destinations):
                backup = backups[destination]
                try:
                    if backup is None:
                        destination.unlink(missing_ok=True)
                    else:
                        try:
                            os.replace(backup, destination)
                        except OSError:
                            # Normally rollback is another same-filesystem
                            # atomic replace.  Retain a copy fallback so a
                            # transient/injected replace failure cannot strand
                            # the newly published file in place.
                            # If publication failed before its rename landed,
                            # the hard-link backup and live destination still
                            # identify the original inode.  Copying one onto
                            # the other would raise SameFileError even though
                            # this destination is already restored.
                            if _same_existing_file(backup, destination):
                                continue
                            if backup.exists():
                                shutil.copy2(backup, destination)
                            elif not destination.exists():
                                raise
                except Exception as rollback_error:
                    rollback_errors.append(
                        f"{destination}: {rollback_error}"
                    )

            if rollback_errors:
                rollback_summary = "; ".join(rollback_errors)
                raise RuntimeError(
                    "File-set publication failed and rollback was "
                    f"incomplete: {rollback_summary}"
                ) from publish_error
            raise


def _publish_chromosome_file_set(staged_dir, output_dir, filename_prefix, chromosomes):
    """Atomically publish a complete managed chromosome-file set.

    Individual files are replaced atomically and the complete set is rolled
    back after a partial failure.  Stale managed files are deliberately
    removed only after the complete new set is in place.
    """
    staged_dir = Path(staged_dir)
    output_dir = Path(output_dir)
    chromosomes = list(chromosomes)
    destinations = [
        output_dir / f"{filename_prefix}{chrom}.txt" for chrom in chromosomes
    ]
    old_files = {
        path.name: path
        for path in output_dir.glob(f"{filename_prefix}*.txt")
        if path.is_file()
    }
    replacements = [
        (staged_dir / f"{filename_prefix}{chrom}.txt", destination)
        for chrom, destination in zip(chromosomes, destinations)
    ]
    _publish_file_set(
        replacements,
        backup_prefix=f".{filename_prefix}backup-",
    )

    new_names = {destination.name for destination in destinations}
    for stale_name, stale_file in old_files.items():
        if stale_name not in new_names:
            stale_file.unlink(missing_ok=True)


class _ChromosomeFilePool:
    """Bound the number of chromosome output files kept open at one time."""

    def __init__(self, output_dir, filename_prefix, header=None, max_open_files=None):
        self.output_dir = Path(output_dir)
        self.filename_prefix = filename_prefix
        self.header = header
        self.max_open_files = (
            MAX_OPEN_CHROMOSOME_FILES
            if max_open_files is None
            else max_open_files
        )
        if self.max_open_files < 1:
            raise ValueError("max_open_files must be at least 1")
        self._handles = OrderedDict()
        self.chromosomes = []
        self._initialized = set()

    def _open(self, chrom):
        handle = self._handles.pop(chrom, None)
        if handle is not None:
            self._handles[chrom] = handle
            return handle

        if len(self._handles) >= self.max_open_files:
            _, oldest_handle = self._handles.popitem(last=False)
            oldest_handle.close()

        first_write = chrom not in self._initialized
        output_file = self.output_dir / f"{self.filename_prefix}{chrom}.txt"
        handle = open(output_file, "w" if first_write else "a")
        try:
            if first_write and self.header is not None:
                handle.write(self.header + "\n")
        except Exception:
            handle.close()
            raise
        if first_write:
            self._initialized.add(chrom)
            self.chromosomes.append(chrom)
        self._handles[chrom] = handle
        return handle

    def write(self, chrom, line):
        self._open(chrom).write(line + "\n")

    def close(self):
        first_error = None
        while self._handles:
            _, handle = self._handles.popitem(last=False)
            try:
                handle.close()
            except Exception as error:
                if first_error is None:
                    first_error = error
        if first_error is not None:
            raise first_error

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


def get_pipeline_base_dir():
    """Get the base directory of the SAS Pipeline based on script location"""
    # Get the directory containing this script
    script_dir = Path(__file__).parent.resolve()
    
    # The script is in src/, so the base directory is the parent
    base_dir = script_dir.parent
    
    return base_dir


class BedMask:
    """Class to handle bed file masking functionality using bedtools intersect -v"""
    
    def __init__(self, bed_files, output_dir):
        self.bed_files = bed_files
        self.output_dir = Path(output_dir)
        self._mask_regions_by_chrom = self._load_mask_regions()
        self._mask_starts_by_chrom = {
            chrom: [start for start, _ in intervals]
            for chrom, intervals in self._mask_regions_by_chrom.items()
        }

    def _load_mask_regions(self):
        """Load mask BEDs once and cache their interval union by chromosome."""
        regions_by_chrom = defaultdict(list)

        for bed_file in self.bed_files:
            if not Path(bed_file).exists():
                print(f"Warning: Bed file {bed_file} does not exist, skipping")
                continue

            with open(bed_file, 'r') as bed_f:
                for line_number, line in enumerate(bed_f, 1):
                    stripped = line.strip()
                    if not stripped or stripped.startswith('#'):
                        continue
                    parts = stripped.split('\t')
                    if len(parts) < 3:
                        continue
                    try:
                        start = int(parts[1])
                        end = int(parts[2])
                    except ValueError:
                        print(
                            f"Warning: Invalid coordinates in bed file {bed_file} "
                            f"at line {line_number}: {stripped}"
                        )
                        continue
                    if end <= start:
                        print(
                            f"Warning: Invalid BED interval in {bed_file} "
                            f"at line {line_number}: {stripped}"
                        )
                        continue
                    regions_by_chrom[parts[0]].append((start, end))

        merged_by_chrom = {}
        for chrom, intervals in regions_by_chrom.items():
            merged = []
            for start, end in sorted(intervals):
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            merged_by_chrom[chrom] = tuple(merged)

        return merged_by_chrom

    def _overlaps_mask(self, chrom, start, end):
        """Return whether BED interval [start, end) overlaps a cached mask."""
        if end <= start:
            return False
        intervals = self._mask_regions_by_chrom.get(chrom)
        if not intervals:
            return False

        # The last mask beginning before query end is the only merged interval
        # that can overlap it.  This keeps each lookup O(log masks_per_chrom).
        index = bisect_left(self._mask_starts_by_chrom[chrom], end) - 1
        return index >= 0 and intervals[index][1] > start

    def mask_bed_regions(self, input_file, output_file, description=""):
        """
        Prefer bedtools intersect -v to subtract user-provided mask beds from
        the BE candidate BED produced by smallerror_short.
        """
        if not self.bed_files:
            shutil.copy(input_file, output_file)
            return
        
        input_file = Path(input_file)
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        # Try fast path with bedtools; fallback to Python implementation on error
        tmp_file = None
        try:
            cmd = [
                "bedtools", "intersect", "-v",
                "-a", str(input_file),
                "-b", *[str(b) for b in self.bed_files]
            ]
            # Use temporary file to filter empty lines
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed') as tmp_out:
                tmp_file = tmp_out.name
                subprocess.run(cmd, check=True, stdout=tmp_out)
            
            # Filter empty lines from bedtools output
            with open(tmp_file, 'r') as fin, open(output_file, 'w') as fout:
                for line in fin:
                    if line.strip():  # Only write non-empty lines
                        fout.write(line)
            
        except Exception as e:
            print(f"Warning: bedtools masking failed for {description}: {e}. Falling back to Python masking.")
            self.apply_bedtools_intersect_v(input_file, output_file, description)
        finally:
            if tmp_file is not None:
                try:
                    os.unlink(tmp_file)
                except FileNotFoundError:
                    pass
                except OSError as e:
                    print(f"Warning: Could not remove temporary BED file {tmp_file}: {e}")
    
    def apply_bedtools_intersect_v(
        self, input_file, output_file, description="", coordinate_mode="bed"
    ):
        """Apply cached BED masking with explicit input coordinate semantics.

        ``bed`` means zero-based, half-open coordinates. ``one_based_closed``
        is used for SE TSV output whose Start/End columns are one-based and
        inclusive; those values are converted to BED before overlap checks.
        """
        if coordinate_mode not in {"bed", "one_based_closed"}:
            raise ValueError(f"Unsupported coordinate mode: {coordinate_mode}")

        if not self.bed_files:
            # If no bed files, just copy the input file
            shutil.copy(input_file, output_file)
            return
        
        # Check if input file exists
        if not Path(input_file).exists():
            print(f"Warning: Input file {input_file} does not exist, skipping bed masking")
            return
        
        if not self._mask_regions_by_chrom:
            print(f"Warning: No valid mask regions found, copying original file without masking")
            shutil.copy(input_file, output_file)
            return

        region_count = sum(len(v) for v in self._mask_regions_by_chrom.values())
        print(f"Loaded {region_count} merged mask regions")
        
        # Process input file
        kept_lines = 0
        removed_lines = 0
        
        with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
            for line in fin:
                # Skip empty lines
                if not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    try:
                        chrom = fields[0]
                        start = int(fields[1])
                        end = int(fields[2])

                        if coordinate_mode == "one_based_closed":
                            start -= 1

                        overlaps = self._overlaps_mask(chrom, start, end)
                        
                        if not overlaps:
                            fout.write(line.rstrip() + '\n')  # Remove trailing whitespace
                            kept_lines += 1
                        else:
                            removed_lines += 1
                    except ValueError:
                        lower_fields = [field.strip().lower() for field in fields[:3]]
                        is_header = (
                            len(lower_fields) == 3
                            and lower_fields[0] in {"chrom", "chromosome", "chr", "seqname"}
                            and lower_fields[1] in {"start", "window_start"}
                            and lower_fields[2] in {"end", "stop", "window_end"}
                        )
                        if not is_header:
                            print(f"Warning: Invalid coordinates in line: {line.strip()}")
                        if line.strip():
                            fout.write(line.rstrip() + '\n')
                            kept_lines += 1
                else:
                    # Skip lines with insufficient fields (likely empty or malformed)
                    # Only keep if it's a comment line
                    if line.strip().startswith('#'):
                        fout.write(line.rstrip() + '\n')
                        kept_lines += 1
        
        print(f"Masking completed: kept {kept_lines} lines, removed {removed_lines} lines")
    
    def get_masked_length(self, chrom, chrom_length=None):
        """
        Calculate the total length of masked regions for a specified chromosome.
        
        Args:
            chrom: Chromosome name
            chrom_length: Optional reference length used to clip mask intervals
            
        Returns:
            int: Total length of masked regions
        """
        total_masked_length = 0
        for start, end in self._mask_regions_by_chrom.get(chrom, ()):
            clipped_start = max(0, start)
            clipped_end = end if chrom_length is None else min(chrom_length, end)
            if clipped_end > clipped_start:
                total_masked_length += clipped_end - clipped_start

        return total_masked_length


class SASPipeline:
    """Main pipeline class for SAS algorithm workflow"""
    
    def __init__(self, config):
        self.config = config
        self.output_dir = Path(config['output_dir'])
        self.skip_filter = config.get('skip_filter', False)
        self.run_be_analysis = config.get('run_be_analysis', True)
        self.run_se_analysis = config.get('run_se_analysis', True)

        # Keep user-facing deliverables separate from large/restartable work
        # products.  This makes the top-level output directory predictable and
        # lets users archive ``final`` without also copying BAMs and depth data.
        self.final_output_dir = self.output_dir / "final"
        self.work_dir = self.output_dir / "work"
        self.alignments_dir = self.work_dir / "alignments"
        self.logs_dir = self.output_dir / "logs"
        self.stage_logs_dir = self.logs_dir / "stages"
        self.metadata_dir = self.output_dir / "metadata"
        self.be_results_dir = self.final_output_dir / "be"
        self.se_results_dir = self.final_output_dir / "se"
        self.be_chromosome_output_dir = self.be_results_dir / "by_chromosome"
        self.se_chromosome_output_dir = self.se_results_dir / "by_chromosome"
        self.chromosome_summary_file = self.final_output_dir / "chromosome_qv_summary.txt"
        self.summary_file = self.final_output_dir / "sas_pipeline_summary.txt"
        self.config_file = self.metadata_dir / "pipeline_config.json"

        # Create the fixed layout first, before setting up logging.
        self.output_dir.mkdir(parents=True, exist_ok=True)
        for directory in (
            self.final_output_dir,
            self.work_dir,
            self.logs_dir,
            self.stage_logs_dir,
            self.metadata_dir,
        ):
            directory.mkdir(parents=True, exist_ok=True)
        
        # Auto-detect pipeline base directory if not provided
        if config.get('scripts_dir') == 'scripts' or config.get('src_dir') == 'src':
            pipeline_base = get_pipeline_base_dir()
            self.scripts_dir = pipeline_base / 'scripts'
            self.src_dir = pipeline_base / 'src'
        else:
            self.scripts_dir = Path(config.get('scripts_dir', 'scripts'))
            self.src_dir = Path(config.get('src_dir', 'src'))
        
        # Setup logging after creating output directory
        self.logger = self._setup_logging()
        
        # Log the directories being used
        self.logger.info(f"Pipeline base directory: {get_pipeline_base_dir()}")
        self.logger.info(f"Scripts directory: {self.scripts_dir}")
        self.logger.info(f"Source directory: {self.src_dir}")
        
        # File paths - convert to absolute paths to avoid working directory issues
        self.short_read_bam = Path(config['short_read_bam']).resolve()
        self.hifi_bam = Path(config['hifi_bam']).resolve()
        self.ont_bam = Path(config['ont_bam']).resolve()
        self.reference_fa = Path(config['reference_fa']).resolve()
        self.reference_fai = Path(config.get('reference_fai')).resolve() if config.get('reference_fai') else None
        self.threads = int(config.get('threads', 8))
        if self.threads < 1:
            raise ValueError("threads must be greater than 0")
        self.pileup_workers = int(
            config.get('pileup_workers') or min(self.threads, 16)
        )
        self.se_workers = int(
            config.get('se_workers') or min(self.threads, 8)
        )
        if self.pileup_workers < 1:
            raise ValueError("pileup_workers must be greater than 0")
        if self.se_workers < 1:
            raise ValueError("se_workers must be greater than 0")
        self.logger.info(
            "Runtime parallelism: %d tool threads, %d pileup workers, "
            "%d SE workers",
            self.threads,
            self.pileup_workers,
            self.se_workers,
        )
        
        # New features
        self.bed_files = config.get('bed_files', [])
        self.bed_mask = BedMask(self.bed_files, self.output_dir) if self.bed_files else None
        
        # Additional bed file for BE analysis (will be split into 50bp windows and merged with low-coverage regions)
        self.additional_bed_file = config.get('additional_bed_file')
        if self.additional_bed_file:
            self.additional_bed_file = Path(self.additional_bed_file).resolve()
            self.logger.info(f"Additional bed file for BE analysis: {self.additional_bed_file}")
        self.ir_bed_file = None
        self.using_expanded_candidates = False
        
        # IR (Identical Region) generation option
        self.generate_ir = config.get('generate_ir', False)
        if self.generate_ir:
            self.logger.info("Optional IR refinement enabled")
            self.logger.info(
                "IR regions will be generated with genmap and evaluated as additional BE candidates"
            )
        
        # Log the resolved paths
        self.logger.info(f"Short-read BAM: {self.short_read_bam}")
        self.logger.info(f"HiFi BAM: {self.hifi_bam}")
        self.logger.info(f"ONT BAM: {self.ont_bam}")
        self.logger.info(f"Reference FASTA: {self.reference_fa}")
        if self.reference_fai:
            self.logger.info(f"Reference FAI: {self.reference_fai}")
        
        if self.bed_files:
            self.logger.info(f"Bed masking enabled with {len(self.bed_files)} files")
            for bed_file in self.bed_files:
                self.logger.info(f"  Bed file: {bed_file}")
        
        # SE analysis parameters
        self.se_window_size = config.get('se_window_size', 2000)
        self.se_min_indel_length = config.get('se_min_indel_length', 50)
        
        # Setup reference fai file
        self._setup_reference_fai()
        
        # Get chromosome information
        self.chromosomes = self._get_chromosomes()
        
        # Output file paths - adjust based on whether filtering is skipped
        if self.skip_filter:
            # Use input BAM files directly (already absolute paths)
            self.short_filtered_bam = self.short_read_bam
            self.hifi_filtered_bam = self.hifi_bam
            self.ont_filtered_bam = self.ont_bam
        else:
            # Use filtered output files
            self.alignments_dir.mkdir(parents=True, exist_ok=True)
            self.short_filtered_bam = self.alignments_dir / "short_read_filtered.bam"
            self.hifi_filtered_bam = self.alignments_dir / "hifi_filtered.bam"
            self.ont_filtered_bam = self.alignments_dir / "ont_filtered.bam"
        
        # BE work products and stable public results.
        self.be_output_dir = self.work_dir / "be"
        if self.run_be_analysis:
            self.be_output_dir.mkdir(parents=True, exist_ok=True)
            self.be_results_dir.mkdir(parents=True, exist_ok=True)
            self.be_chromosome_output_dir.mkdir(parents=True, exist_ok=True)
        self.smallerror_bed = self.be_output_dir / "smallerror.bed"
        has_extra_candidates = bool(self.additional_bed_file or self.generate_ir)
        self.smallerror_bed_windows = self.be_output_dir / "smallerror_50bp_windows.bed" if has_extra_candidates else None
        self.additional_bed_windows = self.be_output_dir / "additional_bed_windows.bed" if self.additional_bed_file else None
        self.ir_bed_windows = self.be_output_dir / "ir_50bp_windows.bed" if self.generate_ir else None
        self.merged_50bp_windows = self.be_output_dir / "candidate_50bp_windows.bed" if has_extra_candidates else None
        self.unmasked_50bp_windows = (
            self.be_output_dir / "candidate_50bp_windows.unmasked.bed"
            if has_extra_candidates and self.bed_mask
            else None
        )
        self.hifi_filtered_regions = self.be_output_dir / "hifi_filtered_regions.txt"
        self.ont_filtered_regions = self.be_output_dir / "ont_filtered_regions.txt"
        self.merged_results = self.be_output_dir / "merged_results.txt"
        self.final_be_results = self.be_results_dir / "final_be_results.txt"
        
        # SE's verbose depth/statistics products stay in work/se.  A compact
        # public subset is published under final/se after analysis succeeds.
        self.se_output_dir = self.work_dir / "se"
        if self.run_se_analysis:
            self.se_output_dir.mkdir(parents=True, exist_ok=True)
            self.se_results_dir.mkdir(parents=True, exist_ok=True)
            self.se_chromosome_output_dir.mkdir(parents=True, exist_ok=True)
        self.se_results_prefix = self.se_output_dir / "structural_errors"
    
    def _get_chromosomes(self):
        """Get list of chromosomes from reference fai file"""
        chromosomes = []
        try:
            with open(self.reference_fai, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 1:
                        chromosomes.append(parts[0])
            self.logger.info(f"Found {len(chromosomes)} chromosomes: {chromosomes}")
        except Exception as e:
            self.logger.warning(f"Could not read chromosomes from fai file: {e}")
            chromosomes = []
        
        return chromosomes
    
    def _get_chromosome_lengths(self):
        """Get chromosome lengths from fai file"""
        chr_lengths = {}
        try:
            with open(self.reference_fai, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        chrom = parts[0]
                        length = int(parts[1])
                        chr_lengths[chrom] = length
        except Exception as e:
            self.logger.warning(f"Could not read chromosome lengths from fai file: {e}")
        
        return chr_lengths
    
    def _setup_logging(self):
        """Setup logging configuration"""
        log_file = getattr(self, "logs_dir", self.output_dir) / "sas_pipeline.log"
        
        # Clear any existing handlers
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        
        handlers = [logging.StreamHandler(sys.stdout)]
        
        # Try to create file handler, but don't fail if it's not possible
        try:
            handlers.append(logging.FileHandler(log_file))
        except (OSError, IOError, PermissionError) as e:
            print(f"Warning: Could not create log file {log_file}: {e}", file=sys.stderr)
            print("Continuing with console logging only.", file=sys.stderr)
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=handlers
        )
        return logging.getLogger(__name__)
    
    def _setup_reference_fai(self):
        """Setup reference fai file, generate if not provided"""
        canonical_fai = Path(str(self.reference_fa) + ".fai")
        if self.reference_fai is None:
            self.reference_fai = canonical_fai
            
        # Check if fai file exists
        if not self.reference_fai.exists():
            self.logger.info(f"Reference fai file not found: {self.reference_fai}")
            self.logger.info("Generating fai file using samtools faidx...")
            
            # Generate fai file using samtools faidx
            cmd = ["samtools", "faidx", str(self.reference_fa)]
            try:
                result = subprocess.run(
                    cmd,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True
                )
                self.logger.info(f"Successfully generated fai file: {self.reference_fai}")
                
            except subprocess.CalledProcessError as e:
                self.logger.error(f"Failed to generate fai file: {e}")
                self.logger.error(f"Command output: {e.stdout}")
                self.logger.error(f"Command error: {e.stderr}")
                raise
            except FileNotFoundError:
                self.logger.error("samtools not found. Please install samtools.")
                raise
        else:
            self.logger.info(f"Using existing fai file: {self.reference_fai}")

        # pysam/samtools consumers discover the index beside the FASTA even
        # when the user supplied an alternate --reference-fai path.
        if self.reference_fai != canonical_fai and not canonical_fai.exists():
            self.logger.info(
                f"Generating canonical FASTA index required by downstream tools: {canonical_fai}"
            )
            subprocess.run(
                ["samtools", "faidx", str(self.reference_fa)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

    @staticmethod
    def _read_fai_lengths(fai_file):
        lengths = OrderedDict()
        with open(fai_file, "r", encoding="utf-8") as input_handle:
            for line_number, line in enumerate(input_handle, 1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    raise ValueError(
                        f"Malformed FASTA index row {line_number}: {line.rstrip()}"
                    )
                chrom, length = fields[0], int(fields[1])
                if not chrom or length <= 0 or chrom in lengths:
                    raise ValueError(
                        f"Invalid FASTA index row {line_number}: {line.rstrip()}"
                    )
                lengths[chrom] = length
        if not lengths:
            raise ValueError(f"FASTA index contains no sequences: {fai_file}")
        return lengths

    def _validate_reference_compatibility(self):
        """Fail before analysis when BAM and reference coordinate spaces differ."""
        if pysam is None:
            raise RuntimeError("pysam is required for BAM/reference validation")

        expected = self._read_fai_lengths(self.reference_fai)
        canonical_fai = Path(str(self.reference_fa) + ".fai")
        if canonical_fai != self.reference_fai:
            canonical = self._read_fai_lengths(canonical_fai)
            if dict(canonical) != dict(expected):
                raise ValueError(
                    "Provided --reference-fai does not match the FASTA's canonical index"
                )

        for label, bam_path in (
            ("short-read", self.short_read_bam),
            ("HiFi", self.hifi_bam),
            ("ONT", self.ont_bam),
        ):
            try:
                with pysam.AlignmentFile(str(bam_path), "rb") as bam_handle:
                    observed = OrderedDict(
                        zip(bam_handle.references, bam_handle.lengths)
                    )
                    sort_order = (
                        bam_handle.header.to_dict().get("HD", {}).get("SO")
                    )
                    if sort_order and sort_order != "coordinate":
                        raise ValueError(
                            f"{label} BAM is not coordinate sorted (SO={sort_order}): {bam_path}"
                        )
                    if self.skip_filter and not bam_handle.has_index():
                        raise ValueError(
                            f"{label} BAM requires a .bai or .csi index with --skip-filter: {bam_path}"
                        )
            except (OSError, ValueError) as error:
                if isinstance(error, ValueError) and str(error).startswith(label):
                    raise
                raise ValueError(f"Cannot validate {label} BAM {bam_path}: {error}") from error

            if dict(observed) != dict(expected):
                missing = [chrom for chrom in expected if chrom not in observed]
                extra = [chrom for chrom in observed if chrom not in expected]
                length_mismatches = [
                    chrom for chrom in expected
                    if chrom in observed and expected[chrom] != observed[chrom]
                ]
                details = []
                if missing:
                    details.append("missing=" + ",".join(missing[:5]))
                if extra:
                    details.append("extra=" + ",".join(extra[:5]))
                if length_mismatches:
                    details.append(
                        "length mismatch=" + ",".join(length_mismatches[:5])
                    )
                raise ValueError(
                    f"{label} BAM reference does not match {self.reference_fai}: "
                    + "; ".join(details)
                )

        self.logger.info("BAM/FASTA reference compatibility validation passed")
    
    def _run_command(self, cmd, description, cwd=None):
        """Run a command without retaining unbounded output in memory."""
        self.logger.info(f"Starting: {description}")
        
        # Convert all command elements to strings for safe logging and execution
        cmd_str = [str(item) for item in cmd]
        self.logger.info(f"Command: {' '.join(cmd_str)}")
        
        try:
            # Child processes stream to files instead of growing pipeline
            # memory. Real pipeline runs retain the full per-stage output;
            # lightweight test/embedded callers without a log directory keep
            # using disposable files.
            with contextlib.ExitStack() as log_stack:
                stage_logs_dir = getattr(self, "stage_logs_dir", None)
                if stage_logs_dir is not None:
                    stage_logs_dir = Path(stage_logs_dir)
                    stage_logs_dir.mkdir(parents=True, exist_ok=True)
                    counter = getattr(self, "_command_log_counter", 0) + 1
                    self._command_log_counter = counter
                    log_stem = re.sub(
                        r"[^A-Za-z0-9._-]+", "_", description.strip()
                    ).strip("_.") or "command"
                    log_stem = f"{counter:03d}_{log_stem[:100]}"
                    stdout_path = stage_logs_dir / f"{log_stem}.stdout.log"
                    stderr_path = stage_logs_dir / f"{log_stem}.stderr.log"
                    stdout_log = log_stack.enter_context(open(stdout_path, "w+b"))
                    stderr_log = log_stack.enter_context(open(stderr_path, "w+b"))
                    self.logger.info(
                        "Stage logs: %s and %s", stdout_path, stderr_path
                    )
                else:
                    stdout_log = log_stack.enter_context(
                        tempfile.TemporaryFile(mode="w+b")
                    )
                    stderr_log = log_stack.enter_context(
                        tempfile.TemporaryFile(mode="w+b")
                    )
                try:
                    result = subprocess.run(
                        cmd_str,
                        check=True,
                        stdout=stdout_log,
                        stderr=stderr_log,
                        cwd=cwd,
                    )
                except subprocess.CalledProcessError as error:
                    stdout_tail, stdout_truncated = _read_log_tail(stdout_log)
                    stderr_tail, stderr_truncated = _read_log_tail(stderr_log)
                    error.stdout = stdout_tail
                    error.stderr = stderr_tail
                    self.logger.error(f"Failed: {description}")
                    self.logger.error(f"Return code: {error.returncode}")
                    if stdout_tail:
                        marker = " (tail; earlier output omitted)" if stdout_truncated else ""
                        self.logger.error(f"STDOUT{marker}: {stdout_tail.rstrip()}")
                    if stderr_tail:
                        marker = " (tail; earlier output omitted)" if stderr_truncated else ""
                        self.logger.error(f"STDERR{marker}: {stderr_tail.rstrip()}")
                    raise

                stdout_tail, stdout_truncated = _read_log_tail(stdout_log)
                stderr_tail, stderr_truncated = _read_log_tail(stderr_log)
                if stdout_tail:
                    marker = " (tail; earlier output omitted)" if stdout_truncated else ""
                    self.logger.info(f"STDOUT{marker}: {stdout_tail.rstrip()}")
                if stderr_tail:
                    marker = " (tail; earlier output omitted)" if stderr_truncated else ""
                    self.logger.warning(f"STDERR{marker}: {stderr_tail.rstrip()}")

                # Preserve the previous return shape for potential callers,
                # while keeping stdout/stderr bounded.
                result.stdout = stdout_tail
                result.stderr = stderr_tail
                self.logger.info(f"Completed: {description}")
                return result
        except FileNotFoundError as e:
            self.logger.error(f"Command not found: {description}")
            self.logger.error(f"Error: {e}")
            # Check if it's genmap
            if 'genmap' in cmd_str[0].lower():
                self.logger.error("genmap is not found in PATH. Please install genmap or add it to your PATH.")
                self.logger.error("You can install genmap using: conda install -c bioconda genmap")
                self.logger.error("Or download from: https://github.com/cpockrandt/genmap")
            raise
    
    def step1_filter_bam_files(self):
        """Step 1: Filter all three BAM files using appropriate modes"""
        self.logger.info("="*50)
        self.logger.info("STEP 1: Filtering BAM files")
        self.logger.info("="*50)
        
        if self.skip_filter:
            self.logger.info("Skipping BAM filtering - input files are already filtered")
            self.logger.info("Using required input BAM files directly:")
            if self.run_be_analysis:
                self.logger.info(f"  Short-read BAM: {self.short_read_bam}")
                self.logger.info(f"  HiFi BAM: {self.hifi_bam}")
            self.logger.info(f"  ONT BAM: {self.ont_bam}")
        else:
            if self.run_be_analysis:
                # Short-read and HiFi BAMs are only consumed by BE analysis.
                cmd = [
                    "bash",
                    str(self.scripts_dir / "filter_bam.sh"),
                    "short-read",
                    "-i", str(self.short_read_bam),
                    "-o", str(self.short_filtered_bam),
                    "-t", str(self.threads),
                    "--quiet",
                ]
                self._run_command(cmd, "Filtering short-read BAM")

                cmd = [
                    "bash",
                    str(self.scripts_dir / "filter_bam.sh"),
                    "long-read",
                    "-i", str(self.hifi_bam),
                    "-o", str(self.hifi_filtered_bam),
                    "-t", str(self.threads),
                    "--quiet",
                ]
                self._run_command(cmd, "Filtering HiFi BAM")
            
            # Filter ONT BAM (split alignment filtering)
            cmd = [
                "bash",
                str(self.scripts_dir / "filter_bam.sh"),
                "long-read",
                "-i", str(self.ont_bam),
                "-o", str(self.ont_filtered_bam),
                "-t", str(self.threads),
                "--quiet",
            ]
            self._run_command(cmd, "Filtering ONT BAM")
            
            self.logger.info("Step 1 completed: Required BAM files filtered successfully")
        
        self.logger.info("Step 1 completed: BAM file preparation finished")

    def _get_additional_candidate_sources(self):
        """Return every non-short-read BE source without dropping provenance."""
        sources = []
        additional_bed = getattr(self, "additional_bed_file", None)
        if additional_bed:
            window_bed = getattr(self, "additional_bed_windows", None)
            if window_bed is None:
                window_bed = self.be_output_dir / "additional_bed_windows.bed"
                self.additional_bed_windows = window_bed
            sources.append(("additional BED", Path(additional_bed), window_bed))

        ir_bed = getattr(self, "ir_bed_file", None)
        if ir_bed:
            window_bed = getattr(self, "ir_bed_windows", None)
            if window_bed is None:
                window_bed = self.be_output_dir / "ir_50bp_windows.bed"
                self.ir_bed_windows = window_bed
            sources.append(("generated IR", Path(ir_bed), window_bed))
        return sources

    def step2_be_analysis(self):
        """Step 2: Run BE (Base Error) analysis"""
        if not self.run_be_analysis:
            self.logger.info("Skipping BE analysis as requested")
            return
            
        self.logger.info("="*50)
        self.logger.info("STEP 2: Base Error (BE) Analysis")
        self.logger.info("="*50)
        
        # Step 2a: Detect small errors from short-read data
        self.logger.info("Step 2a: Detecting small error regions from short-read data")
        output_prefix = str(self.be_output_dir / "smallerror_short")
        cmd = [
            "bash",
            str(self.scripts_dir / "smallerror_short.sh"),
            "-i", str(self.short_filtered_bam),
            "-o", output_prefix,
            "-t", str(self.threads)
        ]
        self._run_command(cmd, "Detecting small error regions from short-read data")
        
        # The output will be {output_prefix}.smallerror.bed
        expected_bed = Path(f"{output_prefix}.smallerror.bed")
        if expected_bed.exists():
            shutil.copy(expected_bed, self.smallerror_bed)
            self.logger.info(f"Small error BED file created: {self.smallerror_bed}")
        else:
            raise FileNotFoundError(f"Expected smallerror BED file not found: {expected_bed}")
        
        # Apply bed masking if enabled - MUST be done before long platform analysis
        if self.bed_mask:
            self.logger.info("Applying bed masking to small error regions (before long platform analysis)")
            masked_smallerror_bed = self.be_output_dir / "smallerror_masked.bed"
            try:
                self.bed_mask.mask_bed_regions(
                    self.smallerror_bed, 
                    masked_smallerror_bed, 
                    "small error regions"
                )
                # Replace the original smallerror.bed with masked version for all subsequent processing
                self.smallerror_bed = masked_smallerror_bed
                self.logger.info(f"Bed masking completed. Using masked BED file for long platform analysis: {self.smallerror_bed}")
            except Exception as e:
                self.logger.warning(f"Bed masking failed for small error regions: {e}")
                self.logger.warning("Using original small error BED file without masking")
        
        # Step 2a-extra: Generate IR bed file if requested
        if self.generate_ir:
            self.logger.info("Step 2a-extra: Generating IR (Identical Region) bed file")
            self.logger.info(
                "IR refinement adds genmap preprocessing and extra BE candidate evaluation; "
                "runtime depends on the assembly and available hardware."
            )
            ir_bed_file = self._generate_ir_bed()
            if ir_bed_file and ir_bed_file.exists():
                # Preserve a user-supplied -a/--additional-bed-file. Generated
                # IR is a second candidate source, not a replacement for it.
                self.ir_bed_file = ir_bed_file
                self.smallerror_bed_windows = self.be_output_dir / "smallerror_50bp_windows.bed"
                self.ir_bed_windows = self.be_output_dir / "ir_50bp_windows.bed"
                self.merged_50bp_windows = self.be_output_dir / "candidate_50bp_windows.bed"
                self.logger.info(f"IR bed file generated: {ir_bed_file}")
            else:
                raise RuntimeError("Failed to generate IR bed file")

        # Step 2a-extra: Process every additional candidate source. A user BED
        # and generated IR may be enabled together and are both retained.
        additional_sources = self._get_additional_candidate_sources()
        for source_name, source_bed, _ in additional_sources:
            if not source_bed.exists():
                raise FileNotFoundError(
                    f"{source_name} file not found: {source_bed}"
                )

        if additional_sources:
            self.logger.info(
                "Step 2a-extra: Processing candidate sources: "
                + ", ".join(source_name for source_name, _, _ in additional_sources)
            )
            self.smallerror_bed_windows = (
                self.smallerror_bed_windows
                or self.be_output_dir / "smallerror_50bp_windows.bed"
            )
            self.merged_50bp_windows = (
                self.merged_50bp_windows
                or self.be_output_dir / "candidate_50bp_windows.bed"
            )
            
            # Split smallerror_bed into 50bp windows
            self.logger.info("Splitting smallerror_bed into 50bp windows")
            self._split_bed_to_windows(
                self.smallerror_bed,
                self.smallerror_bed_windows,
                window_size=50
            )
            
            candidate_window_files = [self.smallerror_bed_windows]
            for source_name, source_bed, window_bed in additional_sources:
                self.logger.info(f"Splitting {source_name} into 50bp windows")
                self._split_bed_to_windows(
                    source_bed,
                    window_bed,
                    window_size=50,
                )
                candidate_window_files.append(window_bed)

            # Combine both candidate sources while preserving every 50 bp
            # analysis unit.  Query batching happens later in
            # smallerror_long.py and must not destroy window boundaries.
            self.logger.info(
                "Combining small-error and additional 50 bp windows "
                "(sorting and removing exact duplicates)"
            )
            combine_output = self.merged_50bp_windows
            if self.bed_mask:
                self.unmasked_50bp_windows = (
                    self.unmasked_50bp_windows
                    or self.be_output_dir / "candidate_50bp_windows.unmasked.bed"
                )
                combine_output = self.unmasked_50bp_windows
            self._combine_50bp_windows(candidate_window_files, combine_output)
            if self.bed_mask:
                self.logger.info(
                    "Removing masked intervals from expanded candidates before pileup"
                )
                self.bed_mask.apply_bedtools_intersect_v(
                    combine_output,
                    self.merged_50bp_windows,
                    "expanded BE candidates",
                )
            self.using_expanded_candidates = True
            
            # Step 2b: Batch-query independent 50 bp windows with HiFi data.
            self.logger.info("Step 2b: Processing candidate 50 bp windows with HiFi data")
            cmd = [
                sys.executable,
                str(self.scripts_dir / "smallerror_long.py"),
                str(self.merged_50bp_windows),
                str(self.hifi_filtered_bam),
                str(self.reference_fa),
                str(self.hifi_filtered_regions),
                "--processes", str(self.pileup_workers)
            ]
            self._run_command(cmd, "Processing candidate 50 bp windows with HiFi data")
            
            # Step 2c: Batch-query the same windows with ONT data.
            self.logger.info("Step 2c: Processing candidate 50 bp windows with ONT data")
            cmd = [
                sys.executable,
                str(self.scripts_dir / "smallerror_long.py"),
                str(self.merged_50bp_windows),
                str(self.ont_filtered_bam),
                str(self.reference_fa),
                str(self.ont_filtered_regions),
                "--processes", str(self.pileup_workers)
            ]
            self._run_command(cmd, "Processing candidate 50 bp windows with ONT data")
        else:
            self.using_expanded_candidates = False
            # Step 2b: Process with HiFi data (using smallerror_long.py) - only when no additional bed file
            # Note: self.smallerror_bed is already masked if bed_mask is enabled
            self.logger.info("Step 2b: Processing candidate regions with HiFi data")
            if self.bed_mask:
                self.logger.info(f"Using masked BED file for HiFi analysis: {self.smallerror_bed}")
            cmd = [
                sys.executable,
                str(self.scripts_dir / "smallerror_long.py"),
                str(self.smallerror_bed),
                str(self.hifi_filtered_bam),
                str(self.reference_fa),
                str(self.hifi_filtered_regions),
                "--processes", str(self.pileup_workers)
            ]
            self._run_command(cmd, "Processing with HiFi data")
            
            # Step 2c: Process with ONT data (using smallerror_long.py) - only when no additional bed file
            # Note: self.smallerror_bed is already masked if bed_mask is enabled
            self.logger.info("Step 2c: Processing candidate regions with ONT data")
            if self.bed_mask:
                self.logger.info(f"Using masked BED file for ONT analysis: {self.smallerror_bed}")
            cmd = [
                sys.executable,
                str(self.scripts_dir / "smallerror_long.py"),
                str(self.smallerror_bed),
                str(self.ont_filtered_bam),
                str(self.reference_fa),
                str(self.ont_filtered_regions),
                "--processes", str(self.pileup_workers)
            ]
            self._run_command(cmd, "Processing with ONT data")
        
        # Step 2d: Sort BED file and merge HiFi and ONT results
        self.logger.info("Step 2d: Sorting BED file and merging HiFi and ONT analysis results")
        # Expanded candidate windows are already sorted in reference order.
        if not self.using_expanded_candidates:
            self._sort_bed_file()
        self._merge_hifi_ont_results()
        
        # Step 2e: Final BE filtering
        self.logger.info("Step 2e: Final BE filtering")
        raw_final_be_results = (
            self.be_output_dir / "final_be_results.unmasked.txt"
            if self.bed_mask
            else self.final_be_results
        )
        cmd = [
            sys.executable,
            str(self.scripts_dir / "filter_be.py"),
            str(self.merged_results),
            str(self.reference_fai),
            "-o", str(raw_final_be_results)
        ]
        self._run_command(cmd, "Final BE filtering")
        
        # Apply bed masking to final BE results if enabled
        if self.bed_mask:
            self.logger.info("Applying bed masking to final BE results")
            try:
                self.bed_mask.apply_bedtools_intersect_v(
                    raw_final_be_results,
                    self.final_be_results,
                    "final BE results"
                )
            except Exception as e:
                self.logger.warning(f"Bed masking failed for final BE results: {e}")
                self.logger.warning("Using original BE results file without masking")
                shutil.copy2(raw_final_be_results, self.final_be_results)
        
        # Step 2f: Split BE results by chromosome
        self.logger.info("Step 2f: Splitting BE results by chromosome")
        self._split_be_results_by_chromosome()
        
        self.logger.info("BE analysis completed successfully")

    def _merge_hifi_ont_results(self):
        """
        Merge HiFi and ONT results using merge_pileup.py script.
        """
        self.logger.info("Merging HiFi and ONT results...")
        
        # Use the sorted, deduplicated 50 bp candidate set when an additional
        # BED/IR source was supplied; every row remains an analysis unit.
        bed_file_for_merge = (
            self.merged_50bp_windows
            if self.using_expanded_candidates
            else self.smallerror_bed
        )
        
        # Use the existing merge_pileup.py script
        cmd = [
            sys.executable,
            str(self.scripts_dir / "merge_pileup.py"),
            str(bed_file_for_merge),  # BED file with coordinates (50bp windows or smallerror_bed)
            str(self.hifi_filtered_regions),  # File1 (HiFi results)
            str(self.ont_filtered_regions),   # File2 (ONT results)
            str(self.reference_fa),           # Reference FASTA file
            str(self.merged_results)          # Output file
        ]
        
        self._run_command(cmd, "Merging HiFi and ONT results")
        
        self.logger.info(f"Merged results saved to: {self.merged_results}")

    def _extract_haplotypes_from_diploid(self, haplotype_dir=None):
        """
        Extract maternal and paternal haplotypes from a diploid reference genome.
        
        Automatically detects the two most common chromosome suffix patterns (e.g., _MATERNAL/_PATERNAL,
        _hap1/_hap2, etc.) and extracts chromosomes with these suffixes as the two haplotypes.
        
        Returns:
            tuple: (mat_fasta_path, pat_fasta_path) or (None, None) if extraction fails
        """
        self.logger.info("Attempting to extract haplotypes from diploid reference genome...")
        
        # Get chromosomes from fai file
        if not self.reference_fai.exists():
            self.logger.warning("Reference FAI file not found, cannot extract haplotypes")
            return None, None
        
        # Read all chromosomes and analyze suffix patterns
        all_chroms = []
        suffix_patterns = {}
        no_suffix_chroms = []  # Chromosomes without underscore (e.g., chrEBV, chrM)
        
        with open(self.reference_fai, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 1:
                    chrom = parts[0]
                    all_chroms.append(chrom)
                    
                    # Extract suffix (part after the last underscore)
                    if '_' in chrom:
                        suffix = chrom.split('_')[-1]
                        if suffix not in suffix_patterns:
                            suffix_patterns[suffix] = []
                        suffix_patterns[suffix].append(chrom)
                    else:
                        # Chromosome without underscore suffix
                        no_suffix_chroms.append(chrom)
        
        self.logger.info(f"Found {len(all_chroms)} chromosomes in total")
        if no_suffix_chroms:
            self.logger.info(f"Chromosomes without suffix (will be ignored): {no_suffix_chroms}")
        self.logger.info(f"Detected suffix patterns: {list(suffix_patterns.keys())}")
        
        # Find the two most common suffix patterns (excluding chromosomes without suffix)
        # Sort by count and take top 2
        sorted_suffixes = sorted(suffix_patterns.items(), key=lambda x: len(x[1]), reverse=True)
        
        if len(sorted_suffixes) < 2:
            self.logger.warning(
                f"Need at least 2 different suffix patterns to extract haplotypes. "
                f"Found {len(sorted_suffixes)} pattern(s): {[s[0] for s in sorted_suffixes]}"
            )
            if no_suffix_chroms:
                self.logger.warning(f"Also found {len(no_suffix_chroms)} chromosomes without suffix that were ignored")
            return None, None
        
        # Get the two most common suffix patterns
        suffix1, chroms1 = sorted_suffixes[0]
        suffix2, chroms2 = sorted_suffixes[1]
        
        self.logger.info(
            f"Using suffix patterns: '{suffix1}' ({len(chroms1)} chromosomes) and "
            f"'{suffix2}' ({len(chroms2)} chromosomes)"
        )
        
        # If there are more than 2 patterns, log a warning
        if len(sorted_suffixes) > 2:
            other_suffixes = [s[0] for s in sorted_suffixes[2:]]
            self.logger.warning(
                f"Found {len(sorted_suffixes)} suffix patterns. Using the two most common ones. "
                f"Ignoring other patterns: {other_suffixes}"
            )
        
        # Assign the two haplotypes (order doesn't matter for IR generation)
        mat_chroms = chroms1
        pat_chroms = chroms2
        
        # Keep generated haplotypes inside the current IR run when supplied so
        # concurrent runs cannot overwrite one another's inputs or indexes.
        if haplotype_dir is None:
            haplotype_dir = self.be_output_dir / "ir_generation" / "extracted_haplotypes"
        haplotype_dir = Path(haplotype_dir)
        haplotype_dir.mkdir(parents=True, exist_ok=True)
        
        # Extract maternal haplotype
        mat_fasta = haplotype_dir / "mat_haplotype.fasta"
        self.logger.info(f"Extracting maternal haplotype to {mat_fasta}...")
        self._extract_chromosomes(self.reference_fa, mat_chroms, mat_fasta)
        
        # Extract paternal haplotype
        pat_fasta = haplotype_dir / "pat_haplotype.fasta"
        self.logger.info(f"Extracting paternal haplotype to {pat_fasta}...")
        self._extract_chromosomes(self.reference_fa, pat_chroms, pat_fasta)
        
        self.logger.info("Successfully extracted haplotypes from diploid reference genome")
        return mat_fasta, pat_fasta
    
    def _extract_chromosomes(self, input_fasta, chromosomes, output_fasta):
        """
        Extract specified chromosomes from a FASTA file.
        
        Args:
            input_fasta: Path to input FASTA file
            chromosomes: List of chromosome names to extract
            output_fasta: Path to output FASTA file
        """
        self.logger.info(f"Extracting {len(chromosomes)} chromosomes from {input_fasta}...")
        
        output_fasta = Path(output_fasta)
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f".{output_fasta.name}.extract-", dir=output_fasta.parent
        ) as temp_dir:
            staged_fasta = Path(temp_dir) / output_fasta.name

            # Send faidx stdout directly to disk; chromosome sequence content
            # is never materialized in Python memory.
            with open(staged_fasta, "w") as fout:
                for chrom in chromosomes:
                    fout.flush()
                    chromosome_offset = fout.tell()
                    cmd = ["samtools", "faidx", str(input_fasta), chrom]
                    with tempfile.TemporaryFile(mode="w+b") as stderr_log:
                        try:
                            subprocess.run(
                                cmd,
                                check=True,
                                stdout=fout,
                                stderr=stderr_log,
                            )
                        except subprocess.CalledProcessError:
                            # Remove any partial stdout for clearer diagnostics;
                            # the staged directory is then discarded and the
                            # previously published FASTA/index remain untouched.
                            fout.seek(chromosome_offset)
                            fout.truncate()
                            stderr_tail, truncated = _read_log_tail(stderr_log)
                            marker = " (tail)" if truncated else ""
                            self.logger.warning(
                                f"Failed to extract chromosome {chrom}{marker}: "
                                f"{stderr_tail.rstrip()}"
                            )
                            raise

            # Index the staged FASTA first so failed indexing cannot publish a
            # partial replacement of the previous haplotype.
            cmd_index = ["samtools", "faidx", str(staged_fasta)]
            self._run_command(cmd_index, f"Indexing extracted FASTA {output_fasta.name}")
            staged_fai = Path(str(staged_fasta) + ".fai")
            if not staged_fai.exists():
                raise FileNotFoundError(
                    f"samtools faidx did not create expected index: {staged_fai}"
                )
            _publish_file_set(
                (
                    (staged_fasta, output_fasta),
                    (staged_fai, Path(str(output_fasta) + ".fai")),
                ),
                backup_prefix=f".{output_fasta.name}.publish-backup-",
            )
        
        self.logger.info(f"Extracted chromosomes saved to {output_fasta}")

    def _generate_ir_bed(self):
        """Generate IR data in an isolated, automatically cleaned run directory."""
        ir_root = self.be_output_dir / "ir_generation"
        ir_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="run-", dir=ir_root) as run_dir:
            return self._generate_ir_bed_in_work_dir(Path(run_dir))

    def _generate_ir_bed_in_work_dir(self, ir_work_dir):
        """
        Generate IR (Identical Region) bed file using genmap on parental haplotypes.
        
        Process:
        1. Split reference into parental haplotypes (mat and pat)
        2. For each haplotype:
           - Run genmap index
           - Run genmap map
           - Process with awk
           - Merge with bedtools
        3. Merge both haplotypes' results into final IR.bed
        
        Returns:
            Path to generated IR.bed file
        """
        self.logger.info("Starting IR (Identical Region) generation...")
        self.logger.info(
            "IR generation runs an additional genmap indexing/mapping phase; "
            "runtime varies with assembly size, storage, and CPU resources."
        )
        
        # Check if genmap is available by trying to run it
        # Use both shutil.which and direct command execution to handle conda environments
        genmap_path = shutil.which("genmap")
        
        # If not found in PATH, try to run genmap directly (it might be in conda environment)
        if not genmap_path:
            try:
                # Try to run genmap --version to check if it's available
                result = subprocess.run(
                    ["genmap", "--version"],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                # If it runs successfully, use "genmap" as the command
                genmap_path = "genmap"
                self.logger.info("genmap found (via direct execution)")
            except (FileNotFoundError, subprocess.TimeoutExpired, subprocess.CalledProcessError):
                # Try alternative: check if genmap is in common conda locations
                import os
                conda_prefix = os.environ.get("CONDA_PREFIX")
                if conda_prefix:
                    possible_paths = [
                        os.path.join(conda_prefix, "bin", "genmap"),
                        os.path.join(conda_prefix, "bin", "GenMap"),
                    ]
                    for path in possible_paths:
                        if os.path.exists(path) and os.access(path, os.X_OK):
                            genmap_path = path
                            self.logger.info(f"Found genmap at: {genmap_path}")
                            break
                
                if not genmap_path:
                    error_msg = (
                        "genmap is not found. Please ensure:\n"
                        "  1. genmap is installed: conda install -c bioconda genmap\n"
                        "  2. The correct conda environment is activated (the one containing genmap)\n"
                        "  3. genmap is in your PATH\n"
                        "Note: If you activated a conda environment with genmap, make sure to activate it "
                        "before running the pipeline."
                    )
                    self.logger.error(error_msg)
                    raise FileNotFoundError("genmap command not found. Please install genmap or activate the correct conda environment.")
        
        # Store genmap path as instance variable for use in _process_haplotype_for_ir
        self._genmap_path = genmap_path
        self.logger.info(f"Using genmap: {genmap_path}")
        
        # Try to infer parental haplotype files from reference name
        # Common patterns: hg002v0.9.mat-MT-EBV.fasta, hg002v0.9.pat-MT-EBV.fasta
        ref_path = Path(self.reference_fa)
        ref_dir = ref_path.parent
        ref_name = ref_path.stem  # without .fasta extension
        ref_base = ref_name.replace('.fasta', '').replace('.fa', '')
        
        # Try to find mat and pat files
        mat_patterns = ['*mat*.fasta', '*mat*.fa', '*mat-MT-EBV.fasta', '*mat-MT-EBV.fa']
        pat_patterns = ['*pat*.fasta', '*pat*.fa', '*pat-MT-EBV.fasta', '*pat-MT-EBV.fa']
        
        mat_fasta = None
        pat_fasta = None
        
        # Search for mat file
        for pattern in mat_patterns:
            matches = list(ref_dir.glob(pattern))
            if matches:
                mat_fasta = matches[0]
                break
        
        # Search for pat file
        for pattern in pat_patterns:
            matches = list(ref_dir.glob(pattern))
            if matches:
                pat_fasta = matches[0]
                break
        
        # If not found, try to construct from reference name
        if not mat_fasta or not pat_fasta:
            # Try common naming patterns
            possible_mat = ref_dir / f"{ref_base.replace('pat', 'mat')}.fasta"
            possible_pat = ref_dir / f"{ref_base.replace('mat', 'pat')}.fasta"
            
            if not mat_fasta and possible_mat.exists():
                mat_fasta = possible_mat
            if not pat_fasta and possible_pat.exists():
                pat_fasta = possible_pat
        
        # If still not found, try to extract from reference name directly
        if not mat_fasta or not pat_fasta:
            # Try to find files with mat/pat in the same directory
            all_fastas = list(ref_dir.glob("*.fasta")) + list(ref_dir.glob("*.fa"))
            for fa in all_fastas:
                if 'mat' in fa.name.lower() and not mat_fasta:
                    mat_fasta = fa
                if 'pat' in fa.name.lower() and not pat_fasta:
                    pat_fasta = fa
        
        # If not found, try to extract from diploid reference genome
        if not mat_fasta or not pat_fasta:
            self.logger.info("Parental haplotype files not found. Attempting to extract from diploid reference genome...")
            mat_fasta, pat_fasta = self._extract_haplotypes_from_diploid(
                ir_work_dir / "extracted_haplotypes"
            )
        
        if not mat_fasta or not pat_fasta:
            raise FileNotFoundError(
                f"Could not find or extract parental haplotype files. "
                f"Please ensure:\n"
                f"  1. Mat and pat FASTA files are in the same directory as the reference, OR\n"
                f"  2. The diploid reference contains chromosomes with '_MATERNAL' and '_PATERNAL' suffixes\n"
                f"Reference: {self.reference_fa}, "
                f"Directory: {ref_dir}"
            )
        
        self.logger.info(f"Using maternal haplotype: {mat_fasta}")
        self.logger.info(f"Using paternal haplotype: {pat_fasta}")
        
        # Process each haplotype
        mat_ir_bed = self._process_haplotype_for_ir(mat_fasta, "mat", ir_work_dir)
        pat_ir_bed = self._process_haplotype_for_ir(pat_fasta, "pat", ir_work_dir)
        
        # Merge both haplotypes' IR beds
        final_ir_bed = getattr(
            self, "be_results_dir", self.be_output_dir
        ) / "IR.bed"
        self.logger.info("Merging maternal and paternal IR beds...")
        
        try:
            # Combine both beds
            combined_bed = ir_work_dir / "combined_ir.tmp"
            with open(combined_bed, 'w') as fout:
                for bed_file in [mat_ir_bed, pat_ir_bed]:
                    if bed_file and bed_file.exists():
                        with open(bed_file, 'r') as fin:
                            for line in fin:
                                if line.strip() and not line.strip().startswith('#'):
                                    fout.write(line)
            
            # Sort and merge
            sorted_bed = ir_work_dir / "sorted_ir.tmp"
            cmd_sort = [
                "bedtools", "sort",
                "-i", str(combined_bed),
                "-faidx", str(self.reference_fai)
            ]
            with open(sorted_bed, 'w') as out_f, tempfile.TemporaryFile(
                mode="w+b"
            ) as stderr_log:
                subprocess.run(cmd_sort, check=True, stdout=out_f, stderr=stderr_log)
            
            # Merge overlapping regions
            cmd_merge = [
                "bedtools", "merge",
                "-i", str(sorted_bed)
            ]
            staged_ir_bed = ir_work_dir / "IR.bed.publish"
            with open(staged_ir_bed, 'w') as out_f, tempfile.TemporaryFile(
                mode="w+b"
            ) as stderr_log:
                subprocess.run(cmd_merge, check=True, stdout=out_f, stderr=stderr_log)
            os.replace(staged_ir_bed, final_ir_bed)
            
            self.logger.info(f"IR bed file generated successfully: {final_ir_bed}")
            return final_ir_bed
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to merge IR beds: {e}")
            raise
        except FileNotFoundError:
            self.logger.error("bedtools not found. Please install bedtools.")
            raise
    
    def _process_haplotype_for_ir(self, fasta_file, haplotype_name, work_dir):
        """
        Process a single haplotype to generate IR bed file.
        
        Args:
            fasta_file: Path to haplotype FASTA file
            haplotype_name: Name of haplotype (e.g., 'mat', 'pat')
            work_dir: Working directory for intermediate files
        
        Returns:
            Path to generated IR bed file for this haplotype
        """
        self.logger.info(f"Processing {haplotype_name} haplotype for IR generation...")
        
        fasta_file = Path(fasta_file)
        if not fasta_file.exists():
            raise FileNotFoundError(f"Haplotype FASTA file not found: {fasta_file}")
        
        # Create haplotype-specific working directory
        haplo_dir = work_dir / haplotype_name
        haplo_dir.mkdir(exist_ok=True)
        
        # Step 1: Generate chromosome sizes file
        chrom_sizes_file = haplo_dir / f"{haplotype_name}.chrom.sizes"
        self.logger.info(f"Generating chromosome sizes for {haplotype_name}...")
        cmd_faidx = ["samtools", "faidx", str(fasta_file)]
        subprocess.run(cmd_faidx, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Extract chromosome sizes from .fai file
        fai_file = Path(str(fasta_file) + ".fai")
        if not fai_file.exists():
            raise FileNotFoundError(f"FASTA index not found: {fai_file}")
        
        with open(fai_file, 'r') as fin, open(chrom_sizes_file, 'w') as fout:
            for line in fin:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chrom = parts[0]
                    length = parts[1]
                    fout.write(f"{chrom}\t{length}\n")
        
        # Step 2: Run genmap index
        index_name = f"{haplotype_name}-index"
        index_dir = haplo_dir / index_name
        
        # genmap index will create the directory itself, so we should not create it beforehand
        # If it already exists (from a previous failed run), remove it first
        if index_dir.exists():
            self.logger.info(f"Removing existing index directory: {index_dir}")
            shutil.rmtree(index_dir)
        
        self.logger.info(f"Running genmap index for {haplotype_name}...")
        # Get genmap path from the instance (set during _generate_ir_bed)
        genmap_cmd = getattr(self, '_genmap_path', 'genmap')
        
        cmd_index = [
            genmap_cmd, "index",
            "-F", str(fasta_file),
            "-I", str(index_dir)
        ]
        self._run_command(cmd_index, f"genmap index for {haplotype_name}")
        
        # Step 3: Run genmap map
        genmap_output_dir = haplo_dir / "genmap_output"
        genmap_output_dir.mkdir(exist_ok=True)
        
        self.logger.info(f"Running genmap map for {haplotype_name}...")
        # Get genmap path from the instance (set during _generate_ir_bed)
        genmap_cmd = getattr(self, '_genmap_path', 'genmap')
        
        cmd_map = [
            genmap_cmd, "map",
            "-K", "150",
            "-E", "0",
            "-I", str(index_dir),
            "-O", str(genmap_output_dir),
            "-t", "-w", "-bg"
        ]
        self._run_command(cmd_map, f"genmap map for {haplotype_name}")
        
        # Find the bedgraph file (genmap map generates .bedgraph file)
        bedgraph_files = list(genmap_output_dir.glob("*.bedgraph"))
        if not bedgraph_files:
            raise FileNotFoundError(f"No bedgraph file found in {genmap_output_dir}")
        bedgraph_file = bedgraph_files[0]
        
        # Step 4: Process with awk
        # awk 'NR==FNR{chr[$1]=$2; next} $4!=1{end=($3+150>chr[$1]?chr[$1]:$3+150); if(end>$2) print $1"\t"$2"\t"end}' chrom.sizes bedgraph | bedtools merge
        self.logger.info(f"Processing {haplotype_name} bedgraph with awk...")
        awk_processed = haplo_dir / f"{haplotype_name}_awk_processed.bed"
        
        # Run awk command
        awk_cmd = [
            "awk",
            'NR==FNR{chr[$1]=$2; next} $4!=1{end=($3+150>chr[$1]?chr[$1]:$3+150); if(end>$2) print $1"\t"$2"\t"end}',
            str(chrom_sizes_file),
            str(bedgraph_file)
        ]
        
        with open(awk_processed, 'w') as fout:
            result = subprocess.run(
                awk_cmd,
                check=True,
                stdout=fout,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )
        
        # Step 5: Merge with bedtools
        ir_bed = haplo_dir / f"{haplotype_name}_IR.bed"
        self.logger.info(f"Merging {haplotype_name} IR regions with bedtools...")
        
        cmd_merge = [
            "bedtools", "merge",
            "-i", str(awk_processed)
        ]
        
        with open(ir_bed, 'w') as fout:
            subprocess.run(cmd_merge, check=True, stdout=fout, stderr=subprocess.PIPE)
        
        self.logger.info(f"{haplotype_name} IR bed generated: {ir_bed}")
        return ir_bed

    def _split_bed_to_windows(self, input_bed, output_bed, window_size=50):
        """
        Split a BED file into windows of specified size (default 50bp).
        Windows are aligned to window_size boundaries (e.g., 0, 50, 100, 150...).
        Ensures windows don't exceed chromosome boundaries.
        
        Example: region 78-128 becomes windows 50-100 and 100-150
        
        Args:
            input_bed: Path to input BED file
            output_bed: Path to output BED file with windows
            window_size: Size of windows in bp (default: 50)
        """
        self.logger.info(f"Splitting BED file into {window_size}bp windows (aligned to {window_size}bp boundaries): {input_bed}")
        
        # Get chromosome lengths
        chr_lengths = self._get_chromosome_lengths()
        
        total_regions = 0
        total_windows = 0
        output_bed = Path(output_bed)
        output_bed.parent.mkdir(parents=True, exist_ok=True)
        staged_output = None

        try:
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                newline="",
                prefix=f".{output_bed.name}.",
                suffix=".tmp",
                dir=output_bed.parent,
                delete=False,
            ) as output_handle:
                staged_output = Path(output_handle.name)
                with open(input_bed, 'r', encoding="utf-8") as f:
                    for line_number, line in enumerate(f, 1):
                        line = line.strip()
                        if not line or line.startswith('#'):
                            continue

                        parts = line.split('\t')
                        if len(parts) < 3:
                            raise ValueError(
                                f"Malformed BED row {line_number} in {input_bed}: {line}"
                            )

                        chrom = parts[0]
                        try:
                            start = int(parts[1])
                            end = int(parts[2])
                        except ValueError as error:
                            raise ValueError(
                                f"Invalid BED coordinates on row {line_number} in "
                                f"{input_bed}: {line}"
                            ) from error
                        if start < 0 or end <= start:
                            raise ValueError(
                                f"Invalid BED interval on row {line_number} in "
                                f"{input_bed}: {line}"
                            )

                        # Get chromosome length
                        if chrom not in chr_lengths:
                            raise ValueError(
                                f"BED chromosome is absent from the reference index on "
                                f"row {line_number}: {chrom}"
                            )
                        chr_len = chr_lengths[chrom]

                        # Ensure end doesn't exceed chromosome length
                        end = min(end, chr_len)

                        if start >= end:
                            raise ValueError(
                                f"BED interval starts beyond chromosome {chrom} "
                                f"length {chr_len} on row {line_number}"
                            )

                        total_regions += 1

                        # Align start to window_size boundary (round down)
                        aligned_start = (start // window_size) * window_size

                        # Align end to window_size boundary (round up)
                        aligned_end = ((end + window_size - 1) // window_size) * window_size

                        # Ensure aligned_end doesn't exceed chromosome length
                        aligned_end = min(aligned_end, chr_len)

                        # Stream aligned windows directly to disk.  Dense IR
                        # regions can contain millions of windows, so retaining
                        # every tuple in a Python list needlessly consumes RAM.
                        current_start = aligned_start
                        while current_start < aligned_end:
                            current_end = min(current_start + window_size, aligned_end, chr_len)
                            if current_end > current_start:
                                output_handle.write(
                                    f"{chrom}\t{current_start}\t{current_end}\n"
                                )
                                total_windows += 1
                            current_start = current_end

            os.replace(staged_output, output_bed)
            staged_output = None

            self.logger.info(f"Split {total_regions} regions into {total_windows} {window_size}bp windows (aligned to {window_size}bp boundaries)")
            self.logger.info(f"Windows saved to: {output_bed}")

        except Exception as e:
            if staged_output is not None:
                staged_output.unlink(missing_ok=True)
            self.logger.error(f"Error splitting BED file: {e}")
            raise
    
    def _merge_bed_files(self, bed_files, output_bed):
        """
        Merge multiple BED files and remove duplicates.
        
        Args:
            bed_files: List of BED file paths to merge
            output_bed: Path to output merged BED file
        """
        self.logger.info(f"Merging {len(bed_files)} BED files...")
        output_bed = Path(output_bed)
        output_bed.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            with tempfile.TemporaryDirectory(
                prefix=".bed-merge-", dir=output_bed.parent
            ) as temp_dir:
                temp_dir = Path(temp_dir)
                combined_bed = temp_dir / "combined.bed"
                sorted_bed = temp_dir / "sorted.bed"
                staged_output = temp_dir / "merged.bed"

                with open(combined_bed, 'w') as fout:
                    for bed_file in bed_files:
                        if Path(bed_file).exists():
                            with open(bed_file, 'r') as fin:
                                for line in fin:
                                    if line.strip() and not line.strip().startswith('#'):
                                        fout.write(line)

                cmd_sort = [
                    "bedtools", "sort", "-i", str(combined_bed),
                    "-faidx", str(self.reference_fai)
                ]
                with open(sorted_bed, 'w') as out_f, tempfile.TemporaryFile(
                    mode="w+b"
                ) as stderr_log:
                    subprocess.run(cmd_sort, check=True, stdout=out_f, stderr=stderr_log)

                cmd_merge = ["bedtools", "merge", "-i", str(sorted_bed)]
                with open(staged_output, 'w') as out_f, tempfile.TemporaryFile(
                    mode="w+b"
                ) as stderr_log:
                    subprocess.run(cmd_merge, check=True, stdout=out_f, stderr=stderr_log)
                os.replace(staged_output, output_bed)
            
            # Count merged regions
            with open(output_bed, "r") as merged_file:
                merged_count = sum(1 for _ in merged_file)
            self.logger.info(f"Merged BED file contains {merged_count} regions: {output_bed}")
            
        except subprocess.CalledProcessError as e:
            self.logger.warning(f"bedtools merge failed: {e}")
            # Fallback to simple concatenation and Python-based merge
            self._merge_bed_files_python(bed_files, output_bed)
        except FileNotFoundError:
            self.logger.warning("bedtools not found, using Python-based merge")
            self._merge_bed_files_python(bed_files, output_bed)
    
    def _combine_50bp_windows(self, bed_files, output_bed):
        """Sort candidate windows and remove exact duplicates.

        Adjacent windows deliberately remain separate: 50 bp is the BE/IR
        analysis unit.  ``smallerror_long.py`` batches those rows only for the
        mpileup query and dispatches the evidence back to each original row.
        """
        self.logger.info(
            f"Combining {len(bed_files)} 50bp window BED files "
            "without merging adjacent analysis units..."
        )
        output_bed = Path(output_bed)
        output_bed.parent.mkdir(parents=True, exist_ok=True)

        try:
            with tempfile.TemporaryDirectory(
                prefix=".50bp-combine-", dir=output_bed.parent
            ) as temp_dir:
                temp_dir = Path(temp_dir)
                combined_bed = temp_dir / "combined.bed"
                sorted_bed = temp_dir / "sorted.bed"
                staged_output = temp_dir / "candidate_windows.bed"
                with open(combined_bed, 'w', encoding="utf-8") as fout:
                    for bed_file in bed_files:
                        if Path(bed_file).exists():
                            with open(bed_file, 'r', encoding="utf-8") as fin:
                                for line in fin:
                                    if line.strip() and not line.strip().startswith('#'):
                                        fout.write(line.rstrip("\r\n") + "\n")

                cmd_sort = [
                    "bedtools", "sort", "-i", str(combined_bed),
                    "-faidx", str(self.reference_fai)
                ]
                with open(sorted_bed, 'w') as out_f, tempfile.TemporaryFile(
                    mode="w+b"
                ) as stderr_log:
                    subprocess.run(cmd_sort, check=True, stdout=out_f, stderr=stderr_log)

                unique_count = 0
                current_locus = None
                seen_at_locus = set()
                with open(sorted_bed, "r", encoding="utf-8") as sorted_handle, open(
                    staged_output, "w", encoding="utf-8"
                ) as output_handle:
                    for line_number, line in enumerate(sorted_handle, 1):
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 3:
                            raise ValueError(
                                f"Malformed sorted BED line {line_number}: {line.rstrip()}"
                            )
                        chrom, start_text, end_text = parts[:3]
                        start, end = int(start_text), int(end_text)
                        if start < 0 or end <= start or end - start > 50:
                            raise ValueError(
                                "Candidate windows must be valid intervals no longer "
                                f"than 50 bp: {chrom}\t{start}\t{end}"
                            )
                        candidate = (chrom, start, end)
                        locus = (chrom, start)
                        if locus != current_locus:
                            current_locus = locus
                            seen_at_locus.clear()
                        if candidate in seen_at_locus:
                            continue
                        output_handle.write(f"{chrom}\t{start}\t{end}\n")
                        seen_at_locus.add(candidate)
                        unique_count += 1
                os.replace(staged_output, output_bed)

            self.logger.info(
                f"Candidate BED contains {unique_count} independent <=50 bp windows: "
                f"{output_bed}"
            )

        except subprocess.CalledProcessError as e:
            self.logger.warning(f"bedtools sort failed: {e}")
            self._combine_50bp_windows_python(bed_files, output_bed)
        except FileNotFoundError:
            self.logger.warning(
                "bedtools not found, using Python-based 50bp window combination"
            )
            self._combine_50bp_windows_python(bed_files, output_bed)

    def _combine_50bp_windows_python(self, bed_files, output_bed):
        """Disk-backed fallback for sorting and exact window deduplication."""
        self.logger.info("Using disk-backed 50bp window combination...")

        output_bed = Path(output_bed)
        output_bed.parent.mkdir(parents=True, exist_ok=True)
        chrom_order = {
            chrom: index for index, chrom in enumerate(self._get_chromosomes())
        }
        with tempfile.TemporaryDirectory(
            prefix=".50bp-combine-fallback-", dir=output_bed.parent
        ) as temp_dir:
            temp_dir = Path(temp_dir)
            connection = sqlite3.connect(temp_dir / "windows.sqlite3")
            try:
                connection.execute("PRAGMA journal_mode=OFF")
                connection.execute("PRAGMA synchronous=OFF")
                connection.execute(
                    "CREATE TABLE regions ("
                    "chrom_order INTEGER, chrom TEXT, start INTEGER, end INTEGER, "
                    "PRIMARY KEY (chrom, start, end))"
                )
                rows = []
                for bed_file in bed_files:
                    if not Path(bed_file).exists():
                        continue
                    with open(bed_file, 'r', encoding="utf-8") as input_handle:
                        for line_number, line in enumerate(input_handle, 1):
                            stripped = line.strip()
                            if not stripped or stripped.startswith('#'):
                                continue
                            parts = stripped.split('\t')
                            if len(parts) < 3:
                                raise ValueError(
                                    f"Malformed BED row {line_number} in {bed_file}"
                                )
                            chrom = parts[0]
                            start = int(parts[1])
                            end = int(parts[2])
                            if start < 0 or end <= start or end - start > 50:
                                raise ValueError(
                                    "Candidate windows must be valid intervals no longer "
                                    f"than 50 bp: {chrom}\t{start}\t{end}"
                                )
                            rows.append((
                                chrom_order.get(chrom, len(chrom_order)),
                                chrom,
                                start,
                                end,
                            ))
                            if len(rows) >= 10_000:
                                connection.executemany(
                                    "INSERT OR IGNORE INTO regions VALUES (?, ?, ?, ?)",
                                    rows,
                                )
                                rows.clear()
                if rows:
                    connection.executemany(
                        "INSERT OR IGNORE INTO regions VALUES (?, ?, ?, ?)", rows
                    )
                connection.commit()
                region_count = connection.execute(
                    "SELECT COUNT(*) FROM regions"
                ).fetchone()[0]
                staged_output = temp_dir / "candidate_windows.bed"
                with open(staged_output, "w", encoding="utf-8") as output_handle:
                    for chrom, start, end in connection.execute(
                        "SELECT chrom, start, end FROM regions "
                        "ORDER BY chrom_order, chrom, start, end"
                    ):
                        output_handle.write(f"{chrom}\t{start}\t{end}\n")
            finally:
                connection.close()
            os.replace(staged_output, output_bed)

        self.logger.info(
            f"Combined {region_count} independent <=50 bp windows: {output_bed}"
        )

    # Backward-compatible private aliases for callers/tests from earlier releases.
    def _merge_50bp_windows(self, bed_files, output_bed):
        self._combine_50bp_windows(bed_files, output_bed)

    def _merge_50bp_windows_python(self, bed_files, output_bed):
        self._combine_50bp_windows_python(bed_files, output_bed)
    
    def _merge_bed_files_python(self, bed_files, output_bed):
        """
        Python-based BED file merging (fallback when bedtools is not available).
        """
        self.logger.info("Using Python-based BED merging...")
        
        regions = set()
        
        for bed_file in bed_files:
            if not Path(bed_file).exists():
                continue
            with open(bed_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split('\t')
                        if len(parts) >= 3:
                            chrom = parts[0]
                            start = int(parts[1])
                            end = int(parts[2])
                            regions.add((chrom, start, end))
        
        # Sort regions
        sorted_regions = sorted(regions, key=lambda x: (x[0], x[1]))
        
        # Merge overlapping regions
        merged = []
        for chrom, start, end in sorted_regions:
            if merged and merged[-1][0] == chrom and merged[-1][2] >= start:
                # Merge with previous region
                merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
            else:
                merged.append((chrom, start, end))
        
        _write_bed_regions_atomic(merged, output_bed)
        
        self.logger.info(f"Merged {len(merged)} regions: {output_bed}")
    
    def _sort_bed_file(self):
        """
        Sort the BED file by chromosome and position for consistent ordering.
        """
        self.logger.info("Sorting BED file for consistent merge ordering...")
        
        # Sort using bedtools sort
        cmd = [
            "bedtools", "sort",
            "-i", str(self.smallerror_bed),
            "-faidx", str(self.reference_fai)
        ]
        
        try:
            with tempfile.TemporaryDirectory(
                prefix=".bed-sort-", dir=self.smallerror_bed.parent
            ) as temp_dir:
                sorted_bed = Path(temp_dir) / "smallerror.sorted.bed"
                with open(sorted_bed, 'w') as out_f, tempfile.TemporaryFile(
                    mode="w+b"
                ) as stderr_log:
                    subprocess.run(
                        cmd,
                        check=True,
                        stdout=out_f,
                        stderr=stderr_log,
                    )
                os.replace(sorted_bed, self.smallerror_bed)
            self.logger.info(f"BED file sorted and saved to: {self.smallerror_bed}")
            
        except subprocess.CalledProcessError as e:
            self.logger.warning(f"bedtools sort failed: {e}")
            self.logger.warning("Using original BED file without sorting")
        except FileNotFoundError:
            self.logger.warning("bedtools not found, using original BED file without sorting")

    def _split_be_results_by_chromosome(self):
        """
        Split BE results by chromosome and save to separate files.
        """
        if not self.final_be_results.exists():
            self.logger.warning("Final BE results file not found, skipping chromosome split")
            return

        chromosome_output_dir = getattr(
            self, "be_chromosome_output_dir", self.be_output_dir
        )
        chromosome_output_dir.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=".be-split-", dir=chromosome_output_dir
        ) as staged_dir:
            header = None
            with open(self.final_be_results, 'r') as f:
                first_row = None
                for raw_line in f:
                    line = raw_line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if _is_chromosome_header(parts):
                        header = line
                    else:
                        first_row = (parts[0], line)
                    break

                with _ChromosomeFilePool(staged_dir, "be_results_", header) as chr_files:
                    if first_row is not None:
                        chr_files.write(*first_row)
                    for raw_line in f:
                        line = raw_line.strip()
                        if not line or line.startswith('#'):
                            continue
                        chrom = line.split('\t', 1)[0]
                        chr_files.write(chrom, line)

            chromosomes = list(chr_files.chromosomes)
            _publish_chromosome_file_set(
                staged_dir,
                chromosome_output_dir,
                "be_results_",
                chromosomes,
            )

        for chrom in chromosomes:
            self.logger.info(
                f"Saved {chrom} BE results to: "
                f"{chromosome_output_dir}/be_results_{chrom}.txt"
            )

        self.logger.info(f"BE results split into {len(chromosomes)} chromosome files")
    
    def step3_se_analysis(self):
        """Step 3: Run SE (Structural Error) analysis"""
        if not self.run_se_analysis:
            self.logger.info("Skipping SE analysis as requested")
            return
            
        self.logger.info("="*50)
        self.logger.info("STEP 3: Structural Error (SE) Analysis")
        self.logger.info("="*50)
        
        # Run integrated SE analysis using filtered ONT BAM
        cmd = [
            sys.executable,
            str(self.src_dir / "integrated_se_analysis.py"),
            str(self.ont_filtered_bam),
            str(self.reference_fa),
            "structural_errors",  # Only pass the prefix, not full path
            "--output-dir", str(self.se_output_dir),
            "--window-size", str(self.se_window_size),
            "--min-indel-length", str(self.se_min_indel_length),
            "--threads", str(self.threads),
            "--workers", str(self.se_workers),
            "--spanning-window-size", "30000",
            "--spanning-step-size", "2000"
        ]
        self._run_command(cmd, "Running SE analysis with filtered ONT data")
        
        # Update se_results_prefix to point to the actual output files
        self.se_results_prefix = self.se_output_dir / "structural_errors"
        
        # Apply bed masking to SE results if enabled
        if self.bed_mask:
            self.logger.info("Applying bed masking to SE results")
            # Check for the correct SE anomalies file name
            se_anomalies_file = self.se_output_dir / "se.anomalous_regions.tsv"
            if not se_anomalies_file.exists():
                # Try alternative name
                se_anomalies_file = self.se_output_dir / "structural_errors.anomalous_regions.tsv"
            
            if se_anomalies_file.exists():
                masked_se_anomalies = self.se_output_dir / "se.anomalous_regions_masked.tsv"
                try:
                    self.bed_mask.apply_bedtools_intersect_v(
                        se_anomalies_file, 
                        masked_se_anomalies, 
                        "SE anomalous regions",
                        coordinate_mode="one_based_closed"
                    )
                    # Update the file reference
                    self.se_anomalies_file = masked_se_anomalies
                    self.logger.info(f"Successfully applied bed masking to SE results: {masked_se_anomalies}")
                except Exception as e:
                    self.logger.warning(f"Bed masking failed for SE results: {e}")
                    self.logger.warning("Using original SE results file without masking")
                    self.se_anomalies_file = se_anomalies_file
            else:
                self.logger.warning(f"SE anomalies file not found at {se_anomalies_file}, skipping bed masking")
                # List available files in the directory for debugging
                available_files = list(self.se_output_dir.glob("*.tsv"))
                if available_files:
                    self.logger.info(f"Available TSV files in {self.se_output_dir}: {[f.name for f in available_files]}")
                self.se_anomalies_file = se_anomalies_file
        else:
            self.se_anomalies_file = self.se_output_dir / "se.anomalous_regions.tsv"

        self._publish_se_results()
        
        # Split SE results by chromosome
        self.logger.info("Splitting SE results by chromosome")
        self._split_se_results_by_chromosome()
        
        self.logger.info("SE analysis completed successfully")

    def _publish_se_results(self):
        """Publish the compact, stable SE deliverables from the work area."""
        results_dir = getattr(self, "se_results_dir", self.se_output_dir)
        results_dir.mkdir(parents=True, exist_ok=True)

        source_anomalies = getattr(
            self,
            "se_anomalies_file",
            self.se_output_dir / "structural_errors.anomalous_regions.tsv",
        )
        source_bed = self.se_output_dir / "structural_errors.result.bed"
        source_summary = self.se_output_dir / "structural_errors.anomaly_summary.tsv"

        if self.bed_mask:
            masked_bed = self.se_output_dir / "structural_errors.result_masked.bed"
            _write_merged_se_bed_from_anomalies(source_anomalies, masked_bed)
            source_bed = masked_bed

            masked_summary = (
                self.se_output_dir / "structural_errors.anomaly_summary_masked.tsv"
            )
            _write_se_anomaly_summary_from_tsv(
                source_anomalies,
                masked_summary,
                getattr(self, "chromosomes", ()),
            )
            source_summary = masked_summary

        sources = (
            (source_bed, "structural_errors.bed"),
            (source_anomalies, "anomalous_regions.tsv"),
            (source_summary, "anomaly_summary.tsv"),
        )
        missing = [str(source) for source, _ in sources if not source.exists()]
        if missing:
            raise FileNotFoundError(
                "SE analysis did not create required result file(s): "
                + ", ".join(missing)
            )

        with tempfile.TemporaryDirectory(
            prefix=".se-publish-", dir=results_dir
        ) as staged_dir_name:
            staged_dir = Path(staged_dir_name)
            replacements = []
            for source, public_name in sources:
                staged_file = staged_dir / public_name
                shutil.copy2(source, staged_file)
                replacements.append((staged_file, results_dir / public_name))
            _publish_file_set(replacements, backup_prefix=".se-results-backup-")

        self.logger.info(f"Published final SE results to: {results_dir}")

    def _split_se_results_by_chromosome(self):
        """
        Split SE results by chromosome and save to separate files.
        """
        # Determine which SE file to use based on bed masking
        if self.bed_mask:
            # Use masked version if bed masking was applied
            se_file = getattr(self, 'se_anomalies_file', self.se_output_dir / "se.anomalous_regions.tsv")
        else:
            # Use original version if no bed masking
            se_file = self.se_output_dir / "se.anomalous_regions.tsv"
        
        if not se_file.exists():
            self.logger.warning("SE results file not found, skipping chromosome split")
            return
        
        chromosome_output_dir = getattr(
            self, "se_chromosome_output_dir", self.se_output_dir
        )
        chromosome_output_dir.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=".se-split-", dir=chromosome_output_dir
        ) as staged_dir:
            header = None
            with open(se_file, 'r') as f:
                first_row = None
                for raw_line in f:
                    line = raw_line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if _is_chromosome_header(parts):
                        header = line
                    else:
                        first_row = (parts[0], line)
                    break

                with _ChromosomeFilePool(staged_dir, "se_results_", header) as chr_files:
                    if first_row is not None:
                        chr_files.write(*first_row)
                    for raw_line in f:
                        line = raw_line.strip()
                        if not line or line.startswith('#'):
                            continue
                        chrom = line.split('\t', 1)[0]
                        chr_files.write(chrom, line)

            chromosomes = list(chr_files.chromosomes)
            _publish_chromosome_file_set(
                staged_dir,
                chromosome_output_dir,
                "se_results_",
                chromosomes,
            )

        for chrom in chromosomes:
            self.logger.info(
                f"Saved {chrom} SE results to: "
                f"{chromosome_output_dir}/se_results_{chrom}.txt"
            )

        self.logger.info(f"SE results split into {len(chromosomes)} chromosome files")
    
    def step4_chromosome_analysis(self):
        """Step 4: Analyze results by chromosome and calculate chromosome-specific QV"""
        self.logger.info("="*50)
        self.logger.info("STEP 4: Chromosome-specific Analysis")
        self.logger.info("="*50)
        
        # Get chromosome lengths
        chr_lengths = self._get_chromosome_lengths()
        
        # Process results by chromosome
        chr_be_results = self._analyze_chromosome_be_results()
        chr_se_results = self._analyze_chromosome_se_results()
        
        # Calculate chromosome-specific QV values
        self._calculate_chromosome_qv(chr_lengths, chr_be_results, chr_se_results)
        
        self.logger.info("Chromosome-specific analysis completed successfully")
    
    def _analyze_chromosome_be_results(self):
        """Analyze BE results by chromosome"""
        self.logger.info("Analyzing BE results by chromosome")
        
        chr_be_results = {}
        
        if self.run_be_analysis and self.final_be_results.exists():
            with open(self.final_be_results, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split('\t')
                        if len(parts) >= 1:
                            chrom = parts[0]
                            if chrom not in chr_be_results:
                                chr_be_results[chrom] = 0
                            chr_be_results[chrom] += 1
        
        # Log results
        for chrom in self.chromosomes:
            be_count = chr_be_results.get(chrom, 0)
            self.logger.info(f"Chromosome {chrom}: {be_count} BE errors")
        
        return chr_be_results
    
    def _analyze_chromosome_se_results(self):
        """Analyze SE results by chromosome"""
        self.logger.info("Analyzing SE results by chromosome")
        
        chr_se_results = {}
        
        if self.run_se_analysis:
            # Determine which SE file to use based on bed masking
            if self.bed_mask:
                # Use masked version if bed masking was applied
                se_anomalies_file = getattr(self, 'se_anomalies_file', self.se_output_dir / "se.anomalous_regions.tsv")
            else:
                # Use original version if no bed masking
                se_anomalies_file = self.se_output_dir / "se.anomalous_regions.tsv"
            
            if se_anomalies_file.exists():
                with open(se_anomalies_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith('#'):
                            parts = line.split('\t')
                            if len(parts) >= 1:
                                chrom = parts[0]
                                if chrom not in chr_se_results:
                                    chr_se_results[chrom] = 0
                                chr_se_results[chrom] += 1
        
        # Log results
        for chrom in self.chromosomes:
            se_count = chr_se_results.get(chrom, 0)
            self.logger.info(f"Chromosome {chrom}: {se_count} SE errors")
        
        return chr_se_results
    
    def _calculate_chromosome_qv(self, chr_lengths, chr_be_results, chr_se_results):
        """Calculate QV values for each chromosome"""
        self.logger.info("Calculating chromosome-specific QV values")
        
        chr_qv_results = {}
        
        for chrom in self.chromosomes:
            total_length = chr_lengths.get(chrom, 0)
            if self.bed_mask and total_length > 0:
                total_length = max(
                    0,
                    total_length
                    - self.bed_mask.get_masked_length(chrom, total_length),
                )
            be_errors = chr_be_results.get(chrom, 0)
            se_errors = chr_se_results.get(chrom, 0)
            total_errors = be_errors + se_errors
            
            # Calculate QV
            if total_length > 0 and total_errors > 0:
                import math
                qv = -10 * math.log10(total_errors / total_length)
                qv_str = f"{qv:.2f}"
            elif total_errors == 0:
                qv_str = "∞"
            else:
                qv_str = "N/A"
            
            chr_qv_results[chrom] = {
                'total_length': total_length,
                'be_errors': be_errors,
                'se_errors': se_errors,
                'total_errors': total_errors,
                'qv': qv_str
            }
            
            self.logger.info(f"Chromosome {chrom}: QV = {qv_str} "
                           f"(BE: {be_errors}, SE: {se_errors}, "
                           f"Length: {total_length:,} bp)")
        
        # Save chromosome-specific summary
        chr_summary_file = getattr(
            self,
            "chromosome_summary_file",
            self.output_dir / "chromosome_qv_summary.txt",
        )
        chr_summary_file.parent.mkdir(parents=True, exist_ok=True)
        with open(chr_summary_file, 'w', encoding='utf-8') as f:
            f.write("Chromosome-specific QV Analysis\n")
            f.write("="*50 + "\n")
            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("Chromosome\tTotal_Length\tBE_Errors\tSE_Errors\tTotal_Errors\tQV\n")
            
            for chrom in self.chromosomes:
                result = chr_qv_results[chrom]
                f.write(f"{chrom}\t{result['total_length']:,}\t"
                       f"{result['be_errors']}\t{result['se_errors']}\t"
                       f"{result['total_errors']}\t{result['qv']}\n")
        
        self.logger.info(f"Chromosome QV summary saved: {chr_summary_file}")
        self.chr_qv_results = chr_qv_results
    
    def step4_generate_summary(self):
        """Step 4: Generate comprehensive summary report"""
        self.logger.info("="*50)
        self.logger.info("STEP 4: Generating Summary Report") 
        self.logger.info("="*50)
        
        self._generate_genome_wide_summary()
    
    def _generate_genome_wide_summary(self):
        """Generate genome-wide summary report (original method)"""
        summary_file = getattr(
            self, "summary_file", self.output_dir / "sas_pipeline_summary.txt"
        )
        summary_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Count results
        def count_lines(file_path):
            try:
                with open(file_path, 'r') as f:
                    return sum(1 for _ in f)
            except:
                return 0
        
        # Get reference genome length from fai file
        def get_reference_lengths():
            try:
                reference_lengths = {}
                with open(self.reference_fai, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            reference_lengths[parts[0]] = int(parts[1])
                return reference_lengths
            except Exception as e:
                self.logger.warning(f"Could not read reference length from {self.reference_fai}: {e}")
                return None
        
        reference_lengths = get_reference_lengths()
        reference_length = (
            sum(reference_lengths.values()) if reference_lengths is not None else None
        )
        
        # Apply bed masking if enabled
        if self.bed_mask and reference_length is not None:
            try:
                total_masked_length = 0
                for chrom in dict.fromkeys(self.chromosomes):
                    chrom_length = reference_lengths.get(chrom)
                    if chrom_length is None:
                        self.logger.warning(
                            f"Skipping masked-length calculation for unknown chromosome: {chrom}"
                        )
                        continue
                    total_masked_length += self.bed_mask.get_masked_length(
                        chrom, chrom_length
                    )
                reference_length = reference_length - total_masked_length
                self.logger.info(f"Applied bed masking: excluded {total_masked_length:,} bp")
            except Exception as e:
                self.logger.warning(f"Error calculating masked length: {e}")
                self.logger.warning("Using original reference length without masking")
        
        # BE analysis results
        be_results = {}
        total_be_errors = 0
        if self.run_be_analysis:
            be_results = {
                'smallerror_count': count_lines(self.smallerror_bed),
                'hifi_filtered_count': count_lines(self.hifi_filtered_regions),
                'ont_filtered_count': count_lines(self.ont_filtered_regions),
                'merged_count': count_lines(self.merged_results) - 1,               # Exclude header
                'final_be_count': count_lines(self.final_be_results)
            }
            total_be_errors = be_results['final_be_count']
        
        # SE analysis results
        se_results = {}
        total_se_errors = 0
        if self.run_se_analysis:
            # Determine which SE file to use based on bed masking
            if self.bed_mask:
                # Use masked version if bed masking was applied
                se_anomalies_file = getattr(self, 'se_anomalies_file', self.se_output_dir / "se.anomalous_regions.tsv")
            else:
                # Use original version if no bed masking
                se_anomalies_file = self.se_output_dir / "se.anomalous_regions.tsv"
            
            se_count = count_lines(se_anomalies_file) - 1 if se_anomalies_file.exists() else 0  # Exclude header
            se_results = {
                'anomalous_regions_count': se_count
            }
            total_se_errors = se_count
        
        # Calculate QV value
        total_errors = total_be_errors + total_se_errors
        qv_value = None
        qv_display = "N/A"
        
        if reference_length is not None and reference_length > 0:
            if total_errors == 0:
                qv_display = "∞ (Infinite - No errors detected)"
                self.logger.info("QV = ∞ (No errors detected)")
            else:
                import math
                qv_value = -10 * math.log10(total_errors / reference_length)
                qv_display = f"{qv_value:.2f}"
                self.logger.info(f"QV = {qv_value:.2f}")
        else:
            self.logger.warning(
                "Cannot calculate QV: unmasked reference length is unavailable or zero"
            )
        
        # Store QV results as instance variables for final display
        self.final_qv_value = qv_value
        self.final_qv_display = qv_display
        self.final_total_errors = total_errors
        self.final_be_errors = total_be_errors
        self.final_se_errors = total_se_errors
        self.final_reference_length = reference_length

        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("SAS Algorithm Pipeline Summary\n")
            f.write("="*50 + "\n")
            f.write(f"Execution Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Output Directory: {self.output_dir}\n")
            f.write(f"External Tool Threads: {self.threads}\n")
            f.write(
                f"Pileup Workers: {getattr(self, 'pileup_workers', min(self.threads, 16))}\n"
            )
            f.write(
                f"SE Workers: {getattr(self, 'se_workers', min(self.threads, 8))}\n"
            )
            f.write(f"Skip Filter: {self.skip_filter}\n")
            f.write(f"Run BE Analysis: {self.run_be_analysis}\n")
            f.write(f"Run SE Analysis: {self.run_se_analysis}\n")
            if self.bed_files:
                f.write(f"Bed Masking: {len(self.bed_files)} files\n")
            if self.additional_bed_file:
                f.write(f"Additional Bed File for BE Analysis: {self.additional_bed_file}\n")
            if getattr(self, 'ir_bed_file', None):
                f.write(f"Generated IR Bed File: {self.ir_bed_file}\n")
            f.write("\n")
            
            f.write("Input Files:\n")
            f.write(f"  Short-read BAM: {self.short_read_bam}\n")
            f.write(f"  HiFi BAM: {self.hifi_bam}\n")
            f.write(f"  ONT BAM: {self.ont_bam}\n")
            f.write(f"  Reference: {self.reference_fa}\n")
            f.write(f"  Reference FAI: {self.reference_fai}\n")
            if reference_length is not None and reference_length > 0:
                f.write(f"  Reference Length: {reference_length:,} bp\n")
            if self.bed_files:
                f.write("  Bed Mask Files:\n")
                for bed_file in self.bed_files:
                    f.write(f"    {bed_file}\n")
            if self.additional_bed_file:
                f.write(f"  Additional Bed File (split into 50bp windows): {self.additional_bed_file}\n")
            f.write("\n")
            
            f.write("Output Files:\n")
            f.write(f"  Short-read filtered BAM: {self.short_filtered_bam}\n")
            f.write(f"  HiFi filtered BAM: {self.hifi_filtered_bam}\n")
            f.write(f"  ONT filtered BAM: {self.ont_filtered_bam}\n")
            
            if self.run_be_analysis:
                f.write(
                    "  BE Analysis Results: "
                    f"{getattr(self, 'be_results_dir', self.be_output_dir)}\n"
                )
                f.write(f"  Final BE results: {self.final_be_results}\n")
            
            if self.run_se_analysis:
                public_se_dir = getattr(self, "se_results_dir", self.se_output_dir)
                f.write(f"  SE Analysis Results: {public_se_dir}\n")
                se_anomalies_file = public_se_dir / "anomalous_regions.tsv"
                f.write(f"  SE anomalous regions: {se_anomalies_file}\n")
            
            chromosome_summary_file = getattr(
                self,
                "chromosome_summary_file",
                self.output_dir / "chromosome_qv_summary.txt",
            )
            f.write(f"  Chromosome QV summary: {chromosome_summary_file}\n")
            
            f.write("\nResults Statistics:\n")
            
            if self.run_be_analysis:
                f.write("  BE (Base-level Errors) Analysis:\n")
                f.write(f"    Initial small error regions: {be_results['smallerror_count']}\n")
                f.write(f"    After HiFi filtering: {be_results['hifi_filtered_count']}\n")
                f.write(f"    After ONT filtering: {be_results['ont_filtered_count']}\n")
                f.write(f"    After merging: {be_results['merged_count']}\n")
                f.write(f"    Final BE count: {be_results['final_be_count']}\n")
            
            if self.run_se_analysis:
                f.write("  SE (Structural Errors) Analysis:\n")
                f.write(f"    Anomalous regions detected: {se_results['anomalous_regions_count']}\n")
            
            # Chromosome-specific breakdown
            if hasattr(self, 'chr_qv_results'):
                f.write("\nChromosome-specific Results:\n")
                f.write("="*50 + "\n")
                
                for chrom in self.chromosomes:
                    result = self.chr_qv_results[chrom]
                    f.write(f"\nChromosome: {chrom}\n")
                    f.write(f"  Length: {result['total_length']:,} bp\n")
                    f.write(f"  BE Errors: {result['be_errors']}\n")
                    f.write(f"  SE Errors: {result['se_errors']}\n")
                    f.write(f"  Total Errors: {result['total_errors']}\n")
                    f.write(f"  QV: {result['qv']}\n")
            
            # Quality Value calculation
            f.write("\nOverall Quality Assessment:\n")
            f.write(f"  Total BE errors: {total_be_errors}\n")
            f.write(f"  Total SE errors: {total_se_errors}\n")
            f.write(f"  Total errors: {total_errors}\n")
            if reference_length is not None and reference_length > 0:
                f.write(f"  Assembly length: {reference_length:,} bp\n")
                f.write(f"  Error rate: {total_errors/reference_length:.2e}\n")
            f.write(f"  Quality Value (QV): {qv_display}\n")
            
            f.write(f"\nQV Formula: QV = -10 × log₁₀((SE + BE) / Assembly Length)\n")
            if total_errors == 0:
                f.write("Note: QV is infinite when no errors are detected.\n")
            
            f.write("\nPipeline Steps Completed:\n")
            if self.skip_filter:
                f.write("  ✓ Step 1: BAM files preparation (skipped filtering - using pre-filtered files)\n")
            else:
                f.write("  ✓ Step 1: BAM files filtered\n")
            
            if self.run_be_analysis:
                f.write("  ✓ Step 2: BE analysis completed\n")
            else:
                f.write("  - Step 2: BE analysis skipped\n")
            
            if self.run_se_analysis:
                f.write("  ✓ Step 3: SE analysis completed\n")
            else:
                f.write("  - Step 3: SE analysis skipped\n")
            
            f.write("  ✓ Step 4: Chromosome-specific analysis completed\n")
            f.write("  ✓ Step 5: Summary report generated\n")
        
        self.logger.info(f"Summary report saved: {summary_file}")
        
        # Log final QV to console for immediate feedback
        if total_errors == 0:
            self.logger.info("🎉 PERFECT ASSEMBLY: QV = ∞ (No errors detected)")
        elif qv_value is not None:
            self.logger.info(f"📊 FINAL QUALITY VALUE: QV = {qv_value:.2f}")
        
        # Also save configuration for reproducibility
        config_file = getattr(
            self, "config_file", self.output_dir / "pipeline_config.json"
        )
        config_file.parent.mkdir(parents=True, exist_ok=True)
        resolved_config = dict(self.config)
        resolved_config.update(
            threads=self.threads,
            pileup_workers=getattr(
                self, 'pileup_workers', min(self.threads, 16)
            ),
            se_workers=getattr(self, 'se_workers', min(self.threads, 8)),
        )
        with open(config_file, 'w', encoding='utf-8') as f:
            json.dump(resolved_config, f, indent=2, default=str)
        
        self.logger.info(f"Configuration saved: {config_file}")
    


    def run_complete_pipeline(self):
        """Run the complete SAS pipeline"""
        start_time = datetime.now()
        
        self.logger.info("="*60)
        self.logger.info("SAS ALGORITHM INTEGRATED PIPELINE STARTED")
        self.logger.info("="*60)
        self.logger.info(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        
        try:
            self._validate_reference_compatibility()
            self.step1_filter_bam_files()
            self.step2_be_analysis()
            self.step3_se_analysis()
            self.step4_chromosome_analysis()
            self.step4_generate_summary()
            
            end_time = datetime.now()
            duration = end_time - start_time
            
            self.logger.info("="*60)
            self.logger.info("SAS PIPELINE COMPLETED SUCCESSFULLY!")
            self.logger.info("="*60)
            self.logger.info(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
            self.logger.info(f"Total duration: {duration}")
            final_output_dir = getattr(self, "final_output_dir", self.output_dir)
            self.logger.info(f"Final results saved in: {final_output_dir}")
            
            # Display final quality assessment
            if hasattr(self, 'final_qv_display'):
                self.logger.info("="*60)
                self.logger.info("FINAL QUALITY ASSESSMENT")
                self.logger.info("="*60)
                self.logger.info(f"📊 Total BE errors: {getattr(self, 'final_be_errors', 0)}")
                self.logger.info(f"📊 Total SE errors: {getattr(self, 'final_se_errors', 0)}")
                self.logger.info(f"📊 Total errors: {getattr(self, 'final_total_errors', 0)}")
                if hasattr(self, 'final_reference_length') and self.final_reference_length:
                    self.logger.info(f"📏 Assembly length: {self.final_reference_length:,} bp")
                    
                if self.final_total_errors == 0:
                    self.logger.info("🎉 ASSEMBLY QUALITY: QV = ∞ (No errors detected!)")
                elif self.final_qv_value is not None:
                    self.logger.info(f"📊 ASSEMBLY QUALITY: QV = {self.final_qv_value:.2f}")
                else:
                    self.logger.info(f"📊 ASSEMBLY QUALITY: QV = {self.final_qv_display}")
                self.logger.info("="*60)
            
            # Display chromosome-specific results
            if hasattr(self, 'chr_qv_results'):
                self.logger.info("="*60)
                self.logger.info("CHROMOSOME-SPECIFIC QUALITY ASSESSMENT")
                self.logger.info("="*60)
                for chrom in self.chromosomes:
                    result = self.chr_qv_results[chrom]
                    self.logger.info(f"📊 {chrom}: QV = {result['qv']} "
                                   f"(BE: {result['be_errors']}, SE: {result['se_errors']}, "
                                   f"Length: {result['total_length']:,} bp)")
                self.logger.info("="*60)
            
            if self.run_be_analysis:
                self.logger.info(f"BE results: {self.final_be_results}")
            if self.run_se_analysis:
                self.logger.info(
                    f"SE results: {getattr(self, 'se_results_dir', self.se_output_dir)}"
                )
            self.logger.info(
                "Chromosome QV summary: "
                f"{getattr(self, 'chromosome_summary_file', self.output_dir / 'chromosome_qv_summary.txt')}"
            )
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            raise


def create_config_from_args(args):
    """Create configuration dictionary from command line arguments"""
    config = {
        'short_read_bam': args.short_read_bam,
        'hifi_bam': args.hifi_bam,
        'ont_bam': args.ont_bam,
        'reference_fa': args.reference_fa,
        'output_dir': args.output_dir,
        'threads': args.threads,
        'pileup_workers': getattr(args, 'pileup_workers', None),
        'se_workers': getattr(args, 'se_workers', None),
        'scripts_dir': args.scripts_dir,
        'src_dir': args.src_dir,
        'skip_filter': args.skip_filter,
        'run_be_analysis': args.run_be_analysis,
        'run_se_analysis': args.run_se_analysis,
        'se_window_size': args.se_window_size,
        'se_min_indel_length': args.se_min_indel_length,
        'bed_files': args.bed_files,
        'additional_bed_file': getattr(args, 'additional_bed_file', None),
        'generate_ir': getattr(args, 'generate_ir', False)
    }
    
    # Add reference_fai if provided
    if hasattr(args, 'reference_fai') and args.reference_fai:
        config['reference_fai'] = args.reference_fai
    
    return config


def validate_input_files(config):
    """Validate that all required input files exist"""
    required_files = [
        'short_read_bam',
        'hifi_bam', 
        'ont_bam',
        'reference_fa'
    ]
    
    for file_key in required_files:
        file_path = config[file_key]
        # Convert to Path object if it's a string
        if isinstance(file_path, str):
            file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"Required input file not found: {file_path}")
    
    # Check reference_fai if provided
    if 'reference_fai' in config and config['reference_fai']:
        ref_fai_path = config['reference_fai']
        if isinstance(ref_fai_path, str):
            ref_fai_path = Path(ref_fai_path)
        if not ref_fai_path.exists():
            raise FileNotFoundError(f"Reference fai file not found: {ref_fai_path}")
    
    # Check bed files if provided
    if 'bed_files' in config and config['bed_files']:
        for bed_file in config['bed_files']:
            bed_path = Path(bed_file)
            if not bed_path.exists():
                raise FileNotFoundError(f"Bed file not found: {bed_path}")
    
    # Check additional bed file if provided
    if 'additional_bed_file' in config and config['additional_bed_file']:
        additional_bed_path = Path(config['additional_bed_file'])
        if not additional_bed_path.exists():
            raise FileNotFoundError(f"Additional bed file not found: {additional_bed_path}")
    
    # Check scripts directory
    scripts_dir = Path(config['scripts_dir'])
    src_dir = Path(config['src_dir'])
    
    # Required scripts for BE analysis
    required_scripts = []
    if config.get('run_be_analysis', True):
        required_scripts.extend([
            'smallerror_short.sh',
            'smallerror_long.py',
            'merge_pileup.py',
            'filter_be.py'
        ])
        
    # BAM filtering runs before either analysis type and is required even for
    # an SE-only run unless the user explicitly supplies pre-filtered BAMs.
    if not config.get('skip_filter', False):
        required_scripts.append('filter_bam.sh')

    # Required scripts for SE analysis
    if config.get('run_se_analysis', True):
        required_src_files = ['integrated_se_analysis.py']
        for src_file in required_src_files:
            src_path = src_dir / src_file
            if not src_path.exists():
                raise FileNotFoundError(f"Required source file not found: {src_path}")
    
    for script in required_scripts:
        script_path = scripts_dir / script
        if not script_path.exists():
            raise FileNotFoundError(f"Required script not found: {script_path}")


def validate_runtime_dependencies(config):
    """Fail before long-running work when required command-line tools are absent."""
    run_be = config.get('run_be_analysis', True)
    run_se = config.get('run_se_analysis', True)
    skip_filter = config.get('skip_filter', False)

    required_tools = {"samtools", "mosdepth"}
    if run_be or run_se:
        required_tools.add("bedtools")
    if run_be or not skip_filter:
        required_tools.add("bash")
    if run_be and config.get('generate_ir', False):
        required_tools.add("genmap")

    missing_tools = sorted(
        tool for tool in required_tools if shutil.which(tool) is None
    )
    if missing_tools:
        raise FileNotFoundError(
            "Required external tools not found in PATH: "
            + ", ".join(missing_tools)
        )


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Run the SAS base-error (BE) and structural-error (SE) pipeline.",
        formatter_class=CLIHelpFormatter,
        epilog="""
Examples:
  sas-pipeline -i illumina.bam -H hifi.bam -N ont.bam \\
    -r reference.fa -o results -t 32

  sas-pipeline -i illumina.bam -H hifi.bam -N ont.bam \\
    -r reference.fa -o results --no-se-analysis -a candidates.bed

  sas-pipeline -i illumina.bam -H hifi.bam -N ont.filtered.bam \\
    -r reference.fa -o results --skip-filter --no-be-analysis
        """,
    )

    required = parser.add_argument_group("required inputs")
    required.add_argument(
        "-i", "--short-read-bam", required=True, metavar="BAM",
        help="Short-read BAM file (Illumina, Element, etc.)",
    )
    required.add_argument(
        "-H", "--hifi-bam", required=True, metavar="BAM",
        help="PacBio HiFi BAM file",
    )
    required.add_argument(
        "-N", "--ont-bam", required=True, metavar="BAM",
        help="Oxford Nanopore BAM file",
    )
    required.add_argument(
        "-r", "--reference-fa", required=True, metavar="FASTA",
        help="Reference genome FASTA file",
    )
    required.add_argument(
        "-o", "--output-dir", required=True, metavar="DIR",
        help="Output root for final results, work files, logs, and metadata",
    )

    general = parser.add_argument_group("general options")
    general.add_argument(
        "-f", "--reference-fai", metavar="FAI",
        help="Reference FASTA index; generated automatically when omitted",
    )
    general.add_argument(
        "-t", "--threads", type=positive_int, default=8, metavar="N",
        help="Threads available to samtools/mosdepth (default: 8)",
    )
    general.add_argument(
        "--pileup-workers", type=positive_int, metavar="N",
        help="Concurrent HiFi/ONT mpileup batches (default: min(threads, 16))",
    )
    general.add_argument(
        "--se-workers", type=positive_int, metavar="N",
        help="Concurrent SE chromosome/spanning workers (default: min(threads, 8))",
    )

    pipeline_base = get_pipeline_base_dir()
    advanced = parser.add_argument_group("advanced path options")
    advanced.add_argument(
        "--scripts-dir", default=str(pipeline_base / "scripts"), metavar="DIR",
        help=f"Directory containing SAS helper scripts (default: {pipeline_base / 'scripts'})",
    )
    advanced.add_argument(
        "--src-dir", default=str(pipeline_base / "src"), metavar="DIR",
        help=f"Directory containing SAS Python source files (default: {pipeline_base / 'src'})",
    )

    control = parser.add_argument_group("pipeline control")
    control.add_argument(
        "--skip-filter", action="store_true",
        help="Use input BAM files as-is without the filtering step",
    )
    control.add_argument(
        "--no-be-analysis", action="store_true",
        help="Skip base-error analysis",
    )
    control.add_argument(
        "--no-se-analysis", action="store_true",
        help="Skip structural-error analysis",
    )

    se_options = parser.add_argument_group("structural-error options")
    se_options.add_argument(
        "--se-window-size", type=positive_int, default=2000, metavar="BP",
        help="SE analysis window size (default: 2000 bp)",
    )
    se_options.add_argument(
        "--se-min-indel-length", type=positive_int, default=50, metavar="BP",
        help="Minimum INDEL length considered by SE analysis (default: 50 bp)",
    )

    be_options = parser.add_argument_group("base-error and masking options")
    be_options.add_argument(
        "-b", "--bed-files", nargs="+", metavar="BED",
        help="One or more BED files whose intervals are masked",
    )
    be_options.add_argument(
        "-a", "--additional-bed-file", metavar="BED",
        help=(
            "Additional BED file to split into independent 50 bp windows and "
            "combine with low-coverage BE candidates"
        ),
    )
    be_options.add_argument(
        "-g", "--generate-ir", action="store_true",
        help=(
            "Run optional identical-region refinement with genmap using "
            "automatically detected parental FASTA files; runtime depends on "
            "assembly size and available hardware"
        ),
    )

    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )

    args = parser.parse_args()
    
    # Convert negative flags to positive
    args.run_be_analysis = not args.no_be_analysis
    args.run_se_analysis = not args.no_se_analysis
    
    # Validate that at least one analysis is enabled
    if not args.run_be_analysis and not args.run_se_analysis:
        print("Error: At least one analysis type (BE or SE) must be enabled", file=sys.stderr)
        sys.exit(1)
    
    # Create configuration
    config = create_config_from_args(args)
    
    # Validate input files
    try:
        validate_input_files(config)
        validate_runtime_dependencies(config)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Run pipeline
    try:
        pipeline = SASPipeline(config)
        pipeline.run_complete_pipeline()
        
    except KeyboardInterrupt:
        print("\nPipeline interrupted by user", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Pipeline failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main() 
