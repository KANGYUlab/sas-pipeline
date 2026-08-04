#!/usr/bin/env python3
"""
Integrated Structural Error Analysis Tool with Spanning Reads Detection

This script integrates clip, indel, depth and spanning reads analysis to identify anomalous regions
in BAM files using sliding windows with customizable thresholds.

Features:
- Analyzes clips, indels, depth and spanning reads in a unified framework
- Identifies anomalous regions based on configurable criteria
- Special handling for chromosome terminal regions (10kb)
- Multi-threaded processing for performance
- Spanning reads analysis with custom window generation and merging

Anomaly Detection Criteria:
1. Non-terminal regions (outside first/last 10kb):
   - Depth anomalies: < 20% or > 200% of mean depth
   - INDEL/Clip anomalies: ReadsWithInsertion/ReadsWithDeletion/ClipCount ratio > 0.5
   - Non-spanning: Regions without spanning reads
   - Note: Uses reads containing INDELs, not INDEL event counts
2. Terminal regions (first/last 10kb):
   - NO depth anomaly detection (depth variations ignored)
   - Anomaly if ReadsWithoutINDELorClip: count < 1 OR ratio <= 10%
   - Normal if ReadsWithoutINDELorClip: count >= 1 AND ratio > 10%

License: MIT

Requirements:
    - pysam
    - numpy
    - pandas
    - mosdepth (external tool)
    - bedtools (for merging)
    - Python 3.8+

Usage:
    python integrated_se_analysis.py input.bam reference.fa output_prefix [options]
"""

from __future__ import annotations

import sys
import os
import argparse
import subprocess
import gzip
import math
import shutil
import traceback
import tempfile
from pathlib import Path
from collections import defaultdict
from typing import Dict, Tuple, Set, List, DefaultDict, Optional
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor
import multiprocessing

try:
    from ._version import __version__
except ImportError:
    # Keep ``python src/integrated_se_analysis.py`` working as documented.
    try:
        from src._version import __version__
    except ImportError:
        from _version import __version__

try:
    import numpy as np
except ModuleNotFoundError:
    np = None

try:
    import pandas as pd
except ModuleNotFoundError:
    pd = None

try:
    import pysam
except ModuleNotFoundError:
    pysam = None


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


# ===== SPANNING READS ANALYSIS FUNCTIONS =====

def call_spanning_analysis(bam_file_path, fai_file, output_prefix, spanning_window_size=30000, spanning_step_size=2000, threads=multiprocessing.cpu_count()):
    """Call spanning.py script to perform spanning reads analysis (outputs only merged+terminal-corrected bed)"""
    script_dir = Path(__file__).parent.parent / "scripts"
    spanning_script = script_dir / "spanning.py"
    if not spanning_script.exists():
        raise FileNotFoundError(f"spanning.py script not found at {spanning_script}")
    cmd = [
        sys.executable,
        str(spanning_script),
        bam_file_path,
        fai_file,
        "--window-size", str(spanning_window_size),
        "--step-size", str(spanning_step_size),
        "--output", output_prefix,
        "--processes", str(threads)
    ]
    try:
        print(f"Running spanning reads analysis...", file=sys.stderr)
        print(f"Command: {' '.join(cmd)}", file=sys.stderr)
        # Keep stdout/stderr available for diagnostics from the child process.
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        if result.stderr:
            print(result.stderr, file=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error running spanning analysis: {e}", file=sys.stderr)
        if e.stderr:
            print(f"Error output: {e.stderr}", file=sys.stderr)
        raise
    except Exception as e:
        print(f"Unexpected error in spanning analysis: {e}", file=sys.stderr)
        raise


SPANNING_OUTPUT_SUFFIXES = (
    ".non_spanning.bed",
    ".non_spanning_regions.merged.bed",
)


def _same_existing_file(first_path, second_path):
    """Return False when a safe same-file comparison cannot be completed."""
    try:
        return os.path.samefile(first_path, second_path)
    except OSError:
        return False


def _publish_spanning_output_pair(replacements, backup_prefix):
    """Publish both spanning BEDs, rolling back after any partial failure."""
    replacements = [
        (Path(source), Path(destination))
        for source, destination in replacements
    ]
    output_dir = replacements[0][1].parent

    with tempfile.TemporaryDirectory(
        prefix=backup_prefix, dir=output_dir
    ) as backup_dir_name:
        backup_dir = Path(backup_dir_name)
        backups = {}
        for _, destination in replacements:
            if not destination.exists():
                backups[destination] = None
                continue
            backup = backup_dir / destination.name
            try:
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
                            # A pre-rename publication failure leaves the live
                            # path and its hard-link backup on the same inode;
                            # that destination is already in its old state.
                            if _same_existing_file(backup, destination):
                                continue
                            if backup.exists():
                                shutil.copy2(backup, destination)
                            elif not destination.exists():
                                raise
                except Exception as rollback_error:
                    rollback_errors.append(f"{destination}: {rollback_error}")

            if rollback_errors:
                raise RuntimeError(
                    "Spanning output publication failed and rollback was "
                    "incomplete: " + "; ".join(rollback_errors)
                ) from publish_error
            raise


def refresh_spanning_outputs(
    bam_file_path,
    fai_file,
    output_prefix,
    spanning_window_size=30000,
    spanning_step_size=2000,
    threads=multiprocessing.cpu_count(),
):
    """Recompute spanning outputs and publish them only after a successful run.

    ``spanning.py`` writes several files directly beside its output prefix.  A
    unique temporary directory prevents an interrupted run from truncating a
    previous result.  Both public outputs are required (an empty BED is a valid
    result), then published as a rollback-safe pair.  The merged BED, which is
    consumed by this module, remains last in the normal replacement order.
    """
    output_prefix = os.path.abspath(os.fspath(output_prefix))
    output_dir = os.path.dirname(output_prefix) or os.curdir
    output_name = os.path.basename(output_prefix)

    with tempfile.TemporaryDirectory(
        prefix=f".{output_name}.spanning-", dir=output_dir
    ) as temp_dir:
        temporary_prefix = os.path.join(temp_dir, output_name)
        call_spanning_analysis(
            bam_file_path,
            fai_file,
            temporary_prefix,
            spanning_window_size,
            spanning_step_size,
            threads,
        )

        temporary_outputs = {
            suffix: temporary_prefix + suffix
            for suffix in SPANNING_OUTPUT_SUFFIXES
        }
        missing = [
            path for path in temporary_outputs.values()
            if not os.path.isfile(path)
        ]
        if missing:
            raise FileNotFoundError(
                "Spanning analysis completed without required output(s): "
                + ", ".join(missing)
            )

        # Publish the merged BED last while keeping the old pair recoverable
        # until both replacements have completed.
        _publish_spanning_output_pair(
            [
                (temporary_outputs[suffix], output_prefix + suffix)
                for suffix in SPANNING_OUTPUT_SUFFIXES
            ],
            backup_prefix=f".{output_name}.spanning-backup-",
        )

    return output_prefix + ".non_spanning_regions.merged.bed"


def assign_merged_regions_to_2k_windows(merged_regions, window_size=2000):
    """Given merged+terminal-corrected regions (0-based), divide them into 2kb windows (chrom, start)."""
    windows = set()
    for chrom, region_start, region_end in merged_regions:
        start = (region_start // window_size) * window_size
        while start < region_end:
            end = min(start + window_size, region_end)
            if end > start:
                windows.add((chrom, start))  # 0-based start
            start += window_size
    return windows

def load_merged_non_spanning_regions_from_bed(bed_file):
    """Load merged+terminal-corrected non-spanning regions from BED file (chrom, start, end, 0-based)."""
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                regions.append((chrom, start, end))
    return regions

# ===== ORIGINAL SE ANALYSIS CLASSES AND FUNCTIONS =====


def new_window_counts():
    """Create counters for one materialized window.

    Most reference windows contain no read-derived event.  Keeping this factory
    separate lets the analysis retain only touched/non-spanning windows and
    synthesize zero windows while writing results.
    """
    return {
        'ClipCount': 0,
        'Insertion': 0,
        'Deletion': 0,
        'ReadsWithInsertion': 0,
        'ReadsWithDeletion': 0,
        'ReadsWithClip': 0,
        'TotalReads': 0,
        'ReadsWithoutINDELorClip': 0
    }


def iter_complete_windows(windows, chrom_lengths, window_size,
                          sort_chromosomes=False):
    """Yield every reference window, synthesizing zeros without storing them.

    Classification follows BAM header order, represented by the insertion
    order of ``chrom_lengths``, just as the former dense counter dictionary
    did.  Callers that historically sorted their output can request the old
    lexicographic ordering explicitly.
    """
    empty_counts = new_window_counts()
    chromosomes = sorted(chrom_lengths) if sort_chromosomes else chrom_lengths
    for chrom in chromosomes:
        chrom_length = chrom_lengths[chrom]
        for start in range(0, chrom_length, window_size):
            yield (chrom, start), windows.get((chrom, start), empty_counts)


class IntegratedSECounters:
    """
    Thread-safe counters for integrated SE analysis.
    
    This class ensures that multiple threads can safely update all SE counts
    when processing different chromosomes simultaneously.
    """
    
    def __init__(self, chrom_lengths=None, window_size=2000):
        self.lock = threading.Lock()
        self.windows: Dict[Tuple[str, int], Dict[str, int]] = {}
        self.chrom_lengths = chrom_lengths or {}
        self.window_size = window_size
    
    def init_chromosome_windows(self, chrom: str, chrom_length: int, window_size: int):
        """Initialize windows for a specific chromosome."""
        with self.lock:
            for start in range(0, chrom_length, window_size):
                if start >= chrom_length:
                    break
                self.windows[(chrom, start)] = new_window_counts()
    
    def update_windows(self, local_windows: Dict[Tuple[str, int], Dict[str, int]]):
        """
        Thread-safe update of window counters.
        """
        with self.lock:
            for window_key, counts in local_windows.items():
                if window_key in self.windows:
                    for count_type, count_value in counts.items():
                        if count_type == 'HasSpanningReads':
                            # For boolean values, use AND logic
                            current = self.windows[window_key].get(count_type, True)
                            self.windows[window_key][count_type] = current and count_value
                        else:
                            current = self.windows[window_key].get(count_type, 0)
                            self.windows[window_key][count_type] = current + count_value

    def store_chromosome_windows(self, chrom, local_windows):
        """Store a completed chromosome without copying its nested counters."""
        with self.lock:
            if any(key[0] != chrom for key in local_windows):
                raise ValueError("Chromosome worker returned windows for another chromosome")
            self.windows.update(local_windows)
    
    def mark_non_spanning_windows(self, non_spanning_windows: Set[Tuple[str, int]]):
        """Mark windows without spanning reads"""
        with self.lock:
            for window_key in non_spanning_windows:
                chrom, start = window_key
                chrom_length = self.chrom_lengths.get(chrom)
                counts = self.windows.get(window_key)
                if counts is not None:
                    counts['HasSpanningReads'] = False
                elif chrom_length is not None and 0 <= start < chrom_length:
                    counts = new_window_counts()
                    counts['HasSpanningReads'] = False
                    self.windows[window_key] = counts


def get_chromosome_lengths_from_bam(bam: pysam.AlignmentFile) -> Dict[str, int]:
    """
    Extract chromosome lengths from BAM file header.
    
    Args:
        bam: Opened BAM file object
        
    Returns:
        Dictionary mapping chromosome names to their lengths
    """
    chrom_lengths = {}
    
    if bam.has_index():
        for i in range(bam.nreferences):
            ref_name = bam.get_reference_name(i)
            ref_len = bam.get_reference_length(ref_name)
            chrom_lengths[ref_name] = ref_len
    else:
        print("Warning: BAM file has no index or missing header information.", file=sys.stderr)
        # Fallback: try to get from reference list
        for chrom_name in bam.references:
            if chrom_name not in chrom_lengths:
                try:
                    chrom_lengths[chrom_name] = bam.get_reference_length(chrom_name)
                except ValueError:
                    print(f"Warning: Cannot get length for chromosome '{chrom_name}'.", file=sys.stderr)
                    chrom_lengths[chrom_name] = 0
    
    return chrom_lengths


def run_mosdepth(
    bam_file: str, 
    output_prefix: str, 
    threads: int = 4,
    window_size: int = 2000
) -> str:
    """
    Run mosdepth to calculate sequencing depth.
    
    Args:
        bam_file: Path to input BAM file
        output_prefix: Prefix for output files
        threads: Number of threads to use
        window_size: Window size for depth calculation
    
    Returns:
        Path to the generated regions.bed.gz file
    
    Raises:
        subprocess.CalledProcessError: If mosdepth fails
        FileNotFoundError: If mosdepth is not found in PATH
    """
    # Limit threads to reasonable number for mosdepth
    max_threads = min(threads, 64)  # mosdepth doesn't benefit much from >64 threads
    if threads != max_threads:
        print(f"Warning: Limiting mosdepth threads from {threads} to {max_threads}")
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Construct mosdepth command
    cmd = [
        'mosdepth',
        '-t', str(max_threads),
        '-n',  # Don't output per-base depth
        '-x',  # Don't output depth for each position
        '-b', str(window_size),  # Window size
        output_prefix,
        bam_file
    ]
    
    print(f"Running mosdepth command: {' '.join(cmd)}")
    
    try:
        # Run mosdepth
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        print("mosdepth completed successfully")
        
        # Return path to regions file
        regions_file = f"{output_prefix}.regions.bed.gz"
        if not os.path.exists(regions_file):
            raise FileNotFoundError(f"Expected output file not found: {regions_file}")
        
        return regions_file
        
    except FileNotFoundError:
        raise FileNotFoundError("mosdepth not found in PATH. Please install mosdepth first.")
    except subprocess.CalledProcessError as e:
        # Provide more detailed error information
        error_msg = f"mosdepth failed with return code {e.returncode}"
        if hasattr(e, 'stderr') and e.stderr:
            error_msg += f"\nSTDERR: {e.stderr}"
        if hasattr(e, 'stdout') and e.stdout:
            error_msg += f"\nSTDOUT: {e.stdout}"
        
        print(f"Error: {error_msg}", file=sys.stderr)
        raise subprocess.CalledProcessError(e.returncode, cmd, error_msg)


def load_depth_data(regions_gz_file: str, window_size: int = 2000) -> Dict[Tuple[str, int], float]:
    """
    Load depth data from mosdepth regions file.
    
    Args:
        regions_gz_file: Path to gzipped regions file
        window_size: Window size used in analysis
    
    Returns:
        Dictionary mapping (chrom, start) to depth value
    """
    depth_data = {}
    
    try:
        with gzip.open(regions_gz_file, 'rt') as f:
            for line_number, line in enumerate(f, 1):
                parts = line.strip().split('\t')
                if not line.strip():
                    continue
                if len(parts) < 4:
                    raise ValueError(
                        f"malformed mosdepth row {line_number}: {line.strip()}"
                    )
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                depth = float(parts[3])
                if (
                    not chrom
                    or start < 0
                    or end <= start
                    or not math.isfinite(depth)
                    or depth < 0
                ):
                    raise ValueError(
                        f"invalid mosdepth row {line_number}: {line.strip()}"
                    )

                # Convert to window start (0-based to match other analyses)
                window_start = (start // window_size) * window_size
                key = (chrom, window_start)
                if key in depth_data:
                    raise ValueError(
                        f"duplicate mosdepth window on row {line_number}: "
                        f"{chrom}:{window_start}"
                    )
                depth_data[key] = depth
    except Exception as error:
        raise RuntimeError(
            f"Failed to load complete mosdepth data from {regions_gz_file}: {error}"
        ) from error

    if not depth_data:
        raise RuntimeError(
            f"Mosdepth regions file contains no depth windows: {regions_gz_file}"
        )

    return depth_data


def validate_depth_windows(depth_data, chrom_lengths, window_size):
    """Verify that mosdepth produced exactly one row per reference window."""
    expected_count = sum(
        (chrom_length + window_size - 1) // window_size
        for chrom_length in chrom_lengths.values()
    )
    for chrom, start in depth_data:
        chrom_length = chrom_lengths.get(chrom)
        if (
            chrom_length is None
            or start < 0
            or start >= chrom_length
            or start % window_size != 0
        ):
            raise RuntimeError(
                f"Mosdepth produced an unexpected window: {chrom}:{start}"
            )
    if len(depth_data) != expected_count:
        raise RuntimeError(
            "Incomplete mosdepth window set: expected "
            f"{expected_count}, observed {len(depth_data)}"
        )


def calculate_affected_windows(start_pos: int, deletion_length: int, window_size: int) -> List[int]:
    """
    Calculate all windows intersecting a half-open deletion interval.
    """
    if deletion_length <= 0:
        return []

    start_window = (start_pos // window_size) * window_size
    # A deletion consumes [start_pos, start_pos + deletion_length).  Convert
    # its last affected base to a window so a deletion ending exactly on a
    # boundary does not also mark the following window.
    end_pos = start_pos + deletion_length - 1
    end_window = (end_pos // window_size) * window_size
    
    affected_windows = []
    current_window = start_window
    while current_window <= end_window:
        affected_windows.append(current_window)
        current_window += window_size
    return affected_windows


def process_chromosome_integrated(
    bam_file: str,
    chrom: str,
    chrom_length: int,
    min_indel_length: int = 50,
    window_size: int = 2000,
    terminal_region_size: int = 10000
) -> Dict[Tuple[str, int], Dict[str, int]]:
    """
    Process all reads for a specific chromosome and collect integrated SE data.
    
    Args:
        bam_file: Path to BAM file
        chrom: Chromosome name to process
        chrom_length: Length of the chromosome
        min_indel_length: Minimum INDEL length to consider (default: 50bp)
        window_size: Size of sliding window
        terminal_region_size: Size of terminal regions requiring special analysis
    
    Returns:
        Dictionary with local window coordinates and integrated SE counts for this chromosome
    
    Note:
        ReadsWithoutINDELorClip counts reads that have:
        - NO significant INDELs (length >= min_indel_length, default 50bp) in the current window
        - NO clip events in the current window
        This is ONLY calculated in terminal regions and used for terminal region anomaly detection.
        A read is counted if it has no clip/indel events within the current window, regardless of read length.
    """
    # Keep only windows touched by reads/events.  Zero windows are synthesized
    # by iter_complete_windows when results are classified or written.
    local_windows = {}

    def get_window(window_key):
        start = window_key[1]
        if start < 0 or start >= chrom_length:
            return None
        counts = local_windows.get(window_key)
        if counts is None:
            counts = new_window_counts()
            local_windows[window_key] = counts
        return counts

    with pysam.AlignmentFile(bam_file, "rb") as bam:
            # Process reads for this specific chromosome
            for read in bam.fetch(chrom):
                if read.is_unmapped or read.cigartuples is None:
                    continue
                if read.reference_start is None:
                    continue
                
                # Track which windows have insertions, deletions, or clips in this read
                read_insertion_windows: Set[Tuple[str, int]] = set()
                read_deletion_windows: Set[Tuple[str, int]] = set()
                read_clip_windows: Set[Tuple[str, int]] = set()
                
                # Global flags for ReadsWithoutINDELorClip calculation
                read_has_significant_indel = False
                read_has_non_terminal_clip = False
                
                current_pos = read.reference_start  # 0-based
                cigartuples = read.cigartuples
                
                # Process each CIGAR operation
                for i, (op, length) in enumerate(cigartuples):
                    
                    # Check for clipping operations (S or H)
                    if op in (4, 5):  # Soft clip or Hard clip
                        # Check if this is a terminal clip that should NOT be recorded
                        is_terminal_clip = False
                        
                        # Check for start terminal clip: first operation AND read starts at chromosome start (position 0)
                        if i == 0 and read.reference_start == 0:
                            is_terminal_clip = True
                        
                        # Check for end terminal clip: last operation AND read ends at chromosome end
                        if i == len(cigartuples) - 1 and read.reference_end is not None and read.reference_end == chrom_length:
                            is_terminal_clip = True
                        
                        # Only record non-terminal clips
                        if not is_terminal_clip:
                            window_start = (current_pos // window_size) * window_size
                            window_key = (chrom, window_start)
                            
                            counts = get_window(window_key)
                            if counts is not None:
                                counts['ClipCount'] += 1
                                read_clip_windows.add(window_key)
                            
                            # Mark read as having non-terminal clip
                            read_has_non_terminal_clip = True
                    
                    # Check for Insertion (1) or Deletion (2)
                    elif op in [1, 2] and length >= min_indel_length:
                        # Mark read as having significant INDEL
                        read_has_significant_indel = True
                        
                        window_start = (current_pos // window_size) * window_size
                        window_key = (chrom, window_start)
                        
                        # Update insertion counts
                        if op == 1:
                            counts = get_window(window_key)
                            if counts is not None:
                                counts['Insertion'] += 1
                                read_insertion_windows.add(window_key)
                        
                        # Update deletion counts
                        elif op == 2:
                            # Calculate all windows affected by this deletion
                            affected_windows = calculate_affected_windows(current_pos, length, window_size)
                            for affected_window in affected_windows:
                                affected_window_key = (chrom, affected_window)
                                counts = get_window(affected_window_key)
                                if counts is not None:
                                    counts['Deletion'] += 1
                                    read_deletion_windows.add(affected_window_key)
                    
                    # Update reference position for operations that consume reference sequence
                    if op in [0, 2, 3, 7, 8]:  # M, D, N, =, X
                        current_pos += length
                
                # Update read counts for each window this read spans
                # reference_end is exclusive.  Reads without a reference span
                # (None or zero length) do not contribute TotalReads.
                read_end = read.reference_end
                if read_end is None or read_end <= read.reference_start:
                    read_windows = ()
                else:
                    read_start_window = (read.reference_start // window_size) * window_size
                    read_end_window = ((read_end - 1) // window_size) * window_size
                    read_windows = range(
                        read_start_window,
                        read_end_window + window_size,
                        window_size,
                    )

                for window_start in read_windows:
                    window_key = (chrom, window_start)
                    counts = get_window(window_key)
                    if counts is not None:
                        counts['TotalReads'] += 1
                        
                        # Count reads without significant INDEL (≥50bp) or non-terminal clip for terminal region analysis
                        # ONLY calculate ReadsWithoutINDELorClip in terminal regions.
                        # Check if this is a terminal region:
                        # - Front: window start < terminal_region_size
                        # - Back: window start >= chromosome_length - terminal_region_size
                        is_terminal_region = (window_start < terminal_region_size or 
                                            window_start >= chrom_length - terminal_region_size)
                        
                        if is_terminal_region:
                            # Check if this read has any clip or indel events in the current window
                            # A read is counted as ReadsWithoutINDELorClip if it has no clip/indel events in this window
                            read_has_clip_in_window = window_key in read_clip_windows
                            read_has_indel_in_window = (window_key in read_insertion_windows or 
                                                      window_key in read_deletion_windows)
                            
                            if not read_has_clip_in_window and not read_has_indel_in_window:
                                counts['ReadsWithoutINDELorClip'] += 1
                
                # Update read counts for windows with insertions
                for window_key in read_insertion_windows:
                    local_windows[window_key]['ReadsWithInsertion'] += 1
                
                # Update read counts for windows with deletions
                for window_key in read_deletion_windows:
                    local_windows[window_key]['ReadsWithDeletion'] += 1
                
                # Update read counts for windows with clips
                for window_key in read_clip_windows:
                    local_windows[window_key]['ReadsWithClip'] += 1
    
    return local_windows


def process_chromosome_and_store(counters, bam_file, chrom, chrom_length,
                                 min_indel_length=50, window_size=2000,
                                 terminal_region_size=10000):
    """Process and immediately integrate one chromosome.

    The future returns only the chromosome name, so completed futures do not
    retain a second outer dictionary containing all chromosome windows.
    """
    local_windows = process_chromosome_integrated(
        bam_file,
        chrom,
        chrom_length,
        min_indel_length,
        window_size,
        terminal_region_size
    )
    counters.store_chromosome_windows(chrom, local_windows)
    return chrom


def classify_error_type(anomaly_reasons):
    if "Low depth" in anomaly_reasons:
        return "low_depth"
    if "High depth" in anomaly_reasons:
        return "high_depth"
    if "insertion reads ratio" in anomaly_reasons or "High insertion reads ratio" in anomaly_reasons or \
       "deletion reads ratio" in anomaly_reasons or "High deletion reads ratio" in anomaly_reasons:
        return "indel"
    if "clip ratio" in anomaly_reasons or "High clip ratio" in anomaly_reasons:
        return "clip"
    if "Terminal region anomaly" in anomaly_reasons:
        return "terminal"
    return "other"


def identify_anomalous_regions(
    windows: Dict[Tuple[str, int], Dict[str, int]], 
    depth_data: Dict[Tuple[str, int], float],
    chrom_lengths: Dict[str, int],
    window_size: int = 2000,
    terminal_region_size: int = 10000
) -> List[Dict]:
    """
    Identify anomalous regions based on the integrated criteria including spanning reads.
    """
    anomalous_regions = []
    # Include zero-depth windows in the genome-wide mean.  Returning early for
    # an all-zero dataset previously discarded independent non-spanning,
    # clipping and INDEL evidence and could falsely report no SEs.
    all_depths = list(depth_data.values())
    mean_depth = float(np.mean(all_depths)) if all_depths else 0.0
    depth_lower_bound = mean_depth * 0.20
    depth_upper_bound = mean_depth * 2.0
    for (chrom, window_start), counts in iter_complete_windows(
            windows, chrom_lengths, window_size):
        chrom_length = chrom_lengths.get(chrom, 0)
        window_end = min(window_start + window_size, chrom_length)
        local_depth = depth_data.get((chrom, window_start), 0.0)
        is_terminal = (window_start < terminal_region_size) or (window_start >= chrom_length - terminal_region_size)
        anomaly_tags = []
        # Check non-spanning first
        has_spanning = counts.get('HasSpanningReads', True)
        if not has_spanning:
            anomaly_tags.append("non_spanning")
        # Other anomaly logic (clip, indel, depth) as before
        if not is_terminal:
            if mean_depth > 0:
                if local_depth < depth_lower_bound:
                    anomaly_tags.append("low_depth")
                elif local_depth > depth_upper_bound:
                    anomaly_tags.append("high_depth")
            total_reads = counts['TotalReads']
            clip_count = counts['ClipCount']
            if total_reads > 0 and (clip_count / total_reads) > 0.5:
                anomaly_tags.append("clip")
            insertion_reads = counts['ReadsWithInsertion']
            if total_reads > 0 and (insertion_reads / total_reads) > 0.5:
                anomaly_tags.append("indel")
            deletion_reads = counts['ReadsWithDeletion']
            if total_reads > 0 and (deletion_reads / total_reads) > 0.5:
                anomaly_tags.append("indel")
        else:
            reads_without_indel_clip = counts['ReadsWithoutINDELorClip']
            total_reads = counts['TotalReads']
            if total_reads > 0:
                ratio = reads_without_indel_clip / total_reads
                if reads_without_indel_clip < 1 or ratio <= 0.1:
                    if counts['Insertion'] > 0 or counts['Deletion'] > 0:
                        anomaly_tags.append("indel")
                    else:
                        anomaly_tags.append("clip")
        # If any anomaly, output this window
        if anomaly_tags:
            anomalous_regions.append({
                'chromosome': chrom,
                'start': window_start + 1,
                'end': window_end,
                'is_terminal': is_terminal,
                'depth': local_depth,
                'clip_count': counts['ClipCount'],
                'insertion_count': counts['Insertion'],
                'deletion_count': counts['Deletion'],
                'reads_without_indel_clip': counts['ReadsWithoutINDELorClip'],
                'total_reads': counts['TotalReads'],
                'has_spanning_reads': counts.get('HasSpanningReads', True),
                'anomaly_tags': ','.join(anomaly_tags)
            })
    return anomalous_regions


def terminal_reads_analysis(
    windows: Dict[Tuple[str, int], Dict[str, int]],
    depth_data: Dict[Tuple[str, int], float],
    chrom_lengths: Dict[str, int],
    terminal_region_size: int = 10000
):
    """
    For terminal regions (ratio<0.1), assign to indel if Insertion or Deletion >0, else to clip.
    Returns: dict chrom->dict{'indel': count, 'clip': count}, total_indel, total_clip
    """
    from collections import defaultdict
    terminal_type_counts = defaultdict(lambda: {'indel': 0, 'clip': 0})
    total_indel = 0
    total_clip = 0
    for (chrom, window_start), counts in windows.items():
        chrom_length = chrom_lengths.get(chrom, 0)
        is_terminal = (window_start < terminal_region_size) or (window_start >= chrom_length - terminal_region_size)
        if not is_terminal:
            continue
        window_total_reads = counts['TotalReads']
        window_reads_without = counts['ReadsWithoutINDELorClip']
        if window_total_reads > 0:
            ratio = window_reads_without / window_total_reads
            if ratio < 0.1:
                if counts['Insertion'] > 0 or counts['Deletion'] > 0:
                    terminal_type_counts[chrom]['indel'] += 1
                    total_indel += 1
                else:
                    terminal_type_counts[chrom]['clip'] += 1
                    total_clip += 1
    return terminal_type_counts, total_indel, total_clip


def save_results(
    windows: Dict[Tuple[str, int], Dict[str, int]], 
    depth_data: Dict[Tuple[str, int], float],
    anomalous_regions: List[Dict],
    chrom_lengths: Dict[str, int],
    output_prefix: str,
    window_size: int,
    terminal_region_size: int
):
    """
    Save all results to output files.
    """
    # Save detailed window statistics
    stats_file = f"{output_prefix}.window_stats.tsv"
    with open(stats_file, 'w') as f:
        f.write("Chromosome\tStart\tEnd\tDepth\tClipCount\tInsertion\tDeletion\t"
                "ReadsWithInsertion\tReadsWithDeletion\tReadsWithClip\t"
                "TotalReads\tReadsWithoutINDELorClip\tHasSpanningReads\n")
        for (chrom, window_start), counts in iter_complete_windows(
                windows, chrom_lengths, window_size,
                sort_chromosomes=True):
            chrom_length = chrom_lengths.get(chrom, window_start + window_size)
            window_end = min(window_start + window_size, chrom_length)
            depth = depth_data.get((chrom, window_start), 0.0)
            is_terminal = (window_start < terminal_region_size) or (window_start >= chrom_length - terminal_region_size)
            if is_terminal:
                reads_without_indel_clip = str(counts['ReadsWithoutINDELorClip'])
            else:
                reads_without_indel_clip = '-'
            has_spanning = counts.get('HasSpanningReads', True)
            f.write(f"{chrom}\t{window_start + 1}\t{window_end}\t{depth:.4f}\t"
                   f"{counts['ClipCount']}\t{counts['Insertion']}\t{counts['Deletion']}\t"
                   f"{counts['ReadsWithInsertion']}\t{counts['ReadsWithDeletion']}\t{counts['ReadsWithClip']}\t{counts['TotalReads']}\t"
                   f"{reads_without_indel_clip}\t{has_spanning}\n")
    print(f"Window statistics saved to: {stats_file}")

    # Save anomalous regions (no AnomalyReasons column)
    anomalies_file = f"{output_prefix}.anomalous_regions.tsv"
    with open(anomalies_file, 'w') as f:
        f.write("Chromosome\tStart\tEnd\tIsTerminal\tDepth\tClipCount\t"
                "InsertionCount\tDeletionCount\tReadsWithoutINDELClip\t"
                "TotalReads\tHasSpanningReads\tAnomalyTags\n")
        for region in anomalous_regions:
            f.write(f"{region['chromosome']}\t{region['start']}\t{region['end']}\t"
                   f"{region['is_terminal']}\t{region['depth']:.4f}\t"
                   f"{region['clip_count']}\t{region['insertion_count']}\t"
                   f"{region['deletion_count']}\t{region['reads_without_indel_clip']}\t"
                   f"{region['total_reads']}\t{region['has_spanning_reads']}\t"
                   f"{region['anomaly_tags']}\n")
    print(f"Anomalous regions saved to: {anomalies_file}")

    # Generate anomalous regions BED file (chr, start-1, end)
    anomalies_bed_file = f"{output_prefix}.anomalous_regions.bed"
    with open(anomalies_file, 'r') as fin, open(anomalies_bed_file, 'w') as fout:
        next(fin)  # skip header
        for line in fin:
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                start = max(0, int(fields[1]) - 1)
                fout.write(f"{fields[0]}\t{start}\t{fields[2]}\n")
    print(f"Anomalous regions BED saved to: {anomalies_bed_file}")

# se.result.bed: First three columns of anomalous_regions.tsv, start-1, then merge

def generate_se_result_bed(anomalous_regions_file, output_bed):
    """Generate se.result.bed by merging se.anomalous_regions.tsv (chr, start-1, end) with bedtools merge."""
    output_bed = os.path.abspath(output_bed)
    output_dir = os.path.dirname(output_bed)
    output_name = os.path.basename(output_bed)
    bed_tmp = None
    output_tmp = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            dir=output_dir,
            prefix=f".{output_name}.input.",
            suffix=".bed",
            delete=False,
        ) as fout:
            bed_tmp = fout.name
            with open(anomalous_regions_file, 'r') as fin:
                next(fin)  # skip header
                for line in fin:
                    fields = line.strip().split('\t')
                    if len(fields) >= 3:
                        start = max(0, int(fields[1]) - 1)
                        fout.write(f"{fields[0]}\t{start}\t{fields[2]}\n")

        with tempfile.NamedTemporaryFile(
            mode="w",
            dir=output_dir,
            prefix=f".{output_name}.output.",
            suffix=".bed",
            delete=False,
        ) as fout:
            output_tmp = fout.name
            cmd = ["bedtools", "merge", "-i", bed_tmp]
            subprocess.run(cmd, stdout=fout, check=True)

        os.replace(output_tmp, output_bed)
        output_tmp = None
    finally:
        for temporary_path in (bed_tmp, output_tmp):
            if temporary_path is None:
                continue
            try:
                os.remove(temporary_path)
            except FileNotFoundError:
                pass
    print(f"Merged SE result BED saved to: {output_bed}")

# save_anomaly_summary: Change spanning column to non-spanning

def save_anomaly_summary(anomalous_regions, chrom_lengths, output_prefix):
    from collections import defaultdict
    summary = defaultdict(lambda: defaultdict(int))
    nr_se = defaultdict(int)
    tag_types = ["low_depth", "high_depth", "indel", "clip", "non_spanning"]
    for region in anomalous_regions:
        chrom = region['chromosome']
        tags = [t.strip() for t in region['anomaly_tags'].split(',')] if region['anomaly_tags'] else []
        for tag in tag_types:
            if tag in tags:
                summary[chrom][tag] += 1
        nr_se[chrom] += 1
    all_chroms = sorted(chrom_lengths.keys())
    summary_file = f"{output_prefix}.anomaly_summary.tsv"
    with open(summary_file, 'w') as f:
        f.write("Chromosome\tlow_depth\thigh_depth\tindel\tclip\tnon-spanning\tnr-SE\n")
        total_row = ["Total", 0, 0, 0, 0, 0, 0]
        for chrom in all_chroms:
            row = [chrom]
            row.append(str(summary[chrom]["low_depth"]))
            row.append(str(summary[chrom]["high_depth"]))
            row.append(str(summary[chrom]["indel"]))
            row.append(str(summary[chrom]["clip"]))
            row.append(str(summary[chrom]["non_spanning"]))
            row.append(str(nr_se[chrom]))
            for i in range(1, 6):
                total_row[i] += int(row[i])
            total_row[6] += int(row[6])
            f.write("\t".join(row) + "\n")
        total_row = [str(x) for x in total_row]
        f.write("\t".join(total_row) + "\n")
    print(f"Anomaly summary saved to: {summary_file}")


def publish_pipeline_compatibility_outputs(output_dir, output_prefix):
    """Expose legacy ``se.*`` names for sas_pipeline's default invocation.

    The pipeline passes ``structural_errors`` as the positional prefix, while
    existing downstream pipeline code and documented output names still read
    ``se.*``.  Hard links preserve those names without duplicating large files.
    Custom prefixes do not receive aliases.
    """
    if output_prefix != "structural_errors":
        return

    source_stem = output_prefix + "."
    for filename in os.listdir(output_dir):
        if not filename.startswith(source_stem):
            continue
        source = os.path.join(output_dir, filename)
        if not os.path.isfile(source):
            continue
        suffix = filename[len(output_prefix):]
        destination = os.path.join(output_dir, "se" + suffix)
        if os.path.lexists(destination):
            os.remove(destination)
        try:
            os.link(source, destination)
        except OSError:
            # Both paths are in the same directory; a relative symlink remains
            # valid if the complete output directory is moved.
            os.symlink(os.path.basename(source), destination)


def main():
    parser = argparse.ArgumentParser(
        description="Run standalone SAS structural-error analysis.",
        formatter_class=CLIHelpFormatter,
        epilog="""
Examples:
  sas-se-analysis input.bam reference.fa sample_name
  sas-se-analysis input.bam reference.fa sample_name -d results -w 1000 -t 8
        """,
    )

    inputs = parser.add_argument_group("required inputs")
    inputs.add_argument("bam_file", metavar="BAM", help="Input BAM file")
    inputs.add_argument("reference_fa", metavar="FASTA", help="Reference FASTA file")
    inputs.add_argument(
        "output_prefix", metavar="PREFIX", help="Prefix for result filenames"
    )

    output = parser.add_argument_group("output options")
    output.add_argument(
        "-d", "--output-dir", default=".", metavar="DIR",
        help="Directory for result files (default: current directory)",
    )

    analysis = parser.add_argument_group("analysis options")
    analysis.add_argument(
        "-w", "--window-size", type=positive_int, default=2000, metavar="BP",
        help="Main analysis window size (default: 2000 bp)",
    )
    analysis.add_argument(
        "-m", "--min-indel-length", type=positive_int, default=50, metavar="BP",
        help="Minimum INDEL length to consider (default: 50 bp)",
    )
    analysis.add_argument(
        "-t", "--threads", type=positive_int, default=4, metavar="N",
        help="Threads available to external tools such as mosdepth (default: 4)",
    )
    analysis.add_argument(
        "--workers", type=positive_int, metavar="N",
        help="Concurrent chromosome/spanning workers (default: min(threads, 8))",
    )
    analysis.add_argument(
        "-T", "--terminal-region-size", type=positive_int, default=10000,
        metavar="BP", help="Chromosome-terminal region size (default: 10000 bp)",
    )
    analysis.add_argument(
        "-W", "--spanning-window-size", type=positive_int, default=30000,
        metavar="BP", help="Spanning-read window size (default: 30000 bp)",
    )
    analysis.add_argument(
        "-S", "--spanning-step-size", type=positive_int, default=2000,
        metavar="BP", help="Spanning-read window step (default: 2000 bp)",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    
    args = parser.parse_args()
    if args.workers is None:
        args.workers = min(args.threads, 8)

    missing_dependencies = [
        name
        for name, module in (("numpy", np), ("pandas", pd), ("pysam", pysam))
        if module is None
    ]
    if missing_dependencies:
        parser.error(
            "missing required Python dependencies: "
            + ", ".join(missing_dependencies)
        )
    
    if not os.path.isfile(args.bam_file):
        parser.error(f"BAM file not found: {args.bam_file}")
    if not os.path.isfile(args.reference_fa):
        parser.error(f"reference FASTA file not found: {args.reference_fa}")
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = args.output_prefix
    full_output_prefix = os.path.join(str(output_dir), output_prefix)
    try:
        print(f"Starting integrated SE analysis for: {args.bam_file}")
        print(f"Reference: {args.reference_fa}")
        print(f"Output directory: {output_dir}")
        print(f"Output prefix: {args.output_prefix}")
        print(f"Window size: {args.window_size} bp")
        print(f"Minimum INDEL length: {args.min_indel_length} bp")
        print(f"Threads: {args.threads}")
        print(f"Concurrent SE workers: {args.workers}")
        print(f"Terminal region size: {args.terminal_region_size} bp")
        print("\n" + "="*50)
        print("STEP 1: Running mosdepth for depth analysis")
        print("="*50)
        regions_gz_file = run_mosdepth(
            args.bam_file, 
            str(full_output_prefix), 
            args.threads, 
            args.window_size
        )
        print("\n" + "="*50)
        print("STEP 2: Loading depth data")
        print("="*50)
        depth_data = load_depth_data(regions_gz_file, args.window_size)
        print(f"Loaded depth data for {len(depth_data)} windows")
        print("\n" + "="*50)
        print("STEP 3: Getting chromosome information")
        print("="*50)
        with pysam.AlignmentFile(args.bam_file, "rb") as bam:
            all_chrom_lengths = get_chromosome_lengths_from_bam(bam)
            chrom_lengths = {chrom: length for chrom, length in all_chrom_lengths.items() if length > 0}
        print(f"Found {len(chrom_lengths)} valid chromosomes")
        validate_depth_windows(depth_data, chrom_lengths, args.window_size)
        print("Depth windows match the BAM reference coordinate space")
        print("\n" + "="*50)
        print("STEP 4: Processing chromosomes for integrated SE analysis")
        print("="*50)
        global_counters = IntegratedSECounters(chrom_lengths, args.window_size)
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            future_to_chrom = {
                executor.submit(
                    process_chromosome_and_store,
                    global_counters,
                    args.bam_file,
                    chrom,
                    chrom_length,
                    args.min_indel_length,
                    args.window_size,
                    args.terminal_region_size
                ): chrom
                for chrom, chrom_length in chrom_lengths.items()
            }
            completed_chroms = 0
            for future in as_completed(future_to_chrom):
                chrom = future_to_chrom[future]
                try:
                    future.result()
                    completed_chroms += 1
                    print(f"Completed chromosome {chrom} ({completed_chroms}/{len(chrom_lengths)})")
                except Exception as e:
                    raise RuntimeError(
                        f"Integrated SE worker failed for chromosome {chrom}"
                    ) from e
        
        print("\n" + "="*50)
        print("STEP 5: Analyzing spanning reads")
        print("="*50)
        # Generate FAI file if it doesn't exist
        fai_file = args.reference_fa + '.fai'
        if not os.path.exists(fai_file):
            print(f"Generating index for reference file: {args.reference_fa}")
            try:
                subprocess.run(['samtools', 'faidx', args.reference_fa], check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error generating fai file: {e}", file=sys.stderr)
                sys.exit(1)
        # Use the merged+terminal-corrected bed for 2kb window assignment
        merged_bed = refresh_spanning_outputs(
            args.bam_file,
            fai_file,
            str(full_output_prefix),
            args.spanning_window_size,
            args.spanning_step_size,
            args.workers
        )
        merged_regions = load_merged_non_spanning_regions_from_bed(merged_bed)
        print(f"Loaded {len(merged_regions)} merged non-spanning regions from {merged_bed}")
        non_spanning_windows = assign_merged_regions_to_2k_windows(merged_regions, args.window_size)
        print(f"Assigned {len(non_spanning_windows)} 2kb non-spanning windows from merged regions (0-based)")
        # Mark non-spanning windows in global counters
        global_counters.mark_non_spanning_windows(non_spanning_windows)
        
        print("\n" + "="*50)
        print("STEP 6: Integrating depth data")
        print("="*50)
        # Depth remains in its single dedicated mapping.  Duplicating it into
        # every materialized counter window only increases memory and is not
        # used by classification or output.
        print("\n" + "="*50)
        print("STEP 7: Identifying anomalous regions")
        print("="*50)
        anomalous_regions = identify_anomalous_regions(
            global_counters.windows,
            depth_data,
            chrom_lengths,
            args.window_size,
            args.terminal_region_size
        )
        print(f"Found {len(anomalous_regions)} anomalous regions")
        print("\n" + "="*50)
        print("STEP 8: Saving results")
        print("="*50)
        save_results(
            global_counters.windows,
            depth_data,
            anomalous_regions,
            chrom_lengths,
            str(full_output_prefix),
            args.window_size,
            args.terminal_region_size
        )


        anomalies_file = f"{full_output_prefix}.anomalous_regions.tsv"
        merge_bed_file = f"{full_output_prefix}.result.bed"
        generate_se_result_bed(anomalies_file, merge_bed_file)

        save_anomaly_summary(
            anomalous_regions,
            chrom_lengths,
            str(full_output_prefix)
        )
        publish_pipeline_compatibility_outputs(
            str(output_dir),
            output_prefix
        )
        print("\n" + "="*50)
        print("ANALYSIS COMPLETE")
        print("="*50)
        total_windows = sum(
            (length + args.window_size - 1) // args.window_size
            for length in chrom_lengths.values()
        )
        print(f"Total windows analyzed: {total_windows}")
        print(f"Anomalous regions found: {len(anomalous_regions)}")
        print(f"Results saved to: {output_dir}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main() 
