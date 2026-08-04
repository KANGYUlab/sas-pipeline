#!/usr/bin/env python3
"""
Processes a BAM file based on regions defined in a BED file to analyze pileup data.

This script takes a BED file specifying genomic regions, a BAM alignment file, and a
reference FASTA file. Nearby regions are grouped into bounded 'samtools mpileup'
queries, then split back into their original BED intervals before pileup cleaning
and statistics are calculated.

The processing is done in parallel using multiple processes to speed up the analysis.
"""
import sys
import argparse
import subprocess
import os
import shutil
import tempfile
import heapq
from collections import OrderedDict
from multiprocessing import Pool

CLI_VERSION = "2.0.0"


class CLIHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Keep examples and explicitly documented defaults readable."""


def positive_int(value):
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be greater than 0")
    return parsed


def nonnegative_int(value):
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer") from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be 0 or greater")
    return parsed

def process_pileup_lines(raw_pileup_lines, region_start, region_end, verbose=False):
    """
    Calculate statistics from the pileup rows belonging to one BED interval.

    This version implements the proper cross-row alignment logic with correct processing order.
    """
    HIGH_DEPTH_THRESHOLD = 1000

    # Parse the region once. IR analysis units are at most 50 bp, so retaining
    # these rows in memory is bounded while avoiding a second disk read.
    total_depth = 0
    line_count = 0
    first_parts = None
    all_lines = []
    pileup_cols = []
    for line in raw_pileup_lines:
        parts = line.strip().split('\t')
        if first_parts is None and parts:
            first_parts = parts
        if len(parts) >= 4 and parts[3].isdigit():
            total_depth += int(parts[3])
        line_count += 1
        all_lines.append(parts)
        pileup_cols.append(list(parts[4]) if len(parts) >= 5 else [])

    average_depth = (total_depth / line_count) if line_count > 0 else 0.0

    # If depth is too high, short‑circuit and mark region as unsupported for this platform.
    if average_depth > HIGH_DEPTH_THRESHOLD:
        print(f"Depth too high ({average_depth:.2f}) for region {region_start}-{region_end}; marking as unsupported.", file=sys.stderr)
        if not first_parts:
            first_parts = ["chrN", str(region_start), str(region_end), ".", "."]
        # Ensure we have at least 5 columns to mimic original structure
        while len(first_parts) < 5:
            first_parts.append(".")
        cleaned_pileup = "DEPTH_TOO_HIGH"
        min_length = "0"
        output_parts = first_parts[:5] + [cleaned_pileup, first_parts[0], str(region_start), str(region_end), min_length, f"{average_depth:.2f}"]
        return ['\t'.join(output_parts)]

    if not pileup_cols or not all_lines:
        return []

    num_rows = len(pileup_cols)
    if num_rows == 0:
        return []
    max_iterations = 100000  # Safety limit to prevent infinite loops
    iteration_count = 0

    # Processing loop with correct order: ^ → $ → +- → base symbols
    while iteration_count < max_iterations:
        changed = False
        iteration_count += 1
        
        # Recalculate max length
        max_len = max(len(col) for col in pileup_cols) if pileup_cols else 0
        if max_len == 0:
            break

        # Phase 1: Handle ^ symbols (left-to-right scan)
        for i in range(max_len):
            for row_idx in range(num_rows):
                if i < len(pileup_cols[row_idx]) and pileup_cols[row_idx][i] == '^':
                    # Delete '^' and its quality char
                    if i + 1 < len(pileup_cols[row_idx]):
                        del pileup_cols[row_idx][i:i+2]  # Delete '^' and quality char
                    else:
                        del pileup_cols[row_idx][i]  # Delete only '^'
                    
                    # If not in the first row, delete corresponding position in subsequent rows
                    if row_idx > 0:
                        # Check if position i still exists after deletion
                        if i < len(pileup_cols[row_idx]):
                            # The character that was originally the second character after ^ is now at position i
                            # Count which number (k-th) this '.', ',' or '*' is (from left to right)
                            position_count = 0
                            for j in range(min(i + 1, len(pileup_cols[row_idx]))):  # Include position i itself
                                if j < len(pileup_cols[row_idx]) and pileup_cols[row_idx][j] in '.,*':
                                    position_count += 1
                            
                            # Delete the character at position i in current row if it exists
                            if i < len(pileup_cols[row_idx]):
                                del pileup_cols[row_idx][i]
                            
                            # Delete the k-th '.', ',' or '*' in all subsequent rows
                            for r in range(row_idx + 1, num_rows):
                                count = 0
                                j = 0
                                while j < len(pileup_cols[r]):
                                    if pileup_cols[r][j] in '.,*':
                                        count += 1
                                        if count == position_count:
                                            del pileup_cols[r][j]
                                            break
                                    j += 1
                    changed = True
                    break
            if changed:
                break
        
        if changed:
            continue

        # Phase 2: Handle $ symbols (right-to-left scan)
        max_len = max(len(col) for col in pileup_cols) if pileup_cols else 0
        for i in range(max_len - 1, -1, -1):
            for row_idx in range(num_rows):
                if i < len(pileup_cols[row_idx]) and pileup_cols[row_idx][i] == '$':
                    if row_idx == num_rows - 1:
                        # If in the last row, only delete '$'
                        del pileup_cols[row_idx][i]
                    else:
                        # Delete '$' first
                        del pileup_cols[row_idx][i]
                        
                        if i > 0 and i - 1 < len(pileup_cols[row_idx]):
                            # The character to the left of '$' is guaranteed to be '.', ',' or '*'
                            # Count which number (k-th from right) this '.', ',' or '*' is
                            position_count = 0
                            for j in range(len(pileup_cols[row_idx]) - 1, max(i - 2, -1), -1):  # From right to left, including position i-1
                                if j >= 0 and j < len(pileup_cols[row_idx]) and pileup_cols[row_idx][j] in '.,*':
                                    position_count += 1
                            
                            # Delete the character to the left of '$' if it exists
                            if i - 1 < len(pileup_cols[row_idx]):
                                del pileup_cols[row_idx][i-1]
                            
                            # Delete the k-th '.', ',' or '*' from right to left in all previous rows
                            for r in range(row_idx):
                                count = 0
                                for j in range(len(pileup_cols[r]) - 1, -1, -1):
                                    if j >= 0 and j < len(pileup_cols[r]) and pileup_cols[r][j] in '.,*':
                                        count += 1
                                        if count == position_count:
                                            del pileup_cols[r][j]
                                            break
                    changed = True
                    break
            if changed:
                break
        
        if changed:
            continue

        # Phase 3: Handle +/- indels
        max_len = max(len(col) for col in pileup_cols) if pileup_cols else 0
        for i in range(max_len):
            has_indel_in_column = any(i < len(pileup_cols[r]) and pileup_cols[r][i] in "+-" for r in range(num_rows))
            
            if has_indel_in_column:
                # Step 1: Find and delete ALL indel descriptors in this column
                for row_idx in range(num_rows):
                    if i < len(pileup_cols[row_idx]) and pileup_cols[row_idx][i] in "+-":
                        j = i + 1
                        num_str = ""
                        while j < len(pileup_cols[row_idx]) and pileup_cols[row_idx][j].isdigit():
                            num_str += pileup_cols[row_idx][j]
                            j += 1
                        
                        if num_str:
                            indel_length = int(num_str)
                            total_indel_len = 1 + len(num_str) + indel_length
                            end_pos = min(i + total_indel_len, len(pileup_cols[row_idx]))
                            if i < len(pileup_cols[row_idx]) and end_pos <= len(pileup_cols[row_idx]):
                                del pileup_cols[row_idx][i:end_pos]
                        else:
                            # Just delete the +/- if no number follows
                            if i < len(pileup_cols[row_idx]):
                                del pileup_cols[row_idx][i]
                
                # Step 2: Delete preceding column across ALL rows (the . or , before indel)
                if i > 0:
                    pos_to_delete = i - 1
                    for r in range(num_rows):
                        if pos_to_delete < len(pileup_cols[r]):
                            del pileup_cols[r][pos_to_delete]
                
                changed = True
                break
        
        if changed:
            continue

        # Phase 4: Handle base mismatches and deletions
        max_len = max(len(col) for col in pileup_cols) if pileup_cols else 0
        for i in range(max_len):
            is_mismatch_column = any(
                i < len(pileup_cols[r]) and (
                    pileup_cols[r][i].lower() in {'a', 't', 'c', 'g', 'n'} or 
                    pileup_cols[r][i] == '*'
                )
                for r in range(num_rows)
            )
            
            if is_mismatch_column:
                # Delete this column from all rows
                for r in range(num_rows):
                    if i < len(pileup_cols[r]):
                        del pileup_cols[r][i]
                changed = True
                break
        
        if changed:
            continue
            
        # If no changes were made in any phase, we're done
        break

    # Generate final output with BED region coordinates
    processed_lines = []
    
    # Calculate the minimum length of all cleaned pileup strings
    valid_lengths = [len(pileup_cols[i]) for i in range(len(pileup_cols)) if i < len(all_lines) and len(all_lines[i]) >= 5]
    min_length_all = min(valid_lengths) if valid_lengths else 0
    
    if verbose:
        print(
            f"Debug: Total lines available: {len(all_lines)}, "
            f"Total pileup columns: {len(pileup_cols)}",
            file=sys.stderr,
        )
    
    # Only process the first line of each region
    if len(all_lines) > 0 and len(pileup_cols) > 0:
        row_idx = 0
        parts = all_lines[row_idx]
        if len(parts) >= 5 and row_idx < len(pileup_cols):
            cleaned_pileup = ''.join(pileup_cols[row_idx])
            first_col = parts[0]
            # Use minimum length of all cleaned pileup strings for column 10
            min_length = str(min_length_all)
            # Use region_start and region_end for columns 8-9 instead of parts[1] and parts[2]
            output_parts = parts[:5] + [cleaned_pileup, first_col, str(region_start), str(region_end), min_length, f"{average_depth:.2f}"]
            processed_lines.append('\t'.join(output_parts))
            if verbose:
                print("Debug: Processed first line only", file=sys.stderr)
        else:
            if verbose:
                print(
                    f"Debug: Skipping first row, len(parts)={len(parts)}, "
                    f"row_idx < len(pileup_cols): {row_idx < len(pileup_cols)}",
                    file=sys.stderr,
                )
    else:
        if verbose:
            print("Debug: No data to process", file=sys.stderr)
            
    if verbose:
        print(f"Debug: Total processed lines: {len(processed_lines)}", file=sys.stderr)
    return processed_lines


def process_pileup_from_file(raw_pileup_path, region_start, region_end, verbose=False):
    """Backward-compatible file wrapper around :func:`process_pileup_lines`."""
    with open(raw_pileup_path, 'r') as pileup_handle:
        return process_pileup_lines(
            pileup_handle,
            region_start,
            region_end,
            verbose=verbose,
        )

DEFAULT_BATCH_GAP = 10_000
DEFAULT_MAX_BATCH_SPAN = 1_000_000
DEFAULT_MAX_REGIONS_PER_BATCH = 1_000
DEFAULT_MAX_OPEN_RAW_FILES = 64
COPY_BUFFER_SIZE = 8 * 1024 * 1024
MAX_MPILEUP_STDERR_BYTES = 64 * 1024


def create_run_temp_dir(temp_parent=None):
    """Create a private temporary directory for one pipeline invocation.

    A caller-supplied path is a parent directory, never a directory owned by
    this function.  Only the unique child returned here may be removed during
    cleanup.  With no parent, ``tempfile`` selects the platform temp directory.
    """
    if temp_parent is not None:
        os.makedirs(temp_parent, exist_ok=True)
    return tempfile.mkdtemp(prefix="sas_pileup_", dir=temp_parent)


def read_stderr_tail(stderr_file, max_bytes=MAX_MPILEUP_STDERR_BYTES):
    """Return at most ``max_bytes`` from the end of a binary stderr file."""
    if max_bytes <= 0:
        raise ValueError("max_bytes must be > 0")

    stderr_file.flush()
    stderr_file.seek(0, os.SEEK_END)
    stderr_size = stderr_file.tell()
    stderr_file.seek(max(0, stderr_size - max_bytes), os.SEEK_SET)
    error_output = stderr_file.read(max_bytes).decode("utf-8", errors="replace")
    if stderr_size > max_bytes:
        return (
            f"[stderr truncated; showing last {max_bytes} bytes]\n"
            f"{error_output}"
        )
    return error_output


def build_region_batches(
    regions,
    max_gap=DEFAULT_BATCH_GAP,
    max_span=DEFAULT_MAX_BATCH_SPAN,
    max_regions=DEFAULT_MAX_REGIONS_PER_BATCH,
):
    """Group nearby regions for mpileup without changing their analysis units.

    Each region is ``(chrom, start, end, input_index, reserved)``; the final
    field is retained for compatibility with earlier helper callers.
    A batch never crosses a chromosome, the gap between consecutive merged
    ranges may not exceed ``max_gap``, and the complete query span may not
    exceed ``max_span``. ``max_regions`` bounds per-worker bookkeeping.
    """
    if max_gap < 0:
        raise ValueError("max_gap must be >= 0")
    if max_span <= 0:
        raise ValueError("max_span must be > 0")
    if max_regions <= 0:
        raise ValueError("max_regions must be > 0")

    sorted_regions = sorted(regions, key=lambda region: (region[0], region[1], region[2], region[3]))
    batches = []
    current_batch = []
    batch_chrom = None
    batch_start = None
    batch_end = None

    for region in sorted_regions:
        chrom, start, end = region[:3]
        if start < 0 or end <= start:
            raise ValueError(f"Invalid BED region: {chrom}\t{start}\t{end}")

        if not current_batch:
            should_split = False
        else:
            gap = max(0, start - batch_end)
            proposed_end = max(batch_end, end)
            should_split = (
                chrom != batch_chrom
                or gap > max_gap
                or proposed_end - batch_start > max_span
                or len(current_batch) >= max_regions
            )

        if should_split:
            batches.append(current_batch)
            current_batch = []

        if not current_batch:
            batch_chrom = chrom
            batch_start = start
            batch_end = end
        else:
            batch_end = max(batch_end, end)

        current_batch.append(region)

    if current_batch:
        batches.append(current_batch)

    return batches


class BoundedAppendWriter:
    """Append text while keeping only a bounded number of files open."""

    def __init__(self, max_open_files=DEFAULT_MAX_OPEN_RAW_FILES):
        if max_open_files <= 0:
            raise ValueError("max_open_files must be > 0")
        self.max_open_files = max_open_files
        self.handles = OrderedDict()

    def write(self, file_path, text):
        handle = self.handles.pop(file_path, None)
        if handle is None:
            if len(self.handles) >= self.max_open_files:
                _, oldest_handle = self.handles.popitem(last=False)
                oldest_handle.close()
            handle = open(file_path, 'a')
        self.handles[file_path] = handle
        handle.write(text)

    def close(self):
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()


def dispatch_pileup_lines(lines, regions, raw_paths, max_open_files=DEFAULT_MAX_OPEN_RAW_FILES):
    """Split one coordinate-sorted mpileup stream back into original BED regions."""
    sorted_regions = sorted(regions, key=lambda region: (region[1], region[2], region[3]))
    next_region = 0
    active_regions = []
    writer = BoundedAppendWriter(max_open_files=max_open_files)

    try:
        for line in lines:
            fields = line.split('\t', 2)
            if len(fields) < 2:
                continue

            try:
                chrom = fields[0]
                position = int(fields[1]) - 1  # mpileup is 1-based; BED is 0-based.
            except ValueError:
                continue

            active_regions = [
                region for region in active_regions
                if region[0] == chrom and region[2] > position
            ]

            while next_region < len(sorted_regions):
                region = sorted_regions[next_region]
                if region[0] != chrom or region[1] > position:
                    break
                next_region += 1
                if region[2] > position:
                    active_regions.append(region)

            output_line = line if line.endswith('\n') else line + '\n'
            for region in active_regions:
                if region[0] == chrom and region[1] <= position < region[2]:
                    writer.write(raw_paths[region[3]], output_line)
    finally:
        writer.close()


def spool_pileup_lines(lines, regions, spool_path):
    """Write one batch stream once and record the byte slice for each region.

    A pileup position is stored only once even when BED intervals overlap.  As
    coordinate-sorted positions belonging to an interval are contiguous in the
    batch stream, a pair of byte offsets is sufficient to reconstruct exactly
    the rows seen by that interval.  This replaces one raw temporary file per
    50 bp window with one file per mpileup batch.
    """
    sorted_regions = sorted(regions, key=lambda region: (region[1], region[2], region[3]))
    next_region = 0
    active_regions = []
    offsets = {region[3]: [None, None] for region in regions}

    with open(spool_path, "wb") as spool_handle:
        for line in lines:
            fields = line.split('\t', 2)
            if len(fields) < 2:
                continue

            try:
                chrom = fields[0]
                position = int(fields[1]) - 1
            except ValueError:
                continue

            active_regions = [
                region for region in active_regions
                if region[0] == chrom and region[2] > position
            ]

            while next_region < len(sorted_regions):
                region = sorted_regions[next_region]
                if region[0] != chrom or region[1] > position:
                    break
                next_region += 1
                if region[2] > position:
                    active_regions.append(region)

            matching_regions = [
                region for region in active_regions
                if region[0] == chrom and region[1] <= position < region[2]
            ]
            if not matching_regions:
                continue

            output_line = line if line.endswith('\n') else line + '\n'
            encoded_line = output_line.encode("utf-8")
            line_start = spool_handle.tell()
            spool_handle.write(encoded_line)
            line_end = spool_handle.tell()

            for region in matching_regions:
                region_offsets = offsets[region[3]]
                if region_offsets[0] is None:
                    region_offsets[0] = line_start
                region_offsets[1] = line_end

    return {
        input_index: tuple(region_offsets)
        for input_index, region_offsets in offsets.items()
        if region_offsets[0] is not None
    }


def read_pileup_slice(spool_handle, offsets):
    """Read and decode one region's bounded slice from a batch spool."""
    start_offset, end_offset = offsets
    spool_handle.seek(start_offset)
    raw_data = spool_handle.read(end_offset - start_offset)
    return raw_data.decode("utf-8").splitlines(keepends=True)


def process_region_batch(args):
    """Run one mpileup query for nearby regions, then analyze each original region."""
    if len(args) == 4:
        batch, bam_file, ref_fasta, batch_id = args
        verbose = False
        output_hint = batch[0][4]
        temp_dir = (os.path.dirname(output_hint) if output_hint else None) or tempfile.gettempdir()
    elif len(args) == 5:
        batch, bam_file, ref_fasta, batch_id, verbose = args
        output_hint = batch[0][4]
        temp_dir = (os.path.dirname(output_hint) if output_hint else None) or tempfile.gettempdir()
    else:
        batch, bam_file, ref_fasta, batch_id, verbose, temp_dir = args
    chrom = batch[0][0]
    batch_start = min(region[1] for region in batch)
    batch_end = max(region[2] for region in batch)
    batch_bed_file = os.path.join(temp_dir, f"batch_{batch_id:09d}.bed")
    raw_batch_file = os.path.join(temp_dir, f"raw_batch_{batch_id:09d}.txt")
    result_file = os.path.join(temp_dir, f"batch_result_{batch_id:09d}.indexed")

    with open(batch_bed_file, 'w') as bed_out:
        for region in sorted(batch, key=lambda item: (item[1], item[2], item[3])):
            bed_out.write(f"{region[0]}\t{region[1]}\t{region[2]}\n")

    region_str = f"{chrom}:{batch_start + 1}-{batch_end}"
    cmd = [
        'samtools', 'mpileup', '-a', '-A', '-B', '-d', '0',
        '-r', region_str, '-l', batch_bed_file,
        '-f', ref_fasta, '-Q', '0', '--ff', '0', '-x', bam_file
    ]
    print(
        f"Processing batch {batch_id}: {region_str}, {len(batch)} original regions",
        file=sys.stderr,
    )

    completed = False
    try:
        with tempfile.TemporaryFile(mode='w+b') as stderr_file:
            with subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=stderr_file,
                universal_newlines=True,
                bufsize=1024 * 1024,
            ) as proc:
                region_offsets = spool_pileup_lines(
                    proc.stdout,
                    batch,
                    raw_batch_file,
                )
                return_code = proc.wait()

            if return_code != 0:
                error_output = read_stderr_tail(stderr_file)
                raise RuntimeError(
                    f"samtools mpileup failed for batch {batch_id} ({region_str}): {error_output}"
                )

        output_count = 0
        with open(raw_batch_file, "rb") as spool_handle, open(
            result_file, "w", encoding="utf-8"
        ) as output_handle:
            for _chrom_name, start, end, input_index, _unused in sorted(
                batch, key=lambda item: item[3]
            ):
                offsets = region_offsets.get(input_index)
                if offsets is None:
                    continue
                formatted_output = process_pileup_lines(
                    read_pileup_slice(spool_handle, offsets),
                    start,
                    end,
                    verbose=verbose,
                )
                for output_line in formatted_output:
                    output_handle.write(f"{input_index}\t{output_line}\n")
                    output_count += 1

        completed = True
        if output_count:
            return [result_file]
        os.remove(result_file)
        return []
    finally:
        if os.path.exists(batch_bed_file):
            os.remove(batch_bed_file)
        if os.path.exists(raw_batch_file):
            os.remove(raw_batch_file)
        if not completed and os.path.exists(result_file):
            os.remove(result_file)


def concatenate_files_streaming(file_list, output_file, buffer_size=COPY_BUFFER_SIZE):
    """Concatenate files in order with fixed memory and atomic replacement."""
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    temp_output = f"{output_file}.tmp.{os.getpid()}"
    try:
        with open(temp_output, 'wb') as output_handle:
            for index, file_path in enumerate(file_list, 1):
                with open(file_path, 'rb') as input_handle:
                    shutil.copyfileobj(input_handle, output_handle, length=buffer_size)
                if index % 1000 == 0:
                    print(f"  - Stream-merged {index}/{len(file_list)} files", file=sys.stderr)
        os.replace(temp_output, output_file)
    except Exception:
        if os.path.exists(temp_output):
            os.remove(temp_output)
        raise


def _merge_indexed_file_group(input_files, output_file, strip_index=False):
    """Merge sorted ``input_index<TAB>payload`` files with bounded memory."""
    handles = []
    heap = []
    try:
        for source_id, file_path in enumerate(input_files):
            handle = open(file_path, "r", encoding="utf-8")
            handles.append(handle)
            line = handle.readline()
            if line:
                index_text, separator, _payload = line.partition("\t")
                if not separator:
                    raise ValueError(f"Malformed indexed result line in {file_path}")
                heapq.heappush(heap, (int(index_text), source_id, line))

        with open(output_file, "w", encoding="utf-8") as output_handle:
            while heap:
                _input_index, source_id, line = heapq.heappop(heap)
                if strip_index:
                    output_handle.write(line.partition("\t")[2])
                else:
                    output_handle.write(line)

                next_line = handles[source_id].readline()
                if next_line:
                    index_text, separator, _payload = next_line.partition("\t")
                    if not separator:
                        raise ValueError(
                            f"Malformed indexed result line in {input_files[source_id]}"
                        )
                    heapq.heappush(
                        heap,
                        (int(index_text), source_id, next_line),
                    )
    finally:
        for handle in handles:
            handle.close()


def merge_indexed_results_streaming(
    file_list,
    output_file,
    max_open_files=DEFAULT_MAX_OPEN_RAW_FILES,
):
    """Merge batch results by original BED order using bounded fan-in.

    One result file is produced per mpileup batch instead of per BED interval.
    Multi-pass merging keeps the number of simultaneously open files bounded.
    """
    if max_open_files < 2:
        raise ValueError("max_open_files must be at least 2")

    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    output_staging_dir = output_dir or os.getcwd()
    merge_work_dir = (
        os.path.dirname(file_list[0]) if file_list else output_staging_dir
    ) or os.getcwd()
    created_intermediates = []
    staged_output = None

    try:
        current_files = list(file_list)
        pass_number = 0
        while len(current_files) > max_open_files:
            next_files = []
            for group_number, group_start in enumerate(
                range(0, len(current_files), max_open_files)
            ):
                group = current_files[group_start:group_start + max_open_files]
                intermediate = os.path.join(
                    merge_work_dir,
                    f".sas-index-merge-{os.getpid()}-{pass_number}-{group_number}.tmp",
                )
                _merge_indexed_file_group(group, intermediate)
                created_intermediates.append(intermediate)
                next_files.append(intermediate)
            current_files = next_files
            pass_number += 1

        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=".sas-result-",
            suffix=".tmp",
            dir=output_staging_dir,
            delete=False,
        ) as staged_handle:
            staged_output = staged_handle.name

        _merge_indexed_file_group(
            current_files,
            staged_output,
            strip_index=True,
        )
        os.replace(staged_output, output_file)
        staged_output = None
    finally:
        if staged_output and os.path.exists(staged_output):
            os.remove(staged_output)
        for intermediate in created_intermediates:
            if os.path.exists(intermediate):
                os.remove(intermediate)


def merge_files_parallel(file_list, output_file, num_threads=8, intermediate_batch_size=1000):
    """Backward-compatible wrapper for the fixed-memory streaming merger."""
    del num_threads, intermediate_batch_size
    concatenate_files_streaming(file_list, output_file)


def process_bed_regions(
    bed_file,
    bam_file,
    ref_fasta,
    output_file,
    temp_dir_path=None,
    processes=64,
    batch_gap=DEFAULT_BATCH_GAP,
    max_batch_span=DEFAULT_MAX_BATCH_SPAN,
    max_regions_per_batch=DEFAULT_MAX_REGIONS_PER_BATCH,
    verbose=False,
):
    """Batch nearby BED regions for mpileup while preserving original outputs."""
    if processes <= 0:
        raise ValueError("processes must be > 0")

    temp_dir = create_run_temp_dir(temp_dir_path)
    print(f"Created temporary run directory: {temp_dir}", file=sys.stderr)

    try:
        regions = []
        with open(bed_file, 'r') as bed_handle:
            for line_number, line in enumerate(bed_handle, 1):
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    if line.strip() and not line.startswith('#'):
                        print(f"Warning: Skipping malformed BED line {line_number}: {line.strip()}", file=sys.stderr)
                    continue
                try:
                    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                    input_index = len(regions)
                    # The fifth tuple field is retained for compatibility with
                    # helper callers. Results are now spooled per batch rather
                    # than allocating one path and one file per interval.
                    regions.append((chrom, start, end, input_index, None))
                except ValueError:
                    print(f"Warning: Skipping malformed BED line {line_number}: {line.strip()}", file=sys.stderr)

        batches = build_region_batches(
            regions,
            max_gap=batch_gap,
            max_span=max_batch_span,
            max_regions=max_regions_per_batch,
        )
        print(
            f"Grouped {len(regions)} original regions into {len(batches)} mpileup batches "
            f"(gap <= {batch_gap:,} bp, span <= {max_batch_span:,} bp, "
            f"regions/batch <= {max_regions_per_batch})",
            file=sys.stderr,
        )

        all_successful_files = []
        if batches:
            worker_count = min(processes, len(batches))
            batch_args = [
                (batch, bam_file, ref_fasta, batch_id, verbose, temp_dir)
                for batch_id, batch in enumerate(batches)
            ]
            with Pool(processes=worker_count) as pool:
                for result_files in pool.imap_unordered(process_region_batch, batch_args):
                    all_successful_files.extend(result_files)

        all_successful_files.sort()
        print(
            f"All batches processed. Order-merging {len(all_successful_files)} batch result files "
            f"into {output_file}...",
            file=sys.stderr,
        )
        merge_indexed_results_streaming(all_successful_files, output_file)
    finally:
        print(f"Cleaning up temporary directory: {temp_dir}", file=sys.stderr)
        shutil.rmtree(temp_dir, ignore_errors=True)

def main():
    """Main function to parse arguments and run the processing."""
    parser = argparse.ArgumentParser(
        description=(
            "Analyze long-read pileup evidence for BED regions using bounded, "
            "parallel samtools mpileup batches."
        ),
        formatter_class=CLIHelpFormatter,
        epilog="""
Examples:
  python3 %(prog)s target.bed input.bam ref.fasta output.txt
  python3 %(prog)s target.bed input.bam ref.fasta output.txt -p 16 -t /scratch
  python3 %(prog)s target.bed input.bam ref.fasta output.txt -g 10000 -S 1000000
""",
    )

    inputs = parser.add_argument_group("required inputs")
    inputs.add_argument("bed_file", metavar="BED", help="Candidate-region BED file")
    inputs.add_argument("bam_file", metavar="BAM", help="Input long-read BAM file")
    inputs.add_argument(
        "ref_fasta", metavar="FASTA", help="Reference genome FASTA file"
    )
    inputs.add_argument("output_file", metavar="OUTPUT", help="Output result file")

    runtime = parser.add_argument_group("runtime options")
    runtime.add_argument(
        "-t", "--temp-dir",
        default=None,
        metavar="DIR",
        help=(
            "Parent directory for a unique temporary run directory. "
            "The parent and its existing contents are preserved. "
            "The platform temporary directory is used when omitted"
        )
    )
    runtime.add_argument(
        "-p", "--processes",
        type=positive_int,
        default=64,
        metavar="N",
        help="Number of parallel worker processes (default: 64)"
    )

    batching = parser.add_argument_group("mpileup batching options")
    batching.add_argument(
        "-g", "--batch-gap",
        type=nonnegative_int,
        default=DEFAULT_BATCH_GAP,
        metavar="BP",
        help=f"Maximum gap between regions placed in one batch (default: {DEFAULT_BATCH_GAP} bp)"
    )
    batching.add_argument(
        "-S", "--max-batch-span",
        type=positive_int,
        default=DEFAULT_MAX_BATCH_SPAN,
        metavar="BP",
        help=f"Maximum genomic span covered by one batch (default: {DEFAULT_MAX_BATCH_SPAN} bp)"
    )
    batching.add_argument(
        "-n", "--max-regions-per-batch",
        type=positive_int,
        default=DEFAULT_MAX_REGIONS_PER_BATCH,
        metavar="N",
        help=f"Maximum original BED regions in one batch (default: {DEFAULT_MAX_REGIONS_PER_BATCH})"
    )
    runtime.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print per-region pileup-cleaning diagnostics"
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {CLI_VERSION}"
    )
    
    args = parser.parse_args()
    
    print(f"Starting analysis with a batch-based parallel model ({args.processes} workers)...")
    process_bed_regions(
        args.bed_file,
        args.bam_file,
        args.ref_fasta,
        args.output_file,
        args.temp_dir,
        args.processes,
        args.batch_gap,
        args.max_batch_span,
        args.max_regions_per_batch,
        args.verbose,
    )
    print(f"Analysis complete. Results saved to {args.output_file}")

if __name__ == '__main__':
    main()
