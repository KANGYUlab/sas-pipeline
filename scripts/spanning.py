#!/usr/bin/env python3
"""
Analyze spanning reads in BAM files
Generate sliding windows from fasta.fai file and check for spanning reads in each window
Handle terminal chromosome regions, merge regions, and output BED files
"""

import argparse
import sys
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

try:
    import pysam
except ModuleNotFoundError:
    pysam = None

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


def read_fasta_index(fai_file):
    """Read fasta.fai file to get chromosome length information"""
    chromosomes = {}
    try:
        with open(fai_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    chrom = fields[0]
                    length = int(fields[1])
                    chromosomes[chrom] = length
    except Exception as e:
        print(f"Error reading fasta index file {fai_file}: {e}", file=sys.stderr)
        sys.exit(1)
    return chromosomes


def generate_sliding_windows(chromosomes, window_size=30000, step_size=2000):
    """Generate sliding windows, ensuring the last window covers the chromosome end, with start divisible by step_size (2k) and length 28k-30k (can be <30k), and start is as close as possible to the end."""
    if window_size <= 0 or step_size <= 0:
        raise ValueError("window_size and step_size must be positive")

    windows = []
    for chrom, length in chromosomes.items():
        if length <= 0:
            continue

        chrom_windows = []

        # A short contig still needs one valid window.  Keeping this case local
        # to the chromosome also prevents the previous chromosome's last
        # window from suppressing it.
        if length < 28000:
            windows.append((chrom, 0, length))
            continue

        # Generate normal sliding windows
        for start in range(0, length - window_size + 1, step_size):
            end = start + window_size
            chrom_windows.append((chrom, start, end))

        # Ensure the last window covers the chromosome end
        min_last_window_length = 28000
        max_last_window_length = 30000
        # Find the largest step_size-multiple start <= length-min_last_window_length
        last_possible_start = ((length - min_last_window_length) // step_size) * step_size
        # The largest aligned start gives the terminal window closest to the
        # chromosome end.  With the default 2 kb step its length is 28-30 kb.
        last_window = None
        if last_possible_start >= 0:
            window_length = length - last_possible_start
            if min_last_window_length <= window_length <= max_last_window_length:
                last_window = (chrom, last_possible_start, length)
        # Add the last window if needed and not duplicate
        if last_window and last_window not in chrom_windows:
            chrom_windows.append(last_window)
        elif not last_window and not any(end == length for _, _, end in chrom_windows):
            # fallback: just add from last_possible_start
            fallback_start = max(0, last_possible_start)
            chrom_windows.append((chrom, fallback_start, length))

        windows.extend(chrom_windows)
    return windows


def find_non_spanning_windows(bam_file, chrom, chrom_windows):
    """Find non-spanning windows with one coordinate-ordered BAM scan.

    For a window ``[start, end)``, a spanning read must satisfy exactly the
    original predicate ``read_start <= start and read_end >= end``.  While
    window starts increase, the maximum end among all reads whose starts have
    already been seen is sufficient to evaluate that predicate.
    """
    ordered_windows = sorted(chrom_windows, key=lambda window: (window[1], window[2]))
    no_spanning_windows = []
    read_iterator = iter(bam_file.fetch(chrom))

    def next_mapped_interval():
        for read in read_iterator:
            if read.is_unmapped:
                continue
            read_start = read.reference_start
            read_end = read.reference_end
            if read_start is not None and read_end is not None:
                return read_start, read_end
        return None

    next_read = next_mapped_interval()
    maximum_read_end = -1
    for window_chrom, start, end in ordered_windows:
        if window_chrom != chrom:
            raise ValueError(f"Window chromosome {window_chrom} does not match {chrom}")

        while next_read is not None and next_read[0] <= start:
            maximum_read_end = max(maximum_read_end, next_read[1])
            next_read = next_mapped_interval()

        if maximum_read_end < end:
            no_spanning_windows.append((chrom, start, end))

    return no_spanning_windows


def process_chromosome(args):
    """Open one BAM handle and scan one chromosome in a worker process."""
    bam_file_path, chrom, chrom_windows = args
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    try:
        return find_non_spanning_windows(bam_file, chrom, chrom_windows)
    finally:
        bam_file.close()


def process_terminal_regions(spanning_regions, chromosomes):
    """Process spanning regions at chromosome start/end"""
    processed_regions = []
    
    for chrom, start, end in spanning_regions:
        chrom_length = chromosomes.get(chrom, 0)
        if chrom_length == 0:
            continue
            
        # Check if at chromosome start or end
        is_chromosome_start = (start == 0)
        is_chromosome_end = (end == chrom_length)
        
        # If at chromosome start, remove last 28k
        if is_chromosome_start:
            region_length = end - start
            if region_length > 28000:
                # Remove last 28k, keep the front part
                new_end = end - 28000
                if new_end > start:
                    processed_regions.append((chrom, start, new_end))
            # If region length <= 28k, skip completely
            
        # If at chromosome end, remove first 28k  
        elif is_chromosome_end:
            region_length = end - start
            if region_length > 28000:
                # Remove first 28k, keep the back part
                new_start = start + 28000
                if new_start < end:
                    processed_regions.append((chrom, new_start, end))
            # If region length <= 28k, skip completely
            
        # If neither at start nor end, keep as-is
        else:
            processed_regions.append((chrom, start, end))
    
    return processed_regions


def merge_regions_with_bedtools(regions, output_prefix):
    """Merge regions using bedtools merge and return merged regions as a list of tuples."""
    if not regions:
        return []
    temp_bed = f"{output_prefix}.temp_spanning.bed"
    merged_bed = f"{output_prefix}.spanning_regions.merged.bed"
    with open(temp_bed, 'w') as f:
        for chrom, start, end in regions:
            f.write(f"{chrom}\t{start}\t{end}\n")
    try:
        cmd = ["bedtools", "merge", "-i", temp_bed]
        with open(merged_bed, 'w') as f:
            subprocess.run(cmd, stdout=f, check=True)
        merged_regions = []
        with open(merged_bed, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                    merged_regions.append((chrom, start, end))
        os.remove(temp_bed)
        print(f"Merged {len(regions)} regions into {len(merged_regions)} regions", file=sys.stderr)
        print(f"Final merged and terminal-corrected regions saved to: {merged_bed}", file=sys.stderr)
        return merged_regions
    except subprocess.CalledProcessError as e:
        print(f"Error running bedtools merge: {e}", file=sys.stderr)
        if os.path.exists(temp_bed):
            os.remove(temp_bed)
        return regions
    except Exception as e:
        print(f"Error merging spanning regions: {e}", file=sys.stderr)
        if os.path.exists(temp_bed):
            os.remove(temp_bed)
        return regions


def correct_terminal_regions(regions, chromosomes):
    """Correct regions at chromosome start/end by cutting 28k from the end or start as needed."""
    corrected = []
    for chrom, start, end in regions:
        chrom_length = chromosomes.get(chrom, 0)
        if chrom_length == 0:
            continue
        # If at chromosome start
        if start == 0:
            new_end = end - 28000
            if new_end > start:
                corrected.append((chrom, start, new_end))
        # If at chromosome end
        elif end == chrom_length:
            new_start = start + 28000
            if new_start < end:
                corrected.append((chrom, new_start, end))
        else:
            corrected.append((chrom, start, end))
    return corrected


def assign_corrected_regions_to_2k_windows(corrected_regions, window_size=2000):
    """Given corrected merged regions, divide them into 2kb windows (chrom, start, end)."""
    windows = []
    for chrom, region_start, region_end in corrected_regions:
        start = (region_start // window_size) * window_size
        while start < region_end:
            end = min(start + window_size, region_end)
            if end > start:
                windows.append((chrom, start, end))
            start += window_size
    return windows


def main():
    parser = argparse.ArgumentParser(
        description="Find BAM windows without reads spanning the complete interval.",
        formatter_class=CLIHelpFormatter,
        epilog="""
Examples:
  python3 %(prog)s input.bam reference.fa.fai -o sample
  python3 %(prog)s input.bam reference.fa.fai -o sample -w 30000 -s 2000 -p 8
        """,
    )
    inputs = parser.add_argument_group("required inputs")
    inputs.add_argument("bam_file", metavar="BAM", help="Input BAM file")
    inputs.add_argument("fai_file", metavar="FAI", help="Reference FASTA index")
    inputs.add_argument(
        "-o", "--output", required=True, metavar="PREFIX",
        help="Prefix for generated BED files",
    )

    analysis = parser.add_argument_group("analysis options")
    analysis.add_argument(
        "-w", "--window-size", type=positive_int, default=30000, metavar="BP",
        help="Sliding-window size (default: 30000 bp)",
    )
    analysis.add_argument(
        "-s", "--step-size", type=positive_int, default=2000, metavar="BP",
        help="Distance between consecutive window starts (default: 2000 bp)",
    )
    analysis.add_argument(
        "-p", "--processes", type=positive_int,
        default=multiprocessing.cpu_count(), metavar="N",
        help=f"Number of chromosome worker processes (default: {multiprocessing.cpu_count()})",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {CLI_VERSION}"
    )
    args = parser.parse_args()
    if pysam is None:
        parser.error("missing required Python dependency: pysam")
    print(f"Reading fasta index from {args.fai_file}...", file=sys.stderr)
    chromosomes = read_fasta_index(args.fai_file)
    print(f"Found {len(chromosomes)} chromosomes", file=sys.stderr)
    print(f"Generating sliding windows (window={args.window_size}, step={args.step_size})...", file=sys.stderr)
    windows = generate_sliding_windows(chromosomes, args.window_size, args.step_size)
    print(f"Generated {len(windows)} windows", file=sys.stderr)
    print(f"Verifying BAM file {args.bam_file}...", file=sys.stderr)
    try:
        test_bam = pysam.AlignmentFile(args.bam_file, "rb")
        test_bam.close()
    except Exception as e:
        print(f"Error accessing BAM file {args.bam_file}: {e}", file=sys.stderr)
        sys.exit(1)
    windows_by_chromosome = {}
    for window in windows:
        windows_by_chromosome.setdefault(window[0], []).append(window)
    chromosome_args = [
        (args.bam_file, chrom, chrom_windows)
        for chrom, chrom_windows in windows_by_chromosome.items()
    ]
    actual_processes = min(max(1, args.processes), len(chromosome_args)) if chromosome_args else 0
    print(
        f"Scanning {len(chromosome_args)} chromosomes with {actual_processes} processes...",
        file=sys.stderr,
    )
    total_windows = len(windows)
    no_spanning_windows = []
    print(f"Analyzing {total_windows} windows with {actual_processes} processes...", file=sys.stderr)
    if len(chromosome_args) == 1:
        no_spanning_windows.extend(process_chromosome(chromosome_args[0]))
    elif chromosome_args:
        with ProcessPoolExecutor(max_workers=actual_processes) as executor:
            future_to_chrom = {
                executor.submit(process_chromosome, chromosome_arg): chromosome_arg[1]
                for chromosome_arg in chromosome_args
            }
            completed_chromosomes = 0
            for future in as_completed(future_to_chrom):
                chrom = future_to_chrom[future]
                try:
                    no_spanning_windows.extend(future.result())
                    completed_chromosomes += 1
                    print(
                        f"Completed chromosome {completed_chromosomes}/{len(chromosome_args)}",
                        file=sys.stderr,
                    )
                except Exception as e:
                    print(f"Chromosome {chrom} generated an exception: {e}", file=sys.stderr)
                    raise
    no_spanning_windows.sort(key=lambda x: (x[0], x[1]))
    print("Merging non-spanning windows...", file=sys.stderr)
    # After collecting all no_spanning_windows, write them to a temp BED file
    temp_bed_file = f"{args.output}.non_spanning_windows.temp.bed"
    with open(temp_bed_file, 'w') as f:
        for chrom, start, end in no_spanning_windows:
            if start < 0:
                print(f"Warning: Adjusting negative start to 0: {chrom}\t{start}\t{end}", file=sys.stderr)
                start = 0
            f.write(f"{chrom}\t{start}\t{end}\n")
    print(f"Current working directory: {os.getcwd()}", file=sys.stderr)
    print(f"Temp BED file: {os.path.abspath(temp_bed_file)}", file=sys.stderr)
    # Use bedtools merge to merge all non-spanning windows
    merged_bed_file = f"{args.output}.non_spanning_regions.merged.bed"
    if no_spanning_windows:
        # no_spanning_windows was already sorted in Python, so bedtools can
        # consume it directly without another shell sort process.
        print(f"About to merge: {os.path.abspath(temp_bed_file)} -> {os.path.abspath(merged_bed_file)}", file=sys.stderr)
        cmd = ["bedtools", "merge", "-i", temp_bed_file]
        try:
            with open(merged_bed_file, 'w') as outbed:
                subprocess.run(cmd, stdout=outbed, check=True)
            print(f"bedtools merge completed: {merged_bed_file}", file=sys.stderr)
            print(f"Merged BED file absolute path: {os.path.abspath(merged_bed_file)}", file=sys.stderr)
        except Exception as e:
            print(f"Error running bedtools merge: {e}", file=sys.stderr)
            os.remove(temp_bed_file)
            sys.exit(1)
    else:
        # All windows have spanning reads.  This is a valid result and both
        # public BED outputs must still be created as empty files.
        open(merged_bed_file, 'w').close()
        print("No non-spanning windows found; creating empty output files.", file=sys.stderr)
    # Read merged regions for terminal correction
    merged_regions = []
    with open(merged_bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                merged_regions.append((chrom, start, end))
    os.remove(temp_bed_file)
    if not merged_regions:
        print(f"Warning: No merged non-spanning regions found!", file=sys.stderr)
    # Apply terminal correction
    print("Applying terminal region correction...", file=sys.stderr)
    final_regions = correct_terminal_regions(merged_regions, chromosomes)
    print(f"{len(final_regions)} regions remain after terminal correction", file=sys.stderr)
    # Overwrite merged_bed_file with terminal-corrected regions
    with open(merged_bed_file, 'w') as outfile:
        for chrom, start, end in final_regions:
            if start >= end:
                continue
            outfile.write(f"{chrom}\t{start}\t{end}\n")
    print(f"Results saved to: {merged_bed_file}", file=sys.stderr)
    print(f"Analysis complete!", file=sys.stderr)
    print(f"Total windows analyzed: {total_windows}", file=sys.stderr)
    print(f"Final non-spanning merged regions: {len(final_regions)}", file=sys.stderr)
    # Assign corrected regions to 2kb windows and output
    print("Assigning corrected regions to 2kb windows...", file=sys.stderr)
    windows_2k = assign_corrected_regions_to_2k_windows(final_regions, window_size=2000)
    print(f"{len(windows_2k)} 2kb non-spanning windows generated", file=sys.stderr)
    nonspanning_bed_file = f"{args.output}.non_spanning.bed"
    with open(nonspanning_bed_file, 'w') as outfile:
        for chrom, start, end in windows_2k:
            if start >= end:
                continue
            outfile.write(f"{chrom}\t{start}\t{end}\n")
    print(f"Results saved to: {nonspanning_bed_file}", file=sys.stderr)


if __name__ == "__main__":
    main() 
