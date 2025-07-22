#!/usr/bin/env python3
"""
Analyze spanning reads in BAM files
Generate sliding windows from fasta.fai file and check for spanning reads in each window
Handle terminal chromosome regions, merge regions, and output BED files
"""

import argparse
import pysam
import sys
import os
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing


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
    windows = []
    for chrom, length in chromosomes.items():
        # Generate normal sliding windows
        for start in range(0, length - window_size + 1, step_size):
            end = start + window_size
            windows.append((chrom, start, end))

        # Ensure the last window covers the chromosome end
        min_last_window_length = 28000
        max_last_window_length = 30000
        # Find the largest step_size-multiple start <= length-min_last_window_length
        last_possible_start = ((length - min_last_window_length) // step_size) * step_size
        # Find the closest start to the end that gives length in [28k, 30k]
        last_window = None
        for potential_start in range(last_possible_start, length, step_size):
            window_length = length - potential_start
            if min_last_window_length <= window_length <= max_last_window_length:
                last_window = (chrom, potential_start, length)
        # Add the last window if needed and not duplicate
        if last_window:
            if not windows or windows[-1][2] < length or windows[-1][1] != last_window[1]:
                windows.append(last_window)
        elif not windows or windows[-1][2] < length:
            # fallback: just add from last_possible_start
            fallback_start = last_possible_start
            if fallback_start < 0:
                fallback_start = 0
            windows.append((chrom, fallback_start, length))
    return windows


def has_spanning_reads(bam_file, chrom, start, end):
    """Check if the specified window has spanning reads"""
    try:
        # Spanning read definition: read spans the entire window region
        for read in bam_file.fetch(chrom, start, end):
            # Skip unmapped reads
            if read.is_unmapped:
                continue
                
            read_start = read.reference_start
            read_end = read.reference_end
            
            # Check if it's a spanning read: read start <= window start and read end >= window end
            if read_start is not None and read_end is not None:
                if read_start <= start and read_end >= end:
                    return True  # Found one spanning read is sufficient
        
        return False
        
    except Exception as e:
        print(f"Error checking spanning reads for {chrom}:{start}-{end}: {e}", file=sys.stderr)
        return False


def process_window_batch(args):
    """Process a batch of windows - for multiprocessing"""
    bam_file_path, window_batch, batch_id = args
    
    # Open BAM file in subprocess
    try:
        bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    except Exception as e:
        print(f"Error opening BAM file in process {batch_id}: {e}", file=sys.stderr)
        return []
    
    no_spanning_windows = []
    
    for chrom, start, end in window_batch:
        if not has_spanning_reads(bam_file, chrom, start, end):
            no_spanning_windows.append((chrom, start, end))
    
    bam_file.close()
    return no_spanning_windows


def split_windows_into_batches(windows, num_processes):
    """Split windows into batches, avoiding over-segmentation"""
    total_windows = len(windows)
    min_batch_size = 50  # Minimum batch size
    
    # If too few windows, reduce process count
    if total_windows < min_batch_size:
        # Too few windows, use single batch
        return [windows]
    
    # Calculate reasonable batch size
    # Ensure each batch has at least min_batch_size windows
    max_effective_processes = max(1, total_windows // min_batch_size)
    effective_processes = min(num_processes, max_effective_processes)
    
    batch_size = max(min_batch_size, total_windows // effective_processes)
    
    batches = []
    for i in range(0, total_windows, batch_size):
        batch = windows[i:i + batch_size]
        batches.append(batch)
    
    return batches


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
    parser = argparse.ArgumentParser(description='Analyze spanning reads in BAM files')
    parser.add_argument('bam_file', help='BAM file path')
    parser.add_argument('fai_file', help='FASTA index file path (.fai)')
    parser.add_argument('-w', '--window-size', type=int, default=30000, 
                       help='Window size (default: 30000)')
    parser.add_argument('-s', '--step-size', type=int, default=2000, 
                       help='Step size (default: 2000)')
    parser.add_argument('-o', '--output', required=True, help='Output prefix')
    parser.add_argument('-p', '--processes', type=int, default=multiprocessing.cpu_count(),
                       help=f'Number of parallel processes (default: {multiprocessing.cpu_count()})')
    args = parser.parse_args()
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
    print(f"Splitting windows into batches for {args.processes} processes...", file=sys.stderr)
    window_batches = split_windows_into_batches(windows, args.processes)
    actual_processes = len(window_batches)
    print(f"Created {actual_processes} batches (optimized from {args.processes} processes)", file=sys.stderr)
    if actual_processes < args.processes:
        avg_windows_per_batch = len(windows) // actual_processes if actual_processes > 0 else 0
        print(f"Optimized: Using {actual_processes} processes with ~{avg_windows_per_batch} windows per batch", file=sys.stderr)
    total_windows = len(windows)
    no_spanning_windows = []
    print(f"Analyzing {total_windows} windows with {actual_processes} processes...", file=sys.stderr)
    batch_args = [(args.bam_file, batch, i) for i, batch in enumerate(window_batches)]
    with ProcessPoolExecutor(max_workers=actual_processes) as executor:
        future_to_batch = {executor.submit(process_window_batch, batch_arg): i 
                          for i, batch_arg in enumerate(batch_args)}
        completed_batches = 0
        for future in as_completed(future_to_batch):
            batch_id = future_to_batch[future]
            try:
                batch_result = future.result()
                no_spanning_windows.extend(batch_result)
                completed_batches += 1
                print(f"Completed batch {completed_batches}/{len(window_batches)}", file=sys.stderr)
            except Exception as e:
                print(f"Batch {batch_id} generated an exception: {e}", file=sys.stderr)
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
    # Check if temp BED file exists and is non-empty
    if not os.path.exists(temp_bed_file) or os.path.getsize(temp_bed_file) == 0:
        print(f"Error: {temp_bed_file} does not exist or is empty! No non-spanning windows found.", file=sys.stderr)
        sys.exit(1)
    print(f"Temp BED file exists and is non-empty.", file=sys.stderr)
    # Sort the temp BED file before merging
    sorted_bed_file = temp_bed_file + ".sorted"
    print(f"About to sort: {os.path.abspath(temp_bed_file)} -> {os.path.abspath(sorted_bed_file)}", file=sys.stderr)
    try:
        subprocess.run(f"sort -k1,1 -k2,2n {temp_bed_file} > {sorted_bed_file}", shell=True, check=True)
        print(f"Sorted BED file created: {sorted_bed_file}", file=sys.stderr)
        print(f"Sorted BED file absolute path: {os.path.abspath(sorted_bed_file)}", file=sys.stderr)
    except Exception as e:
        print(f"Error sorting BED file: {e}", file=sys.stderr)
        os.remove(temp_bed_file)
        sys.exit(1)
    # Use bedtools merge to merge all non-spanning windows
    merged_bed_file = f"{args.output}.non_spanning_regions.merged.bed"
    print(f"About to merge: {os.path.abspath(sorted_bed_file)} -> {os.path.abspath(merged_bed_file)}", file=sys.stderr)
    cmd = ["bedtools", "merge", "-i", sorted_bed_file]
    try:
        with open(merged_bed_file, 'w') as outbed:
            subprocess.run(cmd, stdout=outbed, check=True)
        print(f"bedtools merge completed: {merged_bed_file}", file=sys.stderr)
        print(f"Merged BED file absolute path: {os.path.abspath(merged_bed_file)}", file=sys.stderr)
    except Exception as e:
        print(f"Error running bedtools merge: {e}", file=sys.stderr)
        os.remove(temp_bed_file)
        os.remove(sorted_bed_file)
        sys.exit(1)
    # Read merged regions for terminal correction
    merged_regions = []
    with open(merged_bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                merged_regions.append((chrom, start, end))
    os.remove(temp_bed_file)
    os.remove(sorted_bed_file)
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
