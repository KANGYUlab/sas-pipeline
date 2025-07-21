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
    """Merge regions using bedtools merge"""
    if not regions:
        return []
        
    # Create temporary BED file
    temp_bed = f"{output_prefix}.temp_spanning.bed"
    merged_bed = f"{output_prefix}.spanning_regions.merged.bed"
    
    # Write BED format
    with open(temp_bed, 'w') as f:
        for chrom, start, end in regions:
            f.write(f"{chrom}\t{start}\t{end}\n")
    
    # Use bedtools merge to merge regions
    try:
        cmd = ["bedtools", "merge", "-i", temp_bed]
        with open(merged_bed, 'w') as f:
            subprocess.run(cmd, stdout=f, check=True)
        
        # Read merged results
        merged_regions = []
        with open(merged_bed, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                    merged_regions.append((chrom, start, end))
        
        # Clean up temporary file
        os.remove(temp_bed)
        
        print(f"Merged {len(regions)} regions into {len(merged_regions)} regions", file=sys.stderr)
        print(f"Merged spanning regions saved to: {merged_bed}", file=sys.stderr)
        
        return merged_regions
        
    except subprocess.CalledProcessError as e:
        print(f"Error running bedtools merge: {e}", file=sys.stderr)
        # Clean up temporary file
        if os.path.exists(temp_bed):
            os.remove(temp_bed)
        return regions
    except Exception as e:
        print(f"Error merging spanning regions: {e}", file=sys.stderr)
        # Clean up temporary file
        if os.path.exists(temp_bed):
            os.remove(temp_bed)
        return regions


def main():
    parser = argparse.ArgumentParser(description='Analyze spanning reads in BAM files')
    parser.add_argument('bam_file', help='BAM file path')
    parser.add_argument('fai_file', help='FASTA index file path (.fai)')
    parser.add_argument('-w', '--window-size', type=int, default=30000, 
                       help='Window size (default: 30000)')
    parser.add_argument('-s', '--step-size', type=int, default=2000, 
                       help='Step size (default: 2000)')
    parser.add_argument('-o', '--output', required=True, help='Output prefix')
    parser.add_argument('--bed', action='store_true', 
                       help='Output BED format (default: table format with header)')
    parser.add_argument('-p', '--processes', type=int, default=multiprocessing.cpu_count(),
                       help=f'Number of parallel processes (default: {multiprocessing.cpu_count()})')
    parser.add_argument('--merge', action='store_true',
                       help='Merge overlapping regions and process chromosome terminals')
    
    args = parser.parse_args()
    
    # Read fasta index file
    print(f"Reading fasta index from {args.fai_file}...", file=sys.stderr)
    chromosomes = read_fasta_index(args.fai_file)
    print(f"Found {len(chromosomes)} chromosomes", file=sys.stderr)
    
    # Generate sliding windows
    print(f"Generating sliding windows (window={args.window_size}, step={args.step_size})...", file=sys.stderr)
    windows = generate_sliding_windows(chromosomes, args.window_size, args.step_size)
    print(f"Generated {len(windows)} windows", file=sys.stderr)
    
    # Verify BAM file
    print(f"Verifying BAM file {args.bam_file}...", file=sys.stderr)
    try:
        test_bam = pysam.AlignmentFile(args.bam_file, "rb")
        test_bam.close()
    except Exception as e:
        print(f"Error accessing BAM file {args.bam_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Split windows into batches for multiprocessing
    print(f"Splitting windows into batches for {args.processes} processes...", file=sys.stderr)
    window_batches = split_windows_into_batches(windows, args.processes)
    actual_processes = len(window_batches)
    print(f"Created {actual_processes} batches (optimized from {args.processes} processes)", file=sys.stderr)
    
    if actual_processes < args.processes:
        avg_windows_per_batch = len(windows) // actual_processes if actual_processes > 0 else 0
        print(f"Optimized: Using {actual_processes} processes with ~{avg_windows_per_batch} windows per batch", file=sys.stderr)
    
    # Multiprocess analysis of windows
    total_windows = len(windows)
    no_spanning_windows = []
    
    print(f"Analyzing {total_windows} windows with {actual_processes} processes...", file=sys.stderr)
    
    # Prepare arguments
    batch_args = [(args.bam_file, batch, i) for i, batch in enumerate(window_batches)]
    
    # Use multiprocessing
    with ProcessPoolExecutor(max_workers=actual_processes) as executor:
        # Submit all tasks
        future_to_batch = {executor.submit(process_window_batch, batch_arg): i 
                          for i, batch_arg in enumerate(batch_args)}
        
        # Collect results
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
    
    # Sort output by chromosome and position
    no_spanning_windows.sort(key=lambda x: (x[0], x[1]))
    
    # Process terminal regions and merge (if requested)
    if args.merge:
        print("Processing terminal regions...", file=sys.stderr)
        processed_regions = process_terminal_regions(no_spanning_windows, chromosomes)
        print(f"Processed {len(processed_regions)} regions after terminal filtering", file=sys.stderr)
        
        print("Merging overlapping regions...", file=sys.stderr)
        final_regions = merge_regions_with_bedtools(processed_regions, args.output)
        no_spanning_windows = final_regions
    
    # Output results
    output_file = f"{args.output}.non_spanning.bed"
    with open(output_file, 'w') as outfile:
        if not args.bed:
            outfile.write("chromosome\tstart\tend\twindow_size\n")
        
        for chrom, start, end in no_spanning_windows:
            if args.bed:
                outfile.write(f"{chrom}\t{start}\t{end}\n")
            else:
                outfile.write(f"{chrom}\t{start}\t{end}\t{end - start}\n")
    
    print(f"Results saved to: {output_file}", file=sys.stderr)
    
    # Output statistics
    no_spanning_count = len(no_spanning_windows)
    print(f"Analysis complete!", file=sys.stderr)
    print(f"Total windows analyzed: {total_windows}", file=sys.stderr)
    print(f"Windows without spanning reads: {no_spanning_count}", file=sys.stderr)
    print(f"Percentage without spanning reads: {no_spanning_count/total_windows*100:.2f}%", file=sys.stderr)


if __name__ == "__main__":
    main() 