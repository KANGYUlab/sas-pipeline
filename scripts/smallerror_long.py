#!/usr/bin/env python3
"""
Processes a BAM file based on regions defined in a BED file to analyze pileup data.

This script takes a BED file specifying genomic regions, a BAM alignment file, and a
reference FASTA file. For each region, it uses 'samtools mpileup' to generate
pileup data, which is then processed to clean the pileup strings (removing indels
and other special characters) and calculate statistics like average depth.

The processing is done in parallel using multiple processes to speed up the analysis.
"""
import re
import sys
import argparse
import subprocess
import signal
import os
import shutil
from multiprocessing import Pool, cpu_count

def process_pileup_from_file(raw_pileup_path, region_start, region_end):
    """
    Reads a raw pileup file from disk, calculates average depth, then processes it.
    This version implements the proper cross-row alignment logic with correct processing order.
    """
    # Pass 1: Calculate average depth.
    total_depth = 0
    line_count = 0
    with open(raw_pileup_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4 and parts[3].isdigit():
                total_depth += int(parts[3])
            line_count += 1
    
    average_depth = (total_depth / line_count) if line_count > 0 else 0.0

    # Pass 2: Load all lines and extract pileup strings for processing
    all_lines = []
    pileup_cols = []
    with open(raw_pileup_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                all_lines.append(parts)
                pileup_cols.append(list(parts[4]))  # Convert to list for easier manipulation
            else:
                all_lines.append(parts)
                pileup_cols.append([])

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
    
    print(f"Debug: Total lines available: {len(all_lines)}, Total pileup columns: {len(pileup_cols)}", file=sys.stderr)
    
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
            print(f"Debug: Processed first line only", file=sys.stderr)
        else:
            print(f"Debug: Skipping first row, len(parts)={len(parts)}, row_idx < len(pileup_cols): {row_idx < len(pileup_cols)}", file=sys.stderr)
    else:
        print(f"Debug: No data to process", file=sys.stderr)
            
    print(f"Debug: Total processed lines: {len(processed_lines)}", file=sys.stderr)
    return processed_lines

def process_region(args):
    """
    Worker function that streams samtools ouz
    then processes the resulting file.
    """
    chrom, start, end, bam_file, ref_fasta, final_temp_file = args
    
    region_str = f"{chrom}:{start+1}-{end}"
    print(f"Processing region: {region_str}", file=sys.stderr)

    # Define paths for intermediate and final temp files
    temp_dir = os.path.dirname(final_temp_file)
    base_name = os.path.basename(final_temp_file)
    raw_pileup_file = os.path.join(temp_dir, f"raw_{base_name}")

    cmd = [
        'samtools', 'mpileup', '-a', '-A', '-B', '-d', '0',
        '-r', region_str, '-f', ref_fasta, '-Q', '0', '--ff', '0', '-x', bam_file
    ]
    
    try:
        # Step 1: Stream samtools output directly to a raw pileup file on disk.
        with open(raw_pileup_file, 'w') as stdout_f, subprocess.Popen(cmd, stdout=stdout_f, stderr=subprocess.PIPE, universal_newlines=True) as proc:
            _, stderr = proc.communicate() # Wait for process to finish
            if proc.returncode != 0:
                print(f"Samtools error for region {region_str}: {stderr}", file=sys.stderr)
                return None

        # Step 2: Process the raw pileup file from disk, passing start and end coordinates.
        formatted_output = process_pileup_from_file(raw_pileup_file, start, end)
        
        # Step 3: Write the final processed output.
        if formatted_output:
            with open(final_temp_file, 'w') as f_out:
                f_out.write('\n'.join(formatted_output) + '\n')
            return final_temp_file
        
    except Exception as e:
        print(f"Unexpected error processing region {region_str}: {e}", file=sys.stderr)
        return None
    finally:
        # Step 4: Clean up the raw pileup file.
        if os.path.exists(raw_pileup_file):
            os.remove(raw_pileup_file)
    
    return None

def process_bed_regions(bed_file, bam_file, ref_fasta, output_file, temp_dir_path=None, processes=64):
    """
    Reads a BED file in chunks and processes them in parallel.
    This version includes zero-padding for temp files to ensure correct final order.
    """
    
    # Use specified temp directory or default
    if temp_dir_path:
        temp_dir = temp_dir_path
    else:
        temp_dir = "temp_pileup_results"
    
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.makedirs(temp_dir)
    print(f"Created temporary directory: {temp_dir}", file=sys.stderr)

    all_successful_files = []
    region_counter = 0

    with open(bed_file, 'r') as bed_f, Pool(processes=processes) as pool:
        print(f"Processing BED file in chunks of {processes} regions using {processes} parallel processes...", file=sys.stderr)
        
        while True:
            chunk_lines = [line for _, line in zip(range(processes), bed_f)]
            if not chunk_lines:
                break

            chunk_args = []
            for line in chunk_lines:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    try:
                        chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                        # Zero-pad filenames to ensure correct alphabetical sorting later
                        temp_output_file = os.path.join(temp_dir, f"region_{region_counter:09d}.txt")
                        chunk_args.append((chrom, start, end, bam_file, ref_fasta, temp_output_file))
                        region_counter += 1
                    except (ValueError, IndexError):
                        print(f"Warning: Skipping malformed BED line: {line.strip()}", file=sys.stderr)

            if chunk_args:
                print(f"  - Submitting a new chunk of {len(chunk_args)} regions for processing.", file=sys.stderr)
                for result_file in pool.imap_unordered(process_region, chunk_args):
                    if result_file:
                        all_successful_files.append(result_file)

    print(f"All chunks processed. Merging {len(all_successful_files)} result files into {output_file}...", file=sys.stderr)
    # Sort the list of successful files before merging to maintain order.
    all_successful_files.sort()
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Use batched system cat command for fastest merge, with fallback to Python
    batch_size = 1000  # Process files in batches to avoid "Argument list too long" error
    
    try:
        # Try using batched system cat commands for fastest merge
        with open(output_file, 'w') as final_out:
            for i in range(0, len(all_successful_files), batch_size):
                batch_files = all_successful_files[i:i+batch_size]
                batch_num = i // batch_size + 1
                total_batches = (len(all_successful_files) + batch_size - 1) // batch_size
                
                print(f"  - Processing batch {batch_num}/{total_batches} ({len(batch_files)} files)...", file=sys.stderr)
                
                cat_cmd = ['cat'] + batch_files
                result = subprocess.run(cat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                
                if result.returncode != 0:
                    raise subprocess.CalledProcessError(result.returncode, cat_cmd, result.stderr)
                
                final_out.write(result.stdout)
                
        print(f"  - Fast merge completed using batched system 'cat' commands.", file=sys.stderr)
        
    except (subprocess.CalledProcessError, FileNotFoundError, OSError) as e:
        # Fallback to Python-based merge with optimized buffering
        print(f"  - Cat command failed ({e}), using Python-based merge with optimized buffering...", file=sys.stderr)
        with open(output_file, 'w', buffering=1024*1024) as final_out:  # 1MB buffer
            for i, temp_file in enumerate(all_successful_files):
                if i % 1000 == 0:  # Progress indicator every 1000 files
                    print(f"    - Merged {i}/{len(all_successful_files)} files...", file=sys.stderr)
                try:
                    with open(temp_file, 'r', buffering=512*1024) as temp_in:  # 512KB buffer
                        # Use shutil.copyfileobj for efficient copying
                        shutil.copyfileobj(temp_in, final_out, 1024*1024)  # 1MB chunks
                except IOError as e:
                    print(f"Warning: Could not read temp file {temp_file}: {e}", file=sys.stderr)

    print(f"Cleaning up temporary directory: {temp_dir}", file=sys.stderr)
    shutil.rmtree(temp_dir)

def main():
    """Main function to parse arguments and run the processing."""
    parser = argparse.ArgumentParser(
        description="Analyze samtools mpileup output for regions in a BED file, processing in parallel chunks.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Example:
  python3 %(prog)s target.bed input.bam ref.fasta output.txt
  python3 %(prog)s target.bed input.bam ref.fasta output.txt --temp-dir /path/to/temp
  python3 %(prog)s target.bed input.bam ref.fasta output.txt --processes 16
"""
    )
    parser.add_argument("bed_file", help="Path to the input BED file defining regions of interest.")
    parser.add_argument("bam_file", help="Path to the input BAM file.")
    parser.add_argument("ref_fasta", help="Path to the reference genome FASTA file.")
    parser.add_argument("output_file", help="Path to the output file to store results.")
    parser.add_argument(
        "--temp-dir", "-t",
        type=str,
        default=None,
        help="Path to the temporary directory for intermediate files. If not specified, uses 'temp_pileup_results' in current directory."
    )
    parser.add_argument(
        "--processes", "-p",
        type=int,
        default=64,
        help="Number of parallel processes to use for processing regions (default: 64)."
    )
    
    args = parser.parse_args()
    
    print(f"Starting analysis with a chunk-based parallel model ({args.processes} regions at a time)...")
    process_bed_regions(
        args.bed_file,
        args.bam_file,
        args.ref_fasta,
        args.output_file,
        args.temp_dir,
        args.processes
    )
    print(f"Analysis complete. Results saved to {args.output_file}")

if __name__ == '__main__':
    main()



