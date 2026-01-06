#!/usr/bin/env python3
"""
Data Filtering Script
Filter data rows based on BED positions and multiple conditions
EXCLUDES rows that meet certain criteria, keeps the rest
"""

import pandas as pd
import argparse
import sys

def load_fai(fai_file):
    """
    Load chromosome length information from fai file
    fai format: chr_name, length, offset, linebases, linewidth
    """
    chr_lengths = {}
    try:
        with open(fai_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chr_name = parts[0]
                    chr_length = int(parts[1])
                    chr_lengths[chr_name] = chr_length
    except FileNotFoundError:
        print(f"Error: Cannot find fai file {fai_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: Failed to read fai file {e}")
        sys.exit(1)
    
    return chr_lengths

def is_in_terminal_region(chr_name, start, end, chr_lengths, terminal_size=10000):
    """
    Check if BED region is within 10kb of chromosome terminals
    """
    if chr_name not in chr_lengths:
        return False
    
    chr_length = chr_lengths[chr_name]
    
    # Check if within 10kb of chromosome start
    if start <= terminal_size:
        return True
    
    # Check if within 10kb of chromosome end
    if end >= (chr_length - terminal_size):
        return True
    
    return False

def check_complex_conditions(row):
    """
    Check complex conditions (original logic)
    """
    # Main condition A: ((col8>=0.9 OR col9>=0.9) AND col10>=0.6) OR col10>=0.9
    cond_a = (((row[7] >= 0.9 or row[8] >= 0.9) and row[9] >= 0.6) or 
              row[9] >= 0.9)
    
    if cond_a:
        # Condition B: (col4>0 AND col4/col5>=0.1) OR (col6>0 AND col6/col7>=0.1)
        return ((row[3] > 0 and row[4] != 0 and row[3]/row[4] >= 0.1) or 
                (row[5] > 0 and row[6] != 0 and row[5]/row[6] >= 0.1))
    else:
        # Condition C: col4>4 OR col6>3
        return row[3] > 4 or row[5] > 3

def check_simple_conditions(row):
    """
    Check simplified conditions (for terminal 10kb regions)
    """
    # col4>0 AND col4/col5>=0.1 OR col6>0 AND col6/col7>=0.1
    return ((row[3] > 0 and row[4] != 0 and row[3]/row[4] >= 0.1) or 
            (row[5] > 0 and row[6] != 0 and row[5]/row[6] >= 0.1))

def filter_data(data_file, fai_file, output_file=None, exclude_mode=False, sep='\t'):
    """
    Main data filtering function
    """
    HIGH_DEPTH_THRESHOLD = 1000.0

    # Load chromosome length information
    chr_lengths = load_fai(fai_file)
    print(f"Loaded length information for {len(chr_lengths)} chromosomes")
    
    # Read data file
    try:
        df = pd.read_csv(data_file, sep=sep, header=None, skiprows=1)
        print(f"Read {len(df)} rows of data")
    except Exception as e:
        print(f"Error: Failed to read data file {e}")
        sys.exit(1)
    
    if len(df.columns) < 10:
        print(f"Error: Insufficient columns in data file, expected at least 10, got {len(df.columns)}")
        sys.exit(1)
    
    # Filter data
    filtered_indices = []
    high_depth_skipped = 0
    
    for idx, row in df.iterrows():
        try:
            chr_name = str(row[0])
            start = int(row[1])
            end = int(row[2])

            # Skip rows where any platform reports excessively high depth
            try:
                platform1_depth = float(row[4])
                platform2_depth = float(row[6])
            except (ValueError, TypeError):
                platform1_depth = float('nan')
                platform2_depth = float('nan')

            if (not pd.isna(platform1_depth) and platform1_depth > HIGH_DEPTH_THRESHOLD) or \
               (not pd.isna(platform2_depth) and platform2_depth > HIGH_DEPTH_THRESHOLD):
                high_depth_skipped += 1
                continue
            
            # Check if in terminal 10kb regions
            if is_in_terminal_region(chr_name, start, end, chr_lengths):
                # For terminal regions: Keep rows that do NOT meet simple conditions
                should_keep = not check_simple_conditions(row)
            else:
                # For non-terminal regions: Keep rows that do NOT meet complex conditions
                should_keep = not check_complex_conditions(row)
            
            # Include row based on exclude_mode
            if exclude_mode:
                if not should_keep:  # Rows that should be excluded (满足条件的)
                    filtered_indices.append(idx)
            else:
                if should_keep:  # Rows that should be kept (不满足条件的)
                    filtered_indices.append(idx)
                    
        except (ValueError, IndexError) as e:
            print(f"Warning: Error processing row {idx+2} (after skipping header): {e}")
            continue
    
    # Get filtered data
    filtered_df = df.iloc[filtered_indices]
    
    print(f"Filtered result contains {len(filtered_df)} rows")
    if high_depth_skipped > 0:
        print(f"Skipped {high_depth_skipped} rows due to depth > {HIGH_DEPTH_THRESHOLD}")
    
    # Output results
    if output_file:
        try:
            # Create output directory if it doesn't exist
            import os
            output_dir = os.path.dirname(output_file)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)
            
            filtered_df.to_csv(output_file, sep=sep, header=False, index=False)
            print(f"Results saved to {output_file}")
        except Exception as e:
            print(f"Error: Failed to save results {e}")
            sys.exit(1)
    else:
        # Output to stdout
        filtered_df.to_csv(sys.stdout, sep=sep, header=False, index=False)

def main():
    parser = argparse.ArgumentParser(description='Filter data based on BED positions and complex conditions')
    parser.add_argument('data_file', help='Input data file')
    parser.add_argument('fai_file', help='FASTA index file (.fai)')
    parser.add_argument('-o', '--output', help='Output file (default: stdout)')
    parser.add_argument('-e', '--exclude', action='store_true', 
                       help='Output rows that do NOT meet conditions (default: output rows that meet conditions)')
    parser.add_argument('-s', '--sep', default='\t', 
                       help='Field separator (default: tab)')
    
    args = parser.parse_args()
    
    filter_data(args.data_file, args.fai_file, args.output, args.exclude, args.sep)

if __name__ == '__main__':
    main() 