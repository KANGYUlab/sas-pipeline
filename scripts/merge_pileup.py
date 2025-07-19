#!/usr/bin/env python3
"""
BED File Merger and Base Composition Calculator

Merge files based on bed coordinates and calculate base composition ratios
from reference genome sequences.

"""

import sys
import argparse
from collections import defaultdict, Counter
from pathlib import Path


class BedFileMerger:
    """Merge files based on bed coordinates and calculate base composition."""
    
    def __init__(self):
        self.reference_data = {}
        self.file1_data = {}
        self.file2_data = {}
        self.ref_sequences = {}
        
    def merge_files(self, bed_file: str, file1: str, file2: str, ref_file: str, output_file: str) -> None:
        """
        Merge files based on bed coordinates and calculate base composition.
        
        Args:
            bed_file (str): BED format file with coordinates (chr, start, end)
            file1 (str): First file to merge 
            file2 (str): Second file to merge
            ref_file (str): Reference genome FASTA file
            output_file (str): Output file path
            
        Output columns:
            1-3: BED coordinates (chr, start, end)
            4-5: File1 last two columns
            6-7: File2 last two columns
            8-10: AG_ratio, CT_ratio, GC_ratio from reference sequence
        """
        print(f"Reading BED file: {bed_file}")
        self._read_bed_file(bed_file)
        
        print(f"Reading first file: {file1}")
        self._read_data_file(file1, self.file1_data, skip_header=False)
        
        print(f"Reading second file: {file2}")
        self._read_data_file(file2, self.file2_data, skip_header=False)
        
        print(f"Loading reference sequences: {ref_file}")
        self._load_reference_sequences(ref_file)
        
        print(f"Merging files and calculating base composition: {output_file}")
        self._write_merged_file(output_file)
        
    def _read_bed_file(self, filename: str) -> None:
        """Read the BED format file."""
        with open(filename, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                    
                columns = line.split('\t')
                if len(columns) < 3:
                    print(f"Warning: BED file line {line_num} has insufficient columns: {len(columns)}")
                    continue
                
                # Key: chr, start, end (first three columns)
                try:
                    chr_name = columns[0]
                    start = int(columns[1])
                    end = int(columns[2])
                    key = (chr_name, start, end)
                    self.reference_data[key] = True
                except ValueError:
                    print(f"Warning: BED file line {line_num} has invalid coordinates")
                    continue
                    
        print(f"BED file loaded: {len(self.reference_data)} regions")
        
    def _read_data_file(self, filename: str, data_dict: dict, skip_header: bool = True) -> None:
        """Read a data file and store relevant columns."""
        processed_lines = 0
        with open(filename, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                # Skip header line if requested
                if skip_header and line_num == 1:
                    print(f"Skipping header line in {filename}")
                    continue
                    
                line = line.strip()
                if not line:
                    print(f"Skipping empty line {line_num} in {filename}")
                    continue
                    
                columns = line.split('\t')
                if len(columns) < 5:
                    print(f"Warning: File {filename} line {line_num} has insufficient columns: {len(columns)}")
                    continue
                
                processed_lines += 1
                if processed_lines <= 3:  # Debug first 3 processed lines
                    print(f"Processing line {line_num} in {filename}: columns={len(columns)}, key=({columns[-5]}, {columns[-4]}, {columns[-3]})")
                
                # Key: last 5th, 4th, 3rd columns (倒数5、4、3列)
                try:
                    chr_name = columns[-5]
                    start = int(columns[-4])
                    end = int(columns[-3])
                    key = (chr_name, start, end)
                    
                    # Store last two columns (倒数第二和倒数第一列)
                    last_two = (columns[-2], columns[-1])
                    data_dict[key] = last_two
                    
                except (ValueError, IndexError):
                    print(f"Warning: File {filename} line {line_num} has invalid format")
                    continue
                    
        print(f"Data file {filename} loaded: {len(data_dict)} records from {processed_lines} processed lines")
        
    def _load_reference_sequences(self, ref_file: str) -> None:
        """Load reference sequences from FASTA file."""
        current_chr = None
        current_seq = []
        
        print("Loading reference genome sequences...")
        
        with open(ref_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>'):
                    # Save previous sequence
                    if current_chr and current_seq:
                        self.ref_sequences[current_chr] = ''.join(current_seq).upper()
                    
                    # Start new sequence
                    current_chr = line[1:].split()[0]  # Get chromosome name
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Save last sequence
            if current_chr and current_seq:
                self.ref_sequences[current_chr] = ''.join(current_seq).upper()
                
        print(f"Reference sequences loaded: {len(self.ref_sequences)} chromosomes")
        
    def _calculate_base_ratios(self, chr_name: str, start: int, end: int) -> tuple:
        """Calculate base composition ratios for a genomic region."""
        if chr_name not in self.ref_sequences:
            print(f"Warning: Chromosome {chr_name} not found in reference")
            return (0.0, 0.0, 0.0)
            
        sequence = self.ref_sequences[chr_name]
        
        # Extract subsequence (convert to 0-based indexing)
        if start >= len(sequence) or end > len(sequence) or start >= end:
            print(f"Warning: Invalid coordinates {chr_name}:{start}-{end}")
            return (0.0, 0.0, 0.0)
            
        subseq = sequence[start:end]
        
        # Count bases
        base_counts = Counter(subseq)
        total = len(subseq)
        
        if total == 0:
            return (0.0, 0.0, 0.0)
        
        a_count = base_counts.get('A', 0)
        t_count = base_counts.get('T', 0)
        c_count = base_counts.get('C', 0)
        g_count = base_counts.get('G', 0)
        
        ag_ratio = (a_count + g_count) / total
        ct_ratio = (c_count + t_count) / total
        gc_ratio = (g_count + c_count) / total
        
        return (ag_ratio, ct_ratio, gc_ratio)
        
    def _write_merged_file(self, output_file: str) -> None:
        """Write the merged output file with base composition ratios."""
        print(f"Writing merged results to: {output_file}")
        
        # Create output directory if it doesn't exist
        output_dir = Path(output_file).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w', encoding='utf-8') as f:
            # Write header
            header = [
                "Chr", "Start", "End", 
                "File1_Col1", "File1_Col2",
                "File2_Col1", "File2_Col2", 
                "AG_ratio", "CT_ratio", "GC_ratio"
            ]
            f.write('\t'.join(header) + '\n')
            
            # Process each region in BED file
            for key in self.reference_data.keys():
                chr_name, start, end = key
                
                # Get data from file1 (or default to -1)
                file1_cols = self.file1_data.get(key, ("-1", "-1"))
                file1_col1, file1_col2 = file1_cols
                
                # Get data from file2 (or default to -1)
                file2_cols = self.file2_data.get(key, ("-1", "-1"))
                file2_col1, file2_col2 = file2_cols
                
                # Calculate base composition ratios
                ag_ratio, ct_ratio, gc_ratio = self._calculate_base_ratios(chr_name, start, end)
                
                # Write merged line
                output_line = (
                    f"{chr_name}\t{start}\t{end}\t"
                    f"{file1_col1}\t{file1_col2}\t"
                    f"{file2_col1}\t{file2_col2}\t"
                    f"{ag_ratio:.4f}\t{ct_ratio:.4f}\t{gc_ratio:.4f}\n"
                )
                f.write(output_line)
                
        print(f"Merged file written successfully: {output_file}")
        print(f"Total regions processed: {len(self.reference_data)}")


def create_parser() -> argparse.ArgumentParser:
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(
        description="Merge files based on BED coordinates and calculate base composition",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python merge_files.py regions.bed file1.txt file2.txt reference.fasta output.txt
  
Input files:
  - BED file: Contains genomic coordinates (chr, start, end)
  - File1/File2: Data files with coordinates in last 5th, 4th, 3rd columns
  - Reference: FASTA format reference genome
  
Output format:
  - Columns 1-3: BED coordinates (chr, start, end)
  - Columns 4-5: File1 last two columns
  - Columns 6-7: File2 last two columns
  - Columns 8-10: AG_ratio, CT_ratio, GC_ratio from reference sequence
        """
    )
    
    parser.add_argument(
        'bed_file',
        help='BED format file with genomic coordinates'
    )
    
    parser.add_argument(
        'file1',
        help='First data file to merge'
    )
    
    parser.add_argument(
        'file2', 
        help='Second data file to merge'
    )
    
    parser.add_argument(
        'ref_file',
        help='Reference genome FASTA file'
    )
    
    parser.add_argument(
        'output_file',
        help='Output file path for merged results'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 2.0.0'
    )
    
    return parser


def main():
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    try:
        # Validate input files exist
        for file_path in [args.bed_file, args.file1, args.file2, args.ref_file]:
            if not Path(file_path).exists():
                raise FileNotFoundError(f"Input file not found: {file_path}")
        
        merger = BedFileMerger()
        merger.merge_files(args.bed_file, args.file1, args.file2, args.ref_file, args.output_file)
        print("File merging and base composition calculation completed successfully!")
        
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nOperation interrupted by user", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main() 
