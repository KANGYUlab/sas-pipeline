#!/usr/bin/env python3
"""
Integrated Structural Error Analysis Tool

This script integrates clip, indel, and depth analysis to identify anomalous regions
in BAM files using sliding windows with customizable thresholds.

Features:
- Analyzes clips, indels, and depth in a unified framework
- Identifies anomalous regions based on configurable criteria
- Special handling for chromosome terminal regions (10kb)
- Multi-threaded processing for performance

Anomaly Detection Criteria:
1. Non-terminal regions (outside first/last 10kb):
   - Depth anomalies: < 20% or > 200% of mean depth
   - INDEL/Clip anomalies: ReadsWithInsertion/ReadsWithDeletion/ClipCount ratio > 0.5
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
    - Python 3.6+

Usage:
    python integrated_sv_analysis.py input.bam reference.fa output_prefix [options]
"""

import sys
import os
import argparse
import subprocess
import gzip
import shutil
import traceback
from pathlib import Path
from collections import defaultdict
from typing import Dict, Tuple, Set, List, DefaultDict, Optional
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd
import pysam


class IntegratedSECounters:
    """
    Thread-safe counters for integrated SE analysis.
    
    This class ensures that multiple threads can safely update all SE counts
    when processing different chromosomes simultaneously.
    """
    
    def __init__(self):
        self.lock = threading.Lock()
        self.windows: Dict[Tuple[str, int], Dict[str, int]] = {}
    
    def init_chromosome_windows(self, chrom: str, chrom_length: int, window_size: int):
        """Initialize windows for a specific chromosome."""
        with self.lock:
            for start in range(0, chrom_length, window_size):
                if start >= chrom_length:
                    break
                self.windows[(chrom, start)] = {
                    'ClipCount': 0,
                    'Insertion': 0,
                    'Deletion': 0,
                    'ReadsWithInsertion': 0,
                    'ReadsWithDeletion': 0,
                    'ReadsWithClip': 0,
                    'TotalReads': 0,
                    'ReadsWithoutINDELorClip': 0,
                    'Depth': 0.0
                }
    
    def update_windows(self, local_windows: Dict[Tuple[str, int], Dict[str, int]]):
        """
        Thread-safe update of window counters.
        """
        with self.lock:
            for window_key, counts in local_windows.items():
                if window_key in self.windows:
                    for count_type, count_value in counts.items():
                        self.windows[window_key][count_type] += count_value


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
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    chrom = parts[0]
                    start = int(parts[1])
                    depth = float(parts[3])
                    
                    # Convert to window start (0-based to match other analyses)
                    window_start = (start // window_size) * window_size
                    depth_data[(chrom, window_start)] = depth
    
    except Exception as e:
        print(f"Error loading depth data: {e}", file=sys.stderr)
    
    return depth_data


def calculate_affected_windows(start_pos: int, deletion_length: int, window_size: int) -> List[int]:
    """
    Calculate all windows affected by a deletion.
    """
    start_window = (start_pos // window_size) * window_size
    end_pos = start_pos + deletion_length
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
    window_size: int = 2000
) -> Dict[Tuple[str, int], Dict[str, int]]:
    """
    Process all reads for a specific chromosome and collect integrated SE data.
    
    Args:
        bam_file: Path to BAM file
        chrom: Chromosome name to process
        chrom_length: Length of the chromosome
        min_indel_length: Minimum INDEL length to consider (default: 50bp)
        window_size: Size of sliding window
    
    Returns:
        Dictionary with local window coordinates and integrated SE counts for this chromosome
    
    Note:
        ReadsWithoutINDELorClip counts reads that have:
        - NO significant INDELs (length >= min_indel_length, default 50bp) in the current window
        - NO clip events in the current window
        This is ONLY calculated in terminal regions (first/last 10kb) and used for terminal region anomaly detection.
        A read is counted if it has no clip/indel events within the current window, regardless of read length.
    """
    # Initialize local windows for this chromosome
    local_windows = {}
    for start in range(0, chrom_length, window_size):
        if start >= chrom_length:
            break
        local_windows[(chrom, start)] = {
            'ClipCount': 0,
            'Insertion': 0,
            'Deletion': 0,
            'ReadsWithInsertion': 0,
            'ReadsWithDeletion': 0,
            'ReadsWithClip': 0,
            'TotalReads': 0,
            'ReadsWithoutINDELorClip': 0,
            'Depth': 0.0
        }
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            # Process reads for this specific chromosome
            for read in bam.fetch(chrom):
                if read.is_unmapped or read.cigartuples is None:
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
                            
                            if window_key in local_windows:
                                local_windows[window_key]['ClipCount'] += 1
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
                            if window_key in local_windows:
                                local_windows[window_key]['Insertion'] += 1
                                read_insertion_windows.add(window_key)
                        
                        # Update deletion counts
                        elif op == 2:
                            # Calculate all windows affected by this deletion
                            affected_windows = calculate_affected_windows(current_pos, length, window_size)
                            for affected_window in affected_windows:
                                affected_window_key = (chrom, affected_window)
                                if affected_window_key in local_windows:
                                    local_windows[affected_window_key]['Deletion'] += 1
                                    read_deletion_windows.add(affected_window_key)
                    
                    # Update reference position for operations that consume reference sequence
                    if op in [0, 2, 3, 7, 8]:  # M, D, N, =, X
                        current_pos += length
                
                # Update read counts for each window this read spans
                read_start_window = (read.reference_start // window_size) * window_size
                read_end_window = (read.reference_end // window_size) * window_size if read.reference_end else read_start_window
                
                for window_start in range(read_start_window, read_end_window + window_size, window_size):
                    window_key = (chrom, window_start)
                    if window_key in local_windows:
                        local_windows[window_key]['TotalReads'] += 1
                        
                        # Count reads without significant INDEL (≥50bp) or non-terminal clip for terminal region analysis
                        # ONLY calculate ReadsWithoutINDELorClip in terminal regions (first/last 10kb)
                        terminal_region_size = 10000
                        
                        # Calculate actual window end for this window
                        actual_window_end = min(window_start + window_size, chrom_length)
                        
                        # Check if this is a terminal region:
                        # - Front 10kb: window start < 10000  
                        # - Back 10kb: window start >= (chromosome_length - 10000)
                        is_terminal_region = (window_start < terminal_region_size or 
                                            window_start >= chrom_length - terminal_region_size)
                        
                        if is_terminal_region:
                            # Check if this read has any clip or indel events in the current window
                            # A read is counted as ReadsWithoutINDELorClip if it has no clip/indel events in this window
                            read_has_clip_in_window = window_key in read_clip_windows
                            read_has_indel_in_window = (window_key in read_insertion_windows or 
                                                      window_key in read_deletion_windows)
                            
                            if not read_has_clip_in_window and not read_has_indel_in_window:
                                local_windows[window_key]['ReadsWithoutINDELorClip'] += 1
                
                # Update read counts for windows with insertions
                for window_key in read_insertion_windows:
                    if window_key in local_windows:
                        local_windows[window_key]['ReadsWithInsertion'] += 1
                
                # Update read counts for windows with deletions
                for window_key in read_deletion_windows:
                    if window_key in local_windows:
                        local_windows[window_key]['ReadsWithDeletion'] += 1
                
                # Update read counts for windows with clips
                for window_key in read_clip_windows:
                    if window_key in local_windows:
                        local_windows[window_key]['ReadsWithClip'] += 1
    
    except Exception as e:
        print(f"Error processing chromosome {chrom}: {e}", file=sys.stderr)
    
    return local_windows


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
    Identify anomalous regions based on the integrated criteria.
    """
    anomalous_regions = []
    
    # First, calculate depth statistics
    all_depths = [depth for depth in depth_data.values() if depth > 0]
    if not all_depths:
        print("Warning: No depth data found")
        return anomalous_regions
    
    mean_depth = np.mean(all_depths)
    
    # Updated thresholds: < 20% or > 200% of mean depth
    depth_lower_bound = mean_depth * 0.20
    depth_upper_bound = mean_depth * 2.0
    
    print(f"Depth statistics:")
    print(f"  Mean depth: {mean_depth:.4f}")
    print(f"  Lower bound (20% of mean): {depth_lower_bound:.4f}")
    print(f"  Upper bound (200% of mean): {depth_upper_bound:.4f}")
    
    # Process each window
    for (chrom, window_start), counts in windows.items():
        chrom_length = chrom_lengths.get(chrom, 0)
        window_end = min(window_start + window_size, chrom_length)  # Use actual chromosome length for terminal windows
        
        # Get depth for this window
        local_depth = depth_data.get((chrom, window_start), 0.0)
        
        # Check if this is a terminal region (first or last 10kb)
        # Use window_start to ensure strict 10kb boundary
        is_terminal = (window_start < terminal_region_size) or \
                     (window_start >= chrom_length - terminal_region_size)
        
        anomaly_reasons = []
        
        # Check anomalies based on region type
        if not is_terminal:
            # NON-TERMINAL REGIONS: Check both depth and INDEL/clip anomalies
            
            # 1. Check depth anomalies
            if local_depth > 0:
                if local_depth < depth_lower_bound:
                    anomaly_reasons.append(f"Low depth ({local_depth:.2f} < {depth_lower_bound:.2f})")
                elif local_depth > depth_upper_bound:
                    anomaly_reasons.append(f"High depth ({local_depth:.2f} > {depth_upper_bound:.2f})")
            
            # 2. Check INDEL and clip anomalies
            total_reads = counts['TotalReads']
            
            # Check clips
            clip_count = counts['ClipCount']
            if total_reads > 0 and (clip_count / total_reads) > 0.5:
                anomaly_reasons.append(f"High clip ratio ({clip_count}/{total_reads} = {clip_count/total_reads:.3f})")
            
            # Check insertions (use reads with insertions, not insertion count)
            insertion_reads = counts['ReadsWithInsertion']
            if total_reads > 0 and (insertion_reads / total_reads) > 0.5:
                anomaly_reasons.append(f"High insertion reads ratio ({insertion_reads}/{total_reads} = {insertion_reads/total_reads:.3f})")
            
            # Check deletions (use reads with deletions, not deletion count)
            deletion_reads = counts['ReadsWithDeletion']
            if total_reads > 0 and (deletion_reads / total_reads) > 0.5:
                anomaly_reasons.append(f"High deletion reads ratio ({deletion_reads}/{total_reads} = {deletion_reads/total_reads:.3f})")
        
        else:
            # TERMINAL REGIONS (first/last 10kb): Only check reads without INDEL/clip, NO depth anomaly detection
            reads_without_indel_clip = counts['ReadsWithoutINDELorClip']
            total_reads = counts['TotalReads']
            
            # Anomaly: reads_without_indel_clip < 1 OR ratio <= 10%
            if total_reads > 0:
                ratio = reads_without_indel_clip / total_reads
                if reads_without_indel_clip < 1 or ratio <= 0.1:
                    anomaly_reasons.append(f"Terminal region anomaly: reads without INDEL/clip ({reads_without_indel_clip}/{total_reads} = {ratio:.3f}) is too low")
        
        # If any anomalies were found, add this region
        if anomaly_reasons:
            error_type = classify_error_type('; '.join(anomaly_reasons))
            anomalous_regions.append({
                'chromosome': chrom,
                'start': window_start + 1,  # Convert to 1-based
                'end': window_end,
                'is_terminal': is_terminal,
                'depth': local_depth,
                'clip_count': counts['ClipCount'],
                'insertion_count': counts['Insertion'],
                'deletion_count': counts['Deletion'],
                'reads_without_indel_clip': counts['ReadsWithoutINDELorClip'],
                'total_reads': counts['TotalReads'],
                'anomaly_reasons': '; '.join(anomaly_reasons),
                'error_type': error_type
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


def save_anomaly_summary(anomalous_regions, chrom_lengths, output_prefix, terminal_type_counts=None):
    from collections import defaultdict
    summary = defaultdict(lambda: defaultdict(int))
    for region in anomalous_regions:
        etype = region['error_type']
        chrom = region['chromosome']
        if etype in ("low_depth", "high_depth", "indel", "clip"):
            summary[chrom][etype] += 1
    if terminal_type_counts is None:
        terminal_type_counts = defaultdict(lambda: {'indel': 0, 'clip': 0})
    all_chroms = sorted(chrom_lengths.keys())
    summary_file = f"{output_prefix}.anomaly_summary.tsv"
    with open(summary_file, 'w') as f:
        f.write("Chromosome\tlow_depth\thigh_depth\tindel\tclip\tterminal_indel\tterminal_clip\n")
        total_row = ["Total", 0, 0, 0, 0, 0, 0]
        for chrom in all_chroms:
            row = [chrom]
            indel_val = summary[chrom]["indel"]
            clip_val = summary[chrom]["clip"]
            t_indel = terminal_type_counts[chrom]["indel"]
            t_clip = terminal_type_counts[chrom]["clip"]
            row.append(str(summary[chrom]["low_depth"]))
            row.append(str(summary[chrom]["high_depth"]))
            row.append(str(indel_val + t_indel))
            row.append(str(clip_val + t_clip))
            row.append(str(t_indel))
            row.append(str(t_clip))
            total_row[1] += summary[chrom]["low_depth"]
            total_row[2] += summary[chrom]["high_depth"]
            total_row[3] += indel_val + t_indel
            total_row[4] += clip_val + t_clip
            total_row[5] += t_indel
            total_row[6] += t_clip
            f.write("\t".join(row) + "\n")
        total_row = [str(x) for x in total_row]
        f.write("\t".join(total_row) + "\n")
    print(f"Anomaly summary saved to: {summary_file}")


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
                "TotalReads\tReadsWithoutINDELorClip\n")
        
        for (chrom, window_start), counts in sorted(windows.items()):
            chrom_length = chrom_lengths.get(chrom, window_start + window_size)
            window_end = min(window_start + window_size, chrom_length)
            depth = depth_data.get((chrom, window_start), 0.0)
            is_terminal = (window_start < terminal_region_size) or (window_start >= chrom_length - terminal_region_size)
            if is_terminal:
                reads_without_indel_clip = str(counts['ReadsWithoutINDELorClip'])
            else:
                reads_without_indel_clip = '-'
            f.write(f"{chrom}\t{window_start + 1}\t{window_end}\t{depth:.4f}\t"
                   f"{counts['ClipCount']}\t{counts['Insertion']}\t{counts['Deletion']}\t"
                   f"{counts['ReadsWithInsertion']}\t{counts['ReadsWithDeletion']}\t"
                   f"{counts['ReadsWithClip']}\t{counts['TotalReads']}\t"
                   f"{reads_without_indel_clip}\n")
    
    print(f"Window statistics saved to: {stats_file}")
    
    # Save anomalous regions
    anomalies_file = f"{output_prefix}.anomalous_regions.tsv"
    with open(anomalies_file, 'w') as f:
        f.write("Chromosome\tStart\tEnd\tIsTerminal\tDepth\tClipCount\t"
                "InsertionCount\tDeletionCount\tReadsWithoutINDELClip\t"
                "TotalReads\tAnomalyReasons\tErrorType\n")
        
        for region in anomalous_regions:
            f.write(f"{region['chromosome']}\t{region['start']}\t{region['end']}\t"
                   f"{region['is_terminal']}\t{region['depth']:.4f}\t"
                   f"{region['clip_count']}\t{region['insertion_count']}\t"
                   f"{region['deletion_count']}\t{region['reads_without_indel_clip']}\t"
                   f"{region['total_reads']}\t{region['anomaly_reasons']}\t{region['error_type']}\n")
    
    print(f"Anomalous regions saved to: {anomalies_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Integrated Structural Error Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage:
    python integrated_sv_analysis.py input.bam reference.fa sample_name
    
    # Custom parameters:
    python integrated_sv_analysis.py input.bam reference.fa sample_name --window-size 1000 --threads 8 --min-indel-length 30
    
    # Specify output directory:
    python integrated_sv_analysis.py input.bam reference.fa sample_name --output-dir /path/to/output
        """
    )
    
    parser.add_argument("bam_file", help="Input BAM file path")
    parser.add_argument("reference_fa", help="Reference FASTA file path")
    parser.add_argument("output_prefix", help="Output prefix for result files")
    parser.add_argument("--output-dir", type=str, default=".",
                       help="Output directory for result files (default: current directory)")
    parser.add_argument("--window-size", type=int, default=2000,
                       help="Window size for analysis (default: 2000)")
    parser.add_argument("--min-indel-length", type=int, default=50,
                       help="Minimum INDEL length to consider (default: 50)")
    parser.add_argument("--threads", type=int, default=4,
                       help="Number of threads to use (default: 4)")
    parser.add_argument("--terminal-region-size", type=int, default=10000,
                       help="Size of terminal regions for special analysis (default: 10000)")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.bam_file):
        print(f"Error: BAM file not found: {args.bam_file}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.reference_fa):
        print(f"Error: Reference FASTA file not found: {args.reference_fa}", file=sys.stderr)
        sys.exit(1)
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    full_output_prefix = output_dir / args.output_prefix
    try:
        print(f"Starting integrated SE analysis for: {args.bam_file}")
        print(f"Reference: {args.reference_fa}")
        print(f"Output directory: {output_dir}")
        print(f"Output prefix: {args.output_prefix}")
        print(f"Window size: {args.window_size} bp")
        print(f"Minimum INDEL length: {args.min_indel_length} bp")
        print(f"Threads: {args.threads}")
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
        print("\n" + "="*50)
        print("STEP 4: Processing chromosomes for integrated SE analysis")
        print("="*50)
        global_counters = IntegratedSECounters()
        for chrom, chrom_length in chrom_lengths.items():
            global_counters.init_chromosome_windows(chrom, chrom_length, args.window_size)
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            future_to_chrom = {
                executor.submit(
                    process_chromosome_integrated,
                    args.bam_file,
                    chrom,
                    chrom_length,
                    args.min_indel_length,
                    args.window_size
                ): chrom
                for chrom, chrom_length in chrom_lengths.items()
            }
            completed_chroms = 0
            for future in as_completed(future_to_chrom):
                chrom = future_to_chrom[future]
                try:
                    local_windows = future.result()
                    global_counters.update_windows(local_windows)
                    completed_chroms += 1
                    print(f"Completed chromosome {chrom} ({completed_chroms}/{len(chrom_lengths)})")
                except Exception as e:
                    print(f"Error processing chromosome {chrom}: {e}", file=sys.stderr)
        print("\n" + "="*50)
        print("STEP 6: Integrating depth data")
        print("="*50)
        for window_key, depth in depth_data.items():
            if window_key in global_counters.windows:
                global_counters.windows[window_key]['Depth'] = depth
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
        terminal_type_counts, _, _ = terminal_reads_analysis(
            global_counters.windows,
            depth_data,
            chrom_lengths,
            args.terminal_region_size
        )
        save_anomaly_summary(
            anomalous_regions,
            chrom_lengths,
            str(full_output_prefix),
            terminal_type_counts
        )
        print("\n" + "="*50)
        print("ANALYSIS COMPLETE")
        print("="*50)
        print(f"Total windows analyzed: {len(global_counters.windows)}")
        print(f"Anomalous regions found: {len(anomalous_regions)}")
        print(f"Results saved to: {output_dir}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main() 
