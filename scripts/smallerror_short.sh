#!/bin/bash

# Small Error Region Detection for Short-read BAM files
# Description: This script identifies 50bp windows containing bases with low coverage.
# It first runs mosdepth to get per-base coverage and calculates a dynamic depth
# threshold (10% of the mean coverage). It then identifies all bases below this
# threshold and uses bedtools intersect to find the 50bp windows that contain
# one or more of these low-coverage bases.

set -euo pipefail

# Function to display usage
usage() {
    echo "Usage: $0 -i INPUT_BAM -o OUTPUT_PREFIX [-t THREADS]"
    echo ""
    echo "Options:"
    echo "  -i INPUT_BAM        Input BAM file path"
    echo "  -o OUTPUT_PREFIX    Output file prefix"
    echo "  -t THREADS          Number of threads (default: 4)"
    echo "  -h                  Display this help message"
    echo ""
    echo "Example:"
    echo "  $0 -i sample.bam -o sample_smallerror -t 8"
    exit 1
}

# Default values
THREADS=4
INPUT_BAM=""
OUTPUT_PREFIX=""

# Parse command line arguments
while getopts "i:o:t:h" opt; do
    case $opt in
        i)
            INPUT_BAM="$OPTARG"
            ;;
        o)
            OUTPUT_PREFIX="$OPTARG"
            ;;
        t)
            THREADS="$OPTARG"
            ;;
        h)
            usage
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
    esac
done

# Check required arguments
if [[ -z "$INPUT_BAM" || -z "$OUTPUT_PREFIX" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if input BAM file exists
if [[ ! -f "$INPUT_BAM" ]]; then
    echo "Error: Input BAM file '$INPUT_BAM' does not exist"
    exit 1
fi

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_PREFIX")
if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

# Check if required tools are available
for tool in mosdepth bedtools samtools; do
    if ! command -v $tool &> /dev/null; then
        echo "Error: $tool is not installed or not in PATH."
        echo "Please ensure samtools, bedtools, and mosdepth are installed."
        exit 1
    fi
done

echo "Starting small error region detection..."
echo "Input BAM: $INPUT_BAM"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Window size: 50bp"
echo "Threads: $THREADS"

# Step 1: Create a genome file from BAM header
echo "Step 1: Creating genome file from BAM header..."
GENOME_FILE="${OUTPUT_PREFIX}.genome"
samtools view -H "$INPUT_BAM" | grep '^@SQ' | cut -f 2,3 | sed 's/SN://' | sed 's/LN://' > "$GENOME_FILE"

if [[ ! -s "$GENOME_FILE" ]]; then
    echo "Error: Failed to create genome file from BAM header. No @SQ lines found?"
    exit 1
fi

# Step 2: Create 50bp windows using bedtools
echo "Step 2: Creating 50bp windows with bedtools..."
WINDOWS_BED="${OUTPUT_PREFIX}.windows.bed"
bedtools makewindows -g "$GENOME_FILE" -w 50 > "$WINDOWS_BED"

# Step 3: Calculate per-base depth using mosdepth
echo "Step 3: Calculating per-base depth with mosdepth..."
# We run without -x to ensure the .per-base.bed.gz file is created.
# We also run without -n to ensure the summary file is created for threshold calculation.
mosdepth --threads "$THREADS" "$OUTPUT_PREFIX" "$INPUT_BAM"

# Check if mosdepth completed successfully
if [[ $? -ne 0 ]]; then
    echo "Error: mosdepth failed"
    exit 1
fi

# Step 4: Calculate depth threshold from mosdepth summary
echo "Step 4: Calculating depth threshold..."
MOSDEPTH_SUMMARY="${OUTPUT_PREFIX}.mosdepth.summary.txt"
if [[ ! -f "$MOSDEPTH_SUMMARY" ]]; then
    echo "Error: Mosdepth summary file '$MOSDEPTH_SUMMARY' not found."
    exit 1
fi

MEAN_DEPTH=$(tail -n 1 "$MOSDEPTH_SUMMARY" | awk '{print $4}')
if [[ -z "$MEAN_DEPTH" || "$MEAN_DEPTH" == "0" ]]; then
    echo "Error: Could not extract a valid mean depth from mosdepth summary."
    exit 1
fi

DEPTH_THRESHOLD=$(awk -v mean="$MEAN_DEPTH" 'BEGIN {printf "%.4f", mean * 0.1}')

echo "  Average genome-wide depth: $MEAN_DEPTH"
echo "  Calculated depth threshold (10% of average): $DEPTH_THRESHOLD"


# Step 5: Extract low coverage bases from per-base depth file
echo "Step 5: Identifying low coverage bases (depth < $DEPTH_THRESHOLD)..."
PER_BASE_BED_GZ="${OUTPUT_PREFIX}.per-base.bed.gz"
LOW_COVERAGE_BASES_BED="${OUTPUT_PREFIX}.low_coverage_bases.bed"

if [[ ! -f "$PER_BASE_BED_GZ" ]]; then
    echo "Error: Expected per-base depth file '$PER_BASE_BED_GZ' not found"
    exit 1
fi

zcat "$PER_BASE_BED_GZ" | awk -v threshold="$DEPTH_THRESHOLD" 'BEGIN{OFS="\t"} $4 < threshold {print $1, $2, $3}' > "$LOW_COVERAGE_BASES_BED"

# Step 6: Intersect low-coverage bases with 50bp windows
echo "Step 6: Locating 50bp windows containing low-coverage bases..."
SMALLERROR_BED="${OUTPUT_PREFIX}.smallerror.bed"
bedtools intersect -a "$WINDOWS_BED" -b "$LOW_COVERAGE_BASES_BED" -u > "$SMALLERROR_BED"

# Step 7: Generate summary statistics
echo "Step 7: Generating summary statistics..."
TOTAL_WINDOWS=$(wc -l < "$WINDOWS_BED")
SMALLERROR_REGIONS=$(wc -l < "$SMALLERROR_BED")

if [[ $TOTAL_WINDOWS -gt 0 ]]; then
    SMALLERROR_PERCENTAGE=$(awk -v small="$SMALLERROR_REGIONS" -v total="$TOTAL_WINDOWS" 'BEGIN {printf "%.2f", (small/total)*100}')
else
    SMALLERROR_PERCENTAGE="0.00"
fi

echo "Summary:"
echo "  Total 50bp windows analyzed: $TOTAL_WINDOWS"
echo "  Small error windows (containing bases with depth < $DEPTH_THRESHOLD): $SMALLERROR_REGIONS"
echo "  Percentage of small error windows: ${SMALLERROR_PERCENTAGE}%"

# Step 8: Create summary file
SUMMARY_FILE="${OUTPUT_PREFIX}.smallerror.summary.txt"
cat > "$SUMMARY_FILE" << EOF
Small Error Region Detection Summary
=====================================
Input BAM: $INPUT_BAM
Output prefix: $OUTPUT_PREFIX
Window size: 50bp
Analysis date: $(date)

Thresholding:
  Mean genome-wide depth: $MEAN_DEPTH
  Depth threshold (10% of mean): $DEPTH_THRESHOLD

Results:
  Total 50bp windows analyzed: $TOTAL_WINDOWS
  Small error windows (containing bases with depth < $DEPTH_THRESHOLD): $SMALLERROR_REGIONS
  Percentage: ${SMALLERROR_PERCENTAGE}%

Output files:
  - ${SMALLERROR_BED} (small error windows)
  - ${SUMMARY_FILE} (this summary)
  - ${PER_BASE_BED_GZ} (per-base depth)
  - ${WINDOWS_BED} (all 50bp windows generated)
  - ${GENOME_FILE} (genome sizes from BAM)
  - ${LOW_COVERAGE_BASES_BED} (temporary file with low-coverage base positions)
EOF

echo ""
echo "Analysis completed successfully!"
echo "Output files:"
echo "  - Small error regions BED: $SMALLERROR_BED"
echo "  - Summary report: $SUMMARY_FILE"
echo "  - Per-base depth file: ${PER_BASE_BED_GZ}"