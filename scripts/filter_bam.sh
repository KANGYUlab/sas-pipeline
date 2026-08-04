#!/bin/bash

# BAM Filter Tool - High-performance BAM filtering for genomics analysis
# Optimized for different sequencing technologies
set -euo pipefail

PROGRAM_NAME=${0##*/}
VERSION="2.0.0"

# Detect a portable default without assuming GNU coreutils is installed.
detect_default_threads() {
    local detected=""

    if command -v getconf >/dev/null 2>&1; then
        detected=$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)
        if [[ "$detected" =~ ^[1-9][0-9]*$ ]]; then
            printf '%s\n' "$detected"
            return
        fi
    fi

    if command -v sysctl >/dev/null 2>&1; then
        detected=$(sysctl -n hw.ncpu 2>/dev/null || true)
        if [[ "$detected" =~ ^[1-9][0-9]*$ ]]; then
            printf '%s\n' "$detected"
            return
        fi
    fi

    printf '1\n'
}

# Default parameters
MODE=""
MAX_SPLITS_BASE=3.0
MAX_SPLITS_KB=0.05
THREADS=$(detect_default_threads)
INPUT_BAM=""
OUTPUT_BAM=""
QUIET=false
PRIMARY_ONLY=true  # Only used in long-read mode

# Function to show usage
usage() {
    cat << EOF
Filter BAM alignments for SAS analysis.

Usage:
  $PROGRAM_NAME MODE -i INPUT_BAM -o OUTPUT_BAM [OPTIONS]

Modes:
  short-read              Keep perfect matches (NM=0) without S/H clipping
  long-read               Filter by SA split count; primary alignments only by default

Required arguments:
  -i, --input INPUT_BAM   Input BAM file
  -o, --output OUTPUT_BAM Filtered BAM output file

General options:
  -t, --threads N         samtools threads (default: $THREADS)
  -q, --quiet             Suppress progress and summary messages
  -h, --help              Show this help and exit
  -v, --version           Show version and exit

Long-read options:
      --max-splits-base N Base split threshold (default: $MAX_SPLITS_BASE)
      --max-splits-kb N   Additional splits allowed per kb (default: $MAX_SPLITS_KB)
      --include-all-reads Include secondary and supplementary alignments

Examples:
  $PROGRAM_NAME short-read -i illumina.bam -o perfect.bam -t 16
  $PROGRAM_NAME long-read -i ont.bam -o filtered.bam --max-splits-base 3
EOF
}

option_value_required() {
    local option="$1"
    local value="${2-}"
    if [[ -z "$value" ]]; then
        echo "Error: $option requires a value" >&2
        usage >&2
        exit 2
    fi
}

# Short-read filtering function (perfect matches, no clipping)
filter_short_read_mode() {
    local input_bam="$1"
    local output_bam="$2"
    local threads="$3"
    local quiet="$4"
    
    [[ "$quiet" == false ]] && echo "=== Short-Read Mode: Perfect matches with no clipping ==="
    [[ "$quiet" == false ]] && echo "Filtering criteria: NM=0 (perfect matches) + no soft/hard clipping"
    [[ "$quiet" == false ]] && echo "Processing all reads (primary, secondary, and supplementary)"
    
    # Build samtools arguments (no primary-only filtering for short reads)
    local samtools_args=(view -h -@ "$threads" -e '[NM]==0')
    
    # Count input reads for statistics
    local total_reads=0
    if [[ "$quiet" == false ]]; then
        echo "Counting input reads..."
        total_reads=$(samtools view -c -@ "$threads" "$input_bam")
        echo "Total reads: $total_reads"
    fi
    
    [[ "$quiet" == false ]] && echo "Applying filters: NM=0 and no soft/hard clipping..."
    
    # Create output directory if it doesn't exist
    local output_dir
    output_dir=$(dirname "$output_bam")
    if [[ ! -d "$output_dir" ]]; then
        echo "Creating output directory: $output_dir"
        mkdir -p "$output_dir"
    fi
    
    # Ultra-fast short-read filtering pipeline
    if ! samtools "${samtools_args[@]}" "$input_bam" | \
         awk 'BEGIN{OFS="\t"} /^@/ {print; next} $6 !~ /[SH]/ {print}' | \
         samtools view -b -@ "$threads" -o "$output_bam" -; then
        echo "Error: Short-read filtering pipeline failed" >&2
        exit 1
    fi
    
    # Create index
    [[ "$quiet" == false ]] && echo "Creating index..."
    if ! samtools index -@ "$threads" "$output_bam"; then
        echo "Error: Failed to create index" >&2
        exit 1
    fi
    
    # Calculate statistics
    if [[ "$quiet" == false ]]; then
        echo "Calculating statistics..."

        # Count reads from output to get kept reads
        local kept_reads
        kept_reads=$(samtools view -c -@ "$threads" "$output_bam")
        local filtered_reads
        filtered_reads=$((total_reads - kept_reads))

        # Perform detailed analysis for filtered reads (this requires extra passes)
        # 1. Count reads with NM > 0, and among them, count those that are also clipped.
        local stats_nm_gt_0
        stats_nm_gt_0=$(samtools view -@ "$threads" -e '[NM]>0' "$input_bam" | awk 'BEGIN{total=0; clipped=0} {total++; if($6~/[SH]/) clipped++} END{print total, clipped}')
        local nm_gt_0_count
        nm_gt_0_count=$(echo "$stats_nm_gt_0" | cut -d' ' -f1)
        local both_nm_and_clip_count
        both_nm_and_clip_count=$(echo "$stats_nm_gt_0" | cut -d' ' -f2)

        # 2. Count reads that are clipped but have NM=0.
        local clip_only_nm_ok_count
        clip_only_nm_ok_count=$(samtools view -@ "$threads" -e '[NM]==0' "$input_bam" | awk 'BEGIN{count=0} $6~/[SH]/{count++} END{print count}')
        
        # 3. Calculate reads filtered ONLY due to NM > 0 (and not clipping).
        local nm_only_count
        nm_only_count=$((nm_gt_0_count - both_nm_and_clip_count))
        
        local percentage="0.0"
        if (( total_reads > 0 )); then
            percentage=$(awk "BEGIN {printf \"%.1f\", $kept_reads/$total_reads*100}")
        fi
        
        echo "Results:"
        echo "  Total reads: $total_reads"
        echo "  Kept reads: $kept_reads ($percentage%)"
        echo "  Total filtered reads: $filtered_reads"
        echo "    - Filtered due to clipping only (NM=0): $clip_only_nm_ok_count"
        echo "    - Filtered due to NM > 0 only (no clip): $nm_only_count"
        echo "    - Filtered due to both NM > 0 and clipping: $both_nm_and_clip_count"
        echo "Short-read filtering completed successfully!"
    fi
}

# Long-read filtering function (split alignment filtering)
filter_long_read_mode() {
    local input_bam="$1"
    local output_bam="$2"
    local threads="$3"
    local primary_only="$4"
    local quiet="$5"
    local max_splits_base="$6"
    local max_splits_kb="$7"
    
    [[ "$quiet" == false ]] && echo "=== Long-Read Mode: Filter by split alignment count ==="
    [[ "$quiet" == false ]] && echo "Split threshold = $max_splits_base + $max_splits_kb * (read_length_kb)"
    
    # Build samtools view command with optional primary-only filter
    local samtools_args=(view -h -@ "$threads")
    if [[ "$primary_only" == true ]]; then
        # -F 0x904 excludes: unmapped(4) + secondary(256) + supplementary(2048)
        samtools_args+=(-F 0x904)
        [[ "$quiet" == false ]] && echo "Filtering to primary alignments only (FLAG & 0x904 == 0)"
    else
        [[ "$quiet" == false ]] && echo "Processing all alignments (primary, secondary, and supplementary)"
    fi
    
    # Count input reads for statistics
    local total_reads=0
    if [[ "$quiet" == false ]]; then
        echo "Counting input reads..."
        if [[ "$primary_only" == true ]]; then
            total_reads=$(samtools view -c -@ "$threads" -F 0x904 "$input_bam")
            echo "Total primary reads: $total_reads"
        else
            total_reads=$(samtools view -c -@ "$threads" "$input_bam")
            echo "Total reads: $total_reads"
        fi
    fi
    
    # Create AWK script for filtering
    local awk_script
    awk_script=$(cat << 'EOF'
BEGIN {
    max_splits_base = ENVIRON["MAX_SPLITS_BASE"]
    max_splits_kb = ENVIRON["MAX_SPLITS_KB"]
}

# Pass through header lines
/^@/ {
    print
    next
}

# Process alignment lines
{
    # Default: keep the read
    keep = 1
    
    # Check for SA tag
    for (i = 12; i <= NF; i++) {
        if ($i ~ /^SA:Z:/) {
            sa_value = substr($i, 6)  # Remove "SA:Z:" prefix
            
            # Count splits by counting semicolons, excluding empty entries
            split_count = 0
            n = split(sa_value, parts, ";")
            for (j = 1; j <= n; j++) {
                if (parts[j] != "") {
                    split_count++
                }
            }
            
            # Calculate threshold based on read length
            read_length = length($10)
            threshold = max_splits_base + (max_splits_kb * read_length / 1000.0)
            
            # Filter if too many splits
            if (split_count > threshold) {
                keep = 0
            }
            break
        }
    }
    
    # Output read if it should be kept
    if (keep) {
        print
    }
}
EOF
)
    
    # Export variables for AWK
    export MAX_SPLITS_BASE="$max_splits_base"
    export MAX_SPLITS_KB="$max_splits_kb"
    
    # Perform filtering using samtools and awk pipeline
    [[ "$quiet" == false ]] && echo "Filtering reads based on split alignment count..."
    
    # Create output directory if it doesn't exist
    local output_dir
    output_dir=$(dirname "$output_bam")
    if [[ ! -d "$output_dir" ]]; then
        echo "Creating output directory: $output_dir"
        mkdir -p "$output_dir"
    fi
    
    # Use a more efficient pipeline with proper error handling
    if ! samtools "${samtools_args[@]}" "$input_bam" | \
         awk "$awk_script" | \
         samtools view -b -@ "$threads" -o "$output_bam" -; then
        echo "Error: Long-read filtering pipeline failed" >&2
        exit 1
    fi
    
    # Create index
    [[ "$quiet" == false ]] && echo "Creating index..."
    if ! samtools index -@ "$threads" "$output_bam"; then
        echo "Error: Failed to create index" >&2
        exit 1
    fi
    
    # Calculate statistics
    if [[ "$quiet" == false ]]; then
        echo "Calculating statistics..."
        local kept_reads
        kept_reads=$(samtools view -c -@ "$threads" "$output_bam")
        local filtered_reads
        filtered_reads=$((total_reads - kept_reads))

        local percentage="0.0"
        if (( total_reads > 0 )); then
            percentage=$(awk "BEGIN {printf \"%.1f\", $kept_reads/$total_reads*100}")
        fi
        
        echo "Results:"
        if [[ "$primary_only" == true ]]; then
            echo "  Total primary reads: $total_reads"
        else
            echo "  Total reads: $total_reads"
        fi
        echo "  Kept reads: $kept_reads ($percentage%)"
        echo "  Filtered reads: $filtered_reads"
        echo "Long-read filtering completed successfully!"
    fi
}

# Parse command line arguments
if [[ $# -lt 1 ]]; then
    echo "Error: Mode is required" >&2
    usage >&2
    exit 2
fi

case "$1" in
    -h|--help)
        usage
        exit 0
        ;;
    -v|--version)
        echo "$PROGRAM_NAME $VERSION"
        exit 0
        ;;
esac

# Get mode
MODE="$1"
shift

# Validate mode
if [[ "$MODE" != "short-read" && "$MODE" != "long-read" ]]; then
    echo "Error: Invalid mode '$MODE'. Must be 'short-read' or 'long-read'" >&2
    usage >&2
    exit 2
fi

# Parse remaining arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            option_value_required "$1" "${2-}"
            INPUT_BAM="$2"
            shift 2
            ;;
        -o|--output)
            option_value_required "$1" "${2-}"
            OUTPUT_BAM="$2"
            shift 2
            ;;
        -t|--threads)
            option_value_required "$1" "${2-}"
            THREADS="$2"
            shift 2
            ;;
        --max-splits-base)
            option_value_required "$1" "${2-}"
            if [[ "$MODE" != "long-read" ]]; then
                echo "Warning: --max-splits-base is only used in long-read mode" >&2
            fi
            MAX_SPLITS_BASE="$2"
            shift 2
            ;;
        --max-splits-kb)
            option_value_required "$1" "${2-}"
            if [[ "$MODE" != "long-read" ]]; then
                echo "Warning: --max-splits-kb is only used in long-read mode" >&2
            fi
            MAX_SPLITS_KB="$2"
            shift 2
            ;;
        --include-all-reads)
            if [[ "$MODE" != "long-read" ]]; then
                echo "Warning: --include-all-reads is only used in long-read mode" >&2
            else
                PRIMARY_ONLY=false
            fi
            shift
            ;;
        -q|--quiet)
            QUIET=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -v|--version)
            echo "$PROGRAM_NAME $VERSION"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

# Validate required parameters
if [[ -z "$INPUT_BAM" || -z "$OUTPUT_BAM" ]]; then
    echo "Error: Input and output BAM files are required" >&2
    usage >&2
    exit 2
fi

if [[ ! "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: Threads must be a positive integer: $THREADS" >&2
    exit 1
fi

numeric_pattern='^([0-9]+([.][0-9]*)?|[.][0-9]+)$'
if [[ ! "$MAX_SPLITS_BASE" =~ $numeric_pattern ]]; then
    echo "Error: --max-splits-base must be a non-negative number: $MAX_SPLITS_BASE" >&2
    exit 1
fi

if [[ ! "$MAX_SPLITS_KB" =~ $numeric_pattern ]]; then
    echo "Error: --max-splits-kb must be a non-negative number: $MAX_SPLITS_KB" >&2
    exit 1
fi

if [[ ! -f "$INPUT_BAM" ]]; then
    echo "Error: Input BAM file does not exist: $INPUT_BAM" >&2
    exit 1
fi

# Check for required tools
for tool in samtools awk; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Error: $tool is not installed or not in PATH" >&2
        exit 1
    fi
done

# Display initial information
[[ "$QUIET" == false ]] && echo "BAM Filter Tool"
[[ "$QUIET" == false ]] && echo "Mode: $MODE"
[[ "$QUIET" == false ]] && echo "Input: $INPUT_BAM"
[[ "$QUIET" == false ]] && echo "Output: $OUTPUT_BAM"
[[ "$QUIET" == false ]] && echo "Threads: $THREADS"

# Show primary-only status only for long-read mode
if [[ "$MODE" == "long-read" ]]; then
    [[ "$QUIET" == false ]] && echo "Primary only: $PRIMARY_ONLY"
fi

# Execute appropriate filtering mode
case "$MODE" in
    short-read)
        filter_short_read_mode "$INPUT_BAM" "$OUTPUT_BAM" "$THREADS" "$QUIET"
        ;;
    long-read)
        filter_long_read_mode "$INPUT_BAM" "$OUTPUT_BAM" "$THREADS" "$PRIMARY_ONLY" "$QUIET" "$MAX_SPLITS_BASE" "$MAX_SPLITS_KB"
        ;;
esac

[[ "$QUIET" == false ]] && echo "All operations completed successfully!"
