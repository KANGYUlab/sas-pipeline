#!/usr/bin/env python3
"""
Data Filtering Script
Filter data rows based on BED positions and multiple conditions
EXCLUDES rows that meet certain criteria, keeps the rest
"""

import argparse
import os
import sys
import tempfile

try:
    import pandas as pd
except ModuleNotFoundError:
    pd = None

try:
    import numpy as np
except ModuleNotFoundError:
    np = None

CLI_VERSION = "2.0.0"
DEFAULT_CHUNK_SIZE = 100_000


class CLIHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Keep examples and explicitly documented defaults readable."""


def parse_separator(value):
    """Accept readable separator names without breaking literal separators."""
    aliases = {
        "tab": "\t",
        r"\t": "\t",
        "comma": ",",
        "space": " ",
        "pipe": "|",
    }
    return aliases.get(value.lower(), value)

def load_fai(fai_file):
    """
    Load chromosome length information from fai file
    fai format: chr_name, length, offset, linebases, linewidth
    """
    chr_lengths = {}
    with open(fai_file, 'r', encoding='utf-8') as input_handle:
        for line_number, line in enumerate(input_handle, 1):
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                raise ValueError(
                    f"Malformed FASTA index row {line_number}: {line.rstrip()}"
                )
            chr_name = parts[0]
            try:
                chr_length = int(parts[1])
            except ValueError as error:
                raise ValueError(
                    f"Invalid FASTA length on row {line_number}: {parts[1]}"
                ) from error
            if not chr_name or chr_length <= 0 or chr_name in chr_lengths:
                raise ValueError(
                    f"Invalid FASTA index row {line_number}: {line.rstrip()}"
                )
            chr_lengths[chr_name] = chr_length
    if not chr_lengths:
        raise ValueError(f"FASTA index contains no sequences: {fai_file}")
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


def _prepare_coordinates(df):
    """Convert BED coordinates while retaining the old row-error behavior."""
    row_count = len(df)
    valid = np.ones(row_count, dtype=bool)

    # read_csv normally infers BED coordinates as numeric columns.  Keep this
    # common path vectorized; mixed/object columns use exact int() conversion so
    # values such as the string "1.0" are rejected just as they were before.
    if (pd.api.types.is_numeric_dtype(df[1].dtype) and
            pd.api.types.is_numeric_dtype(df[2].dtype)):
        start_numeric = pd.to_numeric(df[1], errors='coerce')
        end_numeric = pd.to_numeric(df[2], errors='coerce')
        invalid = start_numeric.isna() | end_numeric.isna()

        for position in np.flatnonzero(invalid.to_numpy()):
            row = df.iloc[position]
            try:
                int(row[1])
                int(row[2])
            except (ValueError, IndexError) as e:
                print(
                    f"Warning: Error processing row {position+2} "
                    f"(after skipping header): {e}",
                    file=sys.stderr,
                )
                valid[position] = False

        # int(infinity) raises OverflowError in the original code.  Trigger the
        # same exception rather than silently turning it into an invalid row.
        finite = np.isfinite(start_numeric.fillna(0).to_numpy(dtype=float)) & \
                 np.isfinite(end_numeric.fillna(0).to_numpy(dtype=float))
        for position in np.flatnonzero(~finite):
            row = df.iloc[position]
            int(row[1])
            int(row[2])

        starts = np.trunc(start_numeric.fillna(0).to_numpy(dtype=float))
        ends = np.trunc(end_numeric.fillna(0).to_numpy(dtype=float))
        return starts, ends, valid

    starts = np.empty(row_count, dtype=object)
    ends = np.empty(row_count, dtype=object)
    starts[:] = 0
    ends[:] = 0
    for position in range(row_count):
        row = df.iloc[position]
        try:
            starts[position] = int(row[1])
            ends[position] = int(row[2])
        except (ValueError, IndexError) as e:
            print(
                f"Warning: Error processing row {position+2} "
                f"(after skipping header): {e}",
                file=sys.stderr,
            )
            valid[position] = False

    return starts, ends, valid


def _terminal_mask(df, starts, ends, valid, chr_lengths, terminal_size=10000):
    """Return a positional mask for regions touching either chromosome end."""
    chromosome_names = df[0].astype(str)
    chromosome_lengths = chromosome_names.map(chr_lengths)
    known_chromosome = chromosome_lengths.notna().to_numpy()
    length_values = chromosome_lengths.fillna(0).to_numpy(dtype=float)
    return (valid & known_chromosome &
            ((starts <= terminal_size) |
             (ends >= (length_values - terminal_size))))


def _numeric_condition_mask(df, terminal, active, high_depth_threshold):
    """Evaluate the filtering predicates with pandas vector operations."""
    col4, col5 = df[3], df[4]
    col6, col7 = df[5], df[6]
    col8, col9, col10 = df[7], df[8], df[9]

    high_depth = active & (
        ((col5 > high_depth_threshold) & col5.notna())
        .fillna(False).to_numpy(dtype=bool) |
        ((col7 > high_depth_threshold) & col7.notna())
        .fillna(False).to_numpy(dtype=bool)
    )

    # Replacing a zero denominator with NaN mirrors boolean short-circuiting
    # and avoids evaluating a division that the scalar implementation skipped.
    platform1_ratio = col4 / col5.where(col5 != 0)
    platform2_ratio = col6 / col7.where(col7 != 0)
    simple_condition = (
        ((col4 > 0) & (col5 != 0) & (platform1_ratio >= 0.1)) |
        ((col6 > 0) & (col7 != 0) & (platform2_ratio >= 0.1))
    ).fillna(False).to_numpy(dtype=bool)

    condition_a = (
        (((col8 >= 0.9) | (col9 >= 0.9)) & (col10 >= 0.6)) |
        (col10 >= 0.9)
    ).fillna(False).to_numpy(dtype=bool)
    count_condition = ((col4 > 4) | (col6 > 3)).fillna(False).to_numpy(dtype=bool)
    complex_condition = ((condition_a & simple_condition) |
                         (~condition_a & count_condition))
    excluded_condition = ((terminal & simple_condition) |
                          (~terminal & complex_condition))
    return excluded_condition, high_depth


def _scalar_condition_mask(df, terminal, active, high_depth_threshold):
    """Compatibility path for mixed/non-numeric condition columns."""
    excluded_condition = np.zeros(len(df), dtype=bool)
    high_depth = np.zeros(len(df), dtype=bool)

    for position in np.flatnonzero(active):
        # iloc produces the same row-shaped Series as the former iterrows loop,
        # retaining its comparisons and exception behavior for malformed data.
        row = df.iloc[position]
        try:
            platform1_depth = float(row[4])
            platform2_depth = float(row[6])
        except (ValueError, TypeError):
            platform1_depth = float('nan')
            platform2_depth = float('nan')

        if ((not pd.isna(platform1_depth) and
             platform1_depth > high_depth_threshold) or
                (not pd.isna(platform2_depth) and
                 platform2_depth > high_depth_threshold)):
            high_depth[position] = True
            continue

        if terminal[position]:
            excluded_condition[position] = check_simple_conditions(row)
        else:
            excluded_condition[position] = check_complex_conditions(row)

    return excluded_condition, high_depth


def filter_dataframe(df, chr_lengths, exclude_mode=False,
                     high_depth_threshold=1000.0):
    """Filter an already-loaded DataFrame and return it with skip statistics."""
    starts, ends, valid = _prepare_coordinates(df)
    terminal = _terminal_mask(df, starts, ends, valid, chr_lengths)

    condition_columns_are_numeric = all(
        pd.api.types.is_numeric_dtype(df[column].dtype)
        for column in range(3, 10)
    )
    if condition_columns_are_numeric:
        excluded_condition, high_depth = _numeric_condition_mask(
            df, terminal, valid, high_depth_threshold
        )
    else:
        excluded_condition, high_depth = _scalar_condition_mask(
            df, terminal, valid, high_depth_threshold
        )

    selected_condition = excluded_condition if exclude_mode else ~excluded_condition
    output_mask = valid & ~high_depth & selected_condition
    filtered_df = df.iloc[np.flatnonzero(output_mask)]
    return filtered_df, int(high_depth.sum())

def filter_data(
    data_file,
    fai_file,
    output_file=None,
    exclude_mode=False,
    sep='\t',
    chunk_size=DEFAULT_CHUNK_SIZE,
):
    """
    Main data filtering function
    """
    HIGH_DEPTH_THRESHOLD = 1000.0

    # Load chromosome length information
    chr_lengths = load_fai(fai_file)
    print(
        f"Loaded length information for {len(chr_lengths)} chromosomes",
        file=sys.stderr,
    )
    
    if chunk_size <= 0:
        raise ValueError("chunk_size must be > 0")

    output_handle = sys.stdout
    temporary_output = None
    close_output = False
    if output_file:
        output_path = os.path.abspath(output_file)
        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        staged = tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            prefix=f".{os.path.basename(output_path)}.",
            suffix=".tmp",
            dir=output_dir or os.getcwd(),
            delete=False,
        )
        output_handle = staged
        temporary_output = staged.name
        close_output = True

    total_rows = 0
    filtered_rows = 0
    high_depth_skipped = 0
    chunks = None
    try:
        try:
            chunks = pd.read_csv(
                data_file,
                sep=sep,
                header=None,
                skiprows=1,
                chunksize=chunk_size,
            )
            for chunk in chunks:
                if len(chunk.columns) < 10:
                    raise ValueError(
                        "Insufficient columns in data file, expected at least "
                        f"10, got {len(chunk.columns)}"
                    )
                total_rows += len(chunk)
                filtered_chunk, chunk_high_depth = filter_dataframe(
                    chunk,
                    chr_lengths,
                    exclude_mode,
                    HIGH_DEPTH_THRESHOLD,
                )
                high_depth_skipped += chunk_high_depth
                filtered_rows += len(filtered_chunk)
                filtered_chunk.to_csv(
                    output_handle,
                    sep=sep,
                    header=False,
                    index=False,
                )
        except pd.errors.EmptyDataError:
            # A header-only merged table is a valid no-candidate result.
            pass
        finally:
            if chunks is not None:
                chunks.close()

        if close_output:
            output_handle.flush()
            output_handle.close()
            close_output = False
            os.replace(temporary_output, output_path)
            temporary_output = None

        print(f"Read {total_rows} rows of data", file=sys.stderr)
        print(f"Filtered result contains {filtered_rows} rows", file=sys.stderr)
        if high_depth_skipped > 0:
            print(
                f"Skipped {high_depth_skipped} rows due to depth > "
                f"{HIGH_DEPTH_THRESHOLD}",
                file=sys.stderr,
            )
        if output_file:
            print(f"Results saved to {output_file}", file=sys.stderr)
    except Exception:
        if close_output:
            output_handle.close()
            close_output = False
        raise
    finally:
        if close_output:
            output_handle.close()
        if temporary_output and os.path.exists(temporary_output):
            os.remove(temporary_output)

def main():
    parser = argparse.ArgumentParser(
        description="Filter merged BE evidence using SAS thresholds and terminal rules.",
        formatter_class=CLIHelpFormatter,
        epilog="""
Examples:
  python3 %(prog)s merged.tsv reference.fa.fai -o filtered.tsv
  python3 %(prog)s merged.csv reference.fa.fai -s comma -e
        """,
    )
    inputs = parser.add_argument_group("required inputs")
    inputs.add_argument("data_file", metavar="DATA", help="Merged BE evidence table")
    inputs.add_argument("fai_file", metavar="FAI", help="Reference FASTA index")

    output = parser.add_argument_group("filter and output options")
    output.add_argument(
        "-o", "--output", metavar="FILE",
        help="Output file; writes to standard output when omitted",
    )
    output.add_argument(
        "--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE, metavar="N",
        help=f"Rows processed per memory-bounded chunk (default: {DEFAULT_CHUNK_SIZE})",
    )
    output.add_argument(
        "-e", "--exclude", action="store_true",
        help="Invert the filter and output rows that do not meet the criteria",
    )
    output.add_argument(
        "-s", "--sep", type=parse_separator, default="\t", metavar="SEP",
        help="Field separator: tab, comma, space, pipe, or a literal character (default: tab)",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {CLI_VERSION}"
    )
    
    args = parser.parse_args()
    missing = [name for name, module in (("numpy", np), ("pandas", pd)) if module is None]
    if missing:
        parser.error("missing required Python dependencies: " + ", ".join(missing))
    
    if args.chunk_size <= 0:
        parser.error("--chunk-size must be greater than 0")
    try:
        filter_data(
            args.data_file,
            args.fai_file,
            args.output,
            args.exclude,
            args.sep,
            args.chunk_size,
        )
    except Exception as error:
        parser.error(str(error))

if __name__ == '__main__':
    main() 
