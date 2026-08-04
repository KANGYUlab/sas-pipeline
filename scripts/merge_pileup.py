#!/usr/bin/env python3
"""
BED File Merger and Base Composition Calculator

Merge files based on bed coordinates and calculate base composition ratios
from reference genome sequences.

"""

import sys
import argparse
import os
import sqlite3
import tempfile
from pathlib import Path
try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    pysam = None
    PYSAM_AVAILABLE = False

REFERENCE_CACHE_BP = 1_000_000


class BedFileMerger:
    """Merge files based on bed coordinates and calculate base composition."""
    
    def __init__(self):
        self.ref_fasta = None  # pysam FastaFile handle
        self._ref_cache_chrom = None
        self._ref_cache_start = 0
        self._ref_cache_sequence = ""
        
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
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # A disk-backed coordinate join keeps memory bounded for dense IR runs
        # with millions of candidate windows. Evidence rows retain the former
        # last-record-wins behavior through INSERT OR REPLACE.
        with tempfile.TemporaryDirectory(
            prefix=".merge-pileup-", dir=output_path.parent
        ) as temp_dir:
            database_path = Path(temp_dir) / "evidence.sqlite3"
            staged_output = Path(temp_dir) / "merged_results.txt"
            connection = sqlite3.connect(database_path)
            try:
                connection.execute("PRAGMA journal_mode=OFF")
                connection.execute("PRAGMA synchronous=OFF")
                connection.execute("PRAGMA temp_store=FILE")
                connection.execute(
                    "CREATE TABLE regions ("
                    "input_order INTEGER, chrom TEXT, start INTEGER, end INTEGER, "
                    "PRIMARY KEY (chrom, start, end))"
                )
                for table_name in ("platform1", "platform2"):
                    connection.execute(
                        f"CREATE TABLE {table_name} ("
                        "chrom TEXT, start INTEGER, end INTEGER, "
                        "value1 TEXT, value2 TEXT, "
                        "PRIMARY KEY (chrom, start, end))"
                    )

                region_count = self._load_bed_into_database(
                    connection, bed_file
                )
                platform1_count = self._load_evidence_into_database(
                    connection, file1, "platform1"
                )
                platform2_count = self._load_evidence_into_database(
                    connection, file2, "platform2"
                )
                connection.commit()
                print(
                    f"Loaded {region_count} candidates, {platform1_count} HiFi records, "
                    f"and {platform2_count} ONT records"
                )

                if not PYSAM_AVAILABLE:
                    raise RuntimeError(
                        "pysam is required for bounded-memory reference sequence access"
                    )
                self.ref_fasta = pysam.FastaFile(ref_file)
                self._write_database_join(connection, staged_output)
            finally:
                if self.ref_fasta:
                    self.ref_fasta.close()
                    self.ref_fasta = None
                self._ref_cache_chrom = None
                self._ref_cache_sequence = ""
                connection.close()

            os.replace(staged_output, output_path)
        print(f"Merged file written successfully: {output_path}")

    @staticmethod
    def _flush_rows(connection, statement, rows):
        if rows:
            connection.executemany(statement, rows)
            rows.clear()

    def _load_bed_into_database(self, connection, filename):
        statement = "INSERT OR IGNORE INTO regions VALUES (?, ?, ?, ?)"
        rows = []
        valid_count = 0
        with open(filename, "r", encoding="utf-8") as input_handle:
            for line_number, line in enumerate(input_handle, 1):
                fields = line.strip().split("\t")
                if not line.strip() or line.startswith("#"):
                    continue
                if len(fields) < 3:
                    raise ValueError(f"Malformed BED row {line_number}: {line.strip()}")
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                if start < 0 or end <= start:
                    raise ValueError(f"Invalid BED row {line_number}: {line.strip()}")
                rows.append((valid_count, chrom, start, end))
                valid_count += 1
                if len(rows) >= 10_000:
                    self._flush_rows(connection, statement, rows)
        self._flush_rows(connection, statement, rows)
        return connection.execute("SELECT COUNT(*) FROM regions").fetchone()[0]

    def _load_evidence_into_database(self, connection, filename, table_name):
        statement = f"INSERT OR REPLACE INTO {table_name} VALUES (?, ?, ?, ?, ?)"
        rows = []
        processed_count = 0
        with open(filename, "r", encoding="utf-8") as input_handle:
            for line_number, line in enumerate(input_handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 5:
                    raise ValueError(
                        f"Malformed evidence row {line_number} in {filename}"
                    )
                try:
                    chrom = fields[-5]
                    start = int(fields[-4])
                    end = int(fields[-3])
                except (ValueError, IndexError) as error:
                    raise ValueError(
                        f"Invalid evidence coordinates on row {line_number} in {filename}"
                    ) from error
                rows.append((chrom, start, end, fields[-2], fields[-1]))
                processed_count += 1
                if len(rows) >= 10_000:
                    self._flush_rows(connection, statement, rows)
        self._flush_rows(connection, statement, rows)
        return processed_count

    def _write_database_join(self, connection, output_file):
        query = """
            SELECT r.chrom, r.start, r.end,
                   COALESCE(p1.value1, '-1'), COALESCE(p1.value2, '-1'),
                   COALESCE(p2.value1, '-1'), COALESCE(p2.value2, '-1')
            FROM regions AS r
            LEFT JOIN platform1 AS p1
              ON p1.chrom=r.chrom AND p1.start=r.start AND p1.end=r.end
            LEFT JOIN platform2 AS p2
              ON p2.chrom=r.chrom AND p2.start=r.start AND p2.end=r.end
            ORDER BY r.chrom, r.start, r.input_order
        """
        processed = 0
        with open(output_file, "w", encoding="utf-8") as output_handle:
            output_handle.write(
                "Chr\tStart\tEnd\tFile1_Col1\tFile1_Col2\t"
                "File2_Col1\tFile2_Col2\tAG_ratio\tCT_ratio\tGC_ratio\n"
            )
            for chrom, start, end, p1v1, p1v2, p2v1, p2v2 in connection.execute(query):
                ag_ratio, ct_ratio, gc_ratio = self._calculate_base_ratios(
                    chrom, start, end
                )
                output_handle.write(
                    f"{chrom}\t{start}\t{end}\t{p1v1}\t{p1v2}\t"
                    f"{p2v1}\t{p2v2}\t{ag_ratio:.4f}\t{ct_ratio:.4f}\t"
                    f"{gc_ratio:.4f}\n"
                )
                processed += 1
                if processed % 100_000 == 0:
                    print(f"Merged {processed} candidate windows")
        
    def _calculate_base_ratios(self, chr_name: str, start: int, end: int) -> tuple:
        """Calculate base composition ratios for a genomic region."""
        if self.ref_fasta is None:
            raise RuntimeError("Reference FASTA is not open")
        if start < 0 or end <= start:
            raise ValueError(f"Invalid reference interval: {chr_name}:{start}-{end}")
        cache_end = self._ref_cache_start + len(self._ref_cache_sequence)
        if (
            self._ref_cache_chrom == chr_name
            and start >= self._ref_cache_start
            and end <= cache_end
        ):
            offset = start - self._ref_cache_start
            subseq = self._ref_cache_sequence[offset:offset + (end - start)]
        else:
            block_start = (start // REFERENCE_CACHE_BP) * REFERENCE_CACHE_BP
            block_end = max(end, block_start + REFERENCE_CACHE_BP)
            try:
                get_length = getattr(self.ref_fasta, "get_reference_length", None)
                if get_length is not None:
                    reference_length = get_length(chr_name)
                    if end > reference_length:
                        raise ValueError("interval exceeds reference length")
                    block_end = min(block_end, reference_length)
                cached_sequence = self.ref_fasta.fetch(
                    chr_name, block_start, block_end
                ).upper()
            except (KeyError, ValueError) as error:
                raise ValueError(
                    f"Candidate interval is absent from the reference FASTA: "
                    f"{chr_name}:{start}-{end}"
                ) from error
            self._ref_cache_chrom = chr_name
            self._ref_cache_start = block_start
            self._ref_cache_sequence = cached_sequence
            offset = start - block_start
            subseq = cached_sequence[offset:offset + (end - start)]

        if len(subseq) != end - start:
            raise ValueError(
                f"Reference interval extends beyond sequence bounds: "
                f"{chr_name}:{start}-{end}"
            )

        total = len(subseq)
        a_count = subseq.count('A')
        t_count = subseq.count('T')
        c_count = subseq.count('C')
        g_count = subseq.count('G')

        ag_ratio = (a_count + g_count) / total
        ct_ratio = (c_count + t_count) / total
        gc_ratio = (g_count + c_count) / total
        return (ag_ratio, ct_ratio, gc_ratio)

def create_parser() -> argparse.ArgumentParser:
    """Create command line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Merge HiFi and ONT pileup evidence by BED coordinates and calculate "
            "reference base-composition ratios."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 %(prog)s regions.bed hifi.txt ont.txt reference.fa merged.txt
  
Output format:
  CHROM START END FILE1_COL1 FILE1_COL2 FILE2_COL1 FILE2_COL2 AG CT GC
        """,
    )

    inputs = parser.add_argument_group("required inputs")
    inputs.add_argument("bed_file", metavar="BED", help="BED coordinate file")
    inputs.add_argument("file1", metavar="HIFI_DATA", help="HiFi pileup result file")
    inputs.add_argument("file2", metavar="ONT_DATA", help="ONT pileup result file")
    inputs.add_argument("ref_file", metavar="FASTA", help="Reference genome FASTA")
    inputs.add_argument("output_file", metavar="OUTPUT", help="Merged output file")
    
    parser.add_argument(
        '--version',
        action='version',
        version="%(prog)s 2.0.0"
    )
    
    return parser


def main():
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    try:
        # Validate input files exist
        for file_path in [args.bed_file, args.file1, args.file2, args.ref_file]:
            if not Path(file_path).is_file():
                parser.error(f"input file not found: {file_path}")
        
        merger = BedFileMerger()
        merger.merge_files(args.bed_file, args.file1, args.file2, args.ref_file, args.output_file)
        print("File merging and base composition calculation completed successfully!")
        
    except KeyboardInterrupt:
        print("\nOperation interrupted by user", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main() 
