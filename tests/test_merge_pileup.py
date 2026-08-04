import importlib.util
from pathlib import Path
import tempfile
import types
import unittest
from unittest import mock


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "merge_pileup.py"
SPEC = importlib.util.spec_from_file_location("merge_pileup", str(SCRIPT_PATH))
merge_pileup = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(merge_pileup)


class FakeFasta:
    sequences = {"chr1": "ACGT" * 100, "chr2": "AAAA" * 100}
    fetch_calls = 0

    def __init__(self, _filename):
        pass

    def fetch(self, chrom, start, end):
        type(self).fetch_calls += 1
        if chrom not in self.sequences:
            raise KeyError(chrom)
        return self.sequences[chrom][start:end]

    def get_reference_length(self, chrom):
        if chrom not in self.sequences:
            raise KeyError(chrom)
        return len(self.sequences[chrom])

    def close(self):
        pass


class StreamingPileupMergeTests(unittest.TestCase):
    def test_disk_backed_join_preserves_evidence_and_defaults(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            bed = temp_path / "candidate.bed"
            hifi = temp_path / "hifi.txt"
            ont = temp_path / "ont.txt"
            reference = temp_path / "reference.fa"
            output = temp_path / "merged.txt"
            bed.write_text(
                "chr2\t0\t4\nchr1\t0\t4\nchr1\t4\t8\n",
                encoding="utf-8",
            )
            hifi.write_text(
                "chr1\t1\tA\t1\t.\tchr1\t0\t4\t3\t10\n",
                encoding="utf-8",
            )
            ont.write_text(
                "chr2\t1\tA\t1\t.\tchr2\t0\t4\t2\t8\n",
                encoding="utf-8",
            )
            reference.write_text(">placeholder\nA\n", encoding="utf-8")

            fake_pysam = types.SimpleNamespace(FastaFile=FakeFasta)
            FakeFasta.fetch_calls = 0
            with mock.patch.object(merge_pileup, "pysam", fake_pysam), mock.patch.object(
                merge_pileup, "PYSAM_AVAILABLE", True
            ):
                merge_pileup.BedFileMerger().merge_files(
                    str(bed), str(hifi), str(ont), str(reference), str(output)
                )

            rows = output.read_text(encoding="utf-8").splitlines()
            self.assertEqual([row.split("\t")[:3] for row in rows[1:]], [
                ["chr1", "0", "4"],
                ["chr1", "4", "8"],
                ["chr2", "0", "4"],
            ])
            self.assertEqual(rows[1].split("\t")[3:7], ["3", "10", "-1", "-1"])
            self.assertEqual(rows[2].split("\t")[3:7], ["-1", "-1", "-1", "-1"])
            self.assertEqual(rows[3].split("\t")[3:7], ["-1", "-1", "2", "8"])
            self.assertEqual(FakeFasta.fetch_calls, 2)
            self.assertEqual(list(temp_path.glob(".merge-pileup-*")), [])

    def test_reference_lookup_failure_does_not_replace_existing_output(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            bed = temp_path / "candidate.bed"
            hifi = temp_path / "hifi.txt"
            ont = temp_path / "ont.txt"
            reference = temp_path / "reference.fa"
            output = temp_path / "merged.txt"
            bed.write_text("unknown\t0\t4\n", encoding="utf-8")
            hifi.write_text("", encoding="utf-8")
            ont.write_text("", encoding="utf-8")
            reference.write_text(">placeholder\nA\n", encoding="utf-8")
            output.write_text("old result\n", encoding="utf-8")

            fake_pysam = types.SimpleNamespace(FastaFile=FakeFasta)
            with mock.patch.object(merge_pileup, "pysam", fake_pysam), mock.patch.object(
                merge_pileup, "PYSAM_AVAILABLE", True
            ):
                with self.assertRaises(ValueError):
                    merge_pileup.BedFileMerger().merge_files(
                        str(bed), str(hifi), str(ont), str(reference), str(output)
                    )

            self.assertEqual(output.read_text(encoding="utf-8"), "old result\n")


if __name__ == "__main__":
    unittest.main()
