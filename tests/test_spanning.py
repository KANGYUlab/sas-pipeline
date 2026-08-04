import contextlib
import importlib.util
import io
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


# spanning.py imports pysam at module load time.  A tiny stub keeps these unit
# tests independent of native pysam/samtools installations.
PYSAM_STUB = types.ModuleType("pysam")
PYSAM_STUB.AlignmentFile = None
sys.modules["pysam"] = PYSAM_STUB

MODULE_PATH = Path(__file__).parents[1] / "scripts" / "spanning.py"
SPEC = importlib.util.spec_from_file_location("spanning_under_test", MODULE_PATH)
spanning = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(spanning)


class FakeRead:
    def __init__(self, start, end, unmapped=False):
        self.reference_start = start
        self.reference_end = end
        self.is_unmapped = unmapped


class FakeBam:
    def __init__(self, reads_by_chromosome):
        self.reads_by_chromosome = reads_by_chromosome
        self.fetch_calls = []
        self.closed = False

    def fetch(self, *args):
        self.fetch_calls.append(args)
        if len(args) != 1:
            raise AssertionError("the sweep must fetch each chromosome only once")
        return iter(self.reads_by_chromosome.get(args[0], []))

    def close(self):
        self.closed = True


class SlidingWindowTests(unittest.TestCase):
    def test_consecutive_short_contigs_each_get_one_nonnegative_window(self):
        chromosomes = {
            "short_a": 1,
            "short_b": 10_000,
            "short_c": 27_999,
            "boundary": 28_000,
        }

        windows = spanning.generate_sliding_windows(chromosomes)

        self.assertEqual(
            windows,
            [
                ("short_a", 0, 1),
                ("short_b", 0, 10_000),
                ("short_c", 0, 27_999),
                ("boundary", 0, 28_000),
            ],
        )
        self.assertTrue(all(start >= 0 and start < end for _, start, end in windows))

    def test_normal_and_terminal_windows_are_preserved(self):
        self.assertEqual(
            spanning.generate_sliding_windows({"chr1": 31_000}),
            [("chr1", 0, 30_000), ("chr1", 2_000, 31_000)],
        )

    def test_invalid_window_parameters_are_rejected(self):
        with self.assertRaises(ValueError):
            spanning.generate_sliding_windows({"chr1": 100}, step_size=0)


class SweepTests(unittest.TestCase):
    def test_read_starting_before_first_window_can_span_it(self):
        windows = [("chr1", 2_000, 32_000)]
        bam = FakeBam({"chr1": [FakeRead(0, 32_000)]})

        observed = spanning.find_non_spanning_windows(bam, "chr1", windows)

        self.assertEqual(observed, [])
        self.assertEqual(bam.fetch_calls, [("chr1",)])

    def test_mapped_read_with_none_coordinates_is_ignored(self):
        windows = [("chr1", 2_000, 32_000)]
        bam = FakeBam(
            {
                "chr1": [
                    FakeRead(None, 40_000),
                    FakeRead(0, None),
                ]
            }
        )

        observed = spanning.find_non_spanning_windows(bam, "chr1", windows)

        self.assertEqual(observed, windows)
        self.assertEqual(bam.fetch_calls, [("chr1",)])

    def test_single_scan_matches_exact_original_predicate(self):
        windows = [
            ("chr1", 0, 30),
            ("chr1", 2, 32),
            ("chr1", 4, 34),
            ("chr1", 6, 50),
        ]
        reads = [
            FakeRead(None, None, unmapped=True),
            FakeRead(0, 31),
            FakeRead(2, 32),
            FakeRead(3, 40),
            FakeRead(7, 100),
        ]
        bam = FakeBam({"chr1": reads})

        observed = spanning.find_non_spanning_windows(bam, "chr1", windows)
        expected = [
            window
            for window in windows
            if not any(
                not read.is_unmapped
                and read.reference_start is not None
                and read.reference_end is not None
                and read.reference_start <= window[1]
                and read.reference_end >= window[2]
                for read in reads
            )
        ]

        self.assertEqual(observed, expected)
        self.assertEqual(observed, [("chr1", 6, 50)])
        self.assertEqual(bam.fetch_calls, [("chr1",)])

    def test_windows_must_belong_to_the_requested_chromosome(self):
        bam = FakeBam({"chr1": []})
        with self.assertRaises(ValueError):
            spanning.find_non_spanning_windows(bam, "chr1", [("chr2", 0, 10)])


class OutputTests(unittest.TestCase):
    def test_zero_non_spanning_windows_creates_both_empty_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fai = temp_path / "reference.fai"
            fai.write_text("chr1\t30000\t0\t0\t0\n")
            output_prefix = temp_path / "sample"

            created_bams = []

            def alignment_file(_path, _mode):
                # One read spans the sole 30 kb window.
                bam = FakeBam({"chr1": [FakeRead(0, 30_000)]})
                created_bams.append(bam)
                return bam

            argv = [
                "spanning.py",
                "sample.bam",
                str(fai),
                "--output",
                str(output_prefix),
                "--processes",
                "1",
            ]
            with (
                mock.patch.object(spanning.pysam, "AlignmentFile", side_effect=alignment_file),
                mock.patch.object(spanning.subprocess, "run") as run,
                mock.patch.object(sys, "argv", argv),
                contextlib.redirect_stderr(io.StringIO()),
            ):
                spanning.main()

            run.assert_not_called()
            merged = Path(f"{output_prefix}.non_spanning_regions.merged.bed")
            final_windows = Path(f"{output_prefix}.non_spanning.bed")
            self.assertTrue(merged.exists())
            self.assertTrue(final_windows.exists())
            self.assertEqual(merged.read_text(), "")
            self.assertEqual(final_windows.read_text(), "")
            self.assertEqual(created_bams[-1].fetch_calls, [("chr1",)])


class TerminalCorrectionTests(unittest.TestCase):
    def test_terminal_correction_and_2k_schema_are_unchanged(self):
        corrected = spanning.correct_terminal_regions(
            [
                ("chr1", 0, 40_000),
                ("chr1", 40_000, 60_000),
                ("chr1", 60_000, 100_000),
            ],
            {"chr1": 100_000},
        )
        self.assertEqual(
            corrected,
            [
                ("chr1", 0, 12_000),
                ("chr1", 40_000, 60_000),
                ("chr1", 88_000, 100_000),
            ],
        )
        windows = spanning.assign_corrected_regions_to_2k_windows(
            [("chr1", 4_000, 7_000)]
        )
        self.assertEqual(windows, [("chr1", 4_000, 6_000), ("chr1", 6_000, 7_000)])


if __name__ == "__main__":
    unittest.main()
