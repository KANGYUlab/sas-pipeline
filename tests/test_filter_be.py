import contextlib
import importlib.util
import io
from pathlib import Path
import random
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "filter_be.py"
SPEC = importlib.util.spec_from_file_location("filter_be", str(SCRIPT_PATH))
filter_be = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(filter_be)


def reference_filter(df, chr_lengths, exclude_mode=False, threshold=1000.0):
    """The pre-vectorization row loop, retained here as a regression oracle."""
    filtered_indices = []
    high_depth_skipped = 0

    for idx, row in df.iterrows():
        try:
            chr_name = str(row[0])
            start = int(row[1])
            end = int(row[2])
            try:
                platform1_depth = float(row[4])
                platform2_depth = float(row[6])
            except (ValueError, TypeError):
                platform1_depth = float('nan')
                platform2_depth = float('nan')

            if ((not pd.isna(platform1_depth) and
                 platform1_depth > threshold) or
                    (not pd.isna(platform2_depth) and
                     platform2_depth > threshold)):
                high_depth_skipped += 1
                continue

            if filter_be.is_in_terminal_region(
                    chr_name, start, end, chr_lengths):
                should_keep = not filter_be.check_simple_conditions(row)
            else:
                should_keep = not filter_be.check_complex_conditions(row)

            if (exclude_mode and not should_keep) or (
                    not exclude_mode and should_keep):
                filtered_indices.append(idx)
        except (ValueError, IndexError):
            continue

    return df.iloc[filtered_indices], high_depth_skipped


def assert_matches_reference(test_case, df, chr_lengths):
    for exclude_mode in (False, True):
        expected, expected_high_depth = reference_filter(
            df, chr_lengths, exclude_mode
        )
        with contextlib.redirect_stdout(io.StringIO()):
            observed, observed_high_depth = filter_be.filter_dataframe(
                df, chr_lengths, exclude_mode
            )
        pd.testing.assert_frame_equal(observed, expected)
        test_case.assertEqual(observed_high_depth, expected_high_depth)


class BoundaryRegressionTests(unittest.TestCase):
    def test_thresholds_terminal_boundaries_and_high_depth(self):
        rows = [
            # Terminal boundary is inclusive; ratio 0.1 is inclusive.
            ["chr1", 10_000, 10_050, 1, 10, 0, 10, 0.0, 0.0, 0.0],
            ["chr1", 10_001, 10_051, 4, 10, 3, 10, 0.9, 0.0, 0.6],
            ["chr1", 20_000, 90_000, 1, 0, 0, 0, 0.9, 0.0, 0.6],
            # Strict count thresholds: 4 and 3 do not meet condition C.
            ["chr1", 20_000, 20_050, 4, 10, 3, 10, 0.89, 0.89, 0.59],
            ["chr1", 20_000, 20_050, 4.01, 10, 3, 10, 0.89, 0.89, 0.59],
            # col10 >= 0.9 makes condition A true on its own.
            ["chr1", 20_000, 20_050, 1, 10, 0, 10, 0.0, 0.0, 0.9],
            # Unknown chromosomes use the non-terminal complex condition.
            ["unknown", 0, 50, 5, 10, 0, 10, 0.0, 0.0, 0.0],
            # Exactly 1000 is retained; only values greater than it are skipped.
            ["chr1", 20_000, 20_050, 0, 1000, 0, 10, 0.0, 0.0, 0.0],
            ["chr1", 20_000, 20_050, 0, 1000.01, 0, 10, 0.0, 0.0, 0.0],
            ["chr1", 20_000, 20_050, 0, 10, 0, 1001, 0.0, 0.0, 0.0],
            # NaN comparisons remain false.
            ["chr1", 20_000, 20_050, np.nan, np.nan, np.nan, np.nan,
             np.nan, np.nan, np.nan],
        ]
        assert_matches_reference(self, pd.DataFrame(rows), {"chr1": 100_000})

    def test_malformed_coordinate_rows_are_skipped_in_input_order(self):
        rows = [
            ["chr1", 20_000, 20_050, 0, 10, 0, 10, 0, 0, 0],
            ["chr1", "bad", 20_050, 0, 10, 0, 10, 0, 0, 0],
            ["chr1", "1.0", 50, 0, 10, 0, 10, 0, 0, 0],
            ["chr1", 30_000, "bad", 0, 10, 0, 10, 0, 0, 0],
        ]
        assert_matches_reference(self, pd.DataFrame(rows), {"chr1": 100_000})

    def test_mixed_condition_columns_use_original_scalar_behavior(self):
        # A high-depth row must short-circuit before its malformed condition.
        high_depth_bad = pd.DataFrame([
            ["chr1", 20_000, 20_050, "bad", 1001, 0, 10, 0, 0, 0],
            ["chr1", 20_000, 20_050, 0, "10", 0, "10", 0, 0, 0],
        ])
        assert_matches_reference(self, high_depth_bad, {"chr1": 100_000})

        malformed_active = pd.DataFrame([
            ["chr1", 20_000, 20_050, "bad", 10, 0, 10, 0, 0, 0],
        ])
        with self.assertRaises(TypeError):
            reference_filter(malformed_active, {"chr1": 100_000})
        with self.assertRaises(TypeError):
            filter_be.filter_dataframe(malformed_active, {"chr1": 100_000})


class RandomizedRegressionTests(unittest.TestCase):
    def test_vectorized_results_match_scalar_reference(self):
        rng = random.Random(20260713)
        special = [0.0, 0.09, 0.1, 0.59, 0.6, 0.89, 0.9,
                   3.0, 4.0, 10.0, 1000.0, 1000.01, np.nan]
        rows = []
        for _ in range(2000):
            chrom = rng.choice(["chr1", "chr2", "unknown"])
            start = rng.choice([0, 9_999, 10_000, 10_001,
                                rng.randint(10_002, 89_998)])
            end = rng.choice([start + 50, 89_999, 90_000, 90_001])
            values = [rng.choice(special) for _column in range(7)]
            rows.append([chrom, start, end] + values)

        df = pd.DataFrame(rows)
        assert_matches_reference(
            self, df, {"chr1": 100_000, "chr2": 120_000}
        )


class FileOutputRegressionTests(unittest.TestCase):
    def test_filter_data_preserves_row_order_and_csv_format(self):
        rows = [
            ["chr1", 20_000, 20_050, 0, 10, 0, 10, 0, 0, 0, "first"],
            ["chr1", 20_100, 20_150, 5, 10, 0, 10, 0, 0, 0, "excluded"],
            ["chr1", 20_200, 20_250, 0, 10, 0, 10, 0, 0, 0, "last"],
        ]
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            data_path = temp_path / "input.tsv"
            fai_path = temp_path / "reference.fai"
            output_path = temp_path / "output.tsv"
            data_path.write_text(
                "header\n" + "".join(
                    "\t".join(str(value) for value in row) + "\n"
                    for row in rows
                )
            )
            fai_path.write_text("chr1\t100000\t0\t0\t0\n")

            with contextlib.redirect_stdout(io.StringIO()):
                filter_be.filter_data(
                    str(data_path), str(fai_path), str(output_path)
                )

            observed = output_path.read_text()
            self.assertEqual(
                observed,
                "chr1\t20000\t20050\t0\t10\t0\t10\t0\t0\t0\tfirst\n"
                "chr1\t20200\t20250\t0\t10\t0\t10\t0\t0\t0\tlast\n",
            )

    def test_one_row_chunks_match_single_chunk_output(self):
        rows = [
            ["chr1", 20_000, 20_050, 0, 10, 0, 10, 0, 0, 0],
            ["chr1", 20_100, 20_150, 5, 10, 0, 10, 0, 0, 0],
            ["chr1", 20_200, 20_250, 0, 10, 0, 10, 0, 0, 0],
        ]
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            data_path = temp_path / "input.tsv"
            fai_path = temp_path / "reference.fai"
            chunked_path = temp_path / "chunked.tsv"
            single_path = temp_path / "single.tsv"
            data_path.write_text(
                "header\n" + "".join(
                    "\t".join(str(value) for value in row) + "\n"
                    for row in rows
                ),
                encoding="utf-8",
            )
            fai_path.write_text("chr1\t100000\t0\t0\t0\n", encoding="utf-8")

            filter_be.filter_data(
                str(data_path), str(fai_path), str(chunked_path), chunk_size=1
            )
            filter_be.filter_data(
                str(data_path), str(fai_path), str(single_path), chunk_size=100
            )
            self.assertEqual(chunked_path.read_bytes(), single_path.read_bytes())

    def test_header_only_input_produces_empty_output(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            data_path = temp_path / "input.tsv"
            fai_path = temp_path / "reference.fai"
            output_path = temp_path / "output.tsv"
            data_path.write_text("header\n", encoding="utf-8")
            fai_path.write_text("chr1\t100000\t0\t0\t0\n", encoding="utf-8")

            filter_be.filter_data(
                str(data_path), str(fai_path), str(output_path), chunk_size=1
            )
            self.assertEqual(output_path.read_text(encoding="utf-8"), "")

    def test_late_chunk_failure_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            data_path = temp_path / "input.tsv"
            fai_path = temp_path / "reference.fai"
            output_path = temp_path / "output.tsv"
            row = "chr1\t20000\t20050\t0\t10\t0\t10\t0\t0\t0\n"
            data_path.write_text("header\n" + row + row, encoding="utf-8")
            fai_path.write_text("chr1\t100000\t0\t0\t0\n", encoding="utf-8")
            output_path.write_text("old result\n", encoding="utf-8")
            real_filter = filter_be.filter_dataframe
            calls = 0

            def fail_second_chunk(*args, **kwargs):
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise RuntimeError("injected late failure")
                return real_filter(*args, **kwargs)

            with mock.patch.object(
                filter_be, "filter_dataframe", side_effect=fail_second_chunk
            ):
                with self.assertRaisesRegex(RuntimeError, "late failure"):
                    filter_be.filter_data(
                        str(data_path), str(fai_path), str(output_path), chunk_size=1
                    )

            self.assertEqual(
                output_path.read_text(encoding="utf-8"), "old result\n"
            )


if __name__ == "__main__":
    unittest.main()
