import importlib.util
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "smallerror_long.py"
SPEC = importlib.util.spec_from_file_location("smallerror_long", str(SCRIPT_PATH))
smallerror_long = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(smallerror_long)


def region(chrom, start, end, index, output_path="unused"):
    return (chrom, start, end, index, output_path)


class RegionBatchingTests(unittest.TestCase):
    def test_batches_use_gap_span_chromosome_and_count_limits(self):
        regions = [
            region("chr1", 0, 50, 0),
            region("chr1", 10_050, 10_100, 1),  # exactly 10 kb after region 0
            region("chr1", 20_101, 20_151, 2),  # 10,001 bp gap: split
            region("chr2", 0, 50, 3),
        ]

        batches = smallerror_long.build_region_batches(
            regions,
            max_gap=10_000,
            max_span=1_000_000,
            max_regions=1_000,
        )

        self.assertEqual([[item[3] for item in batch] for batch in batches], [[0, 1], [2], [3]])

        span_limited = smallerror_long.build_region_batches(
            [region("chr1", 0, 50, 0), region("chr1", 999_999, 1_000_050, 1)],
            max_gap=2_000_000,
            max_span=1_000_000,
            max_regions=1_000,
        )
        self.assertEqual(len(span_limited), 2)

        count_limited = smallerror_long.build_region_batches(
            [
                region("chr1", 0, 50, 0),
                region("chr1", 100, 150, 1),
                region("chr1", 200, 250, 2),
            ],
            max_gap=10_000,
            max_span=1_000_000,
            max_regions=2,
        )
        self.assertEqual([len(batch) for batch in count_limited], [2, 1])

    def test_dispatches_overlaps_and_filters_gap_positions(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            raw_paths = {
                0: os.path.join(temp_dir, "region0.txt"),
                1: os.path.join(temp_dir, "region1.txt"),
            }
            for raw_path in raw_paths.values():
                Path(raw_path).write_text("")

            regions = [
                region("chr1", 0, 2, 0),
                region("chr1", 1, 3, 1),
            ]
            lines = [
                "chr1\t1\tA\t1\t.\tI\n",
                "chr1\t2\tA\t1\t.\tI\n",
                "chr1\t3\tA\t1\t.\tI\n",
                "chr1\t4\tA\t1\t.\tI\n",  # outside both regions
            ]

            smallerror_long.dispatch_pileup_lines(lines, regions, raw_paths, max_open_files=1)

            self.assertEqual(
                Path(raw_paths[0]).read_text(),
                "".join(lines[:2]),
            )
            self.assertEqual(
                Path(raw_paths[1]).read_text(),
                "".join(lines[1:3]),
            )

    def test_single_batch_spool_reconstructs_overlapping_regions(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            spool_path = os.path.join(temp_dir, "batch.txt")
            regions = [
                region("chr1", 0, 2, 0),
                region("chr1", 1, 3, 1),
            ]
            lines = [
                "chr1\t1\tA\t1\t.\tI\n",
                "chr1\t2\tA\t1\t.\tI\n",
                "chr1\t3\tA\t1\t.\tI\n",
                "chr1\t4\tA\t1\t.\tI\n",
            ]

            offsets = smallerror_long.spool_pileup_lines(lines, regions, spool_path)
            with open(spool_path, "rb") as spool_handle:
                first = smallerror_long.read_pileup_slice(spool_handle, offsets[0])
                second = smallerror_long.read_pileup_slice(spool_handle, offsets[1])

            self.assertEqual(first, lines[:2])
            self.assertEqual(second, lines[1:3])
            self.assertEqual(Path(spool_path).read_text(), "".join(lines[:3]))


class StreamingMergeTests(unittest.TestCase):
    def test_streaming_merge_preserves_order_and_supports_empty_output(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            part_paths = []
            for index, content in enumerate((b"first\n", b"second\n", b"third\n")):
                part_path = os.path.join(temp_dir, f"part{index}.txt")
                Path(part_path).write_bytes(content)
                part_paths.append(part_path)

            output_path = os.path.join(temp_dir, "merged.txt")
            smallerror_long.concatenate_files_streaming(part_paths, output_path, buffer_size=3)
            self.assertEqual(Path(output_path).read_bytes(), b"first\nsecond\nthird\n")

            empty_output_path = os.path.join(temp_dir, "empty.txt")
            smallerror_long.concatenate_files_streaming([], empty_output_path)
            self.assertEqual(Path(empty_output_path).read_bytes(), b"")

    def test_indexed_batch_merge_restores_original_bed_order_with_bounded_fan_in(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            parts = []
            for file_number, content in enumerate((
                "0\tzero\n4\tfour\n",
                "1\tone\n3\tthree\n",
                "2\ttwo\n",
            )):
                part = Path(temp_dir) / f"batch-{file_number}.txt"
                part.write_text(content, encoding="utf-8")
                parts.append(str(part))

            output = Path(temp_dir) / "result.txt"
            smallerror_long.merge_indexed_results_streaming(
                parts,
                str(output),
                max_open_files=2,
            )

            self.assertEqual(
                output.read_text(encoding="utf-8"),
                "zero\none\ntwo\nthree\nfour\n",
            )
            self.assertEqual(list(Path(temp_dir).glob(".sas-index-merge-*")), [])


class PileupLoggingTests(unittest.TestCase):
    def test_per_region_debug_logging_is_opt_in(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            raw_path = Path(temp_dir) / "raw.txt"
            raw_path.write_text("chr1\t1\tA\t1\t.\tI\n")

            quiet_stderr = io.StringIO()
            with mock.patch.object(smallerror_long.sys, "stderr", quiet_stderr):
                quiet_output = smallerror_long.process_pileup_from_file(
                    str(raw_path), 0, 1
                )

            verbose_stderr = io.StringIO()
            with mock.patch.object(smallerror_long.sys, "stderr", verbose_stderr):
                verbose_output = smallerror_long.process_pileup_from_file(
                    str(raw_path), 0, 1, verbose=True
                )

            self.assertEqual(quiet_output, verbose_output)
            self.assertNotIn("Debug:", quiet_stderr.getvalue())
            self.assertIn("Debug: Total lines available", verbose_stderr.getvalue())
            self.assertIn("Debug: Total processed lines", verbose_stderr.getvalue())


class TempDirectorySafetyTests(unittest.TestCase):
    def test_temp_parent_contents_survive_pipeline_cleanup(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_parent = Path(temp_dir) / "caller-owned-temp"
            temp_parent.mkdir()
            sentinel = temp_parent / "do-not-delete.txt"
            sentinel.write_text("keep me")
            bed_path = Path(temp_dir) / "empty.bed"
            bed_path.write_text("")

            smallerror_long.process_bed_regions(
                str(bed_path),
                "input.bam",
                "reference.fa",
                str(Path(temp_dir) / "output.txt"),
                temp_dir_path=str(temp_parent),
                processes=1,
            )

            self.assertEqual(sentinel.read_text(), "keep me")
            self.assertEqual(list(temp_parent.iterdir()), [sentinel])

    def test_two_pipeline_runs_get_distinct_private_directories(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_parent = Path(temp_dir) / "shared-parent"
            bed_path = Path(temp_dir) / "empty.bed"
            bed_path.write_text("")
            created_run_dirs = []
            original_create = smallerror_long.create_run_temp_dir

            def record_created_dir(parent):
                run_dir = original_create(parent)
                created_run_dirs.append(run_dir)
                return run_dir

            with mock.patch.object(
                smallerror_long,
                "create_run_temp_dir",
                side_effect=record_created_dir,
            ):
                for run_number in range(2):
                    smallerror_long.process_bed_regions(
                        str(bed_path),
                        "input.bam",
                        "reference.fa",
                        str(Path(temp_dir) / f"output-{run_number}.txt"),
                        temp_dir_path=str(temp_parent),
                        processes=1,
                    )

            self.assertEqual(len(created_run_dirs), 2)
            self.assertNotEqual(created_run_dirs[0], created_run_dirs[1])
            self.assertTrue(all(Path(path).parent == temp_parent for path in created_run_dirs))
            self.assertTrue(all(not Path(path).exists() for path in created_run_dirs))


@unittest.skipUnless(os.name == "posix", "fake samtools executables require POSIX")
class BatchProcessTests(unittest.TestCase):
    def test_batch_query_is_split_back_to_original_regions(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            fake_samtools = Path(temp_dir) / "samtools"
            fake_samtools.write_text(
                f"#!{sys.executable}\n"
                "for position in range(1, 5):\n"
                "    print(f'chr1\\t{position}\\tA\\t1\\t.\\tI')\n"
            )
            fake_samtools.chmod(0o755)

            outputs = [
                os.path.join(temp_dir, "region_000000000.txt"),
                os.path.join(temp_dir, "region_000000001.txt"),
            ]
            batch = [
                region("chr1", 0, 2, 0, outputs[0]),
                region("chr1", 3, 4, 1, outputs[1]),
            ]

            old_path = os.environ.get("PATH", "")
            os.environ["PATH"] = temp_dir + os.pathsep + old_path
            try:
                result_files = smallerror_long.process_region_batch(
                    (batch, "input.bam", "reference.fa", 0)
                )
            finally:
                os.environ["PATH"] = old_path

            self.assertEqual(len(result_files), 1)
            indexed_lines = Path(result_files[0]).read_text().strip().splitlines()
            self.assertEqual([line.split("\t", 1)[0] for line in indexed_lines], ["0", "1"])
            first_fields = indexed_lines[0].split("\t", 1)[1].split("\t")
            second_fields = indexed_lines[1].split("\t", 1)[1].split("\t")
            self.assertEqual(first_fields[1], "1")
            self.assertEqual(first_fields[7:9], ["0", "2"])
            self.assertEqual(second_fields[1], "4")
            self.assertEqual(second_fields[7:9], ["3", "4"])

    def test_mpileup_failure_message_contains_only_bounded_stderr_tail(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            fake_samtools = Path(temp_dir) / "samtools"
            fake_samtools.write_text(
                f"#!{sys.executable}\n"
                "import sys\n"
                "sys.stderr.write('BEGIN-OF-STDERR\\n')\n"
                f"sys.stderr.write('x' * {smallerror_long.MAX_MPILEUP_STDERR_BYTES + 1024})\n"
                "sys.stderr.write('\\nEND-OF-STDERR\\n')\n"
                "sys.exit(7)\n"
            )
            fake_samtools.chmod(0o755)
            output_path = os.path.join(temp_dir, "region_000000000.txt")
            batch = [region("chr1", 0, 2, 0, output_path)]

            old_path = os.environ.get("PATH", "")
            os.environ["PATH"] = temp_dir + os.pathsep + old_path
            try:
                with self.assertRaises(RuntimeError) as raised:
                    smallerror_long.process_region_batch(
                        (batch, "input.bam", "reference.fa", 0)
                    )
            finally:
                os.environ["PATH"] = old_path

            message = str(raised.exception)
            self.assertIn("stderr truncated", message)
            self.assertIn("END-OF-STDERR", message)
            self.assertNotIn("BEGIN-OF-STDERR", message)
            self.assertLess(
                len(message.encode("utf-8")),
                smallerror_long.MAX_MPILEUP_STDERR_BYTES + 512,
            )

    def test_cli_runs_multiple_batches_with_spawn_workers(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            fake_samtools = Path(temp_dir) / "samtools"
            fake_samtools.write_text(
                f"#!{sys.executable}\n"
                "import sys\n"
                "region = sys.argv[sys.argv.index('-r') + 1]\n"
                "chrom, coordinates = region.rsplit(':', 1)\n"
                "start, end = (int(value) for value in coordinates.split('-', 1))\n"
                "for position in range(start, end + 1):\n"
                "    print(f'{chrom}\\t{position}\\tA\\t1\\t.\\tI')\n"
            )
            fake_samtools.chmod(0o755)

            bed_path = Path(temp_dir) / "regions.bed"
            bed_path.write_text("chr1\t0\t2\nchr1\t3\t4\nchr1\t20000\t20001\n")
            output_path = Path(temp_dir) / "output.txt"
            env = os.environ.copy()
            env["PATH"] = temp_dir + os.pathsep + env.get("PATH", "")

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    str(bed_path),
                    "input.bam",
                    "reference.fa",
                    str(output_path),
                    "--processes",
                    "2",
                ],
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=True,
            )

            self.assertIn("Grouped 3 original regions into 2 mpileup batches", result.stderr)
            output_lines = output_path.read_text().splitlines()
            self.assertEqual(len(output_lines), 3)
            self.assertEqual([line.split("\t")[7:9] for line in output_lines], [
                ["0", "2"],
                ["3", "4"],
                ["20000", "20001"],
            ])


if __name__ == "__main__":
    unittest.main()
