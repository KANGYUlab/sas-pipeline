import contextlib
import gzip
import importlib.util
import io
import os
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "src" / "integrated_se_analysis.py"


def load_module_without_optional_dependencies():
    fake_numpy = types.ModuleType("numpy")
    fake_numpy.mean = lambda values: sum(values) / len(values)
    fake_pandas = types.ModuleType("pandas")
    fake_pysam = types.ModuleType("pysam")
    fake_pysam.AlignmentFile = object

    replacements = {
        "numpy": fake_numpy,
        "pandas": fake_pandas,
        "pysam": fake_pysam,
    }
    previous = {name: sys.modules.get(name) for name in replacements}
    sys.modules.update(replacements)
    try:
        spec = importlib.util.spec_from_file_location(
            "integrated_se_analysis_under_test", str(SCRIPT_PATH)
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        for name, old_value in previous.items():
            if old_value is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = old_value


se = load_module_without_optional_dependencies()


class VersionTests(unittest.TestCase):
    def test_version_flag_reports_package_version(self):
        with mock.patch.object(sys, "argv", ["sas-se-analysis", "--version"]):
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                with self.assertRaises(SystemExit) as exit_context:
                    se.main()

        self.assertEqual(exit_context.exception.code, 0)
        self.assertEqual(stdout.getvalue().strip(), "sas-se-analysis 2.0.0")


class SpanningPublicationTests(unittest.TestCase):
    def test_pre_landing_failure_handles_hardlink_and_copy_backups(self):
        for force_copy_backup in (False, True):
            with self.subTest(force_copy_backup=force_copy_backup):
                with tempfile.TemporaryDirectory() as temp_dir:
                    temp_dir = Path(temp_dir)
                    first_source = temp_dir / "first.new"
                    second_source = temp_dir / "second.new"
                    first_destination = temp_dir / "first.bed"
                    second_destination = temp_dir / "second.bed"
                    first_source.write_text("new first\n")
                    second_source.write_text("new second\n")
                    first_destination.write_text("old first\n")
                    second_destination.write_text("old second\n")

                    real_replace = se.os.replace
                    replace_calls = 0

                    def fail_publish_and_second_rollback(source, destination):
                        nonlocal replace_calls
                        replace_calls += 1
                        if replace_calls == 2:
                            raise OSError("publish failed before landing")
                        if replace_calls == 3:
                            raise OSError("rollback replace failed")
                        return real_replace(source, destination)

                    link_patch = (
                        mock.patch.object(
                            se.os,
                            "link",
                            side_effect=OSError("force copied backup"),
                        )
                        if force_copy_backup
                        else contextlib.nullcontext()
                    )
                    with (
                        link_patch,
                        mock.patch.object(
                            se.os,
                            "replace",
                            side_effect=fail_publish_and_second_rollback,
                        ),
                    ):
                        with self.assertRaisesRegex(
                            OSError, "publish failed before landing"
                        ):
                            se._publish_spanning_output_pair(
                                [
                                    (first_source, first_destination),
                                    (second_source, second_destination),
                                ],
                                ".spanning-backup-",
                            )

                    self.assertEqual(first_destination.read_text(), "old first\n")
                    self.assertEqual(second_destination.read_text(), "old second\n")
                    self.assertEqual(
                        list(temp_dir.glob(".spanning-backup-*")), []
                    )


class EmptyAlignmentFile:
    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return False

    def fetch(self, chrom):
        return iter(())


def alignment_file_with_reads(reads):
    class AlignmentFileWithReads:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            return False

        def fetch(self, chrom):
            return iter(reads)

    return AlignmentFileWithReads


class FakeRead:
    def __init__(self, start, end, cigartuples):
        self.is_unmapped = False
        self.reference_start = start
        self.reference_end = end
        self.cigartuples = cigartuples


class BoundaryCountingTests(unittest.TestCase):
    def test_deletion_ending_on_window_boundary_does_not_mark_next_window(self):
        self.assertEqual(se.calculate_affected_windows(1500, 500, 2000), [0])
        self.assertEqual(
            se.calculate_affected_windows(1500, 501, 2000),
            [0, 2000],
        )
        self.assertEqual(se.calculate_affected_windows(2000, 0, 2000), [])

    def test_reference_end_is_exclusive_for_total_read_windows(self):
        read = FakeRead(0, 2000, [(0, 2000)])
        alignment_file = alignment_file_with_reads([read])

        with mock.patch.object(se.pysam, "AlignmentFile", alignment_file):
            windows = se.process_chromosome_integrated(
                "input.bam",
                "chr1",
                6000,
                window_size=2000,
                terminal_region_size=1,
            )

        self.assertEqual(windows[("chr1", 0)]["TotalReads"], 1)
        self.assertNotIn(("chr1", 2000), windows)

    def test_reads_without_reference_span_do_not_create_total_read_windows(self):
        reads = [
            FakeRead(2000, None, [(1, 50)]),
            FakeRead(4000, 4000, [(1, 50)]),
        ]
        alignment_file = alignment_file_with_reads(reads)

        with mock.patch.object(se.pysam, "AlignmentFile", alignment_file):
            windows = se.process_chromosome_integrated(
                "input.bam", "chr1", 8000, window_size=2000
            )

        self.assertEqual(windows[("chr1", 2000)]["TotalReads"], 0)
        self.assertEqual(windows[("chr1", 4000)]["TotalReads"], 0)

    def test_custom_terminal_region_size_controls_reads_without_count(self):
        read = FakeRead(8000, 9000, [(0, 1000)])
        alignment_file = alignment_file_with_reads([read])

        with mock.patch.object(se.pysam, "AlignmentFile", alignment_file):
            windows = se.process_chromosome_integrated(
                "input.bam",
                "chr1",
                20000,
                window_size=1000,
                terminal_region_size=1000,
            )

        self.assertEqual(windows[("chr1", 8000)]["TotalReads"], 1)
        self.assertEqual(
            windows[("chr1", 8000)]["ReadsWithoutINDELorClip"],
            0,
        )

    def test_store_helper_passes_terminal_region_size_to_worker(self):
        counters = mock.Mock()
        with mock.patch.object(
            se, "process_chromosome_integrated", return_value={}
        ) as process:
            se.process_chromosome_and_store(
                counters,
                "input.bam",
                "chr1",
                20000,
                50,
                2000,
                1234,
            )

        process.assert_called_once_with(
            "input.bam", "chr1", 20000, 50, 2000, 1234
        )
        counters.store_chromosome_windows.assert_called_once_with("chr1", {})


class SparseWindowTests(unittest.TestCase):
    def test_all_zero_depth_keeps_independent_non_spanning_evidence(self):
        counts = se.new_window_counts()
        counts["HasSpanningReads"] = False

        anomalies = se.identify_anomalous_regions(
            {("chr1", 0): counts},
            {("chr1", 0): 0.0},
            {"chr1": 2000},
            2000,
            0,
        )

        self.assertEqual(len(anomalies), 1)
        self.assertEqual(anomalies[0]["anomaly_tags"], "non_spanning")

    def test_zero_depth_window_is_low_depth_when_genome_mean_is_positive(self):
        anomalies = se.identify_anomalous_regions(
            {},
            {("chr1", 0): 0.0, ("chr1", 2000): 10.0},
            {"chr1": 4000},
            2000,
            0,
        )

        self.assertEqual(len(anomalies), 1)
        self.assertEqual(anomalies[0]["start"], 1)
        self.assertEqual(anomalies[0]["anomaly_tags"], "low_depth")

    def test_depth_loader_rejects_empty_or_malformed_mosdepth_output(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            empty_path = Path(temp_dir) / "empty.regions.bed.gz"
            malformed_path = Path(temp_dir) / "bad.regions.bed.gz"
            duplicate_path = Path(temp_dir) / "duplicate.regions.bed.gz"
            with gzip.open(empty_path, "wt"):
                pass
            with gzip.open(malformed_path, "wt") as output:
                output.write("chr1\tbad\t2000\t10\n")
            with gzip.open(duplicate_path, "wt") as output:
                output.write("chr1\t0\t2000\t10\nchr1\t0\t2000\t11\n")

            with self.assertRaises(RuntimeError):
                se.load_depth_data(str(empty_path))
            with self.assertRaises(RuntimeError):
                se.load_depth_data(str(malformed_path))
            with self.assertRaises(RuntimeError):
                se.load_depth_data(str(duplicate_path))

    def test_depth_window_validation_rejects_truncation_and_unknown_contigs(self):
        with self.assertRaisesRegex(RuntimeError, "expected 2, observed 1"):
            se.validate_depth_windows(
                {("chr1", 0): 10.0}, {"chr1": 4000}, 2000
            )
        with self.assertRaisesRegex(RuntimeError, "unexpected window"):
            se.validate_depth_windows(
                {("unknown", 0): 10.0, ("chr1", 2000): 10.0},
                {"chr1": 4000},
                2000,
            )

    def test_empty_chromosome_retains_no_dense_counter_dictionary(self):
        with mock.patch.object(se.pysam, "AlignmentFile", EmptyAlignmentFile):
            windows = se.process_chromosome_integrated(
                "input.bam", "chr1", 1_000_000, window_size=2000
            )
        self.assertEqual(windows, {})

    def test_zero_windows_are_synthesized_for_output(self):
        counts = se.new_window_counts()
        counts["TotalReads"] = 2
        windows = {("chr1", 2000): counts}

        with tempfile.TemporaryDirectory() as temp_dir:
            prefix = os.path.join(temp_dir, "sample")
            se.save_results(
                windows,
                {("chr1", 2000): 7.5},
                [],
                {"chr1": 5000},
                prefix,
                2000,
                1000,
            )
            lines = Path(prefix + ".window_stats.tsv").read_text().splitlines()

        self.assertEqual(len(lines), 4)  # header plus all three reference windows
        self.assertEqual([line.split("\t")[1] for line in lines[1:]], ["1", "2001", "4001"])
        self.assertEqual(lines[2].split("\t")[10], "2")

    def test_sparse_and_dense_inputs_produce_same_anomalies(self):
        chrom_lengths = {"chr1": 6000}
        depth = {
            ("chr1", 0): 10.0,
            ("chr1", 2000): 10.0,
            ("chr1", 4000): 10.0,
        }
        event = se.new_window_counts()
        event["ClipCount"] = 1
        event["ReadsWithClip"] = 1
        event["TotalReads"] = 1
        non_spanning = se.new_window_counts()
        non_spanning["HasSpanningReads"] = False
        sparse = {
            ("chr1", 2000): event,
            ("chr1", 4000): non_spanning,
        }
        dense = {("chr1", start): se.new_window_counts() for start in (0, 2000, 4000)}
        dense.update(sparse)

        sparse_result = se.identify_anomalous_regions(
            sparse, depth, chrom_lengths, 2000, 0
        )
        dense_result = se.identify_anomalous_regions(
            dense, depth, chrom_lengths, 2000, 0
        )
        self.assertEqual(sparse_result, dense_result)

        with tempfile.TemporaryDirectory() as temp_dir:
            sparse_prefix = os.path.join(temp_dir, "sparse")
            dense_prefix = os.path.join(temp_dir, "dense")
            se.save_results(
                sparse, depth, sparse_result, chrom_lengths,
                sparse_prefix, 2000, 0
            )
            se.save_results(
                dense, depth, dense_result, chrom_lengths,
                dense_prefix, 2000, 0
            )
            self.assertEqual(
                Path(sparse_prefix + ".window_stats.tsv").read_bytes(),
                Path(dense_prefix + ".window_stats.tsv").read_bytes(),
            )
            self.assertEqual(
                Path(sparse_prefix + ".anomalous_regions.tsv").read_bytes(),
                Path(dense_prefix + ".anomalous_regions.tsv").read_bytes(),
            )

    def test_anomalies_follow_bam_chromosome_order(self):
        chrom_lengths = {"chr2": 2000, "chr10": 2000}
        depth = {
            ("chr2", 0): 10.0,
            ("chr10", 0): 10.0,
        }
        windows = {}
        for chrom in chrom_lengths:
            counts = se.new_window_counts()
            counts["HasSpanningReads"] = False
            windows[(chrom, 0)] = counts

        anomalies = se.identify_anomalous_regions(
            windows, depth, chrom_lengths, 2000, 0
        )

        self.assertEqual(
            [region["chromosome"] for region in anomalies],
            ["chr2", "chr10"],
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            prefix = os.path.join(temp_dir, "ordered")
            se.save_results(
                windows, depth, anomalies, chrom_lengths,
                prefix, 2000, 0
            )
            rows = Path(prefix + ".anomalous_regions.tsv").read_text().splitlines()

        self.assertEqual(
            [row.split("\t", 1)[0] for row in rows[1:]],
            ["chr2", "chr10"],
        )


class FailureAndBedTests(unittest.TestCase):
    def test_chromosome_open_failure_propagates(self):
        class BrokenAlignmentFile:
            def __init__(self, *args, **kwargs):
                raise RuntimeError("broken BAM")

        with mock.patch.object(se.pysam, "AlignmentFile", BrokenAlignmentFile):
            with self.assertRaisesRegex(RuntimeError, "broken BAM"):
                se.process_chromosome_integrated("bad.bam", "chr1", 1000)

    def test_empty_non_spanning_bed_is_valid(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bed_file = Path(temp_dir) / "empty.bed"
            bed_file.write_text("", encoding="utf-8")
            self.assertEqual(
                se.load_merged_non_spanning_regions_from_bed(str(bed_file)), []
            )

    def test_refresh_spanning_outputs_atomically_replaces_stale_files(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            prefix = os.path.join(temp_dir, "sample")
            merged = Path(prefix + ".non_spanning_regions.merged.bed")
            windows = Path(prefix + ".non_spanning.bed")
            merged.write_text("chrOld\t0\t10\n")
            windows.write_text("chrOld\t0\t10\n")

            def fake_spanning(_bam, _fai, temporary_prefix, *_args):
                self.assertNotEqual(temporary_prefix, prefix)
                Path(temporary_prefix + ".non_spanning.bed").write_text("")
                Path(
                    temporary_prefix + ".non_spanning_regions.merged.bed"
                ).write_text("chr1\t20\t40\n")

            with mock.patch.object(
                    se, "call_spanning_analysis", side_effect=fake_spanning):
                result = se.refresh_spanning_outputs(
                    "input.bam", "reference.fa.fai", prefix, threads=1
                )

            self.assertEqual(result, str(merged))
            self.assertEqual(merged.read_text(), "chr1\t20\t40\n")
            self.assertEqual(windows.read_text(), "")
            self.assertEqual(
                {path.name for path in Path(temp_dir).iterdir()},
                {merged.name, windows.name},
            )

    def test_refresh_spanning_second_publish_failure_rolls_back_pair(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            prefix = os.path.join(temp_dir, "sample")
            windows = Path(prefix + ".non_spanning.bed")
            merged = Path(prefix + ".non_spanning_regions.merged.bed")
            windows.write_text("chrOld\t0\t10\n")
            merged.write_text("chrOld\t10\t20\n")
            old_mtimes = (
                windows.stat().st_mtime_ns,
                merged.stat().st_mtime_ns,
            )

            def fake_spanning(_bam, _fai, temporary_prefix, *_args):
                Path(temporary_prefix + ".non_spanning.bed").write_text(
                    "chr1\t20\t40\n"
                )
                Path(
                    temporary_prefix + ".non_spanning_regions.merged.bed"
                ).write_text("chr1\t20\t60\n")

            real_replace = se.os.replace
            replace_calls = 0

            def fail_after_second_replace(source, destination):
                nonlocal replace_calls
                replace_calls += 1
                real_replace(source, destination)
                if replace_calls == 2:
                    raise OSError("injected second replace failure")

            with (
                mock.patch.object(
                    se, "call_spanning_analysis", side_effect=fake_spanning
                ),
                mock.patch.object(
                    se.os, "replace", side_effect=fail_after_second_replace
                ),
            ):
                with self.assertRaisesRegex(
                    OSError, "injected second replace failure"
                ):
                    se.refresh_spanning_outputs(
                        "input.bam", "reference.fa.fai", prefix, threads=1
                    )

            self.assertEqual(windows.read_text(), "chrOld\t0\t10\n")
            self.assertEqual(merged.read_text(), "chrOld\t10\t20\n")
            self.assertEqual(
                (windows.stat().st_mtime_ns, merged.stat().st_mtime_ns),
                old_mtimes,
            )
            self.assertEqual(
                {path.name for path in Path(temp_dir).iterdir()},
                {merged.name, windows.name},
            )

    def test_result_bed_failure_preserves_old_output_and_removes_temp_files(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            anomalies = Path(temp_dir) / "sample.anomalous_regions.tsv"
            anomalies.write_text(
                "Chromosome\tStart\tEnd\n"
                "chr1\t1\t2000\n"
            )
            output = Path(temp_dir) / "sample.result.bed"
            output.write_text("chrOld\t0\t10\n")
            bedtools_error = se.subprocess.CalledProcessError(1, "bedtools")
            temporary_paths = []

            def fail_after_partial_output(cmd, stdout, check):
                temporary_paths.extend([Path(cmd[-1]), Path(stdout.name)])
                stdout.write("chrPartial\t0\t100\n")
                stdout.flush()
                raise bedtools_error

            with mock.patch.object(
                    se.subprocess, "run", side_effect=fail_after_partial_output):
                with self.assertRaises(se.subprocess.CalledProcessError):
                    se.generate_se_result_bed(str(anomalies), str(output))

            self.assertEqual(output.read_text(), "chrOld\t0\t10\n")
            self.assertEqual(len(temporary_paths), 2)
            self.assertNotEqual(temporary_paths[0], temporary_paths[1])
            self.assertTrue(all(path.parent == output.parent for path in temporary_paths))
            self.assertTrue(all(not path.exists() for path in temporary_paths))
            self.assertEqual(
                {path.name for path in Path(temp_dir).iterdir()},
                {anomalies.name, output.name},
            )

    def test_default_pipeline_prefix_publishes_hard_link_aliases(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            source = Path(temp_dir) / "structural_errors.anomalous_regions.tsv"
            source.write_text("header\n")
            se.publish_pipeline_compatibility_outputs(temp_dir, "structural_errors")
            alias = Path(temp_dir) / "se.anomalous_regions.tsv"
            self.assertEqual(alias.read_text(), "header\n")
            self.assertEqual(os.stat(str(source)).st_ino, os.stat(str(alias)).st_ino)


class MainFlowTests(unittest.TestCase):
    def _input_files(self, temp_dir, prefix):
        bam = Path(temp_dir) / "input.bam"
        reference = Path(temp_dir) / "reference.fa"
        bam.write_bytes(b"")
        reference.write_text(">chr1\nA\n")
        Path(str(reference) + ".fai").write_text("chr1\t4000\t0\t0\t0\n")
        Path(temp_dir, prefix + ".non_spanning_regions.merged.bed").write_text("")
        return bam, reference

    def test_spanning_helper_uses_current_python_interpreter(self):
        completed = se.subprocess.CompletedProcess([], 0, stdout="", stderr="")
        with mock.patch.object(se.subprocess, "run", return_value=completed) as run:
            se.call_spanning_analysis(
                "input.bam", "reference.fa.fai", "output", threads=1
            )

        self.assertEqual(run.call_args.args[0][0], sys.executable)

    def test_main_honors_custom_output_prefix(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bam, reference = self._input_files(temp_dir, "custom")
            argv = [
                str(SCRIPT_PATH), str(bam), str(reference), "custom",
                "--output-dir", temp_dir, "--threads", "1",
                "--terminal-region-size", "1234",
            ]
            generated = []
            process_chromosome = mock.Mock(return_value={})

            def fake_generate(anomalies_file, output_bed):
                generated.append((anomalies_file, output_bed))
                Path(output_bed).write_text("")

            spanning_calls = []

            def fake_spanning(_bam, _fai, temporary_prefix, *_args):
                spanning_calls.append(temporary_prefix)
                Path(temporary_prefix + ".non_spanning.bed").write_text("")
                Path(
                    temporary_prefix + ".non_spanning_regions.merged.bed"
                ).write_text("")

            patches = (
                mock.patch.object(sys, "argv", argv),
                mock.patch.object(se, "run_mosdepth", return_value="unused.gz"),
                mock.patch.object(
                    se,
                    "load_depth_data",
                    return_value={("chr1", 0): 0.0, ("chr1", 2000): 0.0},
                ),
                mock.patch.object(se, "get_chromosome_lengths_from_bam", return_value={"chr1": 4000}),
                mock.patch.object(
                    se, "process_chromosome_integrated", process_chromosome
                ),
                mock.patch.object(se, "call_spanning_analysis", side_effect=fake_spanning),
                mock.patch.object(se, "generate_se_result_bed", side_effect=fake_generate),
                mock.patch.object(se.pysam, "AlignmentFile", EmptyAlignmentFile),
            )
            with contextlib.ExitStack() as stack:
                for patcher in patches:
                    stack.enter_context(patcher)
                se.main()

            custom_prefix = os.path.join(str(Path(temp_dir).resolve()), "custom")
            self.assertTrue(os.path.exists(custom_prefix + ".window_stats.tsv"))
            self.assertTrue(os.path.exists(custom_prefix + ".anomalous_regions.tsv"))
            self.assertEqual(generated, [
                (custom_prefix + ".anomalous_regions.tsv", custom_prefix + ".result.bed")
            ])
            self.assertEqual(len(spanning_calls), 1)
            self.assertNotEqual(spanning_calls[0], custom_prefix)
            self.assertFalse(os.path.exists(os.path.join(temp_dir, "se.anomalous_regions.tsv")))
            self.assertEqual(process_chromosome.call_args.args[-1], 1234)

    def test_stale_spanning_bed_is_not_loaded_when_recompute_fails(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bam, reference = self._input_files(temp_dir, "failed-spanning")
            stale_bed = Path(
                temp_dir, "failed-spanning.non_spanning_regions.merged.bed"
            )
            stale_bed.write_text("chrStale\t0\t4000\n")
            argv = [
                str(SCRIPT_PATH), str(bam), str(reference), "failed-spanning",
                "--output-dir", temp_dir, "--threads", "1"
            ]
            load_merged = mock.Mock()
            spanning_error = se.subprocess.CalledProcessError(1, "spanning.py")
            patches = (
                mock.patch.object(sys, "argv", argv),
                mock.patch.object(se, "run_mosdepth", return_value="unused.gz"),
                mock.patch.object(
                    se,
                    "load_depth_data",
                    return_value={("chr1", 0): 0.0, ("chr1", 2000): 0.0},
                ),
                mock.patch.object(
                    se, "get_chromosome_lengths_from_bam",
                    return_value={"chr1": 4000}
                ),
                mock.patch.object(se, "process_chromosome_integrated", return_value={}),
                mock.patch.object(
                    se, "call_spanning_analysis", side_effect=spanning_error
                ),
                mock.patch.object(
                    se, "load_merged_non_spanning_regions_from_bed", load_merged
                ),
                mock.patch.object(se.pysam, "AlignmentFile", EmptyAlignmentFile),
            )
            with contextlib.ExitStack() as stack:
                for patcher in patches:
                    stack.enter_context(patcher)
                with contextlib.redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit) as raised:
                        se.main()

            self.assertEqual(raised.exception.code, 1)
            load_merged.assert_not_called()
            # The last complete result remains intact because recomputation is
            # isolated under a temporary prefix until all outputs exist.
            self.assertEqual(stale_bed.read_text(), "chrStale\t0\t4000\n")

    def test_worker_failure_exits_before_result_generation(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bam, reference = self._input_files(temp_dir, "failed")
            argv = [
                str(SCRIPT_PATH), str(bam), str(reference), "failed",
                "--output-dir", temp_dir, "--threads", "1"
            ]
            generate = mock.Mock()
            patches = (
                mock.patch.object(sys, "argv", argv),
                mock.patch.object(se, "run_mosdepth", return_value="unused.gz"),
                mock.patch.object(
                    se,
                    "load_depth_data",
                    return_value={("chr1", 0): 0.0, ("chr1", 2000): 0.0},
                ),
                mock.patch.object(se, "get_chromosome_lengths_from_bam", return_value={"chr1": 4000}),
                mock.patch.object(se, "process_chromosome_integrated", side_effect=RuntimeError("worker failed")),
                mock.patch.object(se, "generate_se_result_bed", generate),
                mock.patch.object(se.pysam, "AlignmentFile", EmptyAlignmentFile),
            )
            with contextlib.ExitStack() as stack:
                for patcher in patches:
                    stack.enter_context(patcher)
                with contextlib.redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit) as raised:
                        se.main()

            self.assertEqual(raised.exception.code, 1)
            generate.assert_not_called()


if __name__ == "__main__":
    unittest.main()
