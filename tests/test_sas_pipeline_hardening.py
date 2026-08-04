import contextlib
import logging
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import sas_pipeline
from sas_pipeline import SASPipeline


class PipelineHardeningTests(unittest.TestCase):
    def _pipeline_stub(self, tmpdir):
        pipeline = SASPipeline.__new__(SASPipeline)
        pipeline.logger = mock.Mock(spec=logging.Logger)
        pipeline.be_output_dir = tmpdir / "be"
        pipeline.se_output_dir = tmpdir / "se"
        pipeline.be_output_dir.mkdir()
        pipeline.se_output_dir.mkdir()
        pipeline.bed_mask = False
        pipeline.reference_fai = tmpdir / "reference.fa.fai"
        pipeline.reference_fai.write_text("chr1\t1000\n", encoding="utf-8")
        return pipeline

    def test_file_set_pre_landing_failure_handles_hardlink_and_copy_backups(self):
        for force_copy_backup in (False, True):
            with self.subTest(force_copy_backup=force_copy_backup):
                with tempfile.TemporaryDirectory() as tmpdir:
                    tmpdir = Path(tmpdir)
                    first_source = tmpdir / "first.new"
                    second_source = tmpdir / "second.new"
                    first_destination = tmpdir / "first.txt"
                    second_destination = tmpdir / "second.txt"
                    first_source.write_text("new first\n", encoding="utf-8")
                    second_source.write_text("new second\n", encoding="utf-8")
                    first_destination.write_text("old first\n", encoding="utf-8")
                    second_destination.write_text("old second\n", encoding="utf-8")

                    real_replace = sas_pipeline.os.replace
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
                        mock.patch(
                            "sas_pipeline.os.link",
                            side_effect=OSError("force copied backup"),
                        )
                        if force_copy_backup
                        else contextlib.nullcontext()
                    )
                    with (
                        link_patch,
                        mock.patch(
                            "sas_pipeline.os.replace",
                            side_effect=fail_publish_and_second_rollback,
                        ),
                    ):
                        with self.assertRaisesRegex(
                            OSError, "publish failed before landing"
                        ):
                            sas_pipeline._publish_file_set(
                                [
                                    (first_source, first_destination),
                                    (second_source, second_destination),
                                ]
                            )

                    self.assertEqual(
                        first_destination.read_text(encoding="utf-8"),
                        "old first\n",
                    )
                    self.assertEqual(
                        second_destination.read_text(encoding="utf-8"),
                        "old second\n",
                    )
                    self.assertEqual(
                        list(tmpdir.glob(".publish-backup-*")), []
                    )

    def test_run_command_uses_disk_logs_and_returns_only_bounded_tails(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            pipeline = self._pipeline_stub(Path(tmpdir))
            payload_size = sas_pipeline.MAX_COMMAND_LOG_BYTES * 2

            def fake_run(command, check, stdout, stderr, cwd):
                self.assertIsNot(stdout, subprocess.PIPE)
                self.assertIsNot(stderr, subprocess.PIPE)
                stdout.write(b"o" * payload_size)
                stderr.write(b"e" * payload_size)
                return subprocess.CompletedProcess(command, 0)

            with mock.patch("sas_pipeline.subprocess.run", side_effect=fake_run):
                result = pipeline._run_command(["tool"], "large output")

            self.assertLessEqual(
                len(result.stdout.encode("utf-8")),
                sas_pipeline.MAX_COMMAND_LOG_BYTES,
            )
            self.assertLessEqual(
                len(result.stderr.encode("utf-8")),
                sas_pipeline.MAX_COMMAND_LOG_BYTES,
            )
            self.assertTrue(
                any(
                    "earlier output omitted" in str(call)
                    for call in pipeline.logger.info.call_args_list
                    + pipeline.logger.warning.call_args_list
                )
            )

    def test_run_command_failure_attaches_only_bounded_tails(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            pipeline = self._pipeline_stub(Path(tmpdir))
            payload_size = sas_pipeline.MAX_COMMAND_LOG_BYTES * 2

            def fake_run(command, check, stdout, stderr, cwd):
                stdout.write(b"o" * payload_size)
                stderr.write(b"e" * payload_size)
                raise subprocess.CalledProcessError(9, command)

            with mock.patch("sas_pipeline.subprocess.run", side_effect=fake_run):
                with self.assertRaises(subprocess.CalledProcessError) as caught:
                    pipeline._run_command(["tool"], "failure")

            self.assertLessEqual(len(caught.exception.stdout), sas_pipeline.MAX_COMMAND_LOG_BYTES)
            self.assertLessEqual(len(caught.exception.stderr), sas_pipeline.MAX_COMMAND_LOG_BYTES)

    def test_run_command_retains_complete_numbered_stage_logs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            pipeline = self._pipeline_stub(Path(tmpdir))
            pipeline.stage_logs_dir = Path(tmpdir) / "logs" / "stages"
            payload_size = sas_pipeline.MAX_COMMAND_LOG_BYTES + 123

            def fake_run(command, check, stdout, stderr, cwd):
                stdout.write(b"o" * payload_size)
                stderr.write(b"e" * payload_size)
                return subprocess.CompletedProcess(command, 0)

            with mock.patch("sas_pipeline.subprocess.run", side_effect=fake_run):
                result = pipeline._run_command(["tool"], "large/output stage")

            stdout_logs = list(
                pipeline.stage_logs_dir.glob("001_large_output_stage.stdout.log")
            )
            stderr_logs = list(
                pipeline.stage_logs_dir.glob("001_large_output_stage.stderr.log")
            )
            self.assertEqual(len(stdout_logs), 1)
            self.assertEqual(len(stderr_logs), 1)
            self.assertEqual(stdout_logs[0].stat().st_size, payload_size)
            self.assertEqual(stderr_logs[0].stat().st_size, payload_size)
            self.assertEqual(len(result.stdout.encode("utf-8")), sas_pipeline.MAX_COMMAND_LOG_BYTES)
            self.assertEqual(len(result.stderr.encode("utf-8")), sas_pipeline.MAX_COMMAND_LOG_BYTES)

    def test_extract_chromosomes_streams_faidx_output_and_publishes_index(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            pipeline = self._pipeline_stub(tmpdir)
            output_fasta = tmpdir / "haplotype.fa"
            observed_stdout = []

            def fake_extract(command, check, stdout, stderr):
                observed_stdout.append(stdout)
                stdout.write(f">{command[-1]}\nACGT\n")
                return subprocess.CompletedProcess(command, 0)

            def fake_index(command, description, cwd=None):
                staged_fasta = Path(command[-1])
                Path(str(staged_fasta) + ".fai").write_text(
                    "chr1\t4\n", encoding="utf-8"
                )
                return subprocess.CompletedProcess(command, 0)

            with (
                mock.patch("sas_pipeline.subprocess.run", side_effect=fake_extract),
                mock.patch.object(pipeline, "_run_command", side_effect=fake_index),
            ):
                pipeline._extract_chromosomes(
                    tmpdir / "reference.fa", ["chr1", "chr2"], output_fasta
                )

            self.assertTrue(all(handle is not subprocess.PIPE for handle in observed_stdout))
            self.assertEqual(
                output_fasta.read_text(encoding="utf-8"),
                ">chr1\nACGT\n>chr2\nACGT\n",
            )
            self.assertEqual(
                Path(str(output_fasta) + ".fai").read_text(encoding="utf-8"),
                "chr1\t4\n",
            )
            self.assertEqual(list(tmpdir.glob(".haplotype.fa.extract-*")), [])

    def test_extract_chromosomes_aborts_on_second_contig_failure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            pipeline = self._pipeline_stub(tmpdir)
            output_fasta = tmpdir / "haplotype.fa"
            output_fai = Path(str(output_fasta) + ".fai")
            output_fasta.write_text(">old\nAAAA\n", encoding="utf-8")
            output_fai.write_text("old\t4\n", encoding="utf-8")
            extract_calls = 0

            def fail_second_contig(command, check, stdout, stderr):
                nonlocal extract_calls
                extract_calls += 1
                stdout.write(f">{command[-1]}\n")
                if extract_calls == 2:
                    stdout.write("PARTIAL")
                    stderr.write(b"injected faidx failure")
                    raise subprocess.CalledProcessError(1, command)
                stdout.write("ACGT\n")
                return subprocess.CompletedProcess(command, 0)

            with (
                mock.patch(
                    "sas_pipeline.subprocess.run", side_effect=fail_second_contig
                ),
                mock.patch.object(pipeline, "_run_command") as index_command,
            ):
                with self.assertRaises(subprocess.CalledProcessError):
                    pipeline._extract_chromosomes(
                        tmpdir / "reference.fa", ["chr1", "chr2"], output_fasta
                    )

            index_command.assert_not_called()
            self.assertEqual(output_fasta.read_text(encoding="utf-8"), ">old\nAAAA\n")
            self.assertEqual(output_fai.read_text(encoding="utf-8"), "old\t4\n")
            self.assertEqual(list(tmpdir.glob(".haplotype.fa.extract-*")), [])

    def test_extract_chromosomes_second_publish_failure_rolls_back_pair(self):
        for old_outputs_exist in (True, False):
            with self.subTest(old_outputs_exist=old_outputs_exist):
                with tempfile.TemporaryDirectory() as tmpdir:
                    tmpdir = Path(tmpdir)
                    pipeline = self._pipeline_stub(tmpdir)
                    output_fasta = tmpdir / "haplotype.fa"
                    output_fai = Path(str(output_fasta) + ".fai")
                    if old_outputs_exist:
                        output_fasta.write_text(
                            ">old\nAAAA\n", encoding="utf-8"
                        )
                        output_fai.write_text("old\t4\n", encoding="utf-8")
                        old_mtimes = (
                            output_fasta.stat().st_mtime_ns,
                            output_fai.stat().st_mtime_ns,
                        )

                    def fake_extract(command, check, stdout, stderr):
                        stdout.write(f">{command[-1]}\nACGT\n")
                        return subprocess.CompletedProcess(command, 0)

                    def fake_index(command, description, cwd=None):
                        staged_fasta = Path(command[-1])
                        Path(str(staged_fasta) + ".fai").write_text(
                            "chr1\t4\n", encoding="utf-8"
                        )
                        return subprocess.CompletedProcess(command, 0)

                    real_replace = sas_pipeline.os.replace
                    replace_calls = 0

                    def fail_after_second_replace(source, destination):
                        nonlocal replace_calls
                        replace_calls += 1
                        real_replace(source, destination)
                        if replace_calls == 2:
                            raise OSError("injected second replace failure")

                    with (
                        mock.patch(
                            "sas_pipeline.subprocess.run",
                            side_effect=fake_extract,
                        ),
                        mock.patch.object(
                            pipeline, "_run_command", side_effect=fake_index
                        ),
                        mock.patch(
                            "sas_pipeline.os.replace",
                            side_effect=fail_after_second_replace,
                        ),
                    ):
                        with self.assertRaisesRegex(
                            OSError, "injected second replace failure"
                        ):
                            pipeline._extract_chromosomes(
                                tmpdir / "reference.fa", ["chr1"], output_fasta
                            )

                    if old_outputs_exist:
                        self.assertEqual(
                            output_fasta.read_text(encoding="utf-8"),
                            ">old\nAAAA\n",
                        )
                        self.assertEqual(
                            output_fai.read_text(encoding="utf-8"), "old\t4\n"
                        )
                        self.assertEqual(
                            (
                                output_fasta.stat().st_mtime_ns,
                                output_fai.stat().st_mtime_ns,
                            ),
                            old_mtimes,
                        )
                    else:
                        self.assertFalse(output_fasta.exists())
                        self.assertFalse(output_fai.exists())
                    self.assertEqual(
                        list(tmpdir.glob(".haplotype.fa.extract-*")), []
                    )
                    self.assertEqual(
                        list(tmpdir.glob(".haplotype.fa.publish-backup-*")), []
                    )

    def test_window_sort_failure_cleans_temp_dir_and_falls_back_atomically(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            pipeline = self._pipeline_stub(tmpdir)
            input_bed = tmpdir / "input.bed"
            output_bed = tmpdir / "merged.bed"
            input_bed.write_text("chr1\t0\t50\nchr1\t50\t100\n", encoding="utf-8")
            output_bed.write_text("old partial data\n", encoding="utf-8")

            def fail_sort(command, check, stdout, stderr):
                raise subprocess.CalledProcessError(1, command)

            with mock.patch("sas_pipeline.subprocess.run", side_effect=fail_sort):
                pipeline._combine_50bp_windows([input_bed], output_bed)

            self.assertEqual(
                output_bed.read_text(encoding="utf-8"),
                "chr1\t0\t50\nchr1\t50\t100\n",
            )
            self.assertEqual(list(tmpdir.glob(".50bp-combine-*")), [])
            self.assertFalse((tmpdir / "combined_50bp_windows.tmp").exists())
            self.assertFalse((tmpdir / "sorted_50bp_windows.tmp").exists())

    def test_ir_run_directory_is_private_and_cleaned_on_failure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            pipeline = self._pipeline_stub(Path(tmpdir))
            observed_run_dirs = []

            def fail_generation(run_dir):
                observed_run_dirs.append(Path(run_dir))
                (Path(run_dir) / "mat" / "mat-index").mkdir(parents=True)
                raise RuntimeError("injected failure")

            with mock.patch.object(
                pipeline,
                "_generate_ir_bed_in_work_dir",
                side_effect=fail_generation,
            ):
                with self.assertRaisesRegex(RuntimeError, "injected failure"):
                    pipeline._generate_ir_bed()

            self.assertEqual(len(observed_run_dirs), 1)
            self.assertFalse(observed_run_dirs[0].exists())

    def test_split_removes_only_stale_managed_files_and_validates_full_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            pipeline = self._pipeline_stub(tmpdir)
            pipeline.final_be_results = pipeline.be_output_dir / "final.txt"
            pipeline.final_be_results.write_text(
                "Chromosome\t7\t8\treal_contig_named_chromosome\n"
                "chr1\t1\t2\trow\n",
                encoding="utf-8",
            )
            (pipeline.be_output_dir / "be_results_chr2.txt").write_text(
                "stale\n", encoding="utf-8"
            )
            unrelated = pipeline.be_output_dir / "be_results_notes.csv"
            unrelated.write_text("keep\n", encoding="utf-8")

            pipeline._split_be_results_by_chromosome()

            self.assertFalse((pipeline.be_output_dir / "be_results_chr2.txt").exists())
            self.assertEqual(unrelated.read_text(encoding="utf-8"), "keep\n")
            self.assertEqual(
                (pipeline.be_output_dir / "be_results_Chromosome.txt").read_text(
                    encoding="utf-8"
                ),
                "Chromosome\t7\t8\treal_contig_named_chromosome\n",
            )

            se_source = pipeline.se_output_dir / "se.anomalous_regions.tsv"
            se_source.write_text(
                "Chromosome\tStart\tEnd\tReason\nchr1\t1\t2\trow\n",
                encoding="utf-8",
            )
            (pipeline.se_output_dir / "se_results_chr2.txt").write_text(
                "stale\n", encoding="utf-8"
            )
            se_unrelated = pipeline.se_output_dir / "se_results_notes.csv"
            se_unrelated.write_text("keep\n", encoding="utf-8")

            pipeline._split_se_results_by_chromosome()

            self.assertFalse((pipeline.se_output_dir / "se_results_chr2.txt").exists())
            self.assertEqual(se_unrelated.read_text(encoding="utf-8"), "keep\n")

    def test_split_publish_failure_restores_complete_be_and_se_sets(self):
        header = "Chromosome\tStart\tEnd\tReason"
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            for result_type in ("be", "se"):
                with self.subTest(result_type=result_type):
                    case_dir = tmpdir / result_type
                    case_dir.mkdir()
                    pipeline = self._pipeline_stub(case_dir)
                    if result_type == "be":
                        prefix = "be_results_"
                        output_dir = pipeline.be_output_dir
                        pipeline.final_be_results = output_dir / "final.txt"
                        source = pipeline.final_be_results
                        split = pipeline._split_be_results_by_chromosome
                    else:
                        prefix = "se_results_"
                        output_dir = pipeline.se_output_dir
                        source = output_dir / "se.anomalous_regions.tsv"
                        split = pipeline._split_se_results_by_chromosome

                    source.write_text(
                        header
                        + "\nchr1\t10\t20\tnew replacement\n"
                        + "chr4\t30\t40\tnew chromosome\n",
                        encoding="utf-8",
                    )
                    old_contents = {
                        f"{prefix}chr1.txt": "old chr1\n",
                        f"{prefix}chr2.txt": "old chr2\n",
                        f"{prefix}stale.txt": "old stale\n",
                    }
                    for name, content in old_contents.items():
                        (output_dir / name).write_text(content, encoding="utf-8")

                    real_replace = sas_pipeline.os.replace
                    replace_calls = 0

                    def fail_after_second_replace(source_path, destination_path):
                        nonlocal replace_calls
                        replace_calls += 1
                        real_replace(source_path, destination_path)
                        if replace_calls == 2:
                            raise OSError("injected second replace failure")

                    with mock.patch(
                        "sas_pipeline.os.replace",
                        side_effect=fail_after_second_replace,
                    ):
                        with self.assertRaisesRegex(
                            OSError, "injected second replace failure"
                        ):
                            split()

                    self.assertEqual(
                        {
                            path.name: path.read_text(encoding="utf-8")
                            for path in output_dir.glob(f"{prefix}*.txt")
                        },
                        old_contents,
                    )
                    self.assertFalse((output_dir / f"{prefix}chr4.txt").exists())
                    self.assertEqual(
                        list(output_dir.glob(f".{prefix}backup-*")), []
                    )


if __name__ == "__main__":
    unittest.main()
