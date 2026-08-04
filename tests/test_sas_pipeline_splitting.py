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
from sas_pipeline import BedMask, SASPipeline


class ChromosomeSplittingTests(unittest.TestCase):
    def _pipeline_stub(self, tmpdir):
        pipeline = SASPipeline.__new__(SASPipeline)
        pipeline.logger = logging.getLogger("test_sas_pipeline_splitting")
        pipeline.be_output_dir = tmpdir / "be"
        pipeline.se_output_dir = tmpdir / "se"
        pipeline.be_output_dir.mkdir()
        pipeline.se_output_dir.mkdir()
        pipeline.bed_mask = False
        return pipeline

    def test_se_chr_records_are_not_mistaken_for_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            pipeline = self._pipeline_stub(tmpdir)
            se_file = pipeline.se_output_dir / "se.anomalous_regions.tsv"
            se_file.write_text(
                "Chromosome\tStart\tEnd\tReason\n"
                "chr1\t1\t10\tfirst\n"
                "chr2\t11\t20\tsecond\n"
                "chr1\t21\t30\tthird\n",
                encoding="utf-8",
            )

            pipeline._split_se_results_by_chromosome()

            self.assertEqual(
                (pipeline.se_output_dir / "se_results_chr1.txt").read_text(
                    encoding="utf-8"
                ),
                "Chromosome\tStart\tEnd\tReason\n"
                "chr1\t1\t10\tfirst\n"
                "chr1\t21\t30\tthird\n",
            )
            self.assertEqual(
                (pipeline.se_output_dir / "se_results_chr2.txt").read_text(
                    encoding="utf-8"
                ),
                "Chromosome\tStart\tEnd\tReason\nchr2\t11\t20\tsecond\n",
            )

    def test_be_and_se_splitters_bound_open_files_and_preserve_order(self):
        contig_count = sas_pipeline.MAX_OPEN_CHROMOSOME_FILES + 6
        header = "Chromosome\tStart\tEnd\tReason"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            for result_type in ("be", "se"):
                with self.subTest(result_type=result_type):
                    case_dir = tmpdir / result_type
                    case_dir.mkdir()
                    pipeline = self._pipeline_stub(case_dir)
                    rows = [
                        f"contig{i}\t{i}\t{i + 1}\tfirst"
                        for i in range(contig_count)
                    ]
                    rows.extend(
                        f"contig{i}\t{i + 100}\t{i + 101}\tsecond"
                        for i in range(contig_count)
                    )

                    if result_type == "be":
                        pipeline.final_be_results = pipeline.be_output_dir / "final.txt"
                        source = pipeline.final_be_results
                        split = pipeline._split_be_results_by_chromosome
                        output_dir = pipeline.be_output_dir
                        prefix = "be_results_"
                    else:
                        source = pipeline.se_output_dir / "se.anomalous_regions.tsv"
                        split = pipeline._split_se_results_by_chromosome
                        output_dir = pipeline.se_output_dir
                        prefix = "se_results_"

                    source.write_text(
                        header + "\n" + "\n".join(rows) + "\n", encoding="utf-8"
                    )

                    observed_open_counts = []
                    original_open = sas_pipeline._ChromosomeFilePool._open

                    def tracked_open(pool, chrom):
                        handle = original_open(pool, chrom)
                        observed_open_counts.append(len(pool._handles))
                        return handle

                    with mock.patch.object(
                        sas_pipeline._ChromosomeFilePool,
                        "_open",
                        autospec=True,
                        side_effect=tracked_open,
                    ):
                        split()

                    self.assertLessEqual(
                        max(observed_open_counts),
                        sas_pipeline.MAX_OPEN_CHROMOSOME_FILES,
                    )
                    self.assertEqual(
                        (output_dir / f"{prefix}contig0.txt").read_text(
                            encoding="utf-8"
                        ),
                        header
                        + "\ncontig0\t0\t1\tfirst\n"
                        + "contig0\t100\t101\tsecond\n",
                    )

    def test_file_pool_closes_all_handles_when_split_fails(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            opened_handles = []
            real_open = open

            def tracking_open(*args, **kwargs):
                handle = real_open(*args, **kwargs)
                opened_handles.append(handle)
                return handle

            pool = sas_pipeline._ChromosomeFilePool(
                Path(tmpdir), "results_", max_open_files=2
            )
            with mock.patch("builtins.open", side_effect=tracking_open):
                with self.assertRaisesRegex(RuntimeError, "injected failure"):
                    with pool:
                        pool.write("chr1", "chr1\t1")
                        pool.write("chr2", "chr2\t2")
                        raise RuntimeError("injected failure")

            self.assertEqual(len(pool._handles), 0)
            self.assertTrue(all(handle.closed for handle in opened_handles))


class BedMaskTemporaryFileTests(unittest.TestCase):
    def test_bedtools_failure_removes_named_temporary_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            mask_file = tmpdir / "mask.bed"
            input_file = tmpdir / "input.bed"
            output_file = tmpdir / "output.bed"
            mask_file.write_text("chr1\t10\t20\n", encoding="utf-8")
            input_file.write_text(
                "chr1\t0\t5\nchr1\t12\t13\n", encoding="utf-8"
            )
            bed_mask = BedMask([mask_file], tmpdir)
            created_temp_files = []
            real_named_temporary_file = tempfile.NamedTemporaryFile

            def recording_named_temporary_file(*args, **kwargs):
                handle = real_named_temporary_file(*args, **kwargs)
                created_temp_files.append(handle.name)
                return handle

            with (
                mock.patch(
                    "tempfile.NamedTemporaryFile",
                    side_effect=recording_named_temporary_file,
                ),
                mock.patch(
                    "sas_pipeline.subprocess.run",
                    side_effect=subprocess.CalledProcessError(1, "bedtools"),
                ),
            ):
                bed_mask.mask_bed_regions(input_file, output_file, "test")

            self.assertEqual(len(created_temp_files), 1)
            self.assertFalse(Path(created_temp_files[0]).exists())
            self.assertEqual(output_file.read_text(encoding="utf-8"), "chr1\t0\t5\n")


if __name__ == "__main__":
    unittest.main()
