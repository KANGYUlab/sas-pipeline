import logging
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from sas_pipeline import BedMask, SASPipeline


class BedMaskTests(unittest.TestCase):
    def test_python_masking_preserves_bed_half_open_boundaries_and_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            mask_file = tmpdir / "mask.bed"
            input_file = tmpdir / "input.tsv"
            output_file = tmpdir / "output.tsv"
            mask_file.write_text("chr1\t10\t25\n", encoding="utf-8")
            input_file.write_text(
                "Chromosome\tStart\tEnd\tValue\n"
                "chr1\t0\t10\tleft_bookend\n"
                "chr1\t25\t30\tright_bookend\n"
                "chr1\t24\t26\toverlap\n"
                "chr2\t10\t25\tother_chrom\n",
                encoding="utf-8",
            )

            bed_mask = BedMask([mask_file], tmpdir)
            bed_mask.apply_bedtools_intersect_v(input_file, output_file)

            self.assertEqual(
                output_file.read_text(encoding="utf-8"),
                "Chromosome\tStart\tEnd\tValue\n"
                "chr1\t0\t10\tleft_bookend\n"
                "chr1\t25\t30\tright_bookend\n"
                "chr2\t10\t25\tother_chrom\n",
            )

    def test_mask_intervals_are_cached_unioned_and_clipped_to_reference(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            mask_file = tmpdir / "mask.bed"
            mask_file.write_text(
                "chr1\t-5\t5\n"
                "chr1\t10\t20\n"
                "chr1\t15\t25\n"
                "chr1\t90\t110\n"
                "chr1\t40\t40\n",
                encoding="utf-8",
            )

            bed_mask = BedMask([mask_file], tmpdir)
            self.assertEqual(bed_mask.get_masked_length("chr1", 100), 30)

            # The mask is loaded once; later on-disk changes do not alter the
            # cached intervals used by repeated masking and length queries.
            mask_file.write_text("chr1\t0\t100\n", encoding="utf-8")
            self.assertEqual(bed_mask.get_masked_length("chr1", 100), 30)

    def test_se_one_based_closed_coordinates_are_converted_explicitly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            mask_file = tmpdir / "mask.bed"
            input_file = tmpdir / "se.tsv"
            output_file = tmpdir / "se.masked.tsv"
            mask_file.write_text("chr1\t10\t11\n", encoding="utf-8")
            input_file.write_text(
                "Chromosome\tStart\tEnd\tReason\n"
                "chr1\t11\t11\tmasked_base\n"
                "chr1\t12\t12\tright_bookend\n",
                encoding="utf-8",
            )

            bed_mask = BedMask([mask_file], tmpdir)
            bed_mask.apply_bedtools_intersect_v(
                input_file,
                output_file,
                coordinate_mode="one_based_closed",
            )

            self.assertEqual(
                output_file.read_text(encoding="utf-8"),
                "Chromosome\tStart\tEnd\tReason\n"
                "chr1\t12\t12\tright_bookend\n",
            )


class Combine50bpWindowsTests(unittest.TestCase):
    def _pipeline_stub(self, tmpdir):
        pipeline = SASPipeline.__new__(SASPipeline)
        pipeline.logger = logging.getLogger("test_sas_pipeline")
        pipeline.be_output_dir = tmpdir
        pipeline.reference_fai = tmpdir / "reference.fa.fai"
        pipeline.reference_fai.write_text(
            "chr2\t1000\t0\t0\t0\nchr1\t1000\t0\t0\t0\n",
            encoding="utf-8",
        )
        return pipeline

    def test_python_fallback_deduplicates_without_merging_adjacent_windows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            first = tmpdir / "first.bed"
            second = tmpdir / "second.bed"
            output = tmpdir / "candidates.bed"
            first.write_text(
                "chr1\t0\t50\nchr1\t101\t150\nchr2\t10\t20\n",
                encoding="utf-8",
            )
            second.write_text(
                "chr1\t50\t100\nchr1\t0\t50\n",
                encoding="utf-8",
            )

            pipeline = self._pipeline_stub(tmpdir)
            pipeline._combine_50bp_windows_python([first, second], output)

            self.assertEqual(
                output.read_text(encoding="utf-8"),
                "chr2\t10\t20\n"
                "chr1\t0\t50\n"
                "chr1\t50\t100\n"
                "chr1\t101\t150\n",
            )

    def test_bedtools_path_sorts_then_deduplicates_without_merge(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_bed = tmpdir / "input.bed"
            output = tmpdir / "candidates.bed"
            input_bed.write_text(
                "chr1\t50\t100\nchr1\t0\t50\nchr1\t0\t50\n",
                encoding="utf-8",
            )
            pipeline = self._pipeline_stub(tmpdir)
            commands = []

            def fake_run(command, check, stdout, stderr):
                commands.append(command)
                stdout.write("chr1\t0\t50\nchr1\t0\t50\nchr1\t50\t100\n")
                return mock.Mock(returncode=0)

            with mock.patch("sas_pipeline.subprocess.run", side_effect=fake_run):
                pipeline._combine_50bp_windows([input_bed], output)

            self.assertEqual([command[1] for command in commands], ["sort"])
            self.assertEqual(
                output.read_text(encoding="utf-8"),
                "chr1\t0\t50\nchr1\t50\t100\n",
            )


if __name__ == "__main__":
    unittest.main()
