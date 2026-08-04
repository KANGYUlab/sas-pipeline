import logging
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from sas_pipeline import (
    BedMask,
    SASPipeline,
    validate_input_files,
    validate_runtime_dependencies,
)


class SasPipelineBookkeepingTests(unittest.TestCase):
    def test_pipeline_uses_separate_final_work_log_and_metadata_directories(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "results"
            config = {
                "short_read_bam": "short.bam",
                "hifi_bam": "hifi.bam",
                "ont_bam": "ont.bam",
                "reference_fa": "reference.fa",
                "output_dir": str(output_dir),
                "scripts_dir": str(Path(temp_dir) / "scripts"),
                "src_dir": str(Path(temp_dir) / "src"),
                "run_be_analysis": True,
                "run_se_analysis": False,
                "skip_filter": False,
                "threads": 64,
            }

            with mock.patch.object(
                SASPipeline, "_setup_logging", return_value=logging.getLogger("layout-test")
            ), mock.patch.object(SASPipeline, "_setup_reference_fai"), mock.patch.object(
                SASPipeline, "_get_chromosomes", return_value=[]
            ):
                pipeline = SASPipeline(config)

            self.assertEqual(
                pipeline.final_be_results,
                output_dir / "final" / "be" / "final_be_results.txt",
            )
            self.assertEqual(
                pipeline.short_filtered_bam,
                output_dir / "work" / "alignments" / "short_read_filtered.bam",
            )
            self.assertEqual(
                pipeline.summary_file,
                output_dir / "final" / "sas_pipeline_summary.txt",
            )
            self.assertEqual(
                pipeline.config_file,
                output_dir / "metadata" / "pipeline_config.json",
            )
            self.assertEqual(pipeline.pileup_workers, 16)
            self.assertEqual(pipeline.se_workers, 8)
            self.assertTrue((output_dir / "logs" / "stages").is_dir())
            self.assertFalse((output_dir / "final" / "se").exists())

    def test_publish_se_results_exposes_only_stable_public_names(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            pipeline = SASPipeline.__new__(SASPipeline)
            pipeline.logger = logging.getLogger("se-publish-test")
            pipeline.se_output_dir = root / "work" / "se"
            pipeline.se_results_dir = root / "final" / "se"
            pipeline.se_output_dir.mkdir(parents=True)
            pipeline.bed_mask = None

            source_files = {
                "structural_errors.result.bed": "chr1\t0\t10\n",
                "se.anomalous_regions.tsv": "Chromosome\tStart\tEnd\nchr1\t1\t10\n",
                "structural_errors.anomaly_summary.tsv": "Chromosome\tnr-SE\nchr1\t1\n",
                "structural_errors.window_stats.tsv": "large diagnostic file\n",
            }
            for name, contents in source_files.items():
                (pipeline.se_output_dir / name).write_text(contents, encoding="utf-8")
            pipeline.se_anomalies_file = (
                pipeline.se_output_dir / "se.anomalous_regions.tsv"
            )

            pipeline._publish_se_results()

            self.assertEqual(
                sorted(path.name for path in pipeline.se_results_dir.iterdir()),
                ["anomalous_regions.tsv", "anomaly_summary.tsv", "structural_errors.bed"],
            )
            self.assertTrue(
                (pipeline.se_output_dir / "structural_errors.window_stats.tsv").exists()
            )

    def test_publish_masked_se_results_rebuilds_bed_and_summary(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            work_dir = root / "work" / "se"
            final_dir = root / "final" / "se"
            work_dir.mkdir(parents=True)

            mask_file = root / "mask.bed"
            mask_file.write_text("chr1\t0\t100\n", encoding="utf-8")
            bed_mask = BedMask([mask_file], root)

            unmasked_anomalies = work_dir / "se.anomalous_regions.tsv"
            unmasked_anomalies.write_text(
                "Chromosome\tStart\tEnd\tAnomalyTags\n"
                "chr1\t1\t100\tlow_depth\n"
                "chr1\t101\t200\tclip\n",
                encoding="utf-8",
            )
            masked_anomalies = work_dir / "se.anomalous_regions_masked.tsv"
            bed_mask.apply_bedtools_intersect_v(
                unmasked_anomalies,
                masked_anomalies,
                coordinate_mode="one_based_closed",
            )

            (work_dir / "structural_errors.result.bed").write_text(
                "chr1\t0\t200\n", encoding="utf-8"
            )
            (work_dir / "structural_errors.anomaly_summary.tsv").write_text(
                "unmasked summary must not be published\n", encoding="utf-8"
            )

            pipeline = SASPipeline.__new__(SASPipeline)
            pipeline.logger = logging.getLogger("masked-se-publish-test")
            pipeline.se_output_dir = work_dir
            pipeline.se_results_dir = final_dir
            pipeline.se_anomalies_file = masked_anomalies
            pipeline.bed_mask = bed_mask
            pipeline.chromosomes = ["chr1"]

            pipeline._publish_se_results()

            self.assertEqual(
                (final_dir / "structural_errors.bed").read_text(encoding="utf-8"),
                "chr1\t100\t200\n",
            )
            self.assertEqual(
                (final_dir / "anomalous_regions.tsv").read_text(encoding="utf-8"),
                "Chromosome\tStart\tEnd\tAnomalyTags\nchr1\t101\t200\tclip\n",
            )
            self.assertEqual(
                (final_dir / "anomaly_summary.tsv").read_text(encoding="utf-8"),
                "Chromosome\tlow_depth\thigh_depth\tindel\tclip\t"
                "non-spanning\tnr-SE\n"
                "chr1\t0\t0\t0\t1\t0\t1\n"
                "Total\t0\t0\t0\t1\t0\t1\n",
            )

    def test_validate_input_files_checks_filter_script_once(self):
        config = {
            "short_read_bam": "short.bam",
            "hifi_bam": "hifi.bam",
            "ont_bam": "ont.bam",
            "reference_fa": "reference.fa",
            "scripts_dir": "scripts",
            "src_dir": "src",
            "run_be_analysis": True,
            "run_se_analysis": False,
            "skip_filter": False,
        }

        with mock.patch.object(Path, "exists", autospec=True, return_value=True) as exists:
            validate_input_files(config)

        checked_paths = [Path(call.args[0]) for call in exists.call_args_list]
        self.assertEqual(
            checked_paths.count(Path("scripts") / "filter_bam.sh"), 1
        )

    def test_se_only_unfiltered_run_still_requires_filter_script(self):
        config = {
            "short_read_bam": "short.bam",
            "hifi_bam": "hifi.bam",
            "ont_bam": "ont.bam",
            "reference_fa": "reference.fa",
            "scripts_dir": "scripts",
            "src_dir": "src",
            "run_be_analysis": False,
            "run_se_analysis": True,
            "skip_filter": False,
        }
        with mock.patch.object(Path, "exists", autospec=True, return_value=True) as exists:
            validate_input_files(config)
        checked_paths = [Path(call.args[0]) for call in exists.call_args_list]
        self.assertEqual(
            checked_paths.count(Path("scripts") / "filter_bam.sh"), 1
        )

    def test_se_only_run_filters_only_ont_bam(self):
        pipeline = SASPipeline.__new__(SASPipeline)
        pipeline.logger = mock.Mock(spec=logging.Logger)
        pipeline.skip_filter = False
        pipeline.run_be_analysis = False
        pipeline.scripts_dir = Path("scripts")
        pipeline.threads = 8
        pipeline.short_read_bam = Path("short.bam")
        pipeline.hifi_bam = Path("hifi.bam")
        pipeline.ont_bam = Path("ont.bam")
        pipeline.short_filtered_bam = Path("short.filtered.bam")
        pipeline.hifi_filtered_bam = Path("hifi.filtered.bam")
        pipeline.ont_filtered_bam = Path("ont.filtered.bam")
        pipeline._run_command = mock.Mock()

        pipeline.step1_filter_bam_files()

        pipeline._run_command.assert_called_once()
        command, description = pipeline._run_command.call_args.args
        self.assertEqual(description, "Filtering ONT BAM")
        self.assertIn("ont.bam", command)
        self.assertNotIn("short.bam", command)
        self.assertNotIn("hifi.bam", command)

    def test_runtime_preflight_is_mode_aware_and_reports_all_missing_tools(self):
        config = {
            "run_be_analysis": True,
            "run_se_analysis": False,
            "skip_filter": True,
            "generate_ir": True,
        }

        with mock.patch(
            "sas_pipeline.shutil.which",
            side_effect=lambda tool: None if tool in {"genmap", "mosdepth"} else tool,
        ):
            with self.assertRaisesRegex(
                FileNotFoundError, "genmap, mosdepth"
            ):
                validate_runtime_dependencies(config)

        config.update(run_be_analysis=False, run_se_analysis=True, generate_ir=False)
        observed = []
        with mock.patch(
            "sas_pipeline.shutil.which",
            side_effect=lambda tool: observed.append(tool) or tool,
        ):
            validate_runtime_dependencies(config)
        self.assertNotIn("bash", observed)
        self.assertNotIn("genmap", observed)

    def test_summary_counts_headerless_hifi_and_ont_rows(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            def write_lines(name, count, header=None):
                path = temp_path / name
                lines = ([] if header is None else [header]) + [
                    f"row-{index}" for index in range(count)
                ]
                path.write_text("\n".join(lines) + "\n", encoding="utf-8")
                return path

            pipeline = SASPipeline.__new__(SASPipeline)
            pipeline.output_dir = temp_path
            pipeline.logger = logging.getLogger("sas-pipeline-bookkeeping-test")
            pipeline.reference_fai = temp_path / "reference.fa.fai"
            pipeline.reference_fai.write_text("chr1\t1000\n", encoding="utf-8")
            pipeline.bed_mask = None
            pipeline.chromosomes = ["chr1"]
            pipeline.run_be_analysis = True
            pipeline.run_se_analysis = False
            pipeline.smallerror_bed = write_lines("smallerror.bed", 1)
            pipeline.hifi_filtered_regions = write_lines("hifi.txt", 2)
            pipeline.ont_filtered_regions = write_lines("ont.txt", 3)
            pipeline.merged_results = write_lines("merged.txt", 4, header="header")
            pipeline.final_be_results = write_lines("final.txt", 5)
            pipeline.threads = 2
            pipeline.skip_filter = True
            pipeline.config = {}
            pipeline.bed_files = []
            pipeline.additional_bed_file = None
            pipeline.short_read_bam = Path("short.bam")
            pipeline.hifi_bam = Path("hifi.bam")
            pipeline.ont_bam = Path("ont.bam")
            pipeline.reference_fa = Path("reference.fa")
            pipeline.short_filtered_bam = Path("short.filtered.bam")
            pipeline.hifi_filtered_bam = Path("hifi.filtered.bam")
            pipeline.ont_filtered_bam = Path("ont.filtered.bam")
            pipeline.be_output_dir = temp_path
            pipeline.se_output_dir = temp_path

            pipeline._generate_genome_wide_summary()

            summary = (temp_path / "sas_pipeline_summary.txt").read_text(
                encoding="utf-8"
            )
            self.assertIn("After HiFi filtering: 2", summary)
            self.assertIn("After ONT filtering: 3", summary)
            self.assertIn("After merging: 4", summary)

    def test_user_additional_bed_and_generated_ir_are_both_candidate_sources(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            pipeline = SASPipeline.__new__(SASPipeline)
            pipeline.be_output_dir = temp_path
            pipeline.additional_bed_file = temp_path / "user.bed"
            pipeline.additional_bed_windows = temp_path / "user.windows.bed"
            pipeline.ir_bed_file = temp_path / "IR.bed"
            pipeline.ir_bed_windows = temp_path / "IR.windows.bed"

            sources = pipeline._get_additional_candidate_sources()

            self.assertEqual(
                [(name, source.name, windows.name) for name, source, windows in sources],
                [
                    ("additional BED", "user.bed", "user.windows.bed"),
                    ("generated IR", "IR.bed", "IR.windows.bed"),
                ],
            )

    def test_window_split_is_reference_aligned_and_rejects_unknown_chromosomes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            pipeline = SASPipeline.__new__(SASPipeline)
            pipeline.logger = mock.Mock(spec=logging.Logger)
            pipeline._get_chromosome_lengths = mock.Mock(
                return_value={"chr1": 125}
            )
            input_bed = temp_path / "input.bed"
            output_bed = temp_path / "windows.bed"
            input_bed.write_text("chr1\t78\t130\n", encoding="utf-8")

            pipeline._split_bed_to_windows(input_bed, output_bed)
            self.assertEqual(
                output_bed.read_text(encoding="utf-8"),
                "chr1\t50\t100\nchr1\t100\t125\n",
            )

            input_bed.write_text("unknown\t0\t50\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "absent from the reference"):
                pipeline._split_bed_to_windows(input_bed, output_bed)
            self.assertEqual(
                output_bed.read_text(encoding="utf-8"),
                "chr1\t50\t100\nchr1\t100\t125\n",
            )

    def test_chromosome_qv_uses_unmasked_chromosome_length(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            mask_file = temp_path / "mask.bed"
            mask_file.write_text("chr1\t0\t200\n", encoding="utf-8")

            pipeline = SASPipeline.__new__(SASPipeline)
            pipeline.logger = logging.getLogger("masked-chromosome-qv-test")
            pipeline.output_dir = temp_path
            pipeline.chromosome_summary_file = temp_path / "chromosome_qv.txt"
            pipeline.chromosomes = ["chr1"]
            pipeline.bed_mask = BedMask([mask_file], temp_path)

            pipeline._calculate_chromosome_qv(
                {"chr1": 1000},
                {"chr1": 1},
                {},
            )

            self.assertEqual(pipeline.chr_qv_results["chr1"]["total_length"], 800)

    def test_reference_preflight_rejects_bam_length_mismatch(self):
        class FakeHeader:
            @staticmethod
            def to_dict():
                return {"HD": {"SO": "coordinate"}}

        class FakeBam:
            references = ("chr1",)
            lengths = (999,)
            header = FakeHeader()

            def __enter__(self):
                return self

            def __exit__(self, *_args):
                return False

            @staticmethod
            def has_index():
                return True

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            reference = temp_path / "reference.fa"
            reference.write_text(">chr1\nA\n", encoding="utf-8")
            fai = Path(str(reference) + ".fai")
            fai.write_text("chr1\t1000\t0\t1\t2\n", encoding="utf-8")

            pipeline = SASPipeline.__new__(SASPipeline)
            pipeline.logger = logging.getLogger("reference-preflight-test")
            pipeline.reference_fa = reference
            pipeline.reference_fai = fai
            pipeline.short_read_bam = temp_path / "short.bam"
            pipeline.hifi_bam = temp_path / "hifi.bam"
            pipeline.ont_bam = temp_path / "ont.bam"
            pipeline.skip_filter = True

            fake_pysam = mock.Mock()
            fake_pysam.AlignmentFile.return_value = FakeBam()
            with mock.patch("sas_pipeline.pysam", fake_pysam):
                with self.assertRaisesRegex(ValueError, "length mismatch=chr1"):
                    pipeline._validate_reference_compatibility()


if __name__ == "__main__":
    unittest.main()
