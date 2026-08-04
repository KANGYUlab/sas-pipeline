import os
from pathlib import Path
import subprocess
import sys
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class CLIHelpTests(unittest.TestCase):
    PYTHON_COMMANDS = {
        "src/sas_pipeline.py": (
            "required inputs:",
            "-a BED, --additional-bed-file BED",
            "--pileup-workers N",
        ),
        "src/integrated_se_analysis.py": (
            "BAM FASTA PREFIX",
            "-d DIR, --output-dir DIR",
            "--workers N",
        ),
        "scripts/smallerror_long.py": ("BED BAM FASTA OUTPUT", "mpileup batching options:"),
        "scripts/merge_pileup.py": ("BED HIFI_DATA ONT_DATA FASTA OUTPUT", "Output format:"),
        "scripts/filter_be.py": ("DATA FAI", "-s SEP, --sep SEP"),
        "scripts/spanning.py": ("BAM FAI", "-o PREFIX, --output PREFIX"),
    }

    SHELL_COMMANDS = {
        "scripts/filter_bam.sh": ("MODE -i INPUT_BAM -o OUTPUT_BAM", "Long-read options:"),
        "scripts/smallerror_short.sh": ("-i INPUT_BAM -o OUTPUT_PREFIX", "--threads N"),
    }

    def run_python(self, relative_path, *args):
        return subprocess.run(
            [sys.executable, str(PROJECT_ROOT / relative_path), *args],
            cwd=PROJECT_ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def run_shell(self, relative_path, *args):
        return subprocess.run(
            ["bash", str(PROJECT_ROOT / relative_path), *args],
            cwd=PROJECT_ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def test_every_python_entry_has_clean_help_and_version(self):
        for path, expected_fragments in self.PYTHON_COMMANDS.items():
            with self.subTest(path=path):
                result = self.run_python(path, "--help")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stderr, "")
                self.assertIn("usage:", result.stdout)
                self.assertIn("-h, --help", result.stdout)
                self.assertNotIn("default: None", result.stdout)
                self.assertNotIn("default: False", result.stdout)
                for fragment in expected_fragments:
                    self.assertIn(fragment, result.stdout)

                version = self.run_python(path, "--version")
                self.assertEqual(version.returncode, 0, version.stderr)
                self.assertTrue(version.stdout.strip().endswith("2.0.0"))

    @unittest.skipUnless(os.name == "posix", "Bash helpers require POSIX")
    def test_every_shell_entry_has_clean_help_and_version(self):
        for path, expected_fragments in self.SHELL_COMMANDS.items():
            with self.subTest(path=path):
                result = self.run_shell(path, "--help")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stderr, "")
                self.assertIn("Usage:", result.stdout)
                self.assertIn("-h, --help", result.stdout)
                for fragment in expected_fragments:
                    self.assertIn(fragment, result.stdout)

                version = self.run_shell(path, "--version")
                self.assertEqual(version.returncode, 0, version.stderr)
                self.assertTrue(version.stdout.strip().endswith("2.0.0"))

    def test_python_numeric_options_reject_zero_before_running_analysis(self):
        cases = (
            (
                "src/sas_pipeline.py",
                "must be greater than 0",
                ("-i", "a", "-H", "b", "-N", "c", "-r", "d", "-o", "e", "-t", "0"),
            ),
            (
                "src/integrated_se_analysis.py",
                "must be greater than 0",
                ("a.bam", "ref.fa", "sample", "-w", "0"),
            ),
            (
                "scripts/smallerror_long.py",
                "must be greater than 0",
                ("a.bed", "a.bam", "ref.fa", "out", "-p", "0"),
            ),
            (
                "scripts/spanning.py",
                "must be greater than 0",
                ("a.bam", "ref.fa.fai", "-o", "out", "-w", "0"),
            ),
        )
        for path, message, arguments in cases:
            with self.subTest(path=path):
                result = self.run_python(path, *arguments)
                self.assertEqual(result.returncode, 2)
                self.assertIn(message, result.stderr)
                self.assertNotIn("Traceback", result.stderr)

    @unittest.skipUnless(os.name == "posix", "Bash helpers require POSIX")
    def test_shell_missing_option_values_have_actionable_errors(self):
        cases = (
            ("scripts/filter_bam.sh", ("short-read", "-i"), "-i requires a value"),
            ("scripts/smallerror_short.sh", ("--output",), "--output requires a value"),
        )
        for path, arguments, message in cases:
            with self.subTest(path=path):
                result = self.run_shell(path, *arguments)
                self.assertEqual(result.returncode, 2)
                self.assertIn(message, result.stderr)
                self.assertNotIn("unbound variable", result.stderr)

    def test_readme_documents_current_commands_and_short_aliases(self):
        readme = (PROJECT_ROOT / "README.md").read_text(encoding="utf-8")
        for command in (
            "sas-pipeline",
            "sas-se-analysis",
            "filter_bam.sh",
            "smallerror_short.sh",
            "smallerror_long.py",
            "merge_pileup.py",
            "filter_be.py",
            "spanning.py",
        ):
            with self.subTest(command=command):
                self.assertIn(command, readme)

        for option in (
            "`-i` | `--short-read-bam`",
            "`-H` | `--hifi-bam`",
            "`-N` | `--ont-bam`",
            "`-r` | `--reference-fa`",
            "`-o` | `--output-dir`",
            "`-a` | `--additional-bed-file`",
            "`-b` | `--bed-files`",
            "`-g` | `--generate-ir`",
        ):
            with self.subTest(option=option):
                self.assertIn(option, readme)

        for output_path in (
            "final_be_results.txt",
            "structural_errors.bed",
            "smallerror_short.per-base.bed.gz",
            "structural_errors.window_stats.tsv",
            "se.anomaly_summary.tsv",
            "logs/",
            "metadata/",
        ):
            with self.subTest(output_path=output_path):
                self.assertIn(output_path, readme)


if __name__ == "__main__":
    unittest.main()
