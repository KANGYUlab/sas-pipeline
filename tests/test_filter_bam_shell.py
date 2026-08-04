import os
from pathlib import Path
import subprocess
import tempfile
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "filter_bam.sh"


@unittest.skipUnless(os.name == "posix", "Bash helper runtime requires POSIX")
class FilterBamShellTests(unittest.TestCase):
    def run_script(self, *args, env=None):
        return subprocess.run(
            ["bash", str(SCRIPT), *args],
            cwd=PROJECT_ROOT,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def test_shell_syntax_is_valid(self):
        result = subprocess.run(
            ["bash", "-n", str(SCRIPT)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_source_uses_argument_arrays_and_not_eval_or_nproc(self):
        source = SCRIPT.read_text()
        self.assertNotIn("eval ", source)
        self.assertNotIn("$(nproc)", source)
        self.assertIn('samtools "${samtools_args[@]}" "$input_bam"', source)
        self.assertIn("getconf _NPROCESSORS_ONLN", source)
        self.assertIn("sysctl -n hw.ncpu", source)

    def test_sysctl_is_used_when_getconf_does_not_return_a_cpu_count(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            getconf = temp_path / "getconf"
            getconf.write_text("#!/bin/sh\nexit 1\n")
            getconf.chmod(0o755)
            sysctl = temp_path / "sysctl"
            sysctl.write_text("#!/bin/sh\nprintf '7\\n'\n")
            sysctl.chmod(0o755)

            env = os.environ.copy()
            env["PATH"] = f"{temp_path}:{env.get('PATH', '')}"
            result = self.run_script("short-read", "--help", env=env)

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("samtools threads (default: 7)", result.stdout)

    def test_threads_must_be_a_positive_integer(self):
        result = self.run_script(
            "short-read", "-i", "missing.bam", "-o", "output.bam", "-t", "0"
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Threads must be a positive integer", result.stderr)

    def test_split_parameters_must_be_non_negative_numbers(self):
        for option, value in (
            ("--max-splits-base", "not-a-number"),
            ("--max-splits-kb", "-0.1"),
        ):
            with self.subTest(option=option):
                result = self.run_script(
                    "long-read",
                    "-i",
                    "missing.bam",
                    "-o",
                    "output.bam",
                    option,
                    value,
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("must be a non-negative number", result.stderr)

    def test_empty_long_read_bam_reports_zero_percent(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            input_bam = temp_path / "empty.bam"
            output_bam = temp_path / "filtered.bam"
            input_bam.touch()

            samtools = temp_path / "samtools"
            samtools.write_text(
                """#!/bin/bash
set -euo pipefail

command_name="$1"
shift

case "$command_name" in
    view)
        for argument in "$@"; do
            if [[ "$argument" == "-c" ]]; then
                printf '0\\n'
                exit 0
            fi
        done

        output=""
        binary_output=false
        while [[ $# -gt 0 ]]; do
            case "$1" in
                -b)
                    binary_output=true
                    shift
                    ;;
                -o)
                    output="$2"
                    shift 2
                    ;;
                *)
                    shift
                    ;;
            esac
        done

        if [[ "$binary_output" == true ]]; then
            : > "$output"
            cat >/dev/null
        else
            printf '@HD\\tVN:1.6\\n'
        fi
        ;;
    index)
        exit 0
        ;;
    *)
        printf 'unexpected samtools command: %s\\n' "$command_name" >&2
        exit 2
        ;;
esac
"""
            )
            samtools.chmod(0o755)

            env = os.environ.copy()
            env["PATH"] = f"{temp_path}:{env.get('PATH', '')}"
            result = self.run_script(
                "long-read",
                "-i",
                str(input_bam),
                "-o",
                str(output_bam),
                "-t",
                "1",
                env=env,
            )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Total primary reads: 0", result.stdout)
        self.assertIn("Kept reads: 0 (0.0%)", result.stdout)
        self.assertIn("Filtered reads: 0", result.stdout)
        self.assertNotIn("nan", result.stdout.lower())


if __name__ == "__main__":
    unittest.main()
