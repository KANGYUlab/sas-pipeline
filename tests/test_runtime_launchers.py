import ast
import logging
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SAS_PIPELINE_PATH = PROJECT_ROOT / "src" / "sas_pipeline.py"
SMALLERROR_SHORT = PROJECT_ROOT / "scripts" / "smallerror_short.sh"
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from sas_pipeline import SASPipeline


def _list_string_constants(node):
    return {
        child.value
        for child in ast.walk(node)
        if isinstance(child, ast.Constant) and isinstance(child.value, str)
    }


@unittest.skipUnless(os.name == "posix", "launcher runtime tests require POSIX")
class RuntimeLauncherTests(unittest.TestCase):
    def test_packaged_helpers_use_explicit_runtime_launchers(self):
        tree = ast.parse(SAS_PIPELINE_PATH.read_text(encoding="utf-8"))
        expected_python_helpers = {
            "smallerror_long.py",
            "merge_pileup.py",
            "filter_be.py",
            "integrated_se_analysis.py",
        }
        expected_shell_helpers = {"filter_bam.sh", "smallerror_short.sh"}
        observed_python = set()
        observed_shell = set()

        for assignment in ast.walk(tree):
            if not isinstance(assignment, ast.Assign):
                continue
            if not any(
                isinstance(target, ast.Name) and target.id == "cmd"
                for target in assignment.targets
            ):
                continue
            node = assignment.value
            if not isinstance(node, ast.List) or not node.elts:
                continue
            strings = _list_string_constants(node)
            python_helpers = expected_python_helpers & strings
            shell_helpers = expected_shell_helpers & strings
            if python_helpers:
                first = node.elts[0]
                self.assertIsInstance(first, ast.Attribute)
                self.assertIsInstance(first.value, ast.Name)
                self.assertEqual((first.value.id, first.attr), ("sys", "executable"))
                observed_python.update(python_helpers)
            if shell_helpers:
                first = node.elts[0]
                self.assertIsInstance(first, ast.Constant)
                self.assertEqual(first.value, "bash")
                observed_shell.update(shell_helpers)

        self.assertEqual(observed_python, expected_python_helpers)
        self.assertEqual(observed_shell, expected_shell_helpers)

    def test_pipeline_can_run_a_non_executable_shell_helper_via_bash(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            helper = Path(temp_dir) / "helper.sh"
            marker = Path(temp_dir) / "marker.txt"
            helper.write_text(
                '#!/bin/bash\nprintf "launched" > "$1"\n', encoding="utf-8"
            )
            helper.chmod(stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH)
            self.assertEqual(helper.stat().st_mode & 0o777, 0o644)

            pipeline = SASPipeline.__new__(SASPipeline)
            pipeline.logger = logging.getLogger("runtime-launcher-test")
            pipeline._run_command(
                ["bash", str(helper), str(marker)], "running packaged shell helper"
            )

            self.assertEqual(marker.read_text(encoding="utf-8"), "launched")

    def test_smallerror_short_help_succeeds_but_usage_errors_fail(self):
        help_result = subprocess.run(
            ["bash", str(SMALLERROR_SHORT), "-h"],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(help_result.returncode, 0)
        self.assertIn("Usage:", help_result.stdout)

        missing_result = subprocess.run(
            ["bash", str(SMALLERROR_SHORT)],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertNotEqual(missing_result.returncode, 0)
        self.assertIn("both --input and --output are required", missing_result.stderr)
        self.assertIn("Usage:", missing_result.stderr)

    def test_smallerror_short_uses_portable_gzip_pipeline(self):
        script_text = SMALLERROR_SHORT.read_text(encoding="utf-8")
        self.assertNotIn("zcat", script_text)
        self.assertIn('gzip -cd -- "$PER_BASE_BED_GZ" | awk', script_text)

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            bin_dir = temp_path / "bin"
            output_prefix = temp_path / "output" / "sample"
            input_bam = temp_path / "input.bam"
            zcat_marker = temp_path / "zcat-called"
            bin_dir.mkdir()
            input_bam.touch()

            helpers = {
                "samtools": """#!/bin/bash
set -euo pipefail
if [[ "$1" == "view" && "$2" == "-H" ]]; then
    printf '@SQ\\tSN:chr1\\tLN:100\\n'
else
    exit 2
fi
""",
                "mosdepth": """#!/bin/bash
set -euo pipefail
prefix="$3"
printf 'chrom length bases mean\\nchr1 100 100 20\\n' > "${prefix}.mosdepth.summary.txt"
printf 'chr1\\t0\\t10\\t0\\nchr1\\t10\\t100\\t20\\n' | gzip -c > "${prefix}.per-base.bed.gz"
""",
                "bedtools": """#!/bin/bash
set -euo pipefail
case "$1" in
    makewindows)
        printf 'chr1\\t0\\t50\\nchr1\\t50\\t100\\n'
        ;;
    intersect)
        low_coverage_file=""
        while (($#)); do
            if [[ "$1" == "-b" ]]; then
                low_coverage_file="$2"
                break
            fi
            shift
        done
        [[ -n "$low_coverage_file" ]]
        grep -q $'chr1\\t0\\t10' "$low_coverage_file"
        printf 'chr1\\t0\\t50\\n'
        ;;
    *)
        exit 2
        ;;
esac
""",
                "zcat": f"""#!/bin/bash
printf called > {zcat_marker!s}
exit 99
""",
            }
            for name, contents in helpers.items():
                helper = bin_dir / name
                helper.write_text(contents, encoding="utf-8")
                helper.chmod(0o755)

            env = os.environ.copy()
            env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"
            result = subprocess.run(
                [
                    "bash",
                    str(SMALLERROR_SHORT),
                    "-i",
                    str(input_bam),
                    "-o",
                    str(output_prefix),
                ],
                text=True,
                capture_output=True,
                check=False,
                env=env,
            )

            self.assertEqual(result.returncode, 0, result.stderr or result.stdout)
            self.assertFalse(zcat_marker.exists())
            self.assertEqual(
                output_prefix.with_suffix(".low_coverage_bases.bed").read_text(
                    encoding="utf-8"
                ),
                "chr1\t0\t10\n",
            )
            self.assertEqual(
                output_prefix.with_suffix(".smallerror.bed").read_text(
                    encoding="utf-8"
                ),
                "chr1\t0\t50\n",
            )


if __name__ == "__main__":
    unittest.main()
