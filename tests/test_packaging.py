import ast
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SETUP_PATH = PROJECT_ROOT / "setup.py"


def setup_keyword_values():
    """Return literal setup() keyword values without executing setuptools."""
    tree = ast.parse(SETUP_PATH.read_text(encoding="utf-8"), filename=str(SETUP_PATH))
    setup_call = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "setup"
    )
    return {
        keyword.arg: ast.literal_eval(keyword.value)
        for keyword in setup_call.keywords
        if keyword.arg in {"entry_points", "python_requires"}
    }


class PackagingStructureTests(unittest.TestCase):
    def test_readme_image_is_included_in_source_distributions(self):
        manifest = (PROJECT_ROOT / "MANIFEST.in").read_text(encoding="utf-8")
        self.assertIn("recursive-include docs/images", manifest)
        self.assertTrue((PROJECT_ROOT / "docs" / "images" / "sas-principle.jpg").is_file())

    def test_console_entry_modules_and_callables_exist(self):
        keywords = setup_keyword_values()
        entries = keywords["entry_points"]["console_scripts"]

        self.assertIn("sas-pipeline=src.sas_pipeline:main", entries)
        self.assertIn("sas-se-analysis=src.integrated_se_analysis:main", entries)
        self.assertTrue((PROJECT_ROOT / "src" / "__init__.py").is_file())

        for entry in entries:
            _, target = entry.split("=", 1)
            module_name, callable_name = target.split(":", 1)
            module_path = PROJECT_ROOT.joinpath(*module_name.split(".")).with_suffix(".py")
            self.assertTrue(module_path.is_file(), entry)
            tree = ast.parse(module_path.read_text(encoding="utf-8"), filename=str(module_path))
            functions = {
                node.name for node in tree.body if isinstance(node, ast.FunctionDef)
            }
            self.assertIn(callable_name, functions, entry)

    def test_supported_python_floor_is_declared(self):
        self.assertEqual(setup_keyword_values()["python_requires"], ">=3.8")

    def test_setup_metadata_uses_the_package_version(self):
        from src import __version__

        tree = ast.parse(SETUP_PATH.read_text(encoding="utf-8"))
        setup_call = next(
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "setup"
        )
        version_value = next(
            keyword.value
            for keyword in setup_call.keywords
            if keyword.arg == "version"
        )
        self.assertIsInstance(version_value, ast.Name)
        self.assertEqual(version_value.id, "__version__")
        self.assertEqual(__version__, "2.0.0")

    def test_se_cli_version_does_not_require_analysis_dependencies(self):
        result = subprocess.run(
            [sys.executable, "src/integrated_se_analysis.py", "--version"],
            cwd=str(PROJECT_ROOT),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        self.assertEqual(result.stdout.strip(), "integrated_se_analysis.py 2.0.0")

    def test_visualization_dependencies_are_not_core_requirements(self):
        requirements = {
            line.strip().split(";", 1)[0]
            for line in (PROJECT_ROOT / "requirements.txt").read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }
        for package in ("matplotlib", "seaborn", "plotly"):
            self.assertFalse(
                any(requirement.lower().startswith(package) for requirement in requirements)
            )

    def test_simulated_install_resolves_entry_points_and_runtime_assets(self):
        """Exercise imports using the same sibling layout produced by build_py."""
        with tempfile.TemporaryDirectory() as temp_dir:
            site_dir = Path(temp_dir) / "site"
            shutil.copytree(PROJECT_ROOT / "src", site_dir / "src")
            shutil.copytree(PROJECT_ROOT / "scripts", site_dir / "scripts")
            shutil.copytree(PROJECT_ROOT / "config", site_dir / "config")

            stubs_dir = Path(temp_dir) / "stubs"
            stubs_dir.mkdir()
            (stubs_dir / "numpy.py").write_text("", encoding="utf-8")
            (stubs_dir / "pandas.py").write_text("", encoding="utf-8")
            (stubs_dir / "pysam.py").write_text(
                "class AlignmentFile:\n    pass\n", encoding="utf-8"
            )

            command = """
from pathlib import Path
from src import integrated_se_analysis, sas_pipeline

base = sas_pipeline.get_pipeline_base_dir()
assert base == Path(sas_pipeline.__file__).resolve().parent.parent
assert callable(sas_pipeline.main)
assert callable(integrated_se_analysis.main)
assert (base / 'scripts' / 'smallerror_long.py').is_file()
assert (base / 'scripts' / 'spanning.py').is_file()
assert (base / 'config' / 'config_template.json').is_file()
assert Path(integrated_se_analysis.__file__).resolve().parent.parent == base
"""
            environment = os.environ.copy()
            environment["PYTHONPATH"] = os.pathsep.join(
                (str(stubs_dir), str(site_dir))
            )
            subprocess.run(
                [sys.executable, "-c", command],
                cwd=temp_dir,
                env=environment,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )

    @unittest.skipUnless(importlib.util.find_spec("setuptools"), "setuptools unavailable")
    def test_build_py_preserves_runtime_layout(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            build_dir = Path(temp_dir) / "lib"
            subprocess.run(
                [
                    sys.executable,
                    "setup.py",
                    "build_py",
                    "--build-lib",
                    str(build_dir),
                ],
                cwd=str(PROJECT_ROOT),
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )

            self.assertTrue((build_dir / "src" / "sas_pipeline.py").is_file())
            self.assertTrue((build_dir / "src" / "integrated_se_analysis.py").is_file())
            for source in list((PROJECT_ROOT / "scripts").glob("*.py")) + list(
                (PROJECT_ROOT / "scripts").glob("*.sh")
            ):
                self.assertEqual(
                    (build_dir / "scripts" / source.name).read_bytes(),
                    source.read_bytes(),
                )
            self.assertEqual(
                (build_dir / "config" / "config_template.json").read_bytes(),
                (PROJECT_ROOT / "config" / "config_template.json").read_bytes(),
            )


if __name__ == "__main__":
    unittest.main()
