#!/usr/bin/env python3
"""
SAS Pipeline Setup Script
"""

from pathlib import Path

from setuptools import find_packages, setup
from setuptools.command.build_py import build_py as _build_py

from src._version import __version__


PROJECT_ROOT = Path(__file__).resolve().parent
RUNTIME_ASSET_PATTERNS = (
    "scripts/*.py",
    "scripts/*.sh",
    "config/*.json",
)


class build_py(_build_py):
    """Place runtime scripts beside ``src`` in installed distributions.

    The production modules intentionally locate helpers relative to their own
    file (``../scripts`` and ``../config``).  Mirroring the repository layout
    in ``build_lib`` preserves that behaviour without duplicating or moving the
    scripts in the source tree.
    """

    def _runtime_assets(self):
        for pattern in RUNTIME_ASSET_PATTERNS:
            yield from sorted(PROJECT_ROOT.glob(pattern))

    def run(self):
        super().run()
        for source in self._runtime_assets():
            destination = Path(self.build_lib) / source.relative_to(PROJECT_ROOT)
            self.mkpath(str(destination.parent))
            self.copy_file(str(source), str(destination))

    def get_outputs(self, include_bytecode=1):
        outputs = super().get_outputs(include_bytecode=include_bytecode)
        outputs.extend(
            str(Path(self.build_lib) / source.relative_to(PROJECT_ROOT))
            for source in self._runtime_assets()
        )
        return outputs

# Read the README file
with open(PROJECT_ROOT / "README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements
with open(PROJECT_ROOT / "requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="sas-pipeline",
    version=__version__,
    author="KANGYUlab",
    author_email="",
    description="SAS: an alignment-based metric for assembly errors",
    license="MIT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KANGYUlab/sas-pipeline",
    packages=find_packages(include=("src", "src.*")),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "pytest-cov>=2.10.0",
            "black>=21.0.0",
            "flake8>=3.8.0",
        ],
        "viz": [
            "matplotlib>=3.1.0",
            "seaborn>=0.9.0",
            "plotly>=4.0.0",
        ],
        "all": [
            "matplotlib>=3.1.0",
            "seaborn>=0.9.0", 
            "plotly>=4.0.0",
            "pytest>=6.0.0",
            "pytest-cov>=2.10.0",
            "black>=21.0.0",
            "flake8>=3.8.0",
        ]
    },
    entry_points={
        "console_scripts": [
            "sas-pipeline=src.sas_pipeline:main",
            "sas-se-analysis=src.integrated_se_analysis:main",
        ],
    },
    cmdclass={"build_py": build_py},
    zip_safe=False,
    keywords="bioinformatics sequencing analysis structural-variants base-errors",
    project_urls={
        "Bug Reports": "https://github.com/KANGYUlab/sas-pipeline/issues",
        "Source": "https://github.com/KANGYUlab/sas-pipeline",
        "Documentation": "https://github.com/KANGYUlab/sas-pipeline/wiki",
    },
)
