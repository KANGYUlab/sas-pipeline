#!/usr/bin/env python3
"""
SAS Pipeline Setup Script
"""

from setuptools import setup, find_packages
import os

# Read the README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="sas-pipeline",
    version="2.0.0",
    author="KANGYUlab",
    author_email="",
    description="SAS: an alignment-based metric for assembly errors",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KANGYUlab/sas-pipeline",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.6",
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
            "sas-se-analysis=src.integrated_sv_analysis:main",
        ],
    },
    include_package_data=True,
    package_data={
        "sas_pipeline": [
            "scripts/*.sh",
            "scripts/*.py",
            "config/*.json",
            "examples/*.md",
        ],
    },
    zip_safe=False,
    keywords="bioinformatics sequencing analysis structural-variants base-errors",
    project_urls={
        "Bug Reports": "https://github.com/KANGYUlab/sas-pipeline/issues",
        "Source": "https://github.com/KANGYUlab/sas-pipeline",
        "Documentation": "https://github.com/KANGYUlab/sas-pipeline/wiki",
    },
) 