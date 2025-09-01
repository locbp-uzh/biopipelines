#!/usr/bin/env python3

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="biopipelines",
    version="0.1.0",
    author="LOCBP - University of Zürich",
    description="Computational protein design pipelines for structural biology",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.uzh.ch/locbp/public/biopipelines",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.9",
    install_requires=[
        "pandas",
        "numpy",
        "biopython",
        "PyYAML",
    ],
    entry_points={
        "console_scripts": [
            "biopipeline=PipelineScripts.pipeline:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.yml", "*.yaml", "*.md"],
    },
)