from pathlib import Path

from setuptools import find_packages, setup


ROOT = Path(__file__).parent
README = ROOT / "README.md"


setup(
    name="IRTools",
    version="2.0",
    description="A computational toolset for detection and analysis of intron retention from RNA-Seq libraries.",
    long_description=README.read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    author="Zhouhao Zeng, Andy Cao",
    author_email="zzhlbj23@gwmail.gwu.edu, acao3@email.gwu.edu",
    url="https://github.com/WeiqunPengLab/IRTools/",
    packages=find_packages(include=["IRTools", "IRTools.*"]),
    package_data={"IRTools": ["data/*.gtf.gz", "data/*.gff3.gz"]},
    include_package_data=True,
    scripts=["bin/IRTools"],
    python_requires=">=3.11",
    install_requires=[
        "numpy>=2.4",
        "scipy>=1.17",
        "pandas>=3.0",
        "networkx>=3.6",
        "bx-python>=0.14",
        "HTSeq>=2.1",
        "pysam>=0.23",
        "statsmodels>=0.14",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
    ],
)
