import sys
import os
from setuptools import setup, find_packages
import glob

_MAJOR = 0
_MINOR = 14
_MICRO = 5
version = "%d.%d.%d" % (_MAJOR, _MINOR, _MICRO)
release = "%d.%d" % (_MAJOR, _MINOR)


metainfo = {
    "authors": {"main": ("Thomas Cokelaer", "thomas.cokelaer@pasteur.fr")},
    "maintainer": {"main": ("Thomas Cokelaer", "thomas.cokelaer@pasteur.fr")},
    "version": version,
    "license": "new BSD",
    "download_url": "https://github.com/sequana/sequana/archive/{0}.tar.gz".format(version),
    "url": "http://github.com/sequana/sequana",
    "description": "A set of standalone application and pipelines dedicated to NGS (new generation sequencing) analysis",
    "platforms": ["Linux", "Unix", "MacOsX", "Windows"],
    "keywords": ["NGS", "snakemake"],
    "classifiers": [
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
    ],
}


packages = find_packages()
packages = [this for this in packages if this.startswith("test.") is False]
packages = [this for this in packages if this not in ["test"]]

# load a common list of requirements
requirements = open("requirements.txt").read().split()

on_rtd = os.environ.get("READTHEDOCS", None) == "True"
if on_rtd:
    # pillow, sphinx, numpydoc are  for the doc only
    extra_packages = ["pillow", "numpydoc", "sphinx"]
    requirements += extra_packages
    # replace gseapy by gseapy<0.11 to make sure RTD works
    # indeed gseapy requires rust, not installed by default on RTD.
    gseapy =  [i for i,x in enumerate(requirements) if x.startswith('gseapy')]
    if gseapy: # for future, if dropped
        requirements[gseapy[0]] = 'gseapy<0.11'



if sys.version_info.major == 2 or on_rtd:
    requirements = [x for x in requirements if x.startswith("snakemake") is False]


setup(
    name="sequana",
    version=version,
    maintainer=metainfo["authors"]["main"][0],
    maintainer_email=metainfo["authors"]["main"][1],
    author=metainfo["authors"]["main"][0],
    author_email=metainfo["authors"]["main"][1],
    long_description=open("README.rst").read(),
    keywords=metainfo["keywords"],
    description=metainfo["description"],
    license=metainfo["license"],
    platforms=metainfo["platforms"],
    url=metainfo["url"],
    download_url=metainfo["download_url"],
    classifiers=metainfo["classifiers"],
    # package installation
    packages=packages,
    # pillow, sphinx-gallery and numpydoc are  for the doc only
    # mock is for the test only qtconsole is required by Sequanix
    install_requires=requirements,
    # specific packages for testing
    extras_require={
        "testing": [
            "pytest",
            "pytest-cov",
            "pytest-xdist",
            "pytest-mock",
            "pytest-timeout",
            "pytest-runner",
            "coveralls",
        ],
        "pipelines": [
            "sequana_quality_control",
            "sequana_fastqc",
            "sequana_rnaseq",
            "sequana_variant_calling",
            "sequana_fastqc",
            "sequana_demultiplex",
            "sequana_coverage",
            "sequana_denovo",
            "sequana_mapper",
        ],
        "doc": [
            "sphinx>=3",
            "sphinx_rtd_theme",
            "sequana_sphinxext",
        ],
    },
    # This is recursive include of data files
    exclude_package_data={"": ["__pycache__"]},
    package_data={
        "": [
            "Snakefile*",
            "*html",
            "README.rst",
            "requirements*txt",
            "config*.yaml",
            "*.css",
            "*.js",
            "snpEff.config*",
            "*.fa",
            "*.rules",
        ],
        "sequana.rules": ["*/*.rules", "*/*/*.rules", "*/*/*/*.rules"],
        "sequana.pipelines": ["*/*"],
        "sequana.resources.data": ["*.*"],  # use *.* for files and not ./adapters
        "sequana.resources.examples": ["*"],
        "sequana.resources.templates": ["*.R"],
        "sequana.resources.images": ["*"],
        "sequana.resources.testing": ["*", "*/*", "*/*/*", "*/*/*/*", "*/*/*/*/*", "*/*/*/*/*/*"],
        "sequana.multiqc": ["*yaml"],
    },
    # these files do not need to be added in MANIFEST.in since there are python
    # packages that will be copied from sequana/ into sequana/
    # Note, however, that e.g. ./pipelines must be added
    zip_safe=False,
    entry_points={
        "console_scripts": [
            "sequana=sequana.scripts.main.main:main",
            "sequana_taxonomy=sequana.scripts.taxonomy:main",
            "sequana_coverage=sequana.scripts.coverage:main",
            "sequana_substractor=sequana.scripts.substractor:main",  # march 2019
        ],
        "sequana.module": [
            "sequana_coverage=sequana.modules_report.coverage:CoverageModule",
            "sequana_variant_calling=sequana.modules_report.variant_calling:VariantCallingModule",
            "sequana_summary=sequana.modules_report.summary:SummaryModule",
        ],
        "multiqc.modules.v1": [
            "sequana_pacbio_qc=sequana.multiqc.pacbio_qc:MultiqcModule",
            "sequana_quality_control=sequana.multiqc.quality_control:MultiqcModule",
            "sequana_coverage=sequana.multiqc.coverage:MultiqcModule",
            "sequana_isoseq=sequana.multiqc.isoseq:MultiqcModule",
            "sequana_isoseq_qc=sequana.multiqc.isoseq_qc:MultiqcModule",
            "sequana_bamtools_stats=sequana.multiqc.bamtools_stats:MultiqcModule",
            "sequana_laa=sequana.multiqc.laa:MultiqcModule",
            "sequana_kraken=sequana.multiqc.kraken:MultiqcModule",
            "sequana_pacbio_amplicon=sequana.multiqc.pacbio_amplicon:MultiqcModule",
        ],
        "multiqc.hooks.v1": [
            "before_config = sequana.multiqc.config:load_config",
        ],
    },
)
