# -*- coding: utf-8 -*-
# License: 3-clause BSD
__revision__ = "$Id: $" # for the SVN Id
import sys
import os
from setuptools import setup, find_packages
import glob

_MAJOR               = 0
_MINOR               = 7
_MICRO               = 2
version              = '%d.%d.%d' % (_MAJOR, _MINOR, _MICRO)
release              = '%d.%d' % (_MAJOR, _MINOR)

#version += ".post3"

metainfo = {
    'authors': {"main": ("yourname", "email@whatever.org")},
    'version': version,
    'license' : 'new BSD',
    'download_url' : ['http://pypi.python.org/pypi/sequana'],
    'url' : ["http://github.com/sequana/"],
    'description': "A set of standalone application and pipelines dedicated to NGS (new generation sequencing) analysis" ,
    'platforms' : ['Linux', 'Unix', 'MacOsX', 'Windows'],
    'keywords' : [''],
    'classifiers' : [
          'Development Status :: 4 - Beta',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Topic :: Software Development :: Libraries :: Python Modules',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Information Analysis',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Physics']
    }


packages = find_packages()
packages = [this for this in packages if this.startswith('test.') is False]
packages = [this for this in packages if this not in ['test']]

# load a common list of requirements
# - mock is for the test only
# - qtconsole is required by Sequanix
requirements = open("requirements.txt").read().split()

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd:
    # pillow, sphinx, numpydoc are  for the doc only
    extra_packages = ["pillow", "numpydoc", "sphinx"]
    requirements += extra_packages


if sys.version_info.major == 2 or on_rtd:
    requirements = [x for x in requirements 
                    if x.startswith("snakemake") is False]



setup(
    name             = "sequana",
    version          = version,
    maintainer       = metainfo['authors']['main'][0],
    maintainer_email = metainfo['authors']['main'][1],
    author           = metainfo['authors']['main'][0],
    author_email     = metainfo['authors']['main'][1],
    long_description = open("README.rst").read(),
    keywords         = metainfo['keywords'],
    description      = metainfo['description'],
    license          = metainfo['license'],
    platforms        = metainfo['platforms'],
    url              = metainfo['url'],
    download_url     = metainfo['download_url'],
    classifiers      = metainfo['classifiers'],

    # package installation
    packages = packages,

    # pillow, sphinx-gallery and numpydoc are  for the doc only
    # mock is for the test only qtconsole is required by Sequanix
    install_requires = requirements,

    # specific packages for testing
    tests_require = open('requirements_dev.txt').read().split(),

    # here below '': pattern means include that pattern in all packages
    # so '' :['README.rst'] will include all README.rst recursively
    # required to use python setup.py install

    # This is recursive include of data files
    exclude_package_data = {"": ["__pycache__"]},
    package_data = {
        '': ['Snakefile*', '*html', 'README.rst', "requirements*txt",
             'config*.yaml', '*.css', "*.js",
             "snpEff.config*", "*.fa", "*.rules"],
        'sequana.rules' : ['*/*.rules', "*/*/*.rules"],
        'sequana.pipelines' : ['*/*'],
        'sequana.resources.data' : ['*.*'],  # use *.* for files and not ./adapters
        'sequana.resources.data.adapters' : ['*'],
        'sequana.resources.images' : ['*'],
        'sequana.resources.testing' : ['*'],
        'sequana.resources.busco' : ['*'],
        'sequana.multiqc' : ['*yaml'],
        },

    # these files do not need to be added in MANIFEST.in since there are python
    # packages that will be copied from sequana/ into sequana/
    # Note, however, that e.g. ./pipelines must be added

    zip_safe=False,
    entry_points = {
        'console_scripts':[
           'sequana_gui=sequana.gui.sequana_gui:main',
           'sequanix=sequana.gui.sequana_gui:main',
           'fastq_head=sequana.scripts.fastq_head:main',
           'fastq_count=sequana.scripts.fastq_count:main',
           'sequana_fastq_head=sequana.scripts.fastq_head:main',
           'sequana_fastq_count=sequana.scripts.fastq_count:main',
           'sequana=sequana.scripts.main:main',
           'sequana_taxonomy=sequana.scripts.taxonomy:main',
           'sequana_coverage=sequana.scripts.coverage:main',
           'sequana_summary=sequana.scripts.summary:main',
           'sequana_mapping=sequana.scripts.mapping:main',
           'sequana_compressor=sequana.scripts.compressor:main',
           'sequana_report=sequana.scripts.reports:main',
           'sequana_vcf_filter=sequana.scripts.vcf_filter:main', # june 2018
           'sequana_bam_splitter=sequana.scripts.bam_splitter:main', # aug 2018
        ],
        'sequana.module':[
            'sequana_coverage=sequana.modules_report.coverage:CoverageModule',
            'sequana_variant_calling=sequana.modules_report.variant_calling:VariantCallingModule',
            'sequana_summary=sequana.modules_report.summary:SummaryModule',
        ],
        "multiqc.modules.v1": [
            "sequana_pacbio_qc=sequana.multiqc.pacbio_qc:MultiqcModule",
            "sequana_quality_control=sequana.multiqc.quality_control:MultiqcModule",
            "sequana_coverage=sequana.multiqc.coverage:MultiqcModule",
            "sequana_isoseq=sequana.multiqc.isoseq:MultiqcModule",
            "sequana_isoseq_qc=sequana.multiqc.isoseq_qc:MultiqcModule",
            "sequana_bamtools_stats=sequana.multiqc.bamtools_stats:MultiqcModule",
            "sequana_kraken=sequana.multiqc.kraken:MultiqcModule",
            "sequana_laa=sequana.multiqc.laa:MultiqcModule"
        ],
        'multiqc.hooks.v1': [
            'before_config = sequana.multiqc.config:load_config',
        ]
    },


)
