# -*- coding: utf-8 -*-
__revision__ = "$Id: $" # for the SVN Id
import sys
import os
from setuptools import setup, find_packages
import glob

_MAJOR               = 0
_MINOR               = 0
_MICRO               = 4
version              = '%d.%d.%d' % (_MAJOR, _MINOR, _MICRO)
release              = '%d.%d' % (_MAJOR, _MINOR)

metainfo = {
    'authors': {"main": ("yourname", "email@whatever.org")},
    'version': version,
    'license' : 'GPL',
    'download_url' : ['http://pypi.python.org/pypi/sequana'],
    'url' : ["http://github.com/sequana/"],
    'description': "Put a short description here" ,
    'platforms' : ['Linux', 'Unix', 'MacOsX', 'Windows'],
    'keywords' : [''],
    'classifiers' : [
          'Development Status :: 4 - Beta',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Topic :: Software Development :: Libraries :: Python Modules',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Information Analysis',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Physics']
    }


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
    packages = ['sequana', 'sequana.funny', 'rules'],
    package_dir = {
        "sequana": "src",
        "sequana.funny": 'src/funny',
        "rules": "rules"
    },
    include_package_data = True,

    install_requires = ["easydev==0.9.17", "reports==0.1.2", "matplotlib", "pandas",
        "cutadapt==1.9.1", "pysam", "pyVCF"],



    # here below '': pattern means include that pattern in all packages
    # so '' :['README.rst'] will include all README.rst recursively
    package_data = {
        'share.data' : ["*"],
        '' : ['*.html', "*.css"],
        },


    zip_safe=False,
    #entry_points = {
    #    'console_scripts':[
     #       'fastq_head=scripts.fastq_head:main',
     #       'fastq_count=scripts.fastq_count:main',
      #  ]
    #},

)
