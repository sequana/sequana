Sequana documentation
##########################################

|version|, |today|


.. image:: https://badge.fury.io/py/sequana.svg
    :target: https://pypi.python.org/pypi/sequana

.. image:: https://travis-ci.org/sequana/sequana.svg?branch=master
    :target: https://travis-ci.org/sequana/sequana

.. image::  https://coveralls.io/repos/github/sequana/sequana/badge.svg?branch=master
   :target: https://coveralls.io/github/sequana/sequana?branch=master 

.. image:: http://readthedocs.org/projects/sequana/badge/?version=latest
    :target: http://sequana.readthedocs.org/en/latest/?badge=latest
    :alt: Documentation Status

What is Sequana ?
=====================

Sequana provides (i) modules and pipelines dedicated to NGS in the form of Snakefiles (Makefile-like syntax with Python syntax) and (ii) tools to help in the creation of such pipelines and associtated HTML reports.


Installation
#################


:: 

    pip install sequana

With conda::

    conda install numpy pandas matplotlib
    pip install sequana

Examples
##########

Given a set of FASTQ files in a local directory, you can run a pipeline as
follows (given its name)::

    sequana --init fix_removal
    # populate the fastq_raw with your FASTQ files (R1 and R2)
    # execute the script sequana.sh
    sh sequana.sh
    # open the report in ./report
    open report/index.html



This will copy the Snakefile locally as well as the corresponding configuration
file.

Pipelines can be found in the directory ./pipelines or ./rules



Issues
##########

Please fill bug report in https://github.com/sequana/sequana



.. toctree::
    :maxdepth: 2

    userguide
    references
    Changelog


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
