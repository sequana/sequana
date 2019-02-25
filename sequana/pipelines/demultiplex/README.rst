:Overview: demultiplex
:Input: bcl
:Output: fastq

Usage
~~~~~~~

::

    sequanix -w sequanix -p demultiplex -i 20190207_FS10000482_2_BPA73010-1129/

Requirements
~~~~~~~~~~~~~~~~~~

This pipeline is used with bcl2fastq 2.20.0

.. include:: ../sequana/pipelines/demultiplex/requirements.txt

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/demultiplex/dag.png


Details
~~~~~~~~~

bcl2fastq
^^^^^^^^^^^^
.. snakemakerule:: bcl2fastq

