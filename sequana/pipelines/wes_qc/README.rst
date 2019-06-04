:Overview: Quality control, trimming (adapter removal) and taxonomic overview
:Input: A set of FastQ files (paired or single-end)
:Output: fastqc, cleaned FastQ files, Krona plot

Usage
~~~~~~~

::

    sequana --pipeline wes_qc --input-directory . --working-directory analysis --no-adapters


Or use :ref:`sequanix_tutorial` interface.

Requirements
~~~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/quality_control/requirements.txt


Details
~~~~~~~~~



Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/quality_control/config.yaml` to be used with the pipeline. Each rule used in the pipeline may have a section in the
configuration file. In the *quality_control* pipeline, we use the *bwa_mem*,
*fastqc*, *cutadapt*, *fastq_stats* and *kraken* rules described here below.


FastQC
^^^^^^^^^^^
.. snakemakerule:: fastqc_dynamic

Cutadapt
^^^^^^^^^
.. snakemakerule:: cutadapt

Kraken
^^^^^^^
.. snakemakerule:: kraken

BWA-mem
^^^^^^^^
.. snakemakerule:: bwa_mem_dynamic

FastQ stats
^^^^^^^^^^^^
.. snakemakerule:: fastq_stats_dynamic
