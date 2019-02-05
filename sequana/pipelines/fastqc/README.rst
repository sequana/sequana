:Overview: fastqc + multiqc 
:Input: A set of FastQ files (paired or single-end)
:Output: summary.html

Usage
~~~~~~~

::

    sequana --pipeline fastqc -i . -o analysis

This creates a directory **analysis**. You just need to execute the pipeline::

    cd analysis
    sh runme.sh


If you are familiar with snakemake, you can retrieve the fastqc.rules and config.yaml files and then execute the pipeline::

    snakemake -s fastqc.rules --cores 4 --stats stats.txt

Or use :ref:`sequanix_tutorial` interface.

Requirements
~~~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/fastqc/requirements.txt

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/fastqc/dag.png


Details
~~~~~~~~~

This pipeline runs fastqc in parallel on the input fastq files (paired or not)
and then execute multiqc. A brief sequana summary report is also produced.


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/fastqc/config.yaml` to be used with the pipeline. Each rule used in the pipeline may have a section in the
configuration file. 


FastQC
^^^^^^^^^^^
.. snakemakerule:: fastqc_dynamic

mutliqc
^^^^^^^^^^^^^^^
.. snakemakerule:: multiqc2

