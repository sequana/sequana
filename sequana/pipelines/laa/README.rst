:Overview: laa pipeline
:Input: A set of FastQ files or BAM files (see details)
:Output: summary.html

Usage
~~~~~~~

::

    sequana --pipeline laa --input-directory . --working-directory analysis


Or use :ref:`sequanix_tutorial` interface.

Requirements
~~~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/laa/requirements.txt

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/laa/dag.png


Details
~~~~~~~~~

This pipeline takes as input CCS reads from a pacbio amplicon analysis. It maps
the reads on a reference that must be provided. A consensus is created from the
coverage (using IGVtools). From the refernce, freebayes is called to obtain VCF
files. Snpeff is used for annotation given a genbank to be provided by the user.
Kraken is used for a quick taxonomy on the input CCS reads. The consensus
obtained are used to build a phylogentic tree based on mafft and raxml. The
output is sent to itol website to retrieve a phylogenetic tree image. Finally,
various HTML reports are created including a multiqc report.


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/fastqc/config.yaml` to be used with the pipeline. Each rule used in the pipeline may have a section in the
configuration file. 



mutliqc
^^^^^^^^^^^^^^^
.. snakemakerule:: multiqc2

