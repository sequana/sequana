.. _wrappers:

Wrappers
##########

As of August 2021, **Sequana** team created the e `sequana wrappers repository <https://github.com/sequana/sequana-wrappers>`_, which is intended to replace the rules. The adavantage is that wrappers can be tested with a continuous integration.  


Wrappers are used within a Snakemake rule. When you call your Snakemake
pipeline, you will need to add::

    --wrapper-prefix git+file:https://github.com/sequana/sequana-wrappers/

We provide documentation for each wrapper. It can be included in this
documentation thanks to a sphinx extension. For example::

    .. sequana_wrapper:: fastqc

Here is a non exhaustive list of documented wrappers. 


bowtie2/align
==============

.. sequana_wrapper:: bowtie2/align

bowtie2/build
==============

.. sequana_wrapper:: bowtie2/build


add_read_group
==============

.. sequana_wrapper:: add_read_group


bam_coverage
============

.. sequana_wrapper:: bam_coverage


bcl2fastq
=========

.. sequana_wrapper:: bcl2fastq


blast
=====

.. sequana_wrapper:: blast


busco
=====

.. sequana_wrapper:: busco


bz2_to_gz
=========

.. sequana_wrapper:: bz2_to_gz


canu
====

.. sequana_wrapper:: canu


consensus
=========

.. sequana_wrapper:: consensus


digital_normalisation
=====================

.. sequana_wrapper:: digital_normalisation


dsrc_to_gz
==========

.. sequana_wrapper:: dsrc_to_gz


falco
=====

.. sequana_wrapper:: falco


fastp
=====

.. sequana_wrapper:: fastp


fastq_stats
===========

.. sequana_wrapper:: fastq_stats


fastqc
======

.. sequana_wrapper:: fastqc


feature_counts
==============

.. sequana_wrapper:: feature_counts


freebayes
=========

.. sequana_wrapper:: freebayes


freebayes_vcf_filter
====================

.. sequana_wrapper:: freebayes_vcf_filter


gz_to_bz2
=========

.. sequana_wrapper:: gz_to_bz2


hmmbuild
========

.. sequana_wrapper:: hmmbuild


hmmscan
=======

.. sequana_wrapper:: hmmscan


index
=====

.. sequana_wrapper:: index


longorfs
========

.. sequana_wrapper:: longorfs


macs3
=====

.. sequana_wrapper:: macs3


makeblastdb
===========

.. sequana_wrapper:: makeblastdb


mark_duplicates
===============

.. sequana_wrapper:: mark_duplicates


minimap2
========

.. sequana_wrapper:: minimap2


multiqc
=======

.. sequana_wrapper:: multiqc


polypolish
==========

.. sequana_wrapper:: polypolish


predict
=======

.. sequana_wrapper:: predict


prokka
======

.. sequana_wrapper:: prokka


quast
=====

.. sequana_wrapper:: quast


rulegraph
=========

.. sequana_wrapper:: rulegraph


sambamba_filter
===============

.. sequana_wrapper:: sambamba_filter


sambamba_markdup
================

.. sequana_wrapper:: sambamba_markdup


samtools_depth
==============

.. sequana_wrapper:: samtools_depth


sequana_coverage
================

.. sequana_wrapper:: sequana_coverage


sequana_taxonomy
================

.. sequana_wrapper:: sequana_taxonomy


snpeff
======

.. sequana_wrapper:: snpeff


snpeff_add_locus_in_fasta
=========================

.. sequana_wrapper:: snpeff_add_locus_in_fasta


sort
====

.. sequana_wrapper:: sort


spades
======

.. sequana_wrapper:: spades


trinity
=======

.. sequana_wrapper:: trinity


trinity_quantify
================

.. sequana_wrapper:: trinity_quantify


unicycler
=========

.. sequana_wrapper:: unicycler

