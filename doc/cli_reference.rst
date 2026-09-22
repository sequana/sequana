.. _cli_reference:

CLI reference
=============

The top-level ``sequana`` command groups ~30 sub-commands (FASTQ/FASTA
utilities, GFF/GTF fixers, enrichment helpers, summaries, …). The reference
below is generated directly from the source.

Each sub-command also accepts ``--help`` from the shell, e.g.::

    sequana fastq --help
    sequana enrichment-kegg --help

Quick lookup by task
====================

.. list-table::
   :widths: 25 30
   :header-rows: 1

   * - Task
     - Command
   * - FASTQ/FASTA utilities
     - ``sequana fastq``, ``sequana fasta``, ``sequana fasta-compare``
   * - Genome annotation (GFF/GTF)
     - ``sequana gff3``, ``sequana gtf``
   * - Variant analysis
     - ``sequana vcf``, ``sequana vcf-filter``
   * - Sequence properties
     - ``sequana seq-length``, ``sequana gc-content``
   * - Coverage & depth
     - ``sequana coverage``, ``sequana repeats``
   * - Taxonomy & enrichment
     - ``sequana taxonomy``, ``sequana enrichment-*``
   * - Reports & summaries
     - ``sequana summary``, ``sequana mapping``

Full reference
==============

.. click:: sequana.scripts.main.main:main
   :prog: sequana
   :nested: full
