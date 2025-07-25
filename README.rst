SEQUANA
############


.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)
   :target: http://bioconda.github.io/recipes/sequana/README.html

.. image:: https://badge.fury.io/py/sequana.svg
    :target: https://pypi.python.org/pypi/sequana

.. image:: https://github.com/sequana/sequana/actions/workflows/main.yml/badge.svg?branch=main
    :target: https://github.com/sequana/sequana/actions/workflows/main.yml

.. image:: https://coveralls.io/repos/github/sequana/sequana/badge.svg?branch=main
    :target: https://coveralls.io/github/sequana/sequana?branch=main

.. image:: http://readthedocs.org/projects/sequana/badge/?version=main
    :target: http://sequana.readthedocs.org/en/latest/?badge=main
    :alt: Documentation Status

.. image:: http://joss.theoj.org/papers/10.21105/joss.00352/status.svg
   :target: http://joss.theoj.org/papers/10.21105/joss.00352
   :alt: JOSS (journal of open source software) DOI

.. image:: https://img.shields.io/badge/python-3.8%20%7C%203.9%20%7C3.10-blue.svg
    :target: https://pypi.python.org/pypi/sequana
    :alt: Python 3.9 | 3.10 | 3.11 | 3.12

.. image:: https://img.shields.io/github/issues/sequana/sequana.svg
    :target: https://github.com/sequana/sequana/issues
    :alt: GitHub Issues

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black


:How to cite: Citations are important for us to carry on developments.
    For Sequana library (including the pipelines), please use

    Cokelaer et al, (2017), 'Sequana': a Set of Snakemake NGS pipelines, Journal of
    Open Source Software, 2(16), 352, `JOSS DOI doi:10.21105/joss.00352 <https://joss.theoj.org/papers/10.21105/joss.00352>`_

    For the **genome coverage** tool (sequana_coverage):  Desvillechabrol et al, 2018:
    detection and characterization of genomic variations using running median and
    mixture models. GigaScience, 7(12), 2018. https://doi.org/10.1093/gigascience/giy110

    For **Sequanix**: Desvillechabrol et al.
    Sequanix: A Dynamic Graphical Interface for Snakemake Workflows
    Bioinformatics, bty034, https://doi.org/10.1093/bioinformatics/bty034
    Also available on bioRxiv (DOI: https://doi.org/10.1101/162701)


🔧 Overview and Installation
============================

**Sequana** is a Python library dedicated to bioinformatics. It is also a project that includes a set of pipelines related to NGS (new generation sequencing) including quality control, variant calling, coverage, taxonomy, transcriptomics. We also ship **Sequanix**, a graphical user interface for Snakemake pipelines.


Pipelines and related projects
==============================

Here is a non exahustive list of tools and pipelines from the project. with users and developers audience.


.. list-table::
    :widths: 15 35 20 15 15
    :header-rows: 1

    * - **name/github**
      - **description**
      - **Latest Pypi version**
      - **Test passing**
      - **apptainers**
    * - `sequana_pipetools <https://github.com/sequana/sequana_pipetools>`_
      - Create and Manage Sequana pipeline
      - .. image:: https://badge.fury.io/py/sequana-pipetools.svg
            :target: https://pypi.python.org/pypi/sequana_pipetools
      - .. image:: https://github.com/sequana/sequana_pipetools/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/sequana_pipetools/actions/workflows/main.yml
      -  Not required
    * - `sequana-wrappers <https://github.com/sequana/sequana-wrappers>`_
      - Set of wrappers to build pipelines
      - Not on pypi
      - .. image:: https://github.com/sequana/sequana-wrappers/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/sequana-wrappers/actions/workflows/main.yml
      - Not required
    * - `demultiplex <https://github.com/sequana/demultiplex>`_
      - Demultiplex your raw data
      - .. image:: https://badge.fury.io/py/sequana-demultiplex.svg
            :target: https://pypi.python.org/pypi/sequana-demultiplex
      - .. image:: https://github.com/sequana/demultiplex/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/demultiplex/actions/workflows/main.yml
      - License restriction
    * - `denovo <https://github.com/sequana/denovo>`_
      - denovo sequencing data
      - .. image:: https://badge.fury.io/py/sequana-denovo.svg
            :target: https://pypi.python.org/pypi/sequana-denovo
      - .. image:: https://github.com/sequana/denovo/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/denovo/actions/workflows/main.yml
      - .. image:: https://github.com/sequana/denovo/actions/workflows/apptainer.yml/badge.svg
            :target: https://github.com/sequana/denovo/actions/workflows/apptainer.yml
    * - `fastqc <https://github.com/sequana/fastqc>`_
      - Get Sequencing Quality control
      - .. image:: https://badge.fury.io/py/sequana-fastqc.svg
            :target: https://pypi.python.org/pypi/sequana-fastqc
      - .. image:: https://github.com/sequana/fastqc/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/fastqc/actions/workflows/main.yml
      - .. image:: https://github.com/sequana/fastqc/actions/workflows/apptainer.yml/badge.svg
            :target: https://github.com/sequana/fastqc/actions/workflows/apptainer.yml
    * - `LORA <https://github.com/sequana/lora>`_
      - Map sequences on target genome
      - .. image:: https://badge.fury.io/py/sequana-lora.svg
            :target: https://pypi.python.org/pypi/sequana-lora
      - .. image:: https://github.com/sequana/lora/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/lora/actions/workflows/main.yml
      - .. image:: https://github.com/sequana/lora/actions/workflows/apptainer.yml/badge.svg
            :target: https://github.com/sequana/lora/actions/workflows/apptainer.yml
    * - `mapper <https://github.com/sequana/mapper>`_
      - Map sequences on target genome
      - .. image:: https://badge.fury.io/py/sequana-mapper.svg
            :target: https://pypi.python.org/pypi/sequana-mapper
      - .. image:: https://github.com/sequana/mapper/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/mapper/actions/workflows/main.yml
      - .. image:: https://github.com/sequana/mapper/actions/workflows/apptainer.yml/badge.svg
            :target: https://github.com/sequana/mapper/actions/workflows/apptainer.yml
    * - `nanomerge <https://github.com/sequana/nanomerge>`_
      - Merge barcoded (or unbarcoded) nanopore fastq and reporting
      - .. image:: https://badge.fury.io/py/sequana-nanomerge.svg
            :target: https://pypi.python.org/pypi/sequana-nanomerge
      - .. image:: https://github.com/sequana/nanomerge/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/nanomerge/actions/workflows/main.yml
      - .. image:: https://github.com/sequana/nanomerge/actions/workflows/apptainer.yml/badge.svg
            :target: https://github.com/sequana/nanomerge/actions/workflows/apptainer.yml
    * - `pacbio_qc <https://github.com/sequana/pacbio_qc>`_
      - Pacbio quality control
      - .. image:: https://badge.fury.io/py/sequana-pacbio-qc.svg
            :target: https://pypi.python.org/pypi/sequana-pacbio-qc
      - .. image:: https://github.com/sequana/pacbio_qc/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/pacbio_qc/actions/workflows/main.yml
      - .. image:: https://github.com/sequana/pacbio_qc/actions/workflows/apptainer.yml/badge.svg
            :target: https://github.com/sequana/pacbio_qc/actions/workflows/apptainer.yml
    * - `ribofinder <https://github.com/sequana/ribofinder>`_
      - Find ribosomal content
      - .. image:: https://badge.fury.io/py/sequana-ribofinder.svg
            :target: https://pypi.python.org/pypi/sequana-ribofinder
      - .. image:: https://github.com/sequana/ribofinder/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/ribofinder/actions/workflows/main.yml
      - .. image:: https://github.com/sequana/ribofinder/actions/workflows/apptainer.yml/badge.svg
            :target: https://github.com/sequana/ribofinder/actions/workflows/apptainer.yml
    * - `rnaseq <https://github.com/sequana/rnaseq>`_
      - RNA-seq analysis
      - .. image:: https://badge.fury.io/py/sequana-rnaseq.svg
            :target: https://pypi.python.org/pypi/sequana-rnaseq
      - .. image:: https://github.com/sequana/rnaseq/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/rnaseq/actions/workflows/main.yml
      - .. image:: https://github.com/sequana/rnaseq/actions/workflows/apptainer.yml/badge.svg
            :target: https://github.com/sequana/rnaseq/actions/workflows/apptainer.yml
    * - `variant_calling <https://github.com/sequana/variant_calling>`_
      - Variant Calling
      - .. image:: https://badge.fury.io/py/sequana-variant-calling.svg
            :target: https://pypi.python.org/pypi/sequana-variant-calling
      - .. image:: https://github.com/sequana/variant_calling/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/variant_calling/actions/workflows/main.yml
      - .. image:: https://github.com/sequana/variant_calling/actions/workflows/apptainer.yml/badge.svg
            :target: https://github.com/sequana/variant_calling/actions/workflows/apptainer.yml
    * - `multicov <https://github.com/sequana/multicov>`_
      - Coverage (mapping)
      - .. image:: https://badge.fury.io/py/sequana-multicov.svg
            :target: https://pypi.python.org/pypi/sequana-multicov
      - .. image:: https://github.com/sequana/multicov/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/multicov/actions/workflows/main.yml
      - .. image:: https://github.com/sequana/coverage/actions/workflows/apptainer.yml/badge.svg
            :target: https://github.com/sequana/coverage/actions/workflows/apptainer.yml
    * - `laa <https://github.com/sequana/laa>`_
      - Long read Amplicon Analysis
      - .. image:: https://badge.fury.io/py/sequana-laa.svg
            :target: https://pypi.python.org/pypi/sequana-laa
      - .. image:: https://github.com/sequana/laa/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/laa/actions/workflows/main.yml
      - .. image:: https://github.com/sequana/laa/actions/workflows/apptainer.yml/badge.svg
            :target: https://github.com/sequana/laa/actions/workflows/apptainer.yml
    * - `revcomp <https://github.com/sequana/revcomp>`_
      - reverse complement of sequence data
      - .. image:: https://badge.fury.io/py/sequana-revcomp.svg
            :target: https://pypi.python.org/pypi/sequana-revcomp
      - .. image:: https://github.com/sequana/revcomp/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/revcomp/actions/workflows/main.yml
      - .. image:: https://github.com/sequana/revcomp/actions/workflows/apptainer.yml/badge.svg
            :target: https://github.com/sequana/revcomp/actions/workflows/apptainer.yml
    * - `downsampling <https://github.com/sequana/downsampling>`_
      - downsample sequencing data
      - .. image:: https://badge.fury.io/py/sequana-downsampling.svg
            :target: https://pypi.python.org/pypi/sequana-downsampling
      - .. image:: https://github.com/sequana/downsampling/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/downsampling/actions/workflows/main.yml
      - Not required
    * - `depletion <https://github.com/sequana/depletion>`_
      - remove/select reads mapping a reference
      - .. image:: https://badge.fury.io/py/sequana-downsampling.svg
            :target: https://pypi.python.org/pypi/sequana-depletion
      - .. image:: https://github.com/sequana/depletion/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/depletion/actions/workflows/main.yml
      -





.. list-table:: Pipelines not yet released
    :widths: 20 40 20 20
    :header-rows: 1

    * - **name/github**
      - **description**
      - **Latest Pypi version**
      - **Test passing**
    * - `trf <https://github.com/sequana/trf>`_
      - Find repeats
      - .. image:: https://badge.fury.io/py/sequana-trf.svg
            :target: https://pypi.python.org/pypi/sequana-trf
      - .. image:: https://github.com/sequana/trf/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/trf/actions/workflows/main.yml
    * - `multitax <https://github.com/sequana/multitax>`_
      - Taxonomy analysis
      - .. image:: https://badge.fury.io/py/sequana-multitax.svg
            :target: https://pypi.python.org/pypi/sequana-multitax
      - .. image:: https://github.com/sequana/multitax/actions/workflows/main.yml/badge.svg
            :target: https://github.com/sequana/multitax/actions/workflows/main.yml

**Please see the** `documentation <http://sequana.readthedocs.org>`_ for an
up-to-date status and documentation.


Contributors
============

Maintaining Sequana would not have been possible without users and contributors.
Each contribution has been an encouragement to pursue this project. Thanks to all:

.. image:: https://contrib.rocks/image?repo=sequana/sequana
    :target: https://github.com/sequana/sequana/graphs/contributors


Changelog :memo:
~~~~~~~~~~~~~~~~~

========= ==========================================================================
Version   Description
========= ==========================================================================
0.19.2    * NEW new modules related to genomic metrics (zdna, imotif, cruciform,etc)
          * NEW new visualisation tools
          * NEW parser for hmmtools
0.19.1    * update pyproject with matploblib, scipy, and all dependencie except
            snakemake (keep <8.X and its pulp dependency <2.8) because 8.X breaks
            the containers and therefore all pipelines. no workaround. we will need
            a snakemake 'frozen library'. This is unfortunate.
          * CHANGES: remove multiqc plugin quality_control (obsolet) and updated
            many plugin due to increased version ofmultiqc.
0.19.0    * UPDATE pyproject layout to be poetry2.0 compatible. drop support py3.8
          * NEW: kozak module, somy score, telomere, biomol, rnafold, restriction
            enzyme
          * UPDATES: sequence module has more metrics (e.g. karlin signature),
            fastq has a histogram for long reads,
0.18.0    * NEW: somy scores (module and standalone)
          * CHANGES: coverage now uses mosdepth for the bam2cov. can also handle
            creation of the 4-column coverage file from BAM on the fly.
          * Drop python 3.8 support (more and more dependencies issues)
0.17.3    * CHANGES: fix https://github.com/sequana/rnaseq/issues/45 to better
            handle multiqc files especially, feature counts plugin. Fix #855 to
            be able to name batch effect column arbitrary to any name. Fix batch
            effect column (and all others) to be strings rather than integer. If
            integer, factor are sorted based on the integer values, not the factor
            itself.
          * NEWS: add size factor comments and table.
0.17.2    * Pin down pulp<2.8 and snakemake to <8.0 (too many changes in snakemake)
0.17.1    * NEW: tsne plot
          * CHANGES: update iem module with additional specs
0.17.0    * CHANGE. remove the substractor utility (use sequana_depletion pipeline
            instead)
          * CHANGE. remove get_max_gc_correlation function from bedtools. not used.
          * CHANGE. Got rid of freebayes_bcf_filter redundant with
            freebayes_vcf_filter; replace scipy fisher test with own implementation.
            Remove useless VCF code.
          * FIXES. rnadiff HTML report
          * IMPROV. speedup kegg enrichment using multiprocess
          * IMPROV. sequana_taxonomy can now download toydb and viruses_masking DBs
            from zenodo
          * NEW function to retrieve the GO terms and genes from panther website
          * NEW keep_reads function in fastq class
          * Major update of the IEM module (renamed IEM class into SampleSheet)
          * NEW addW find-integrated-genes standalone
0.16.9    * Major fix on PCA and add batch effect plots in RNAdiff analysis
          * count matrix and DESeq2 output files' headers fixed with missing index
            (no impact on analysis but only for those willing to use the CSV files
            in excel)
          * Taxonomy revisited to save taxonomy.dat in gzipped CSV format.
0.16.8    * update IEM for more testing
          * better handling of error in RNADiff
          * Add new methods for ribodesigner
0.16.7    * Stable release (fix doc), deprecated.
0.16.6    * Refactor IEM to make it more robust with more tests.
0.16.5    * refactor to use pyproject instead of setuptools
          * remove pkg_resources (future deprecation)
          * remove unused requirements (cookiecutter, adjusttext, docutuils, mock,
            psutil, pykwalify)
          * cleanup resources (e.g. moving canvas/bar.py into viz)
0.16.4    * hot fixes on RNAdiff reports and enrichments
0.16.3    * Remove all rules (see https://github.com/sequana/sequana-wrappers)
            instead
          * add precommit for developers and applied to all modules and doc
          * Fix wrong import for sequana standalone (regression)
0.16.2    * save coverage PNG image (regression)
          * Update taxonomy/coverage standalone (regression) and more tests
0.16.1    * hotfix missing module
0.16.0    * add mpileup module
          * homogenization enrichment + fixup rnadiff
          * Complete refactoring of sequana coverage module.
            Allow sequana_coverage to handle small eukaryotes in a more memory
            efficient way.
          * use click for the sequana_taxonomy and sequana_coverage and
            sequana rnadiff command
          * Small fixup on homer, idr and phantom modules (for chipseq pipeline)
0.15.4    * add plot for rnaseq/rnadiff
0.15.3    * add sequana.viz.plotly module. use tqdm in bamtools module
          * KEGG API changed. We update sequana to use headless server and keep
            the feature of annotated and colored pathway.
          * Various improvements on KEGG enrichment including saving pathways,
            addition --comparison option in sequana sub-command, plotly plots, etc
0.15.2    * ribodesigner can now accept an input fasta with no GFF assuming the
            fasta already contains the rRNA sequences
          * Fix IEM module when dealing with double indexing
          * Fix anchors in HTML reports (rnadiff module)
          * refactorise compare module to take several rnadiff results as input
          * enrichment improvements (export KEGG and GO as csv files
0.15.1    * Fix creation of images directory in modules report
          * add missing test related to gff
          * Fix #804
0.15.0    * add logo in reports
          * RNADiff reports can now use shrinkage or not (optional)
          * remove useless rules now in sequana-wrappers
          * update main README to add LORA in list of pipelines
          * Log2FC values are now **shrinked log2FC** values in volcano plot
            and report table. "NotShrinked" columns for Log2FC and Log2FCSE
            prior shrinkage are displayed in report table.
0.14.6    * add fasta_and_gff_annotation module to correct fasta and gff given a
            vcf file.
          * add macs3 module to read output of macs3 peak detector.
          * add idr module to read results of idr analysis
          * add phantom module to compute phantom peaks
          * add homer module to read annotation files from annotatePeaks
0.14.5    * move start_pipeline standalone in
            https://github.com/sequana/sequana_pipetools
          * update snpeff module to allows build command to have options
0.14.4    * hotfix bug on kegg colorised pathways
          * Fix the hover_name in rnadiff volcano plot to include the
            index/attribute.
          * pin snakemake to be >=7.16
0.14.3    * new fisher metric in variant calling
          * ability to use several feature in rnaseq/rnadiff
          * pin several libaries due to regression during installs
0.14.2    * Update ribodesigner
0.14.1    * Kegg enrichment: add gene list 'all' and fix incomplete annotation case
          * New uniprot module for GO term enrichment and enrichment
            refactorisation (transparent for users)
0.14.0    * pinned click>=8.1.0 due to API change (autocomplete)
          * moved tests around to decrease packaging from 16 to 4Mb
          * ribodesigner: new plots, clustering and notebook
0.13.X    * Remove useless standalones or moved to main **sequana** command
          * Move sequana_lane_merging into a subcommand (sequana lane_merging)
          * General cleanup of documentation, test and links to pipelines
          * add new ribodesigner subcommand
0.12.X    * remove some rules now in https://github.com/sequana/sequana-wrappers
          * refactorisation of VCF tools/modules to use vcfpy instead of pyVCF
          * complete change log before 0.12.4 in the github /doc/Changelog.txt
========= ==========================================================================

Any :question: Feel free to [open an issue](https://github.com/sequana/sequana/issues)
