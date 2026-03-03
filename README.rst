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

.. image:: https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C3.12-blue.svg
    :target: https://pypi.python.org/pypi/sequana
    :alt: Python 3.10 | 3.11 | 3.12

.. image:: https://img.shields.io/github/issues/sequana/sequana.svg
    :target: https://github.com/sequana/sequana/issues
    :alt: GitHub Issues

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black

10.5281/zenodo.853159

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.853159.svg
   :target: https://doi.org/10.5281/zenodo.853159
   :alt: DOI

.. image:: https://static.pepy.tech/badge/sequana
   :target: https://pepy.tech/project/sequana
   :alt: downloads

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

Here is a non-exhaustive list of tools and pipelines from the project, with users and developers audience.


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
0.20.0    * Refactoring FASTA and GFF3 modules/classes (less memory, handles large
            eukaryotes). Update Kozak module.
0.19.6    * add gff command; rename and improve the add_CDS function in GFF class;
            fix regression in kraken analysis
0.19.5    * add embl2fasta, fix CDS parents in GFF file
0.19.4    * improved TRF module; add fastq_split and html_report commands
0.19.3    * fix plotly issue in RNAdiff plot (#872)
0.19.2    * new modules for genomic metrics (zdna, imotif, cruciform, etc.);
            new visualisation tools; parser for hmmtools
0.19.1    * update pyproject dependencies; fix multiqc plugins
0.19.0    * pyproject updated for Poetry 2.0; drop Python 3.8 support;
            new modules: kozak, somy score, telomere, biomol, rnafold,
            restriction enzyme; more sequence metrics
0.18.0    * new somy scores module and standalone; coverage now uses mosdepth
            for bam2cov; drop Python 3.8 support
0.17.3    * fix multiqc/feature counts handling in RNAdiff; add size factor table
0.17.2    * pin pulp<2.8 and snakemake<8.0
0.17.1    * new tsne plot; update IEM module with additional specs
0.17.0    * remove substractor utility; speedup KEGG enrichment; IEM class
            renamed to SampleSheet; new find-integrated-genes standalone;
            fixes in rnadiff HTML report and VCF code cleanup
0.16.9    * fix PCA and add batch effect plots in RNAdiff; taxonomy saves
            taxonomy.dat as gzipped CSV
0.16.8    * update IEM; better error handling in RNAdiff; new ribodesigner methods
0.16.6    * refactor IEM to be more robust with more tests
0.16.5    * refactor to use pyproject; remove pkg_resources; cleanup resources
0.16.3    * remove all rules (moved to sequana-wrappers); add precommit
0.16.0    * add mpileup module; complete refactoring of sequana coverage module
            (small eukaryotes, memory efficient); use click for CLI tools
0.15.3    * add sequana.viz.plotly module; update KEGG to use headless server;
            improvements on KEGG enrichment (pathways, plotly, CSV export)
0.15.2    * ribodesigner accepts fasta with no GFF; fix IEM double indexing;
            refactorisation of compare module; enrichment improvements
0.15.0    * add logo in reports; RNAdiff shrinkage support; shrunken log2FC
            values in volcano plot and report table
0.14.6    * add fasta_and_gff_annotation, macs3, idr, phantom and homer modules
0.14.4    * fix kegg colorised pathways; pin snakemake>=7.16
0.14.3    * new fisher metric in variant calling; support several features in
            rnaseq/rnadiff
0.14.1    * KEGG enrichment improvements; new uniprot module for GO term enrichment
0.14.0    * pin click>=8.1.0; ribodesigner new plots, clustering and notebook
0.13.X    * cleanup standalones; new ribodesigner subcommand; lane_merging moved
            to sequana sub-command
0.12.X    * rules moved to sequana-wrappers; VCF tools refactored to use vcfpy;
            complete changelog before 0.12.4 in the github /doc/Changelog.txt
========= ==========================================================================

Any :question: Feel free to [open an issue](https://github.com/sequana/sequana/issues)
