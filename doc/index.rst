Sequana documentation
##########################################

|release|

.. include:: ../README.rst


What is Sequana?
================

**Sequana** is a versatile bioinformatics tool that provides:

#. A Python library dedicated to NGS analysis (file-format wrappers, plots,
   coverage, taxonomy, enrichment, RNA-seq, …).
#. A family of :ref:`pipelines <pipelines>` based on
   `Snakemake <https://snakemake.readthedocs.io>`_, each living in its own
   GitHub repository.
#. A set of :ref:`standalone applications <applications>`:
    #. :ref:`sequana_coverage <standalone_sequana_coverage>` — genome coverage
       analysis with confidence intervals.
    #. :ref:`sequana_taxonomy <standalone_sequana_taxonomy>` — quick taxonomy
       on a FastQ file (Kraken + Krona).
    #. ~30 other utility sub-commands grouped under the top-level ``sequana``
       CLI (see :ref:`applications`).

Pipelines cover NGS quality control, variant calling, coverage analysis,
taxonomy, de-novo assembly, :ref:`RNA-seq <pipeline_rnaseq>`,
:ref:`variant calling <pipeline_vc>` and more — see the :ref:`pipelines`
catalogue.

**Sequana** can be used either as a Python library (developers, library users)
or via its pipelines and standalones (end users). To join the project, please
let us know on `github <https://github.com/sequana/sequana/issues/306>`__.

Quick start
===========

Install Sequana and run a quality-control pipeline in 2 minutes::

    pip install sequana_fastqc
    sequana_fastqc --input-directory . --execute

Or use the Python library to read and analyze genomic data::

    from sequana import FastA, SequanaCoverage
    fasta = FastA("genome.fasta")
    cov = SequanaCoverage("coverage.bed")
    cov[0].plot_coverage()

Full guides below.

Choosing your entry point
==========================

Not sure where to start? This table helps:

.. list-table::
   :widths: 20 50
   :header-rows: 1

   * - What you want to do
     - Use this
   * - Run a pre-built NGS pipeline (fastqc, variant calling, RNA-seq, etc.)
     - :ref:`pipelines` → :ref:`tutorial`
   * - Use Python to read/analyze sequences, coverage, or variants
     - :ref:`userguide` → :ref:`references`
   * - Find a quick helper tool (rename files, filter GFF, count bases)
     - ``sequana --help`` → :ref:`cli_reference`
   * - See working code examples
     - :doc:`auto_examples/index` or :doc:`notebooks`


.. |bam| image::
    ./auto_examples/images/sphx_glr_plot_bam_001.png
    :target: auto_examples/plot_bam.html

.. |coverage| image::
    ./auto_examples/images/sphx_glr_plot_coverage_001.png
    :target: auto_examples/plot_coverage.html

.. |fastqc| image::
    ./auto_examples/images/sphx_glr_plot_fastqc_hist_001.png
    :target: auto_examples/plot_fastqc_hist.html

.. |kraken| image::
    ./auto_examples/images/sphx_glr_plot_kraken_001.png
    :target: auto_examples/plot_kraken.html

.. |pacbio| image::
    ./auto_examples/images/sphx_glr_plot_qc_pacbio_002.png
    :target: auto_examples/plot_qc_pacbio.html


.. raw:: html

   <div class="body">
   <div id="index-grid" class="section group">
   <div class="col span_1_of_3">
        <h3><a href="installation.html">Installation</a></h3>
        <p>pip install sequana</p>
        <h3><a href="userguide.html">Library User Guide</a></h3>
        <p>Use Sequana from Python to read, analyze, and visualize genomic data</p>
        <h3><a href="tutorial.html">Tutorial</a></h3>
        <p>Step-by-step guides for common workflows and pipelines</p>
        <h3><a href="auto_examples/index.html">Code Examples</a></h3>
        <p>Runnable Python examples from the gallery</p>
        <h3><a href="pipelines.html">NGS Pipelines</a></h3>
        <p>Browse the Snakemake pipelines catalogue</p>
        <h3><a href="cli_reference.html">CLI Reference</a></h3>
        <p>The sequana command-line tools (30+ sub-commands)</p>
    </div>
    <div class="col span_2_of_3">
    <div class="jcarousel-wrapper">
    <div class="jcarousel">

* |coverage|
* |fastqc|
* |kraken|
* |bam|
* |pacbio|

.. raw:: html

            </div>
        <a href="#" class="jcarousel-control-prev">&lsaquo;</a>
        <a href="#" class="jcarousel-control-next">&rsaquo;</a>
        <p class="jcarousel-pagination">
        </p>
        </div>
        </div>
        </div>
   </div>
   <div style="clear: left"></div>


.. _quick_start:


.. toctree::
    :caption: Getting started
    :maxdepth: 2

    installation

.. toctree::
    :caption: Library user guide
    :maxdepth: 2

    userguide
    auto_examples/index

.. toctree::
    :caption: Pipelines
    :maxdepth: 2

    pipeline_user_guide
    pipelines
    tutorial
    tutorial_phylo

.. toctree::
    :caption: Standalone applications
    :maxdepth: 2

    applications
    cli_reference

.. toctree::
    :caption: More
    :maxdepth: 2

    case_examples
    faqs
    notebooks
    developers
    sequanix
    wrappers
    references
    glossary


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
