Sequana documentation
##########################################

|version|, |today|

.. include:: ../README.rst



What is Sequana ?
=====================

**Sequana** is a versatile tool that provides

#. A Python library dedicated to NGS analysis (e.g., tools to visualise standard NGS formats).
#. A set of :ref:`pipelines <Pipelines>` dedicated to NGS in the form of Snakefiles
   (Makefile-like with Python syntax based on snakemake framework).
#. Original tools to help in the creation of such pipelines including HTML reports.
#. :ref:`Standalone applications<applications>`:
    #. :ref:`sequana_coverage<standalone_sequana_coverage>` ease the
       extraction of genomic regions of interest and genome coverage information
    #. :ref:`sequana_taxonomy<standalone_sequana_taxonomy>` performs a quick
       taxonomy of your FastQ. This requires dedicated databases to be downloaded.
    #. :ref:`Sequanix`, a GUI for Snakemake workflows (hence Sequana pipelines as well)

The sequana pipelines are various. Since March 2020, they have their own independent life within dedicated github repositories. You may find pipelines for NGS quality control (e.g. adapters removal,
phix removal, trimming of bad quality bases), variant calling, characterisation
of the genome coverage, taxonomic classification, de-novo assembly, 
:ref:`Variant calling <pipeline_vc>`, :ref:`RNA-seq <pipeline_rnaseq>`, etc. See the :ref:`pipelines`
section for more information.


**Sequana** can be used by developers to create new pipelines and by users in the
form of applications ready for production. Moreover, **Sequanix** can be used to
set the parameters of pipelines and execute them easily with a graphical user
interface.

To join the project, please let us know on `github <https://github.com/sequana/sequana/issues/306>`__.


.. Here we are building the carrousel? Note that html and pdf version look for
   images in different folders...

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

.. |sequanix| image::
    _static/sequanix.png
    :target: applications.html#sequanix

.. |pacbio| image::
    ./auto_examples/images/sphx_glr_plot_qc_pacbio_002.png
    :target: auto_examples/plot_qc_pacbio.html



.. raw:: html

   <div class="body">
   <div id="index-grid" class="section group">
   <div class="col span_1_of_3">
        <h3><a href="installation.html">Installation</a></h3>
        <p>conda install sequana</p>
        <h3><a href="auto_examples/index.html">Examples</a></h3>
        <p>Visit our example gallery to use the Python library</p>
        <h3><a href="pipelines.html">NGS pipelines</a></h3>
        <p>Learn about available Snakemake pipelines</p>
        <h3><a href="applications.html">Standalone applications</a></h3>
        <p>Standalone applications including Sequanix (GUI for snakemake)
        and the sequana_coverage tool.</p>
    </div>
    <div class="col span_2_of_3">
    <div class="jcarousel-wrapper">
    <div class="jcarousel">

* |coverage|
* |fastqc|
* |kraken|
* |bam|
* |sequanix|
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


User guide and reference
###########################


.. toctree::
    :numbered:
    :maxdepth: 2

    installation.rst
    userguide
    tutorial
    pipelines
    auto_examples/index
    case_examples
    applications
    sequanix.rst
    developers
    wrappers
    references
    references_enrich
    references_stats
    references_viz
    faqs
    glossary




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

