Sequana documentation
##########################################

|version|, |today|


.. raw:: html

    <div style="width:80%"><p>


    <a href="http://bioconda.github.io/recipes/sequana/README.html">
    <img src="https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square"></a>

    <a href="https://pypi.python.org/pypi/sequana">
    <img src="https://badge.fury.io/py/sequana.svg"></a>

    <a href="https://travis-ci.org/sequana/sequana">
    <img src="https://travis-ci.org/sequana/sequana.svg?branch=master"></a>

    <a href="https://coveralls.io/github/sequana/sequana?branch=master">
    <img src="https://coveralls.io/repos/github/sequana/sequana/badge.svg?branch=master"></a>

    <a href="http://sequana.readthedocs.org/en/master/?badge=master">
    <img src="http://readthedocs.org/projects/sequana/badge/?version=master"></a>

    <a href="http://joss.theoj.org/papers/10.21105/joss.00352">
    <img src="http://joss.theoj.org/papers/10.21105/joss.00352/status.svg"></a>

    <a href="https://gitter.im/sequana/sequana?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge">
    <img src="https://badges.gitter.im/sequana/sequana.svg">


    </p>
    </div>


:Python version: Python 3.5, 3.6 although most modules are Python2.7 compatible.
:Source: See  `http://github.com/sequana/sequana <https://github.com/sequana/sequana/>`_.
:Issues: Please fill a report on `github <https://github.com/sequana/sequana/issues>`_
:How to cite: For Sequana in general including the pipelines, please use

    Cokelaer et al, (2017), 'Sequana': a Set of Snakemake NGS pipelines, Journal of
    Open Source Software, 2(16), 352, `JOSS DOI doi:10.21105/joss.00352 <http://www.doi2bib.org/bib/10.21105%2Fjoss.00352>`_


    For the **genome coverage** tool (sequana_coverage), please cite:

    Dimitri Desvillechabrol, Christiane Bouchier, Sean Kennedy, Thomas Cokelaer
    http://biorxiv.org/content/early/2016/12/08/092478

    For **Sequanix** (GUI for Snakemake pipeline), please cite:

    Dimitri Desvillechabrol, Rachel Legendre, Claire Rioualen,
    Christiane Bouchier, Jacques van Helden, Sean Kennedy, Thomas Cokelaer
    Sequanix: A Dynamic Graphical Interface for Snakemake Workflows
    Bioinformatics, bty034, https://doi.org/10.1093/bioinformatics/bty034
    Also available on bioRxiv(DOI: https://doi.org/10.1101/162701)


What is Sequana ?
=====================

**Sequana** is a versatile tool that provides

#. A Python library dedicated to NGS analysis (e.g., tools to visualise standard NGS formats).
#. A set of :ref:`pipelines <Pipelines>` dedicated to NGS in the form of Snakefiles
   (Makefile-like with Python syntax based on snakemake framework) with more
   than 80 re-usable rules (see :ref:`rules`).
#. Original tools to help in the creation of such pipelines including HTML reports.
#. :ref:`Standalone applications<applications>`:
    #. :ref:`sequana_coverage<standalone_sequana_coverage>` ease the
       extraction of genomic regions of interest and genome coverage information
    #. :ref:`sequana_taxonomy<standalone_sequana_taxonomy>` performs a quick
       taxonomy of your FastQ. This requires dedicated databases to be downloaded.
    #. :ref:`Sequanix`, a GUI for Snakemake workflows (hence Sequana pipelines as well)

Currently, the available pipelines cover quality control (e.g. adapters removal,
phix removal, trimming of bad quality bases), variant calling, characterisation
of the genome coverage, taxonomic classification, de-novo assembly, RNA-seq. See the :ref:`pipelines`
section for more information.

**Sequana** can be used by developers to create new pipelines and by users in the
form of applications ready for production. Moreover, **Sequanix** can be used to
set the parameters of pipelines and execute them easily with a graphical user
interface.

To join the project, please let us know on `github <https://github.com/sequana/sequana/issues/306>`_.


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
    rules
    references
    faqs
    Changelog
    glossary


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

