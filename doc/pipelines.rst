
.. _pipelines:

Pipelines
##############

Sequana ships many pipelines related to NGS. Some are very simple (fastqc,
demultiplex, downsampling), others are real-life NGS pipelines used in
production.

.. warning:: Since v0.8.0, pipelines are now independent from **Sequana**. 
    They must be installed separetely and their dependencies must also be
    installed by the user/developer.

Quick Start
===========

If **Sequana** is installed, installing a pipeline is straightforward. For
example, to install the variant calling pipeline::

    pip install sequana_variant_calling --upgrade

Since version 0.8.1, you can check whether you have the required dependencies.
If not, an error message will appear anyway::

    sequana_variant_calling --deps

Once those dependencies are available, you can run the pipeline::

    sequana_variant_calling --help



Overview
========

In **Sequana** parlance, a pipeline is an application based on Snakemake that consists of a Snakefile and a configuration file. Although since v0.8.0, we augmented a pipeline with other optional files such as a schema to check the config file, a logo, a dag image representing the pipeline, a requirements file with external dependencies and so on. 

All pipelines are based on Snakemake. For a tutorial, you can have a look at the Snakemake page or
online-tutorials (e.g. http://slowkow.com/notes/snakemake-tutorial/).


.. note:: **Pipeline naming convention**

    A pipeline is named **sequana_pipelines_name** where name is to be replaced by the
    pipeline name. The name can contain underscores. For instance, the
    **variant_calling** pipeline is called **sequana_pipelines_variant_calling**.
    Actually, we have aliases and pipeline have usually a short name where
    **_pipelines** is dropped. So you can refer to a pipeline as
    **sequana_pipelines_variant_calling** or **sequana_variant_calling**. The reason 
    for having the long and short versions is to avoid conflict name with Sequana 
    standalones. For instance, the **sequana_coverage** tool 
    exists. It is a standalone that study the coverage on a unique sample. We 
    also have a pipeline to analyse several samples in parallel. Therefore the 
    sequana_pipeline_coverage pipeline has no alias.

    Future version will use the short version only. 


Installation
============

Given its name, and provided you have installed Sequana, you can install the
pipeline **name** using::

    pip install sequana_name

where name is replaced by the pipeline name. For instance::

    pip install sequana_fastqc

Since, you want to be up-to-date, add the --upgrade argument::

    pip install sequana_fastqc --upgrade

Usage
======

Each pipeline is different. We recomment to look this complementary section
:ref:`facilitator`. Generally speaking, the --help argument should be sufficient
to run most of the pipelines::

    sequana_name --help

The input arguments --input-directory, --input-pattern and --input-readtag will
help you selecting the input data for the pipeline. Then, you will have to
introspect the help and the documentation of the pipelines. Each pipeline has
its own repositoty and living documentaion, which are available in the link here
below.

List pipelines
==============
This is a non-exhaustive list of pipelines

.. toctree::
    :maxdepth: 1

    pipeline_demultiplex.rst
    pipeline_fastqc.rst
    pipeline_quality_control.rst
    pipeline_rnaseq.rst

Please see the https://github.com/sequana organisation to get the full list.
