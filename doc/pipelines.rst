
.. _pipelines:

Pipelines
##############

In **Sequana** parlance, a pipeline is an application based on Snakemake that consists of a Snakefile and a configuration file.

For snakemake tutorial, you can have a look at the Snakemake page or
online-tutorials (e.g. http://slowkow.com/notes/snakemake-tutorial/)

Pipelines can be initialised and run via a command line interface called
:ref:`sequana_app` but we would recommend to use :ref:`sequanix` instead.

The following sections are dedicated to each pipeline.

.. toctree::
    :maxdepth: 1

    pipeline_coverage.rst
    pipeline_denovo_assembly.rst
    pipeline_quality_control.rst
    pipeline_rnaseq.rst
    pipeline_smallrnaseq.rst
    pipeline_variant_calling.rst
    pipeline_compressor.rst
    pipeline_pacbio_qc.rst
    pipeline_pacbio_denovo.rst
    pipeline_fastqc.rst


