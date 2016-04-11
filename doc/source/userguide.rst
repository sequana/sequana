User Guide
############

.. contents::

Snakemake tools
================

We provide so utilities in the :mod:`sequana.snakemake` to be used within
Snakefile. Although Python code can be added the print command output may
interfer with the Snakefile interpretation. The :func:`sequana.snakemake.message` is provided to prevent that issue (see example below).


Another feature of Sequana is to provide a set of modules that are independent
and generic enough. For instance, Snakemake has a nice option to create the dag
file out of the snakefile. However, one has to then use **dot** command and
possibly other commands. In Sequana, we also have a dot parser to annotate the
dot (useful in HTML document for clickable image). The module that takes care of
that is called ... dag and can be included in Your Snakefile as follows::

    from sequana import snakemake as sm

    sm.message("Include dag rule in this pipeline" )

    include: sm.rules['dag']

Save this example in a file called Snakefile and type::

    snakemake

If the analysis is successful, you can open the report in ::

    report/index.html

If you change the config file or wish to restart, you will need to force the
analysis to overewrite existing files::

    snakemake --forceall

Or simply delete all files that have been created::

    snakemake cleanup

.. todo:: explain the cleanup module.


Sequana Toolkit
====================

We provide a set of tools to perform post-analysis on standard files (e.g., BAM,
FastQ, VCF). As an example, we will play here below with a smaple of BAM file:

.. plot::
    :include-source:

    from sequana import sequana_data, BAM
    b = BAM(sequana_data("test.bam", "testing"))
    b.plot_bar_mapq()


Sequana Reports
==================

In snakemake, we provide also with different snakefile and standard formats,
the ability to create or generate reports. Here is a BAM report file::

    from sequana import BAM, sequana_data, BAMReport
    b = BAM(sequana_data("test.bam", "testing"))

    r = BAMReport()
    r.set_data(b)
    r.create_report()

that results can be shown in `report/bam.html <_static/report/bam.html>`_

