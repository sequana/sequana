Overview
############

.. contents::

**Sequana** provides standalone applications (e.g., **sequana_coverage**,
**sequana_taxonomy**) and pipelines in the form of Snakefiles. Although the standalone
applications are usually simpler, they may not have all features or parameters
offered by the pipelines. Since version 0.8.0, most pipelines have been moved to
different repository with one repository per pipeline. This was done to make
pipelines independe and the Sequana more modular and effective for deployment in
production mode. 

The :ref:`Tutorial`, :ref:`pipelines`, :ref:`case_examples` 
sections provide many examples on their usage. Check also the Gallery section
for code snippets.

This section will not describe all available standalones and pipelines.
We will focus on one example (coverage) to show how one can use
the **Sequana** library, or standalone application, or pipeline to get
information about the coverage of a set of mapped reads onto a reference.


**Sequana** is a Python library
===============================

Example 1 : running median on coverage
----------------------------------------

**Sequana** is a Python library. It contains many functionalities, which are
fully documented and available in the :ref:`references` section. We can first
look at the coverage contained within a BED file using the library. First, we
need some data. **Sequana** provides some test examples, which can be accessed
using :func:`~sequana.sequana_data` function. The test case is a virus (about
18,000 bases)::

    from sequana import sequana_data
    filename = sequana_data('JB409847.bed')


We can then use the :class:`~sequana.bedtools.GenomeCov` class to read the
file::

    from sequana import GenomeCov
    gc = GenomeCov(filename)

Select a chromosome (first one) and compute the running median::

    chrom = gc[0] 
    chrom.running_median(n=5001, circular=True)
    chrom.compute_zscore()

and finally plot the coverage together with confidence interval (3 sigma)::

    chrom.plot_coverage()


.. plot::

    from sequana import sequana_data
    filename = sequana_data('JB409847.bed')
    from sequana import GenomeCov
    gc = GenomeCov(filename)

    chrom = gc[0]
    chrom.running_median(n=5001, circular=True)
    chrom.compute_zscore()
    chrom.plot_coverage()

.. seealso:: notebook available in the `github repository
   <https://github.com/sequana/sequana/blob/master/notebooks/coverage.ipynb>`_


As you can see, **Sequana** is a standard Python library where developers can
select functions, classes, modules to help them design new tools and pipelines.


Example2: read a fastq file
------------------------------

Let us use the :class:`FastQC` class to get the distribution of the bases ACGT
across all reads of a FastQ file.


.. plot::
    :include-source:

    from sequana import FastQC
    from sequana import sequana_data
    filename = sequana_data("test.fastq")

    fastqc = FastQC(filename)
    print(fastqc.fastq)
    for x in 'ACGT': 
        fastqc.get_actg_content()[x].hist(alpha=0.5, label=x, histtype='step', lw=3, bins=10)

    from pylab import legend
    legend()



Many more functionalities are available. The reference guide should help you.

**Sequana** provides standalone applications
============================================

The Python example about the coverage is actually quite useful. We 
therefore decided to provide a standalone
application. There are other standalone applications listed in
:ref:`applications` section.

The one related to the coverage example shown above is named
**sequana_coverage**. If you have a BED file, type::

    sequana_coverage  -i <BEDFILENAME> 

If your organism has a circular DNA, add ``-o``. You can play with the window
size for the running median using ``-w``.

Using the BED file and reference mentionned in the previous section you should
obtain the same figure as above.

An additional feature is the report using  ``--show-html`` option.

.. _facilitator:

**Sequana**, a pipeline construction facilitator
==================================================

In **Sequana**, in addition to the library and standalone applications, we also
provide a set of pipelines (see :ref:`pipelines` section). Originally, pipeline
were provided with Sequana, inside the same source repository. Since version
0.8.0, pipeline have their own repository. For instance, 
:ref:`pipeline_vc` is available on
https://github.com/sequana/variant_calling.
We will not describe all pipelines here below since new ones may appear now and
then. Instead, let us explain the way pipelines can be designed and run.

Installation of a pipeline
--------------------------

With the new design implemented in v0.8.0, pipelines are independent Python
packages posted on Pypi. You can now install a pipeline (e.g., variant calling)
as follows in your virtual environment::

    pip install sequana_variant_calling --upgrade

The --upgrade option is to make sure you install the newest version.

To check if the installation is successful, just type::

    sequana_ariant_calling --help


Usage
-----

A very simple and useful pipeline for this explanation is the
**sequana_fastqc** pipeline. Install it as follows::

    pip install sequana_fastqc

and check the help message::

    sequana_fastqc --help

You will see 4 sections some of which are common to all **Sequana** pipelines.

The generic section allows use to print the help with --help, to set the level
of information printed to the screen (--level), the version (--version).
Pipelines can be run locally or on a SLURM clusters. This can be set with the
--run-mode option. Note, however, that this option is set automatically to
slurm-mode if slurm commands are found (e.g. sbatch). 

The *slurm* section can be used to set slurm options for Snakemake. If you do
not know what it means, let it be the defaults values. Just note that memory
usage is set to 4Gb by default and number of cores is limited to 4 per job.

The *snakemake* section allows you to set to maximum number of jobs to be used,
which is set to 4 (if run-mode is set to local) and 40 (if run-mode is set to
slurm). 

The --working-directory is set to the name of the pipeline and is the important
parameter. It tells sequana where to store the pipeline files (e.g., snakemake,
configuration files). You can change it to your will but if it exists already,
the pipeline zill not be set up and you will need to use the --force option to
overwrite existing files.

The next section is about your input data. Most of the pipelines expect to find
Illumina data with single or paired-end data sets. The directory where to find
the data is defined by the --input-directory parameter. You can refine the
search by providing an input pattern, which is set to `*fastq.gz` by default.
Since, Illumina data may be paired, we have a mechanism to check and discovered
paired data for each sample. By default, the paired data are differentiate
thanks to a pattern _R1_ or _R2_ to be found in the filenames. The common
pattern set with --input-readtag is set to _R[12]_ but can be easily changed.
For instance if your files do not contain the R or if the _R1 is to be found at
the end of the file, just change it accordingly.


So, let us now perform the fastqc of a bunch of samples. You could type::

    sequana_fastqc --input-directory my_data_directory --working-directory test1

This will copy the snakefile, the configuration files and useful files to run
the analysis. Follow the instructions that is::

    cd test1

In this directory, you can find  The configuration file called
**config.yaml**. This pipeline is very simple but you can see the parameters related to
your input data::

    input_directory: /home/login_example/data_example/my_data_directory
    input_readtag: _R[12]_
    input_pattern: '*fastq.gz'

So you can edit this file to correct it or change other parameters. If you are
happy with those choices, it is now time to run the pipeline. If you know
snakemake, you can just use it. For example::

    snakemake -s fastqc.rules

or just type::

    sh fastq.sh

Wait and see. Once done. If every went well, you can keep the configuration
files and pipeline-related files, or delete them using::

    make clean

.. seealso:: :ref:`pipelines` section for more information.



Using **Sequanix** standalone
---------------------------------

An even easier way is to use our graphical interface named **Sequanix**. A
snapshot can be found in the :ref:`sequanix` section and a tutorial in
:ref:`sequanix_tutorial`.

Note, however, that the Sequanix interface is slightly different. The content of
the working directory may differ slightly for the time being. The advantage of
using Sequanix is that complex configuration pipeline can be tuned easily
through its graphical interface.



**Sequana** Reports
=====================


Pipelines and standalone make use of internal reporting. Since they are part of
the **Sequana** library, they can also be used with your own code. For instance,
if you have a BAM file, you can use the following code to create a basic
report::

    from sequana import BAM, sequana_data
    from sequana.modules_report.bamqc import BAMQCModule
    filename = sequana_data("test.bam", "testing")
    r = BAMQCModule(filename, "bam.html")

that results can be shown in `bam.html <_static/bam.html>`_

