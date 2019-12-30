:Overview: **compressor** can be used to compress/uncompress FastQ files recursively.
   It was designed to uncompress gzipped files and to compress them back into a
   bzip2 format. It was then extended to dsrc and non-compressed files.
   Supported formats are gz and bz2 and dsrc
   (http://sun.aei.polsl.pl/dsrc/download.html).
:Input: Any number of FastQ files (compressed or not)
:Output: The input files (compressed or not)
:requirements: gzip and bzip and their parallel versions (pigz and pbzip2) as
    well as dsrc (DNA compression tool).

Usage
~~~~~~~

A standalone named **sequana_compressor** is provided. Please see::

    sequana_compressor --help

to get detailled information about the arguments. The following example
converts all file ending in fastq.gz into a new compression format (bz2).
Note that "fastq" before the extension is required::

    sequana_compressor --source fastq.gz --target fastq.bz2

If you want to add recursivity, add the ``--recursive`` argument. On a
distributed system (e.g. slurm), you should use the ``--snakemake-cluster``.
For example on SLURM, add::

    --snakemake-cluster "sbatch --qos normal"



**compressor** allows one to go from one format in (fastq, fastq.gz, fastq.bz2,
fastq.dsrc) to any format in the same list. So this is a fully connected
network as shown below:


.. image:: _static/compressor_codecs.png
   :width: 60%

Here is another example::

    sequana_compressor --source fastq.gz --target fastq.bz2 --threads 8
    --snakemake-cluster "sbatch --qos normal"  --recursive --jobs 20

The number of jobs is set to 4 by default and limited to 20 to have a
reasonnable IO access. You can use more using the ``--bypass`` argument.
If nodes have 8 CPUs, use threads=8, this means 20 nodes will be used. 



.. warning:: During the conversion, a .snakemake is created in each processed
   directory. If you interrupt the process, snakemake locks the directory. If
   you get an error message about locked directories, relaunch your previous
   command with ``--unlock`` to unlock the directories and start again

Requirements
~~~~~~~~~~~~~~~~~~
Parallel version of gzip and bzip, as well as dsrc:

.. include:: ../sequana/pipelines/compressor/requirements.txt

Config file
~~~~~~~~~~~~~~

In principle you should not use any config file if you use the standalone.
Note, however, the format of the underlying config file (for pbzip2 and pigz,
the number of threads is automatically set to the number of available threads).

::

    compressor:
        source: fastq.gz
        target: fastq.bz2
        threads: 4
        recursive: True
        verbose: True


DAG
~~~~~~~~~~~

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/compressor/dag.png




Rules used by the pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Depending on the value of the target and source, only one rule is included in the
pipeline. For example if your source is **fastq.gz** and the target **fastq.bz2**, the
**gz_to_bz2** rule is included. Its documentation is here below:

.. snakemakerule:: gz_to_bz2

Others similar rules that convert from one compressed format to another
compressed formats are:

.. snakemakerule:: gz_to_dsrc
.. snakemakerule:: bz2_to_dsrc
.. snakemakerule:: bz2_to_gz
.. snakemakerule:: dsrc_to_gz
.. snakemakerule:: dsrc_to_bz2


