:Overview: ChIPSeq: Differential expressed peaks analysis
:Input: FastQ raw data from Illumina Sequencer (either paired or not)
:Output: BAM, BED and HTML files



Usage
~~~~~~~~~

Example::

   sequana --pipeline chipseq -i data/ -o analysis --no-adapter -t "_R[12]"
   cd analysis
   sbatch snakemake -s chipseq.rules --stats stats.txt -p -j 12 --nolock --cluster-config cluster_config.json --cluster "sbatch --mem={cluster.ram} --cpus-per-task={threads}"

Or use :ref:`sequanix_tutorial` interface.

Requirements
~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/chipseq/requirements.txt

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/chipseq/dag.png


Details
~~~~~~~~~

ChIPuana is a snakemake-based workflow for the analysis of epigenomic data (ChIPseq) from the raw fastq files to the
differential analysis of transcription factor binding or histone modification marking. It streamlines critical steps
like the quality assessment of the immunoprecipitation using the cross correlation and the replicate comparison for
both narrow and broad peaks. For the differential analysis ChIPuana provides linear and non linear methods for
normalisation between samples as well as conservative and stringent models for estimating the variance and testing the
significance of the observed differences. We show examples of how various settings can allow users to improve the
discriminative power of their comparisons depending on the dynamics of the epigenomic factor under study.




Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/chipseq/config.yaml` to be used with the pipeline.
Each rule used in the pipeline may have a section in the configuration file.
Here are the rules and their developer and user documentation.



FastQC
^^^^^^^^^^^

FastQC is used to check quality of sequenced reads.

.. snakemakerule:: fastqc_dynamic

Cutadapt
^^^^^^^^^

Cutadapt is used to trim and filter sequences.

.. snakemakerule:: cutadapt

Mapping with bowtie2
^^^^^^^^^^^^^^^^^^^^^

Bowtie2 is used to aligned read against reference genome (and spike-in genome if needed)

.. warning:: with paired-end data use `--dovetail --no-mixed --no-discordant` parameters


.. snakemakerule:: bowtie2_mapping_dynamic

Deduplication
^^^^^^^^^^^^

Picard tools Markduplicate is used to remove duplicates from BAM files.

.. snakemakerule:: mark_duplicates_dynamic


Remove blacklist
^^^^^^^^^^^^^^^^

Bedtools is used to remove reads that aligned in blacklisted zones of the genome (optional step)

.. snakemakerule:: remove_blacklist

Peak Calling
^^^^^^^^^^^^

MACS2 is used to call peak on IP Bam files.

.. warning:: For paired-end reads, it is very recommended to put `--keep-dup 2500`

.. snakemakerule::  macs2_dynamic



IDR step
^^^^^^^^^
According to ENCODE guidelines (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/), reproducible peaks are selected
using IDR pipeline.


.. snakemakerule:: compute_idr

Summarise results
^^^^^^^^^^^^^^^^^^^

MultiQC (https://multiqc.info/docs/) is used to summarise pipeline steps in an HTML file.

.. snakemakerule:: multiqc

.. warning:: a custom config file is used to help you during data interpretation