:Overview: Quality control, trimming (adapter removal) and taxonomic overview
:Input: A set of FastQ files (paired or single-end)
:Output: fastqc, cleaned FastQ files, Krona plot

Usage
~~~~~~~

::

    sequana --pipeline quality_control --input-directory . --working-directory analysis --no-adapters


Or use :ref:`sequanix_tutorial` interface.

Requirements
~~~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/quality_control/requirements.txt

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/quality_control/dag.png


Details
~~~~~~~~~

This pipeline is used to check the quality of the input FastQ files. It is used
to remove the Phix (coliphage) that may be present in the data. Low quality reads
are trimmed using a dedicated tool such as **cutadapt** or **atropos**. If one specifies
the quality trimming option in the config file, then we trim low-quality ends from 
reads BEFORE adapter removal.

The quality trimming algorithm from cutadapt/atropos is the same as in BWA. That is: substract the
cutoff (e.g. 30) from all qualities; compute partial sums from the end of the
sequence; cut the sequence at the index at which the sum is minimal.

::

    # Original qualities
    42, 40, 26, 27, 8, 7, 11, 4, 2, 3
    # Subtracting the threshold gives:
    32, 30, 16, 17, -2, -3, 1, -6, -8, -7
    # Partial sum from the end. Stop early if the sum is greater than zero:
    (70), (38), 8, -8, -25, -23, -20, -21, -15, -7

Minimum is -25, we keep the bases 1,2,3,4::

    42, 40, 26, 27

Another important point is that all searches for adapter
sequences are error tolerant (allowing errors such as
mismatches, insertions and deletions). The level of error tolerance
is 10% by default.


In the cutadapt configuration section, you can provide the type of adapters to
remove such as PCRFree, Nextera, TruSeq, etc. Note, however, that this removes
all adapters. It is useless and may take more time than expected. Instead, you
can provide an *experimental design* file that associate each sample with
its adapter. See the reference for more information :class:`~sequana.expdesign.ExpDesignAdapter`
Typically, you will provide a Sample Sheet (generated with IEM tool from
Illumina), which looks like::

    [Header]
    IEMFileVersion,4
    Investigator Name,lauma
    Experiment Name,181207_PE-MPdbl_300v3_std
    Date,12/6/2018
    Workflow,GenerateFASTQ
    Application,FASTQ Only
    Assay,Nextera XT
    Description,
    Chemistry,Amplicon

    [Reads]
    300
    300

    [Settings]
    ReverseComplement,0
    Adapter,CTGTCTCTTATACACATCT

    [Data]
    Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
    AL_LB_C_1,AL_LB_C_1,TP-third,A01,N701,TAAGGCGA,S517,GCGTAAGA,,
    AL_LB_C_2,AL_LB_C_2,TP-third,A02,N702,CGTACTAG,S517,GCGTAAGA,,
    AL_LB_C_3,AL_LB_C_3,TP-third,A03,N703,AGGCAGAA,S517,GCGTAAGA,,
    AL_LB_C_4,AL_LB_C_4,TP-third,A04,N704,TCCTGAGC,S517,GCGTAAGA,,




Another optional step is the taxonomy analysis. This is performed with Kraken
using a dedicated database. We provide a couple of database or tools to download
them. See sequana_taxonomy tool to download databases. The minikraken database
is provided by the author of Kraken while sequana_db1 is a 8Gb database as
described here https://github.com/sequana/data/tree/master/sequana_db1 ::

    sequana_taxonomy --download minikraken
    sequana_taxonomy --download sequana_db1

For the second database, you will need **synapseclient**::

    pip install synapseclient

and an account on synapse website (https://www.synapse.org/).

FAQS
~~~~~~~~~~~~~~

If you use the quality_control pipeline for the first time, you will need a DB
of taxons, which is downloaded automatically. 

However, if you are on a server that do not allow network connection, it will
not be downloaded and an error may occur. If you have a network connection on
the head / master where you start sequanix/sequana, then you will need to
download the DB yourself manually. To do so, open an IPython or Python shell and
type::

    from biokit import Taxonomy
    tax = Taxonomy(verbose=True)
    tax._load_flat_file()


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/quality_control/config.yaml` to be used with the pipeline. Each rule used in the pipeline may have a section in the
configuration file. In the *quality_control* pipeline, we use the *bwa_mem*,
*fastqc*, *cutadapt*, *fastq_stats* and *kraken* rules described here below.


FastQC
^^^^^^^^^^^
.. snakemakerule:: fastqc_dynamic

Cutadapt
^^^^^^^^^
.. snakemakerule:: cutadapt

Kraken
^^^^^^^
.. snakemakerule:: kraken

BWA-mem 
^^^^^^^^
.. snakemakerule:: bwa_mem_dynamic

FastQ stats
^^^^^^^^^^^^
.. snakemakerule:: fastq_stats_dynamic

