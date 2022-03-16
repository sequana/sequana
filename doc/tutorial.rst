.. _tutorial:

Tutorial
==========

Here below you can find tutorials on how to use some of the sequana pipelines
and standalones. This gives you a flavor of what can be achieved using
**Sequana**. Dedicated sections on pipelines and applications are found in other
pages and can be a complement to this quick overview. 

.. contents::
   :depth: 2


The standalone Sequana
----------------------

New since version 0.9.0. We are a single entry point for a set of tools used in
pipelines or as standalone applications. You can type::

    sequana --help 


to get the list of applications. Would you need completion, this is possible
using e.g. for bash users::

    eval "$(_SEQUANA_COMPLETE=source_bash sequana)"

The fastqc pipeline
--------------------

The following example will show how to run the fastqc pipeline 
(https://github.com/sequana/fastqc) on a pair of
FastQ files. The data comes from a sequencing (using HiSeq technology) of a
Measles virus. For testing purposes, you can download :download:`R1
<../sequana/resources/data/Hm2_GTGAAA_L005_R1_001.fastq.gz>` and
:download:`R2 <../sequana/resources/data/Hm2_GTGAAA_L005_R2_001.fastq.gz>`)
files that contain only 1500 reads. Copy them in a local directory.

Those files are from an HiSeq2500 run. The adapters are PCRFree. There is
only one sample for which the index is GTGAAA. You should have 10% of adapters.

Then, initiate the pipeline::

    sequana_fastqc --input-directory . 
    cd fastqc
    sh fastq.sh

Open the summary.html file that is generated for you.



Quality Control pipelines
--------------------------

The quality_control pipeline  (https://github.com/sequana/quality_control)
is not maintained anymore and has been split into several smaller pipelines.

For book-keeping, we keep this section though.


The following example will show how to run the quality control pipeline
(https://github.com/sequana/quality_control) on a pair of
FastQ files. The data comes from a sequencing (using HiSeq technology) of a
Measles virus. For testing purposes, you can download :download:`R1
<../sequana/resources/data/Hm2_GTGAAA_L005_R1_001.fastq.gz>` and
:download:`R2 <../sequana/resources/data/Hm2_GTGAAA_L005_R2_001.fastq.gz>`)
files that contain only 1500 reads. Copy them in a local directory.

Those files are from an HiSeq2500 run. The adapters are PCRFree. There is
only one sample for which the index is GTGAAA. You should have 10% of adapters.

Make sure you have installed the pipeline::

    pip install sequana_quality_control --upgrade

This example show hows to initialise and run the quality control
pipeline on a pair of FastQ files. Copy the two data files (link above) into the
local directory where you will initiate the pipeline.

First, run the sequana standalone application to initialise the pipeline
**quality_control**::

    sequana_quality_control  --cutadapt-adapter-choice TruSeq

Since your data is in the current directory, no need to provide the
--input-directory argument for now since its default value is *.*. 

This command fill the required configuration file(s) and copy it along the
pipeline itself inside the default working directory (quality_control)

The pipeline does 3 things:

1. remove the Phix if present
2. apply cutadapt to trim the bases with quality below 30 and removes adapters 
   (here TruSeq)
3. taxonomy if Kraken databases are provided.

in particular
the config file and the pipeline itself. This example should work out of
the box but you may want to look at the
configuration file **config.yaml**. For instance, you may want to change the
reference to the *phix* (by default we use *phix174.fa*, which is provided in
Sequana) or
change the adapter_removal section to your needs (cutadapt parameters, in
particular the forward and reverse complement list of adapters; None by
default).

By default, the output directory is called **quality_control** and can be overwritten
with the ``--working-directory`` parameter. Then, run the pipeline and wait for
completion.::

    cd quality_control
    snakemake -s quality_control.rules --stats stats.txt -p -j 4 --forceall

The -p option shows the commands, -j 4 means use 4 threads when possible.
Alternatively, there is also a **runme.sh** script.

.. note:: you can also use the shell script **sh quality_control.sh** instead of
   the snakemake command.

You should now have a directory with a HTML report corresponding to the sample::

    open index.html




Taxonomy (standalone)
----------------------

To perform a quick taxonomy of your reads, you can use :ref:`standalone_sequana_taxonomy`
either from Python or as a standalone.

Here we show how to use the Python approach (see :ref:`standalones`) for the
other approach.

Download a toy kraken database designed for this problem (contains only 100
FASTA files mixing measles viruses and others viruses)::


    from sequana import KrakenDownload, sequana_config_path
    kd = KrakenDownload()
    kd.download("toydb")
    database_path = sequana_config_path + "/kraken_toydb"

Then, you may use the following code to perform the analysis (using :mod:`sequana.kraken`)::

    from sequana import KrakenPipeline
    kp = KrakenPipeline(["R1.fastq.gz", "R2.fastq.gz"], database="~/.config/sequana/kraken_toydb")
    kp.run()

Alternatively, you can use the standalone application::

    sequana_taxonomy  --file1 Test_R1.cutadapt.fastq.gz
        --file2 Test_R2.cutadapt.fastq.gz --database  <database_path>



Open the local HTML file taxonomy/kraken.html. An example is available
in  `Krona example <_static/krona.html>`_


Variant calling pipeline
--------------------------

The following example will show how to initialise and run the variant calling
pipeline on a pair of FastQ files.
For testing purposes, you can download :download:`R1
<../sequana/resources/data/Hm2_GTGAAA_L005_R1_001.fastq.gz>` and
:download:`R2 <../sequana/resources/data/Hm2_GTGAAA_L005_R2_001.fastq.gz>`)
files that contain only 1500 reads. Copy them in a local directory.

Note that this does the variant calling + snpEff + coverage.
See more information in the :ref:`pipeline_vc` section.

Make sure you have installed the pipeline::

    pip install sequana_variant_calling --upgrade

The variant calling requires input files. Since you want to map your reads onto
a reference, you must have a reference. Besides, you may want to annotate your
results with a specific annotation file. So, let us download those files first.

Get the genbank reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can use `BioServices <https://bioservices.readthedocs.io/en/master/>`_ to
download those files.


Assuming the reference is **K01711.1** (Measles virus), we first need to fetch
the genbank file from NCBI::

    from bioservices import EUtils
    eu = EUtils()
    data = eu.EFetch(db="nuccore",id="K01711.1", rettype="gbwithparts", retmode="text")
    with open("measles.gbk", "w") as fout:
        fout.write(data.decode())

Get the FASTA reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We will also get the FASTA from ENA::

    from bioservices import ENA
    ena = ENA()
    data = ena.get_data('K01711', 'fasta')
    with open("measles.fa", "w") as fout:
        fout.write(data.decode())


Assuming the genbank and reference have the same name, you can simply
type::

    from sequana.snpeff import download_fasta_and_genbank
    download_fasta_and_genbank("K01711", "measles")

.. Get a snpEff config file and update it
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Then you need to initialise a config file for snpEff tool::
       from sequana import snpeff
       v = snpeff.SnpEff("measles.gbk")




Run the pipeline
~~~~~~~~~~~~~~~~~~~~


::

    sequana_variant_calling --input-directory . --reference measles.fa --annotation measles.gbk 
    cd variant_calling
    sh variant_calling.sh

Wait and see. If the run is succesful, you can just type ::

    make clean

to remove some temporary files. Finally, open the file **index.html** and
explore summary HTML report pages (multiqc page). Then, you can go to individual
HTML report page for each sample. The individual report page are in
**report_SAMPLENAME/summary.html**.

About the configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We strongly recommend to look at the configuration file **config.yaml** and to
check or change the parameters according to your needs. In principle, the
reference and annotation file have been set up for you when initiating the
pipeline. 

For example, you should see those lines at the top of the config file::

    annotation_file: measles.gbk
    reference_file: measles.fa

.. warning:: In the configuration file, in the mark_duplicates section,
    some output files are huge and requires temporary directory on cluster.

.. warning:: in the configuration file (coverage section), 
    you may need to decrease the window size for short genomes.


De novo
-------

The denovo_assembly pipeline can be initialised in the same way::

    sequana_denovo --input-directory . --working-directory denovo_test

Go to the **denovo_test** directory and edit the config file. 

.. warning:: this is very time and computationally expensive. The
   **digital_normalisation** section is one that controls the memory footprint.
   In particular, you can check change max-tablesize to a small value for
   test-purposes (set the value to 3e6)




RNA-seq
-------------------


See more information in the :ref:`pipeline_rnaseq` section.
The following example will show you how to initialise and run the RNAseq pipeline on a couple of FastQ files (in single-end mode).
The data comes from a sequencing (using HiSeq2500 technology) of a saccharomyces cerevisiae strain.
For testing purposes, you can download :download:`Fastq1
<../sequana/resources/data/WT_ATCACG_L001_R1_001.fastq.gz>` and
:download:`Fastq2 <../sequana/resources/data/KO_ATCACG_L001_R1_001.fastq.gz>`)
files that contain only 100,000 reads. Copy them in a local directory.


Initialise the pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Call **sequana** standalone as follows::

    sequana_rnaseq --working-directory EXAMPLE

This command download the pipeline and its configuration file. The configuration
file is prefilled with adapter information and input data files found in the
input directory provided. You can change the configuration afterwards.

Go to the project directory and execute the script
::

    cd EXAMPLE
    sh rnaseq.sh


Get the fasta and GFF reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Assuming the reference is **Saccer3** (Saccharomyces cerevisiae), we first need to fetch
the fasta and the GFF files from SGD before to run the pipeline::

    mkdir Saccer3
    cd Saccer3
    wget http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz
    tar -xvzf chromFa.tar.gz
    cat *.fa > Saccer3.fa
    wget http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff -O Saccer3.gff
    rm -f chr*
    cd ..

.. warning:: All files (fasta, GFF, GTF...) used in RNA-seq pipeline must have 
    the same prefix (Saccer3 in the example) and must be placed in a new directory, 
    named as the prefix or not.

.. warning:: For the counting step, the RNA-seq pipeline take only GFF files. GTF and SAF files will be integrated soon.

Initiate the pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:: 

    sequana_rnaseq --genome-directory Saccer3  --aligner bowtie2


Run the pipeline
~~~~~~~~~~~~~~~~~~~~

On local::

    snakemake -s rnaseq.rules --stats stats.txt -p -j 12 --nolock

on SGE cluster::

    snakemake -s rnaseq.rules --stats stats.txt -p -j 12 --nolock --cluster-config cluster_config.json
    --cluster "qsub -l mem_total={cluster.ram} -pe thread {threads} -cwd -e logs -o logs -V -b y "

on slurm cluster ::

    sbatch snakemake -s rnaseq.rules --stats stats.txt -p -j 12 --nolock --cluster-config cluster_config.json
    --cluster "sbatch --mem={cluster.ram} --cpus-per-task={threads} "



Singularity and Sequanix
----------------------------

.. warning:: FOR LINUX USERS ONLY IF YOU WANT TO USE SEQUANIX. YOU CAN STILL USE
   THE SEQUANA STANDALONE

Here we will use a singularity container to run Sequanix and the quality pipeline to analyse
local data sets stored in your /home/user/data directory.

First, Install singularity (http://singularity.lbl.gov/). Check also the
:ref:`Installation` for information.

Second, download this specific container::

    singularity pull --name sequana.img shub://sequana/sequana

This is about 1.5Go of data. Once downloaded, you can play with the container in
**shell** or **exec** mode. 

**shell** mode means that you enter in the container where you have an
isolated environement. Because the isolated environment is protected, only the
directory from where you start singularity, and optional bound directories are
writable. So, if you want to read/write data in a specific directory, you must
use the -B option (see section bind path here below)::

    singularity shell -B /home/user/data/:/data sequana.img

Once in the container, you should see a prompt like this::

    Singularity: Invoking an interactive shell within container...
    Singularity sequana-sequana-release_0_5_2.img:~/Work/github/sequana/singularity>

Just move to the *data* directory::

    cd data

You should see your input files. You can now analyse your data following the
quality pipeline tutorial (top of the page), or use Sequanix::

    sequanix -i . -w analysis -p quality_tutorial

In **exec** mode, this is even simpler::

    singularity exec sequana.img sequanix

or with pre-filled parameters:: 

    sequanix -i . -w analysis -p quality_tutorial

A Sequanix window should appear. You can now follow the Sequanix tutorial
:ref:`sequanix`


binding path (Mounting)
~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have data on a non standard path or want to mount a path so that the
container can see it, use the binding method (see also above). 

Imagine that your data on the host machine is located on /projets/1/data and
that the file to analyse is called virus.bed, you can use the sequana_coverage
tool as follows to analyse your data::

    singularity exec -B /projets/1/data/:/data sequana.simg sequana_coverage --input /data/virus.bed

Here we bind the /projects/1/data directory (host) on the /data directory
available in the container. Other directories available within the container are
/mounting and /scratch.








