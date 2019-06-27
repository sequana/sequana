.. _tutorial:

Tutorial
==========

The following example will show how to run the quality control on a pair of
FastQ files. The data comes from a sequencing (using HiSeq technology) of a
Measles virus. For testing purposes, you can download :download:`R1
<../sequana/resources/data/Hm2_GTGAAA_L005_R1_001.fastq.gz>` and
:download:`R2 <../sequana/resources/data/Hm2_GTGAAA_L005_R2_001.fastq.gz>`)
files that contain only 1500 reads. Copy them in a local directory.

Those files are from an HiSeq2500 run. The adapters are PCRFree. There is
only one sample for which the index is GTGAAA. You should have 10% of adapters.


Quality pipeline
---------------------

**Sequana** comes with standalone applications and pipelines in the form of
Snakefile (`snakemake <https://bitbucket.org/snakemake/snakemake/wiki/Home>`_)

The following example will show how to initialise and run the quality control
pipeline
on a pair of FastQ files.
The data comes from a sequencing (using HiSeq technology) of a
Measles virus. For testing purposes, you can download :download:`R1
<../sequana/resources/data/Hm2_GTGAAA_L005_R1_001.fastq.gz>` and
:download:`R2 <../sequana/resources/data/Hm2_GTGAAA_L005_R2_001.fastq.gz>`)
files that contain only 1500 reads. Copy them in a local directory.

First, run the sequana standalone application to initialise the pipeline
**quality_control**::

    sequana --pipeline quality_control --working-directory TEST --adapters PCRFree

This command downloads the required configuration file(s) in particular
the config file and the pipeline itself. This example should work out of
the box but you may want to look at the
configuration file **config.yaml**. For instance, you may want to change the
reference to the *phix* (by default we use *phix174.fa*, which is provided in
Sequana) or
change the adapter_removal section to your needs (cutadapt parameters, in
particular the forward and reverse complement list of adapters; None by
default).

By default, the output directory is called **analysis** and can be overwritten
with the ``--working-directory`` parameter. Then, run the pipeline and wait for
completion.::

    cd TEST
    snakemake -s quality_control.rules --stats stats.txt -p -j 4 --forceall

The -p option shows the commands, -j 4 means use 4 threads when possible.
Alternatively, there is also a **runme.sh** script.

You should now have a directory with a HTML report corresponding to the sample::

    open Hm2_GTGAAA_L005/report_qc_Hm2_GTGAAA_L005/summary.html


See :ref:`quick_start`


Taxonomy
-------------------------------

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

Then, you may use a Sequana pipeline (see :ref:`pipeline_taxon` and
 :mod:`sequana.kraken`) or this standalone application::

    sequana_taxonomy  --file1 Test_R1.cutadapt.fastq.gz
        --file2 Test_R2.cutadapt.fastq.gz --database  <database_path>

where <database_path> must be replaced with the proper path.


Open the local HTML file krona.html. An example is available
in  `Krona example <_static/krona.html>`_


Variant calling
-------------------

The following example will show how to initialise and run the variant calling
pipeline on a pair of FastQ files.
For testing purposes, you can download :download:`R1
<../sequana/resources/data/Hm2_GTGAAA_L005_R1_001.fastq.gz>` and
:download:`R2 <../sequana/resources/data/Hm2_GTGAAA_L005_R2_001.fastq.gz>`)
files that contain only 1500 reads. Copy them in a local directory.

Note that this does the variant calling + snpEff + coverage.
See more information in the :ref:`pipeline_variant_calling` section.



Initialise the pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Call **sequana** standalone as follows::

    sequana --pipeline variant_calling --input-directory . --working-directory TUTORIAL

Or use Sequanix. 

Go to the project directory
::

    cd TUTORIAL


Get the genbank reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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


New in v0.10
~~~~~~~~~~~~~~~~

Assuming the genbank and reference have the same name, you can simply
type::

    from sequana.snpeff import download_fasta_and_genbank
    download_fasta_and_genbank("K01711", "measles")

Get a snpEff config file and update it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Then you need to initialise a config file for snpEff tool::

    from sequana import snpeff
    v = snpeff.SnpEff("measles.gbk")

Update the snpeff config file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Edit the config file **config.yaml** and add the filename *measles.gbk* in the
snpEff section::

    # snpEff parameter
    snpeff:
        do: yes
        reference: "measles.gbk"

and bwa_ref section::

    # Bwa parameter for reference mapping
    bwa_mem_ref:
      reference: "measles.fa"

.. warning:: In the configuration file, in the mark_duplicates section,
    some output files are huge and requires temporary directory on cluster.


.. warning:: in the configuration file -- coverage section -- note that for short genomes, 
    you may need to decrease the window size.

.. warning:: the mark_duplicates may be changed in the close future to use
   another tool. 


Run the pipeline
~~~~~~~~~~~~~~~~~~~~


::

    snakemake -s variant_calling.rules --stats stats.txt -p -j 4 --forceall


De novo
-------------

The denovo_assembly pipeline can be initialised in the same way::

    sequana --pipeline denovo_assembly --input-directory . --working-directory denovo_test

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

    sequana --pipeline rnaseq --input-directory . --working-directory EXAMPLE
        --adapter-fwd GATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter-rev GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC

This command download the pipeline and its configuration file. The configuration
file is prefilled with adapter information and input data files found in the
input directory provided. You can change the configuration afterwards.

An alternative is to use :ref:`sequanix`.


Go to the project directory
::

    cd EXAMPLE


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

Edit the config file(s)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Edit the config file **config.yaml** and fill the genome section::

    genome:
      do: yes
      genome_directory: Saccer3
      name: Saccer3 #path to index name
      fasta_file: Saccer3/Saccer3.fa
      gff_file: Saccer3/Saccer3.gff
      rRNA_file:
      rRNA_feature: "rRNA"


.. warning:: Note that fastq_screen is off by default. Indeed, Sequana does not provide 
   a fastq_screen database so far. Therefore, if you want to run fastq_screen, please 
   see the manual (https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
   and add the config file in the tool section.


Finally, also edit the **multi_config.yaml** file and replace::

    custom_logo: "Institut_Pasteur.png"

with yours or as follows (empty, not an empty string like "" ) ::

    custom_logo:

.. note:: there are other places with hard-coded path but the corresponding
   sections are not used by default. If you decide to use them (e.g.
    fastq_screen), you will need to edit the configuration file accordingly.





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








