Notes
=======

Jan 2018. ruamel.yaml installed is 0.12.13. In the requirements, right now we
have ruamel.yaml<=0.15 as recommended in the ruamel.yaml doc. l
ocally I (TC) use version 0.15 without trouble but in the singularity, somehow,
it is version 0.12.13 and this causes an error in the sequanix interface. 



ubuntu image 
=============

sudo singularity build ubuntu.img Singularity.ubuntu

This singularity install development tools (for qt and others, git, gcc, etc).
Then, we install a miniconda version and python 3.6. There is no specific
environement, just the base environement.

    conda install python=3.6


R image
=======

Based on ubuntu docker, we install cran package from the core, dev, recommended
channels as well as development librairies such as libcurl4-openssl-dev libssl-dev 
libxml2-dev  libcairo2-dev  libxt-dev.

We set a CRAN mirror and install the devtools package.


sartools image
==============

Bootstrap from local R image

install RNADiffAnalysis and latex libraries.



 conda install pigz pbzip2 dsrc cutadapt atropos
    conda install fastqc                                   # for fasqc pipeline
    conda install bcftools deeptools samtools bamtools bwa # general
    conda install multiqc                                  # general
    conda install freebayes sambamba snpEff                # variant calling
    conda install bowtie bowtie2 STAR                      # rnaseq
    conda install subread   # rnaseq (featureCounts)
    conda install igvtools   # rnaseq
    conda install fastq-screen   # rnaseq

    conda install picard 

sequana image
=============
inherit from: ubuntu
adds: sequana and its dependencies (conda)

sequana_tools image
===================

inherit from: sequana
adds: tools for the rnaseq pipeline and possibly other pipelines



