Docker **sequana_conda_ngs**
====================================

Overview
---------

This is a Docker image on top on sequana/sequana_ngs that installs
conda packages (NGS specific) for **Sequana**.


Description
----------------

Se the sequana_install.sh for detailled list of packages.

Here are some of them:


- snakemake
- bioservices
- pyvcf
- spades
- bwa
- freebates
- kraken
- bowtie2




For developers:
------------------

Build the image::

    git clone https://github.com/sequana/sequana
    cd sequana/docker/sequana_conda_ngs
    sh build.sh

Push on hub.docker.com::

    sh push.sh

Run the local image (not a pulled one in this example)::

    sudo docker run -it sequana_conda_ngs
