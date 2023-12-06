**sequana_taxonomy** docker
===============================

IN PROGRESS

This docker creates a single entry point for the sequana_coverage **standalone**. It is a simple layer on top of the **sequana** image, which makes available the **sequana_taxonomy** standalone only. For the full library, please see the main docker image.

#. Make sure you are logged in hub.docker.com::

    docker login

#. Get the docker image and rename it::

    docker pull sequana/sequana_taxonomy

#. Assuming you have a bed file named test.bed, type::

    docker run --user 1000 -v $(pwd):/home/sequana/data -w /home/sequana/data -it sequana/sequana_coverage --file1 test.fastq.gz

and browse the report in ./report

You may want to create an alias as follows::

    alias sequana_cov = "docker run --user 1000 -v $(pwd):/home/sequana/data -w /home/sequana/data -it sequana/sequana_coverage"

For developers
-----------------

Build the image::

    git clone https://github.com/sequana/sequana
    cd sequana/docker/sequana/sequana_taxonomy
    sh build.sh

Push it on hub.docker.com::

    sh push.sh
