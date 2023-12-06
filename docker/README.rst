Docker containers for **Sequana**
====================================

We do not provide Docker containers anymore. However, since sequana is posted on  bioconda, one can get some Dokcer containers. For example version 0.4.1 is available as explained here below. For full list please checkout
https://quay.io/repository/biocontainers/sequana

Example: sequana_coverage
--------------------------

To pull a Sequana container (here version 0.4.1), use this type of command::

    docker pull quay.io/biocontainers/sequana:0.4.1--py35_0

Checkout the `quay.io <https://quay.io/repository/biocontainers/sequana>`_
website. After pulling the image above, you can use it as follows::

    docker run -v $PWD:/home/default -it quay.io/biocontainers/sequana:0.4.1--py35_0

.. warning:: once in the docker shell, go to /home/default. Here, this directory
    is linked to your real directory where you type "docker run..." so what you
    modify here is directly reflected in your directory !


Assuming you have a BED file JB409847 in your directory,  otherwise uncomment
the commented line here below::

    cd /home/default
    export MPLBACKEND="agg"
    # wget https://tinyurl.com/y9j69t3k -O JB409847.bed
    sequana_coverage --input JB409847.bed
    exit

Back on your local directory, you should now see a ./report directory with the
results of the analysis.
