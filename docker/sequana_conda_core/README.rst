Docker **sequana_conda_core**
====================================

Overview
---------

This is a Docker image on top on sequana/sequana_core that installs
conda and common scientific Python packages.

Description
----------------
This sets conda and the following packages:

- numpy
- matplotlib
- pandas
- scipy
- graphviz
- ipython

This installs other libraries such as pyqt5.

The environmental variables for matplotlib are also setup such as the back end
set to Agg



For developers:
------------------

Build the image::

    git clone https://github.com/sequana/sequana
    cd sequana/docker/sequana_conda_core
    sh build.sh

Push on hub.docker.com::

    sh push.sh

Run the local image (not a pulled one in this example)::

    sudo docker run -it sequana_conda_core
