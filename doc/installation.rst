.. _installation:

Installation
##########################################

If you are a developer, you would want to install **Sequana** from source.
There are lots of dependencies that require compilation and may be time
consuming. We therefore recommend the **Anaconda** solution. Sequana is indeed
available on **bioconda**. Note, however, that
releases of Sequana are also available on Pypi so you could also use **pip**. 

If you just want to test **Sequana** or **Sequanix** or one of the Sequana
standalone, we also provide **Singularity** containers. This is a great 
solution for reproducibility as well. Containers are
available on https://singularity-hub.org/collections/114/. 

Latest recommended installation for developers
==============================================

::

    conda install -c anaconda qt
    conda install -c anaconda pyqt=5.6.0
    python setup.py install


Overview of installation methods
====================================

We support 3 types of installations:

#. Singularity (tested with version 2.4.2; see below for installation) . Strictly speaking, there is no compilation. This method is for testing and production. It downloads an image / container that is ready-to-use (here the latest available release)::

      singularity pull --name sequana.img shub://sequana/sequana

   and can be used as follows (for example)::

      singularity exec sequana.img sequanix --help

   See :ref:`Singularity <singularity_details>` section to install a specific release and more details.

#. Bioconda. **Sequana** is available on conda/bioconda as a pre-compiled package::

       conda install sequana

#. From source. If you prefer to install everything yourself, the source code is available on
   github (http://github.com/sequana/sequana) and releases are posted on Pypi::

        pip install sequana

These three methods are detailled hereafter.

.. _installation_conda:


From bioconda (Recommended )
===================================

If you have not installed **Sequana**, be aware that it relies on many dependencies
that needs to be compiled (i.e., it is time consumming and requires proper C compilator).
Besides, many pipelines rely on third-party software such as BWA or samtools that are not
Python libraries. Yet, using **conda**, this process is simplified.

Install conda executable
----------------------------

In practice, we do use `Anaconda <https://conda.readthedocs.io/>`_ . We recommend to
install **conda** executable via the manual installer (`download <https//continuum.io/downloads>`_). 
You may have the choice between Python 2 and 3. We recommend to choose a Python version 3.

Add bioconda channels
------------------------

When you want to install a new package, you have to use this type of syntax::

    conda install ipython

where **ipython** is the package you wish to install. Note that by default,
**conda** looks on the official Anaconda website (channel). However, there are
many channels available. We will use the **bioconda** channel. To use it, type
these commands (once for all)::

    conda config --add channels r
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

.. warning:: **it is important to add them in this order**, as mentionned on bioconda webpage
    (https://bioconda.github.io/).

If you have already set the channels, please check that the order is correct.
With the following command::

    conda config --get channels

You should see::

    --add channels 'r'   # lowest priority
    --add channels 'defaults'
    --add channels 'conda-forge'
    --add channels 'bioconda'   # highest priority

Create an environement
-------------------------

Once **conda** is installed and the channels set, open a new shell.
Although this is not required strictly speaking, we would
recommend to create an environment dedicated to Sequana. This environment can
later be removed without affecting your system or conda installation. A
**conda** environment is nothing else than a directory and can be created as
follows::

    conda create --name sequana_env python=3.5

Then, since you may have several environments, you must activate the **sequana**
environment itself (each time you open a new shell)::

    source activate sequana_env


Installation
-------------------

Sequana is on `bioconda <https://bioconda.github.io/>`_. You can follow these `instructions <http://bioconda.github.io/recipes/sequana/README.html>`_ or type::

    conda install sequana


This does not provide all dependencies needed by the different pipelines. So,
you may need to install extra packages as listed in this requirement file that
can be used with conda::

    conda install --file https://raw.githubusercontent.com/sequana/sequana/master/requirements_pipelines.txt

Additional tools may need to be installed. Thos large packages are stored in
another requirements to keep the main distribution lighter::

    conda install --file https://raw.githubusercontent.com/sequana/sequana/master/requirements_pipelines_extra.txt


From Pypi website (released source code)
==========================================
If you do not want to use **conda**, we provide releases on the Python Package Index website (pip tool)::

    pip install sequana
    pip install PyQt5


.. warning:: we do not support this methods but it should work. The main
    issues being that you will need to install the dependencies yourself. See
    hereafter for some of the tool used by the pipelines


From GitHub Source code
===========================

Finally, if you are a developer and wish to use the latest code, you 
can install **sequana** from source::

    git clone git@github.com:sequana/sequana.git
    cd sequana
    python setup.py install

This should install most of the required dependencies. However, you may need to
install more packages depending on the pipeline used. See hereafter.

.. _singularity_details:

Singularity
================

We provide Singularity images on https://singularity-hub.org/collections/114/ .
They contain Sequana standalones and some of the pipelines dependencies
as well as Sequanix. Note, however, that Sequanix relies on PyQt (graphical
environment) and would work for Linux users only for the time being. The main
reason being that under Mac and windows a virtualbox is used by Singularity
preventing a X connection. This should be solved in the near future.

First, install singularity (http://singularity.lbl.gov/). You must use at least
version 2.4. We tested this recipe with version 2.4.2 (Dec 2017)::

    VERSION=2.4.2
    wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-$VERSION.tar.gz
    tar xvf singularity-$VERSION.tar.gz
    cd singularity-$VERSION
    ./configure --prefix=/usr/local
    make
    sudo make install

Second, download a Sequana image. For instance, for the latest master version::

    singularity pull --name sequana.img shub://sequana/sequana:latest

or for the release 0.6.3::

    singularity pull --name sequana_0_6_3.img shub://sequana/sequana:0_6_3

The term latest in Singularity Hub will pull, across all of your branches and
tags, the most recent image, so if you come back in a year and get the latest (or ommit tha tag), you may not get the same container ! So, it is best using a specific tag. 

Do not interrupt the download (1.5Go). Once downloaded,
you can use, for instance, the sequana_coverage executable::

    singularity exec sequana.img sequana_coverage --help

or sequanix::

    singularity exec sequana.img sequanix

Would you miss a dependency, just enter into the singularity container and install the missing dependencies. You will need writable permission::

    sudo singularity shell -w sequana.img

Then, inside the container, install or fix the problem and type exit to save the
container.

.. note:: method tested with success on Fedora 23, ubuntu and Centos 6.
.. seealso:: Notes for developers about :ref:`dev_singularity` especially to get
   specific versions.


.. note:: you may need to install squashfs-tools (e.g. yum install squashfs-tools )

Notes about dependencies
===========================

When installing **Sequana** with conda and from the source, it should install
all the Python dependencies and you should be ready to go to use the Sequana
Python library.

However, note that most of the pipelines rely on extra dependencies that are not
necesseraly Python-based. For instance **bwa** is in C, others may be in R or
perl. 

The list of requirements is available in the source code::

    https://raw.githubusercontent.com/sequana/sequana/master/requirements_pipelines.txt

and conda may be used to install those dependencies automatically::

    conda install --file https://raw.githubusercontent.com/sequana/sequana/master/requirements_pipelines.txt

Otherwise you need to proceed to the installation of those dependencies by
yourself.

.. note:: atropos is an alternative to cutadapt with additional options but same
   type of functionalties and arguments. We use version 1.0.23 and above though.


.. note:: the denovo_assembly pipelines uses Quast tool, which we ported to
    python 3.5 and was pulled on Quast official github page. This is not
    yet in bioconda but one can get it from the quast github (sept 2016). This is
    to be installed manually by users (due to licensing restrictions)

.. note:: For GATK (variant caller), please go to
   https://software.broadinstitute.org/gatk/download/auth?package=GATK and
   download the file GenomeAnalysisTK-3.7.tar.bz2 ; then type::

       gatk-register GenomeAnalysisTK-3.7.tar.bz2



.. include:: ../docker/README.rst






