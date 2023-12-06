.. _installation:

Installation
##########################################

Here below are the instructions to install Sequana. There are different ways (source, bioconda, singularity, conda environment, pip). Let us summarize the different methods for you.

If you want the latest version of Sequana, you should install it from source (see :ref:`github_method`). Otherwise, you can install a release of **Sequana** from the Pypi website (using **pip**). Note that for pipelines, which are now independent Python packages, we also use Pypi releases. However, third-party dependencies (not Python) should be installed manually. Most of them are provided through **Anaconda** channels.  See the :ref:`installation_conda` Section for details on how to set up Conda.

For instance, if you want to use the sequana_fastqc pipelinem you must install **fastqc** yourself, which is not a
Python package.

If you just want to test **Sequana** or **Sequanix** (see note here below) or one of the Sequana
standalone, we also provide **Singularity** containers as explained in the
:ref:`singularity_details` section.


.. topic:: Design choice

    Since version 0.8.0, we decided to move the pipelines outside of the main
    sequana library. This choice was made to face the increase of pipelines
    available in the Sequana project. Indeed, each pipeline comes with its own
    dependencies, which are not neccesseraly Python. The full installation of
    Sequana started to be cumbersome even for experienced users. We dealt with this
    issue using bioconda. Yet, even with such solutions it started to be
    difficult to manage easy installation. So, as usual, divide and conquer:
    each pipeline has now its own life cycle outside of Sequana. For example,
    the variant calling pipeline is hosted on
    https://github.com/sequana/variant_calling. This way, you can install
    Sequana quite easily using pip.

.. topic:: Sequanix

    Sequanix has now its own repository here: https://github.com/sequana/sequanix and should
    be installed independently.


Latest recommended installation method
======================================

Sequana is maintained under Python 3.8 and above  (Dev 2023).

We strongly recommend to use a virtual environment so that (i)
you can install all requirements without root permissions and (ii) you do
not interfer with your system.

We will use `conda <https://docs.conda.io/en/latest>`_ for that. Before starting
you should install and set the channels as explained in the  :ref:`installation_conda` section. Then, create an environment. Here we set Python to 3.8 but could be 3.9 or 3.10::

    conda create --name sequana_env python=3.8
    source activate sequana_env

pip installation
----------------

For the latest release of Sequana::

    pip install sequana --upgrade

This will install all Python dependencies such as Pandas, Numpy, etc.


pipelines
----------
Sequana pipelines are now easily installable using **pip**::

    pip install sequana_rnaseq
    pip install sequana_fastqc
    pip install sequana_demultiplex
    pip install sequana_pacbio_qc
    # etc

The dependencies of this pipeline must be dealt with by the developer or users.
Each pipeline has its own repository on github (https://github.com/sequana/sequana_PIPELINENAME)
where more details about specific dependencies are provided. Note, however, that all pipelines
can be launched with apptainers so theoritecally, only Python is required (and apptainer).

A set of predefined pipelines can be installed using::

    pip install sequana[pipelines]


Other solutions (overview)
========================================


#. Bioconda. **Sequana** is available on conda/bioconda as a pre-compiled package::

        # Note that its version may be behind the pypi releases
        conda install sequana

#. From source. If you prefer to install everything yourself, the source code is available on
   github (http://github.com/sequana/sequana)::

        git clone https://github.com/sequana/sequana
        cd sequana
        pip install sequana .

# Singularity/Apptainer

We provide container with sequana shipped inside (no pipelines) within the damona project
(https:://github.com/cokelaer/damo/Apptainer

We provide container with sequana shipped inside (no pipelines) within the damona project
(https:://github.com/cokelaer/damon/Apptainer

We provide container with sequana shipped inside (no pipelines) within the damona project
(https:://github.com/cokelaer/damon/Apptainer

We provide container with sequana shipped inside (no pipelines) within the damona project
(https:://github.com/cokelaer/damona))



These three methods are detailled hereafter.

.. _installation_conda:

From bioconda
==============

If you have not installed **Sequana**, be aware that many dependencies need to
be compiled (i.e., time consumming and requires proper C compilator).
Besides, many pipelines rely on third-party software such as BWA or samtools that are not
Python libraries. We therefore recommend to use **conda** that provides pre-compiled
software for you.

Install conda executable
----------------------------

.. warning:: this is currently broken on bioconda. We advise you to install sequana
   with Python (pip) for the latest versions.


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

As of May 2020, the recommended order is now::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

Create an environment
-------------------------

Once **conda** is installed and the channels set, open a new shell.
Although this is not required strictly speaking, we would
recommend to create an environment dedicated to Sequana. This environment can
later be removed without affecting your system or conda installation. A
**conda** environment is nothing else than a directory and can be created as
follows::

    conda create --name sequana_env 'python=3.8'

Then, since you may have several environments, you must activate the **sequana**
environment itself (each time you open a new shell)::

    source activate sequana_env


Installation
-------------------

Sequana is on `bioconda <https://bioconda.github.io/>`_. You can follow these `instructions <http://bioconda.github.io/recipes/sequana/README.html>`_ or type::

    conda install sequana

.. _github_method:

From GitHub Source code
===========================

Finally, if you are a developer and wish to use the latest code, you
can install **sequana** in develop mode as follows::

    conda create --name sequana 'python=3.8'
    source activate sequana
    git clone git@https://github.com:sequana/sequana.git
    cd sequana
    pip install -e .

    # to perform testing and documentation:
    pip install -e .[doc,testing]


This should install most of the required dependencies. However, you may need to
install more packages depending on the pipeline used (related to Qt for
instance).

.. _singularity_details:

Singularity/Apptainer
======================

We maintain a version of sequana within the Damona project.

You can download e.g version 0.16.2 and use it as follows::

    wget https://zenodo.org/record/10258126/files/sequana_0.16.2.img
    singularity sequana_0.16.2.img sequana --help
