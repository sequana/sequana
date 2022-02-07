.. _installation:

Installation
##########################################

Here below are the instructions to install Sequana. There are different ways (source, bioconda, singularity, conda environment, pip). Let us summarize the different methods for you.

If you want the latest version of Sequana, you should install it from source (see :ref:`github_method`). Otherwise, you can install a release of **Sequana** from the Pypi website (using **pip**). Note that for pipelines, which are now independent Python packages, we also use Pypi releases. However, third-party dependencies (not Python) should be installed manually. Most of them are provided through **Anaconda** channels.  See the :ref:`installation_conda` Section for details on how to set up Conda. 

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
    Sequana quite easily using pip, or bioconda, or virtual environment as shown
    here below.

.. topic:: Sequanix

    Sequanix has now its own repository here: https://github.com/sequana/sequanix and should 
    be installed independently.


Latest recommended installation method
======================================

Sequana is maintained under Python 3.7 and above  (Feb 2022).

.. warning:: For Sequanix, for a while we adviced to use Python 3.7.3 due to a PyQt library issue
    preventing Sequanix to work.

Lots of dependencies have been dropped in version 0.8.0 so that you could simply
use **pip** to install Sequana.


In any case we strongly recommend to use a virtual environment so that (i))
you can install all requirements without root permissions and (ii) you do
not interfer with your system.

We will use `conda <https://docs.conda.io/en/latest>`_ for that. Before starting
you should install and set the channels as explained in the  :ref:`installation_conda` section. Then, create an environment:
::

    conda create --name sequana_env python=3.7.3
    source activate sequana_env

pip installation
----------------

For the latest release of Sequana::

    pip install sequana --upgrade

This will install all Python dependencies such as Pandas, Numpy, etc. It will take about 5-10 minutes to install this version.

.. note:: If you want to use Sequanix, which rely on PyQt5, please install PyQt5 using conda::

        conda install -c anaconda qt pyqt>5

    Using pip may lead to compatibility issues with your underlying Qt library,
    which must be available to install PyQt. PyQt5, v5.9.2 is known to work.
    v5.15.0 fails (PyQt5.QtWebEngineWidgets)

    With Python 3.6.12, this is compatible::

        pip install "pyqt5<=5.10"


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
where more details about specific dependencies are provided. 


Other solutions (overview)
========================================

#. Singularity (tested with version 2.4.2; see below for installation) . Strictly speaking, there is no compilation. This method is for testing and production. It downloads an image / container that is ready-to-use (here the latest available release)::

      # NOTE THAT THIS IS AN OLD RELEASE 0.6.5
      singularity pull --name sequana.img shub://sequana/sequana

   and can be used as follows (for example)::

      singularity exec sequana.img sequanix --help

   See :ref:`Singularity <singularity_details>` section to install a specific release and more details.

#. Bioconda. **Sequana** is available on conda/bioconda as a pre-compiled package::

        # Note that its version may be behind the pypi releases
        conda install sequana

#. From source. If you prefer to install everything yourself, the source code is available on
   github (http://github.com/sequana/sequana) and releases are posted on Pypi::

        pip install sequana

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

Create an environement
-------------------------

Once **conda** is installed and the channels set, open a new shell.
Although this is not required strictly speaking, we would
recommend to create an environment dedicated to Sequana. This environment can
later be removed without affecting your system or conda installation. A
**conda** environment is nothing else than a directory and can be created as
follows::

    conda create --name sequana_env python=3.7.3

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
can install **sequana** from source::

    conda create --name sequana python=3.7.3
    source activate sequana
    git clone git@https://github.com:sequana/sequana.git
    cd sequana
    python setup.py install

    # to use sequanix interface:
    conda install -c anaconda qt pyqt>5

    # to perform testing and documentation:
    pip install -r requirements_dev.txt


This should install most of the required dependencies. However, you may need to
install more packages depending on the pipeline used (related to Qt for
instance).

.. _singularity_details:

Singularity
============
.. warning:: this is now up-to-date. Come back later or contribute to this
   section.

We provide Singularity images on https://singularity-hub.org/collections/114/ .
They contain Sequana standalones and some of the pipelines dependencies
as well as Sequanix. Note, however, that Sequanix relies on PyQt (graphical
environment) and would work for Linux users only for the time being. The main
reason being that under Mac and windows a virtualbox is used by Singularity
preventing a X connection. 

First, install singularity (http://singularity.lbl.gov/). You must use at least
version 3.5. We suggest users to look at the l=singularity installation page
itself to install the tool.
 
Once done, you can either build an image yourself or download a Sequana image. 
For instance, for the latest master version::

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

.. note:: you may need to install squashfs-tools (e.g. yum install squashfs-tools )


.. .. include:: ../docker/README.rst






