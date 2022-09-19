FAQS
======

Conda related
---------------

Create a conda environment on IP cluster::

    module load conda
    conda create --name py37 python=3.7 
    conda activate py37

add channel from where to download packages::

    conda config --add channels r bioconda
    conda install sequana


What are the dependencies
-----------------------------

There are two kind of dependencies. First, the Python libraries such as
matplotlib or Pandas. Second, the external tools such as BWA (alignment) or
Kraken (taxonomy). The first kind of tools can be installed using Anaconda and the
default conda channel. For instance::

    conda install pandas

The second kind of tools can also be installed using another conda channel
called **bioconda**. For instance::

    conda install bwa

Since version 0.12, most pipelines have been moved outside of sequana. Sequana itself only requires::

    conda install kraken2 cd-hit krona


Installation issues
-----------------------


As explained in the previous section, most of the dependencies can be installed
via Conda. If not, pip is recommended. Yet there are still a few dependencies
that needs manual installation. 

quast
~~~~~~~~~

http://quast.bioinf.spbau.ru/manual.html#sec1

::

    wget https://downloads.sourceforge.net/project/quast/quast-4.2.tar.gz
    tar -xzf quast-4.2.tar.gz
    cd quast-4.2

Alternatively, get the source code from their GitHub (takes a while)::

    git clone https://github.com/ablab/quast
    cd quast
    python setup.py install

graphviz
~~~~~~~~~~~~~~~~~~

graphviz provides an executable called **dot**. If you type **dot** in a shell
and get this error message::

    Warning: Could not load
    ...lib/graphviz/libgvplugin_gd.so.6" - file not found

This may be solved by re-installation graphviz using the main anaconda channel
(instead of bioconda)::

    conda install --override-channels -c anaconda graphviz=2.38.0 

:Update April 2017: replace anaconda with conda-forge


matplotlib
~~~~~~~~~~~~~~~~~

If you get errors related to the X connection, you may need to change the
backend of matplotlib. To do so, go in your home directory and in this directory

    cd /home/user/.config/matplotlib/

Check if the file **matplotlibrc** exits, if not, type::

    echo "backend: Agg" > matplotlibrc

or edit the file and make sure the line starting with "backend" uses the Agg
backend::

    backend: Agg

Save, exit the shell, start a new shell.


pysam / samtools / bzip2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We have experienced few issues with pysam and samtools. Here are some solutions.


::

    from pysam.libchtslib import *
    ...ImportError: libhts.so.1: cannot open shared object file: No such file or directory


This may be solved by removing conda installation and using pip instead::

     conda remove pysam
     pip install pysam

Another error know for pysam version 0.11.2.2 raises this error::

    ImportError: libbz2.so.1.0: cannot open shared object file: No such file or
    directory

Downgrading to version 0.11.2.1 and upgrading to working version solves the problem::

    conda install pysam=0.11.2.1

but one reason was also related to the order of the channel in the .condarc
file. You may get bzip2 from the default channel and not from
conda-forge (reference: https://github.com/bioconda/bioconda-recipes/issues/5188)
::

    conda install --override-channels -c conda-forge bzip2

pysam may not compile due to a missing dependency on lzma. Under fedora,
type::

    yum install liblzma liblzma-devel




qt and pyqt
~~~~~~~~~~~~~~~~~~

Qt Version
^^^^^^^^^^

With PyQt 5.12.3 and python3.7, we got lots of errors::

    SystemError: <built-in function connectSlotsByName> returned a result with an error set

This seems to be a PyQt bug according to several github projets based on pyqt.
It may be fixed a version above. Dowgrading e.g. to pyqt 5.9.2 does not solve
the problem.




Qt compatibility across platform
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

    from PyQt5.QtWebKitWidgets import QWebView
    ...ImportError: libQt5WebKitWidgets.so.5: cannot open shared object file: No such file or directory

This may be solved by re-installation qt using the main anaconda channel
(instead of bioconda)::

    conda install --override-channels -c anaconda qt

and possibly::

    pip install PyQtWebEngine

If we believe this issue: https://github.com/conda-forge/pyqt-feedstock/issues/19


libselinux
~~~~~~~~~~~~~~~~~

If you get this error (using **conda install sequana**)::

    ImportError: libselinux.so.1: cannot open shared object file: No such file or directory

it looks like you need to install libselinux on your environment as reported 
`here <https://github.com/sequana/sequana/issues/438>`_.


pytz installation
~~~~~~~~~~~~~~~~~~~~


If you get this error::

    ImportError: C extension: No module named 'pytz.tzinfo' not built. If you
    want to import pandas from the source directory, you may need to run 'python
    setup.py build_ext --inplace --force' to build the C extensions first.

try this::

    pip uninstall pytz
    pip install --pre pytz

reference: https://github.com/sequana/sequana/issues/499



Expected input format
----------------------------

Most of the pipelines and standalone expect FastQ files with the extension
**fastq.gz** meaning that files are gzipped.


Besides, the filename convention is as follows::

    PREFIX_R1_.fastq.gz

that is **_R1_** and **_R2_** indicates the paired or single-ended files and
the PREFIX is used to create directories or reports; it must be present.

.. versionadded:: 0.2
    more flexible tags are now possible in sequana pipelines and sequanix using
    e.g. _R[12] in the **input_readtag** in the configuration file of the
    pipelines.


Sequanix related
----------------------

For question related to Sequanix, we have a dedicated section in
:ref:`sequanix_faqs`.


QXcbConnection issue
----------------------
If you get this error::

    QXcbConnection: Could not connect to display localhost:10.0

this is an issue with your Qt backend. You need to change it to Agg.




Variant Calling pipeline
----------------------------

If snpeff fails with this type of errors::

    java.lang.RuntimeException: Error reading file 'null'
    java.lang.RuntimeException: Cannot find sequence for 'LN831026.gbk'

this may be because your genbank does not contain the sequences.

Another type of errors is that the sequence and genbank are not synchrone. We
would recommend to use the code here to download the Fasta and genbank:

http://sequana.readthedocs.io/en/main/tutorial.html#new-in-v0-10


Quality Control pipeline
---------------------------

Please see the tutorial, user guide or pipelines section and look for the quality control.

Then, if you do not find your solution, please open an issue on github: https://github.com/sequana/sequana/issues


Singularity
-----------------

If you use the singularity container and get this kind of error::

    singularity shell sequana-sequana-master.img
    ERROR  : Base home directory does not exist within the container: /pasteur
    ABORT  : Retval = 255

it means the container does not know about the Base home directory.

If you have sudo access, add the missing path as follows::

    sudo singularity shell --writable sequana-sequana-master.img
    mkdir /pasteur
    exit

If you do not have sudo permissions, copy the image on a computer where you have
such permission, use the same code as above and copy back the new image on the
computer where you had the issue. 

Finally, try to use the container again using this code::

    singularity shell sequana-sequana-master.img


I got a error "main thread is not in the main loop"
---------------------------------------------------

::

    Traceback (most recent call last):
      File
    ".../lib/python3.5/tkinter/__init__.py",
    line 627, in after_cancel
        data = self.tk.call('after', 'info', id)
    RuntimeError: main thread is not in main loop

This is related to the backend used by matplotlib. This can be ignored. We do
not have any solution for now, except finding an alternated backend for
matplotlib. This can be done using a special file called matplotlibrc with this
content::

    backend: tkagg

where you can replace tkagg with e.g. qt5agg

Installation issue on Mac
--------------------------

On a MacOSx conda environment (PYthon3.9), I could not build **datrie** with this kinf of error message::

    error: command ‘llvm-ar’ failed: No such file or directory
    ERROR: Failed building wheel for datrie
    Failed to build datrie
    Failed to build datrie
    ERROR: Could not build wheels for datrie, which is required to install pyproject.toml-based projects

The solution was to set the AR variable::

    export AR=/usr/bin/ar
