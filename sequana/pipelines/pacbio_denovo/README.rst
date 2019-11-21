:Overview: denovo with canu dedicated to pacbio raw data set followed by quality assessment
:Input: ***.fasta** files  to be found in a directory. This is not the BAM file.
    To convert the BAM to fasta, you may use **bioconvert** at https://github.com/bioconvert

:Output: assembly
:Status: in progress

Usage
~~~~~~~

For now, please use :ref:`sequanix_tutorial` interface. 

Requirements
~~~~~~~~~~~~~~~~~~

- canu
- busco

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/pacbio_denovo/dag.png


Details
~~~~~~~~~

- Canu step:
- Busco step:

you would need to download the datasets required by BUSCO. The data is stored in
your HOME in the .config/sequana/busco directory. To download and uncompress the
data automatically, you can use::

    from sequana import busco
    busco.BuscoDownload().download()

Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documented configuration file :download:`../sequana/pipelines/pacbio_denovo/config.yaml` to be used with the pipeline. Each rule used in the pipeline may have a section in the
configuration file. In the *pacbio_denovo* pipeline, we use the *canu*, *busco* rules described here below.


Canu
^^^^^^^^^^^
.. snakemakerule:: canu

Busco
^^^^^^^^^
.. snakemakerule:: busco


