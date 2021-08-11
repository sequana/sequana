.. _wrappers:

Wrappers
##########

As of August 2021, **Sequana** team created the e `sequana wrappers repository <https://github.com/sequana/sequana-wrappers>`_, which is intended to replace the rules. The adavantage is that wrappers can be tested with a continuous integration.  


Wrappers are used within a Snakemake rule. When you call your Snakemake
pipeline, you will need to add::

    --wrapper-prefix git+file:https://github.com/sequana/sequana-wrappers/

We provide documentation for each wrapper. It can be included in this
documentation thanks to a sphinx extension. For example::

    .. sequana_wrapper:: fastqc

Here is a non exhaustive list of documented wrappers. 


.. sequana_wrapper:: fastqc
.. sequana_wrapper:: rulegraph
