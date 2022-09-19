.. _rules:

Rules
##########

As of August 2017, **Sequana** has about 80 different rules. 
The list is available from the `source code <https://github.com/sequana/sequana/tree/main/sequana/rules>`_. We design our rules following some strict conventions as explained in the :ref:`developers` section.

Rules are documented and we developed a Sphinx extension to automatically add
their docstring in this documentation. For example, the documentation of the
rule **fastq_sampling** looks like:

.. snakemakerule:: fastq_sampling


In order to use a Sequana rule in your pipeline, add this code::

    from sequana_pipetools import snaketools as sm
    include: sm.modules["fastq_sampling"]

This takes care of the physical location of the rule. 
Of course, you will then need to look at the documentation and define the
required variables in your pipeline. For instance, in the example above, given
the documentation, you will need to define those two variables::

    __fastq_sampling_input_fastq
    __fastq_sampling_output_fastq

and have a configuration file with::

    fastq_sampling:
        N: 1000

Many rules are used inside the Sequana pipelines but not all. For instance, the
**codecs** rules (e.g. gz_to_bzip) are used in standalones.

Please see the :ref:`pipelines` section for other rule documentation (e.g. bwa,
fastqc, ...).
