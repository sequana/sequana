:Overview: Denovo Assembly from FASTQ files
:Input: FASTQ file(s) from Illumina Sequencing instrument
:Output: FASTA, VCF and HTML files

Usage
~~~~~~~~~

Command line interface
########################

Example::

    sequana --pipeline denovo_assembly \
            --input-directory data/ \
            --input-readtag _[12].fastq \
            --working-dir analysis
    cd analysis
    snakemake -s denovo_assembly.rules --stats stats.txt

Requirements
~~~~~~~~~~~~~~~~

.. include:: ../sequana/pipelines/denovo_assembly/requirements.txt

.. image:: https://raw.githubusercontent.com/sequana/sequana/master/sequana/pipelines/denovo_assembly/denovo_dag.png

Details
~~~~~~~~~

Snakemake *de-novo* assembly pipeline dedicates to small genome like bacteria.
It bases on `SPAdes <http://cab.spbu.ru/software/spades/>`_.
The assembler corrects reads then assemble them using different size of kmer.
If the careful option is set, SPAdes corrects mismatches and short INDELs in
the contigs using BWA.

The sequencing depth can be normalised with `khmer <https://github.com/dib-lab/khmer>`_.
Digital normalisation converts the existing high coverage regions into a Gaussian
distributions centered around a lower sequencing depth. To put it another way,
genome regions covered at 200x will be covered at 20x after normalisation. Thus,
some reads from high coverage regions are discarded to reduce the quantity of data.
Although the coverage is drastically reduce, the assembly will be as good or better
than assembling the unnormalised data. Furthermore, SPAdes with normalised data
is notably speeder and cost less memory than without digital normalisation.
Above all, khmer does this in fixed, low memory and without any reference
sequence needed.

The pipeline assess the assembly with several tools and approach. The first one
is `Quast <http://quast.sourceforge.net/>`_, a tools for genome assemblies
evaluation and comparison. It provides a HTML report with useful metrics like
N50, number of mismatch and so on. Furthermore, it creates a viewer of contigs
called `Icarus <http://quast.sourceforge.net/icarus.html>`_.
The second approach is to characterise coverage with sequana coverage and
to detect mismatchs and short INDELs with
`Freebayes <https://github.com/ekg/freebayes>`_.
The last approach but not the least is `BUSCO <http://busco.ezlab.org/>`_, that
provides quantitative measures for the assessment of genome assembly based on
expectations of gene content from near-universal single-copy orthologs selected
from `OrthoDB <http://www.orthodb.org/>`_.

Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a documenteted configuration file :download:`../sequana/pipelines/denovo_assembly/config.yaml` to be used with the pipeline. Each rule used in the pipeline may have a section in the configuration file. Here are the rules and their developer and user documentation.

*De-novo* assembly
######################

Digital normalisation
^^^^^^^^^^^^^^^^^^^^^^^^
.. snakemakerule:: digital_normalisation

SPAdes
^^^^^^^^^^
.. snakemakerule:: spades

Format contigs
^^^^^^^^^^^^^^^^
.. snakemakerule:: format_contigs

Quality assessment metrics
##############################

QUAST
^^^^^^^^^
.. snakemakerule:: quast

BUSCO
^^^^^^
TODO

Genome annotation
###################

Prokka
^^^^^^^^
.. snakemakerule:: prokka

Re-mapping
############

BWA
^^^^
.. snakemakerule:: bwa_mem_dynamic

Sambamba markdup
^^^^^^^^^^^^^^^^^^
.. snakemakerule:: sambamba_markdup

Sambamba filter
^^^^^^^^^^^^^^^^^^
.. snakemakerule:: sambamba_filter

Mismatch detection
###################

Freebayes
^^^^^^^^^^
.. snakemakerule:: freebayes

Freebayes filter
^^^^^^^^^^^^^^^^^^
.. snakemakerule:: freebayes_vcf_filter

Coverage analysis
###################

Samtools depth
^^^^^^^^^^^^^^^^
.. snakemakerule:: samtools_depth

Sequana coverage
^^^^^^^^^^^^^^^^^^^
.. snakemakerule:: sequana_coverage
