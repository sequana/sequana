Effect of the trimming on SNPs detection
===========================================

:Description: Effect of trimming (or not trimming) on the the SNPs detection.

In this case example, we will take a paired-end data set, and apply the quality
pipeline using trimming quality (removing bases with quality below 30). Then, we
will run the variant calling pipeline to perform the mapping on a reference
and detect SNPs.

We will repeat this analysis without trimming low quality reads at all.

We will finally compare the two sets of SNPs showing that the trimming quality
is not important in this example. Meaning that the mapping tool used (freebayes)
is able to cope with low quality reads.


The data
------------

We will use a paired-end data set (MiSeq 250bp). It contains 250,000 reads
(X2). The organism sequenced is *Bordetella*. As a reference, we use
the ENA accession *CP010347.1*. The data will be posted later but the original
data were generated at Pole Biomics (Institut Pasteur) and named **Tohama-R0_S4_L001_R1_001**
from which we used only the first 250,000 reads.

Here is a boxplot of the base quality across the reads showing that the quality
is quite high and falls below 30 after 200 bases.

.. image:: _static/case_examples/fastqc_raw.png
   :width: 80%



Quality pipeline
------------------

Assuming DATA (fastq.gz files) are in <DIR1> directory, type this command to create the
**quality** pipeline and config file automatically::

    sequana_quality_control --input-directory <DIR1> --working-directory trimming

Then go to the project and execute the pipeline::

    cd trimming
    snakemake -s quality_control.rules -p -j 4 --forceall

.. note:: you can also use Sequanix to help in the configuration design.


The final cleaned reads are in trimming/Tohama-R0_S4_L001/report/outputs (refered
to <DIR2> hereafter)
and named after the project: (`trimming_R1_.cutadapt.fastq.gz` and
`trimming_R2_.cutadapt.fastq.gz`). These two files should be used later as
the input of the variant_calling pipeline, as shown hereafter.


There is no adapters in the data so in the config file, the adapter sections are
empty (no forward or reverse adapters). Note, however, that bad quality bases
below 30 (default) are removed. In order to set the quality to another values,
use **sequana** with the `--quality` option

.. seealso:: See the :ref:`tutorial` and :ref:`quick_start` sections for more details.

Quality pipeline (No trimming)
-----------------------------------

Repeat the previous two steps. In the first step, change the adapter section
(cutadapt) to set the quality to zero (this prevents the trimming of bad quality
bases)::

    cutadapt:
        quality: 0,0

Change the project name e.g. *no_trimming* as a tag to the project in the first
step and *variant_no_trimming*.


Variant analysis
------------------

The output of the **quality** pipeline will be the input of the **variant calling**
pipeline::

    sequana_variant_calling --input-directory <DIR2> --working-directory variant_trimming


Here you need to make sure that the **config.yaml** configuration file
has the correct reference. See the :ref:`tutorial` section (variant section).

::

    reference = "CP010347"
    from bioservices import EUtils
    eu = EUtils()
    data = eu.EFetch(db="nuccore",id=reference, rettype="gbwithparts",
        retmode="text")
    with open("data.gbk", "w") as fout:
        fout.write(data.decode())
    from bioservices import ENA
    ena = ENA()
    data = ena.get_data(reference', 'fasta')
    with open("data.fa", "w") as fout:
        fout.write(data.decode())
    from sequana import snpeff
    v = snpeff.SnpEff("data.gbk")

Edit the config.yaml to change those sections::

    # snpEff parameter
    snpeff:
        do: yes
        reference: "data.gbk"

    # Bwa parameter for reference mapping
    bwa_ref:
        reference: "data.fa"


Run the analysis::

    cd variant_trimming
    snakemake -s variant_calling.rules --stats report/stats.txt -p -j 4 --forceall

Once done, you should have VCF files in **variant/report/** named **cutadapt.ann.vcf**





SNPs results comparison
----------------------------

You should now have two VCF files. Here below we plot the read depth versus
strand balance. The color will indicates the overall freebayes score (normalised
by the largest score). A *good* candidate should have large score and balance
value around 0.5 The y-axis shows the read depth.

::

    from pylab import *

    from sequana import freebayes_vcf_filter

    vcf1 = freebayes_vcf_filter.VCF("variant/report/cutadapt.ann.vcf")
    vcf2 = freebayes_vcf_filter.VCF("variant_no_trimming/report/variant_no_trimming.ann.vcf")

    df1 = vcf1.filter_vcf()._vcf_to_df()
    df2 = vcf2.filter_vcf()._vcf_to_df()

    subplot(1,2,1)
    scatter(list(df1.strand_balance.values), list(df1.depth.values),
    c=list(df1.freebayes_score.values/1240))
    xlabel("strand balance")
    ylabel("Depth")
    grid()
    ylim([0,90])

    subplot(1,2,2)
    scatter(list(df2.strand_balance.values), list(df2.depth.values),
        c=list(df2.freebayes_score.values/1240))
    grid()
    ylim([0,90])
    title("Trimming quality (left) vs no trimming (right)
    \n")


.. image:: _static/case_examples/trim_vs_notrim_on_snp.png
   :width: 80%


In this figure the LHS (trimming) 294 SNPs were found while in the RHS (no
trimming)  309  were found. The additional SNPs all have low coverage below 20.
A third of them have low balance strand.

There is one SNP found in the trim case not found in no_trim. However, it is
marginal with strand balance of 0.12, depth of 11, frequence of 0.73 and one of the lowest score



Conclusions
--------------

The detection of SNPs does not suffer from not trimming low quality bases. 
Actually, some new SNPs are found. However, the are usually not significant
(low depth, low score or unbalanced). Interestingly, the distribution of the
SNPs in the depth vs strand balance plane seems to be more centered on strand
balance=0.5. We also notice that the depth is 10% better which means that the
low quality bases have contributed to the improvments of the depth and freebayes
sorfes. It could be interesting to extend the analysis to more data, lower
quality, or higher quality threshold. Note also that because there are more low
quality bases, there much more false alarms; However setting a freebayes score
threshold around 5 removes most of them.






