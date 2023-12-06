cd /home/sequana
source /home/sequana/miniconda3/bin/activate

conda install pigz -y
conda install snakemake -y
conda install bioservices
conda install cutadapt pysam pyvcf snpeff -y
conda install biokit spades khmer -y
conda install bwa bcftools samtools==1.3.0 bedtools picard freebayes fastqc -y
conda install kraken krona bowtie2 -y
conda install bowtie bowtie2 fastq-screen subread tophat -y
conda install pbr
