Notes
=======

The phiX174.fa file provided here is the one used within the Biomics platform and
differs slightly from the ENA sequence: http://www.ebi.ac.uk/ena/data/view/J02482
and http://www.ncbi.nlm.nih.gov/nuccore/9626372?report=fasta (NC_001422.1) with 
5 different bases.

Paired-end
===========

Those 2 files are paired-end and contains 2500 reads amongst which about 10 can
map to the phiX174 reference also provided in this directory:

- Hm2_GTGAAA_L005_R1_001.fastq.gz
- Hm2_GTGAAA_L005_R2_001.fastq.gz

Others
========

- adapters_others.fa
- adapters_netflex_pcr_free_1_rev.fa
- adapters_netflex_pcr_free_1_fwd.fa
- test_10000.fastq.gz


BAM
=========
- measles.fa.sorted.bam created with sequana_mapping of Hm*gz with reference measles.fa

RNAseq
========

- data/KO_ATCACG_R1_test.fastq.gz : a small sample of 500 reads for testing

SIRV 
=======

The SIRV.fa file is used for testing but can also be used for ISOSEQ analysis to
identify spikes that were possibly injected. This file is made of 68 spikes (Lot 001603 ) 
from lexogen.
