Validated on real data
----------------------

**Feb 2019**

TruSeqCD double indexing for Hiseq
TruSeqCD double indexing for Iseq100 on 3 indices only

Those 2 are different. See illumina specs




--------------------------------------------------------------

NEXTERA OK
PCRFree renamed into NEXTFlex48_DNA
TruSeqCD_DNA --> NEW

adapters_NEBNext2_fwd.fa       adapters_Nextera_fwd.fa   adapters_TruSeq_fwd.fa
TruSeqRNA_fwd.fa
adapters_NEBNext_fwd.fa        adapters_Small_fwd.fa     Nextera_fwd_bbmap.fa
TruSeqUD_fwd.fa
adapters_NEBNextOligos_fwd.fa  adapters_SMARTer_fwd.fa   NEXTFlex48_DNA_fwd.fa
adapters_Nextera_2_fwd.fa      NEXTFlex96_DNA_fwd.fa


TODO:
-----
NEXTERA_no_transposase.fa.gz
This should be used instead of NEXTERA when trimming Nextera Long-Mate Pair
libraries.



The adapters are stored as FASTA files compatible with cutadapt.
The files are exemplars of what can be used and are not universal. 

They are provided in 3 different versions: standard, reverse and reverse
complement; these 3 files contain the same information.



Summary of what's included
--------------------------

- Nextera i5 Index 2 including all S/N indexes in adapters_Nextera
- Nextera i7 Index 1 including only N indexes in adapters_Nextera
- TruSeq Single Indexes in adapters_TruSeq
- TruSeq Small RNA in adapters_Small






Nextera DNA Indexes
--------------------
contains only Index 2 (i5) for now.

Transposase
~~~~~~~~~~~~~~~

The transposase adapters are used for Nextera tagmentation and are removed. The
two transposases are:

Read 1::

    5′ TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG

Read 2::

    5′ GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG

Indexes
~~~~~~~~~
The Nextera file contains i5 Bases for Sample Sheet MiSeq, HiSeq 2000/2500
as quoted in page 13 of http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences_1000000002694-01.pdf

known as S/E/N

N50X Nextera DNA, S50X Nextera XT and E50X Nextera Enrichment and Nextera Rapid Capture.

The series S/E/N have actually the same index.


PCRFree
----------
Validé 18 Janvier:

Bar code are 6-bases long
http://www.biooscientific.com/Portals/0/IEM/Bioo-Scientific-PCR-Free-Barcode-Indices-v1-1-15.pdf

This is actually NEXTFlex-PCRFree DNA-seq kit for illumina

D'apres la doc, on a DNA adapter1::

>Universal_Adapter|name:universal
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
> et barcoded
GATCGGAAGAGCACACGTCTGAACTCCAGTCACXXXXXXATCTCGTATGCCGTCTTCTGCTTG

Tested on project 1142

SB1785,,,,NF_15,ATGTCA,,
SB1793,,,,NF_16,CCGTCC,,
SB3207,,,,NF_17,GTAGAG,,
SB3209,,,,NF_26,ATGAGC,,
SB3213,,,,NF_27,ATTCCT,,
SB3218,,,,NF_37,CGGAAT,,
SB3219,,,,NF_41,GCGCTA,,
SB1783,,,,NF_42,TAATCG,,
SB3177,,,,NF_43,TACAGC,,


Ceci est bon pour le kit de 48 adapters 6-basse long.


Ajout du 96 barcodes cas également.

16S
-------
http://support.illumina.com/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf

The dual indexing strategy uses two 8 base indices, Index1 (i7) adjacent to the P7 sequence, and Index2 (i5) adjacent to the P5 sequence. Dual indexing is enabled by adding a unique Index 1 (i7) and Index 2 (i5) to each  sample.

The  96 sample Nextera XT Index Kit (FC‐131- 1002) use12 different Index1 (i7) adapters (N701-N712) and 8 different Index2 (i5) adapters (S501-S508).

The 24 sample Nextera XT Index Kit (FC‐131-1001) uses 6 different Index 1 (i7) adapters (N701-N706) and 4 different Index 2 (i5) adapters (S501-S504). 

In the Index adapter name, the N or S refers
to Nextera XT sample preparation, 7 or 5 refers to Index1 (i7) or Index2 (i5), respectively.

The 01-12 refers to the Index number.

============ ========= ============ =======================
Index 1 (i7) Sequence, Index 2 (i5) Sequence
============ ========= ============ =======================
N701         TAAGGCGA  S501         TAGATCGC
N702         CGTACTAG  S502         CTCTCTAT
N703         AGGCAGAA  S503         TATCCTCT
N704         TCCTGAGC  S504         AGAGTAGA
N705         GGACTCCT  S505         GTAAGGAG
N706         TAGGCATG  S506         ACTGCATA
N707         CTCTCTAC  S507         AAGGAGTA
N708         CAGAGAGG  S508         CTAAGCCT
N709         GCTACGCT
N710         CGAGGCTG
N711         AAGAGGCA
N712         GTAGAGGA
============ ========= ============ =======================


TruSeq CD Indexes 
------------------------------------
:name: TruSeqCD_DNA

Combinatorial dual (CD) index adapters (formerly TruSeq HT). Extracted from IEM
using sequana.iem.IEM class.



TruSeq Single Indexes
----------------------
**current name is TruSeq**. Not to be confused with double indexing.


from TruSeq **LT** Kits and TruSeq v1/v2 Kits for DNA and RNA.

The following sequences are used for adapter trimming.

Read 1::

    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

Read 2::

    AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

TruSeq UniversalAdapter::

    5′ AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT

Index adapter sequences are six bases. The index numbering is not sequential, 
so indexes 17, 24, and 26 are skipped. Additionally, the bases
preceding each index adapter sequence are the same, but the two bases following
the index adapter sequence can vary !

Note that all the indexed adapters should be 5’-Phosphorylated. For unknown
reasons adapters 13-27 have an extra 2 bases (these are not used for the
indexing). Illumina also reserve certain numbers e.g. 17, 24 and 26. 

We included the UA, Read1 and Read2 and 23 index adapters in TruSeq_DNA_SI 


This commands works actually pretty well removing 11% of adapters symetrically::

    time atropos trim -T 4 -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -B
    AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -pe1 Hm2_GTGAAA_L005_R1_001.fastq.gz -pe2
    Hm2_GTGAAA_L005_R2_001.fastq.gz -o test1.fastq.gz -p test2.fastq.gz --progress
    bar --max-reads 200000 -m 20 -q 30 -O 6 --trim-n 

using -a/-A would work as well

TruSeq LT and HQ
-------------------

Read1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
Read2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT


ScriptSeq and TruSeq DNA Methylation:
-------------------------------------

    Read 1: AGATCGGAAGAGCACACGTCTGAAC
    Read 2: AGATCGGAAGAGCGTCGTGTAGGGA




TruSeq Small RNA
-----------------
from Illumina TruSeq small RNA Kits.

Based on document 10000002694 v09

!! in the document above, six-bases long index adapters are reverse
complemented. !!

Reference are eg RPI5 for index 5.

::

    TGGAATTCTCGGGTGCCAAGG

Others
--------
AmpliSeq, Nextera, Nextera DNA Flex, Nextera DNA, Nextera XT, Nextera
Enrichment, Nextera Rapid Capture Enrichment, TruSight Enrichment, TruSight
Rapid Capture Enrichment, TruSight HLA::

    CTGTCTCTTATACACATCT



SMARTer
-------
from SMARTer Stranded RNA-Seq Kit


NEBNext
------------

Contains the single index 
Universal found in manualE7335_index_primers_set1 of illumina.


Why do we have 

  4 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG

in the file

and 

5´-CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTG-
GAGTTCAGACGTGTGCTCTTCCGATC-s-T-3´

in manuals




Nextera Mate Pair:
--------------------
    CTGTCTCTTATACACATCT+AGATGTGTATAAGAGACAG











NExtFlex


TODO TOCHECK
==============

**PCRFree**, TruSeqCD + TruSeqUD. Faut il include le premier A dans les adapteurs. 
En effet d'après IEM, c'est inclus mais mais sur qu'il le faut. Par exemple dans
TruSeq, il n'y ait pas

A+GATCGGAAGAGCACACGTCTGAACTCCAGTCA



Validations
===========

**NEXTERA** validated on PG-P16-1 data using thr sample sheet HV7FYBCX2.csv
Results show that R1 and R2 have index N701 and N505 as well as transposase 1
for second iand first read. universal is not found. 

**TruSeqCD_DNA** example is BHNKWWBCX2. double indexing
