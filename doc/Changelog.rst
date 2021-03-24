Changelog
=============

.. contents::
    :depth: 2




0.10.0
-------
MAJOR release merging the rnadiff branch into the dev branch so as to
finalise the way differential analysis are handle in Sequana with respect to
RNA-seq analysis

Module concerned: rnadiff, compare, featurecounts and enrichment

* CHANGES:

  * removed expdesign module. Improved iem module to cope with old formats. From
    now on, we wil use IEM formats only for the illumina sample sheets. No more
    home-made formats will be supported.

* NEWS:

   * module kegg and command in sequana (sequana taxonomy) to search for a given
     name in all KEGG orgnanism so as to easily retrieve the taxon and organism name


0.9.8
-----

* BUGS:

  * fastqc.FastQC class can now hndle empty samples (no reads... yes it may
    happen)
  * fastq random read selection  corrected to handle paired data

* IMPROVEMENTS:

  * CLI sequana summary can parse GFF and print summary information
  * add custom section for Makefiles in pipelines
  * add an unlock.sh for all pipelines
  * featurecounts: allow input to use the --extraAttributes option; more tests

* NEWS:
  * gff 1. add save_gff_filtered method to select interesting features, 2. add
    gff2gtf conversion
  * fasta module has GC content method
  * new rule macs3_dynamic


0.9.7
-----

* rnadiff HTML report can now show the locus_tag in volcano plot
* rules star_index uses floor instead of round for small genomes to agree with
  star index recommandation
* add MultiFeatureCounts() class to read several feature counts 

0.9.6
-----

* BUGS:

    * Fix star_index rule (regression bug) reported in https://github.com/sequana/sequana_rnaseq/issues/7
    

* IMPROVMENTS/CHANGES:

    * add plot() method in Repeats class
    * featurecounts: sort the output dataframe columns 
    * sequana CLI: summary/enrichment sub command allows to switch 
      between log/linear scale + option to compute the GO levels and 
      empty plots are now shown 





0.9.5
-----

* rules fastqc: allows handling of bam files. remove sanity checks on empty files. 
* rules falco: new rules added (similar to fastqc)
* module snaketools: PipelineManagerGeneric set paired field to "undefined"
  since this is suppose to deal with bam, or other non-fastq files
* module fastqc now plots the pacbio phred score on the full y-axis range
* for developers: In SequanaManager teardown, option to skip check_fastq_files


0.9.4
-----

* faster gff module to be used by sequana_rnaseq>=0.9.19

0.9.3
-----

general
~~~~~~~~~

* finalise the rnadiff HTML report with plotly visualisation
* finalise the enrichmentf HTML report
* template for rnadiff (rna differential analysis)
* method in FastQ to extract some specific reads
* add new application in main sequana tool:
  - biomart: a tool to produce a CSV for a given dataset!mart/set of attributes
  - summary: a tool to produce HTML for bam, fasta/fastq/bam, rnadiff results, enrichment for a rnadiff result
  - fastq: can now merge several files
  - gtf_fixer that replaces the script sequana_gtf_fixer
* add quick_fix for sample sheet errors in iem.IEM class (used in
  sequana_demultiplex pipeline)
* Fix featurecounts package to include tolerance parameter
* Fix bug for category in Dendogram

for developers
~~~~~~~~~~~~~~~

* rules: star_index now sets automatically the IndexNbases parameter
* removed sequana_summary from standalone now in the general sequana CLI
* improved find_motif module with plots, pep8, better API
* js: reverse papaparse

0.9.2
-----


* add salmon module for the rnaseq pipeline
* fix bug in multiqc plugin for bamtools used in sequana_mapper
* some rules updates: (1) add_read_group has now the option
  "VALIDATION_STRINGENCY=SILENT" that seem to be required in newer version of
  picard. surprising that this error did not appear earlier since the issue seems
  old (https://github.com/Cibiv/NextGenMap/issues/3) (2)  mark_duplicates and
  features_counts rules needed o change some parameters, now keywords in snakemake
* add new viz module called dendogram
* add new report module for rnadiff. 
* add new report module for KeggPathway and Panther
* Add solution to convert ensembl ID to gene name fr the KeggPathwayEnrichment
* Sequana main script has now a 'summary' command. Integrated the RNADiff and
  BAMQC modules for now.

0.9.1
-----
* Fix html syntax in the templates
* Fix demultiplex.barplot_samples layout
* Created a sequana main script with currentl one command, fastq, which replaces
  sequana_fastq_head, sequana_fastq_count, etc

0.9.0
-----


* Fix requirements (missing stattmodels)
* Fix bug kraken sequential when no reads are classified
* Fix FileFactory to make it much faster (speed up pipeline initialisation). This should also fix a recurrent issue with read tag (see https://github.com/sequana/sequana/issues/526
* NEWS:

    * implemented rnadiff module
    * implemented go term enrichment
    * implemented feature counts 
    * implemented plot.viz.venn
    * implemented compare module (for rnadiff comparison)


0.8.6
-----

* CHANGES:

    * module kraken: better label in plot2 method to avoid long label to be cut
    * module compare: fix label in comparative plots 
    * module modules_report/bases now include title in the embedded image 


0.8.5
------

* CHANGES:

    * coverage multiqc plugin to remove unique chromomosome name from the sample name, if it is unique
    * bamtools_stats multiqc plugin to remove unique chromomosome name from the sample name, if it is unique
    * fix regression bug in report summary page for the taxonomy tools

* NEWs:

    * Contigs now acceps Spades fasta files
    * added a new viz module for PCA plotting
    * added a new stats module for data normalisation
    * added a new enrichment module to include KeggPathway enrichment



0.8.4
-----

* MultiKrakenResults set nan tp zeros
* fasta module: new filter() method to keep or exclude some entries
* CHANGES: 

    * pipeline_common: PipelineManager renamed into SequanaManager to avoid same
      name as in snaketools.PipelineManager
    * snaketools: PipelineManagerGeneric is now the base class for
      PipelineManager. factorise common methods into PipelineManagerGeneric
    * SnakeMakeStats now skip the figure if stats is empty

* NEWs:

    * snaketools: new PipelineManagerDirectory class (used by
      sequana_demultiplex) when no input files are used, but only an input
      directory. we can now import exiting project to fill the config file automatically
    * module fastqc to plot static images of fastqc contents. The parsing method
      was adapted from multiqc.module.fastqc itself as a quick hack to
      incorporate in the sequana_fastqc pipeline 
    

0.8.3
-----

* Pipeline related:

    * add plot_stats function in PipelineManager
    * add setup/teardown function in PipelineManager

* CHANGES:

    * new dependency (bx) to include new features in bamtools module (infer
      strandness, insert size etc)
    * sequana_taxonomy now handles both kraken1 and kraken2 databases
      the taxonomy databases is now downloaded from ncbi directly and build
      within sequana. The standalone can download the taxonomy data itself with
      --download-taxonomy option. New pipeline available on 
      github/sequana/sequana_taxonom
    * rename KrakenHierarchical into KrakenSequential
    * refactoring of kraken.py and kraken_builder to include kraken2
    * taxonomy module can now aggregate several taxonomic databases. We now
      retrieve the input data from NCBI instead of EBI to be up-to-date whenever
      we want (e.g. daily basis)

* BUGs:

    * add package *packaging* in requirements (for travis and pip installation)

* NEWs

    * 7 new modules in sequana.viz dedicated to visualisation, in particular
      dendogram, corrplot, hist2D from biokit project (merging with sequana) 
      together with a bunch of notebooks
    * NEW bed module to include a BED class to read 12-column BED files.
    * gff3 module: a new method to create utility files for the rnaseq pipeline
    * New rule for the rnaseq pipeline: salmon
    * kraken2 database can now be used within sequana_taxonomy

* General cleaning: lane_merging tools with tests, deprecated warnings
* Move from 70 to 85% of testing coverage


0.8.2
-----

* a quick fix in the requirements to build proper recipes on bioconda, now
  available on bioconda channel (1 April 2020).


0.8.1
-----

* add missing dependency *packaging* in requirements.txt
* add new sphinx extension to include external pipeline documentation
* for pipelines, we created a common epilog and prolog statement and a
  init_pipeline to insert before parsing the options. This allows to have the
  --deps argument to print all dependencies of a pipeline

0.8.0 (24 March 2020)
---------------------

This is an stable release made to includes lots of new features and pipelines
and bug fixes made.

We removed all pipelines from sequana. Pipelines have now their own repositories
on github to ease the developpement of sequana and those pipelines. The
rationale being that we do not need to update sequana when a pipeline changes
and a pipeline can have its own biocontainer and life cycle. We kept just one
for testing purposes (compressor). Other pipelines have now their own repositories:

- sequana_coverage
- sequana_demultiplex
- sequana_fastqc
- sequana_revcomp
- sequana_rnaseq
- sequana_variant_calling
- sequana_denovo
- sequana_pacbio_amplicon
- etc

New pipelines have been added such as sequana_downsampling and sequana_mapper.


* NEWS:

    * snpeff now includes GFF input file input as an option (in addition to 
      the GBK)
    * new module trf for tandem repeat finder output.
    * new scripts:

        * sequana_start_pipelines to initiate a new pipeline from the cookiecutter recipes
        * sequana_gtf_fixer to help fixing GTF issues
    * new module *gff3* to read GFF files
    * Module can now encapsulate logos
    * Module version implemented
* BUG:

    * snpeff_add_locus_tag: if contig name and length in GFF header not in the
      same order as in the fasta, a new fasta with wrong header was created.
      This caused trouble in the variant_calling pipeline
    * kraken: Fix kraken plot (matplotlib version) when 100% of the reads are
      classified
    * Header of igvtools count output may vary. Make the consensus.get_bases more
      robust to automatically identify number of lines to skip.
    * Fix the kraken multiqc report
    * Fix bug in gui/browser to fix import of QWebPage on travis
    * bowtie2 dynamic rule now uses templating correctly (RNASeq pipeline)
    * Fix issue in snaketools for input_readtag set to _[12] for paired data
      The paired attribute wass wrongly set to unpaired. Besides, we make it
      more robust for those who tag their paired data with _1 and _2 instead of
      _R1_/_R2_
    * Repeats: for multi fasta with similar header, we were expecting the chrom
      name to be unique but underlying tool uses regular expression. So, this was
      buggy when chrom name were starting with same string. e.g chr1 anc chr11.
    * multiqc section of sequana_coverage: duplicate chrom names across multiple
      samples were shown as a single entry in the report. 
    * draft version of multiqc for sequna_quality_control now available
* MAJOR CHANGES/FIXES:

   * The main script 'sequana' is redundant with the new framework of
     pipelines. It has been removed in this version
   * sequana_coverage now handles low coverage correctly in the 
     HTML reports.Fix the ylimits of the coverage plot for low coverage.
   * cutadapt rules was failing due to a stricter optional/positional argument
     handling. Fixed the rule accordingly.
   * sequana_lane_merging is now ready for production. changes made: copy of
     the script in the local directory, not the data directory. 
   * RNASeq pipeline: removed sartools, kraken. Fixed bamCoverage rule. Simplify
     usage related to indexing and mapping. Fixed igvtools rule. Fixed the
     reoderSam rule (wrong executable). Fixed a incorrect parameter name in
     bamCoverage rule. Fixed incorrect Snakemake syntax in the fastq_screen
     rule and RNAseQC. Fixed another deprecated rule: fastq_screen_report.
   * New pipeline_common module to be used by all pipelines 
* MINOR CHANGES/FIXES

    * snaketools:

          * pipelines discovery updated in ModuleFinderSingleton. Finally
            fixed the lost of comments in the config when saved. 
          * Removed onweb() method. 
          * Fixed the loss of comments when saving yaml file after an update
            of the key/value. 
          * remove check_sequana_fields.
          * more tests and cleanup 
    * demultiplex: fix a Pandas deprecated warning (add sort argument in pd.concat)
    * python dependencies not in conda are not harcoded inside the setup.py
      (itolapi). add cython into the list of requirements.
    * Fix deprecated bamCoverage rule to use newest deeptools version.
    * The check_config_with_schema function now performs the validation
      correctly
    * Fix stdout of the fastqc, unpigz, bowtie1 and bowtie2 rules
    * Atropos 2.0 changed its API. fastq module compat with atropos 1.0 and 2.0 


0.7.2
----------

* NEWS:


    * New script: sequana_fastq_summary included in fastqc pipeline
    * New script: sequana_substractor to remove reads that mapped against a reference(s)
    * added a new module to upload/export phylogenetic tree on itol website.
      Used in the laa pipeline
    * added backspace2fusion code to merge lanes in Illumina raw data
    * added new pipeline called fastqc to simply run fastqc + multiqc in parallel
    * added laa pacbio pipeline
    * multiqc modules: bamtools_stats and kraken module for the laa pipeline
    * added test file and test for SIRVRerence class (partial fix of issue #504)
    * added Makefile class in snaketools to help building pipeline
    * added MultiKrakenResults class
    * sequanix and snaketools now handle the presence of a multiqc_config 
      file in the pipeline module
    * add laa multiqc

* BUGS:

    * in quality_control when using the design file in cutadapt rule
    * Fix multiqc report for pacbio_qc pipeline


* CHANGES:

    * adapters added: TruSeqCD, TruSeqUD, etc
    * adapters removed: rubicon
    * remove clean_ngs rule and code related to this software, not used in sequana

* CHANGES for developers:

    * adapters are now named NAME_fwd.fa instead of adapters_NAME_fwd. This
      should not affect the user interface. Also, the index sequence stored in the
      adapter files are now identical in the forward/reverse/revcomp versions
      to simplify the code. We also added a script in ./resources/data/adapters
      to create the rev and revcomp version automatically.
    * add missing xlrd dependencies in requirements


0.7.1
---------

* NEWS:

    * added metropolis hastings module
    * added a sniffer module for BAM/SAM/CRAM
    * added a SMA/CRAM reader

* CHANGES:

    * refactoring of bamtools. added SAM and CRAM classes. remove the
      plot_acgt_content method. Instead of inheriting from pysam.Alignement, 
      we store the data as an attribute.

* FIXES:

    * cutadapt rules and expdesign can now handle sample names with several
      underscores
    * Issue 515: sequanix should now be able to handle list in YAML files
    * Issues 520: level info in sequanix was always set to INFO at start time
    * Issue 519: fix issues in sequanix due to different API in new ruamel.yaml version
    * Issue #522: fix bam_splitter tool


0.7.0
------

* BUGS:

    * add /1 and /2 in quality control pipeline https://github.com/sequana/sequana/issues/508
    * Fix test failure due to freebayes version 1 and 1.2 https://github.com/sequana/sequana/issues/512
    * Fix reading of SampleSheet for MiSeq: https://github.com/sequana/sequana/issues/511
    * Add Exp Design checked in quality control pipeline: https://github.com/sequana/sequana/issues/500

* CHANGES:

    * sequana_vcf_filter: finalised version with INDEL removal, filters on DP4
      and AF1 fields
    * rename PacbioBAM into PacbioSubreads

0.6.X
-----
0.6.5
~~~~~~~~~~~

* CHANGES:

    * sequana_coverage. Major refactoring of bedtools module to handle large
      data sets (human), and provide ability to focus on CNVs using an
      additional naive clustering (merge_rois_into_cnvs method) and binning. 
      We can also analyse data chunk by chunk (to avoid filling the memory). 
      added a plot_rois function
    * sequana_coverage standalone: add the --chunksize, --cnv-clustering and
      --binning options.

* NEWS:

    * add cnvnator class
    * coverage pipeline added in the pipelines

* BUGS:

    * Fix silent warning (regex) in snpeff module
    * double indexing adapters issue for Nextera fixed: https://github.com/sequana/sequana/issues/501


0.6.4
~~~~~~~~~~

* BUGS:

    * Fix issue https://github.com/sequana/sequana/issues/380 is_sorted property
      of the BAM class.
    * Fix --no-report option in sequana_coverage and add --clustering (double
      threshold option)
    * pacbio_qc pipeline is now able to also read old pacbio format

* NEWS:

    * SARTools rule added and used in the RNAseq pipeline
    * add summary module to store summary in json formats.
    * simple vcf_filter standalone

* CHANGES:

    * pin kraken version to 1.1 (newest on bioconda)
    * MAJOR REFACTORING of bedtools and sequana_coverage standalone. In
      particular, change default window size to 20,001 or a fifth of genome 
      length (for small genome); speed up code; add plot_roi function, uses
      multiqc for summary page; add log2 ratio column. See
      https://github.com/sequana/sequana/issues/495 for details. 
      Scan large files by chunk. Add a snakemake that can be used in sequanix.
    * remove the sequana_report standalone, which was not finalised and won't be
      used in the future. We will use multiqc instead.

0.6.3.post1
~~~~~~~~~~~~~

- a bug fix in the sequanix GUI and singularity. a statement related to the 
  ruamel.yaml package causes trouble if version is not 0.15. A temporary fix
  consisted in adding a try/except (the statement is just a warning.filter and
  has no impact on analysis)

0.6.3
~~~~~~~~~~~

* BUGS:

    * Fix bug in the copy of the fastqc data sets in the quality control
    * atropos bug in the reports (not full). Bug reported to atropos
      github. https://github.com/jdidion/atropos/issues/57. Need to use version
      1.1.16
    * kraken report table were not sorted by percentage (as expected). Also,
      if the case of poor databases with few entries, the output may contain
      lots of classified sequences with Taxon 1, which was not reported
      correctly in the krona plot.

* NEWS:

    * example of a schema.yaml implemented for the quality control.
    * sequanix: reads schema.yaml automatically for sequana pipelines
      and can import one for generic cases. An option in the preference
      was added to switch on/off the validation of the config file with this
      schema. Can also import schema file for the generic case.

* CHANGES:

    * Taxonomy file is downloaded for Kraken only when Kraken is used, not in
      the main __init__  file anymore.


0.6.2
~~~~~~~~~~~~

* BUGS:

    * Fix regression bug (https://github.com/sequana/sequana/issues/484)
    * Fix missing N_final column in table of the quality_control multi-summary
      page
    * Remove phix174.fa requirements in RNAseq pipeline config file
    * Fix path starting with tilde (https://github.com/sequana/sequana/issues/486)

* NEWS:

    * add isoseq Class
    * add vcf_filter module back to help in filtering VCF files created with
      mpileup for instance
    * add sequana_vcf_filter standalone
    * add cigar module to help deciphering CIGAR strings

0.6.1
~~~~~~~~~~

* BUGS:

   * pipeline quality control: fix https://github.com/sequana/sequana/issues/477
   * Fix empty dependency list in HTML report if sequana installed with conda

0.6.0
~~~~~~~~~~~~~

 * BUGS:

   * add missing file for the RNAseq pipeline in the setup.py
   * Fix RTD building
   * Fix reag_tag filtering https://github.com/sequana/sequana/issues/480 
   * Set singularity hub (v2.4)

Prior 0.5.X
-----------
0.5.2
~~~~~~~~~~~~~~~

* BUGS:

    * cutadapt rule: remove the '--progress bar' for now because of a bug in atropos
      (reported) that fails in the progress bar code

* Updates:

    * pipeline pacbio_qc: finalise output tree structure.
    * pipeline quality_control: add sanity check (thread must be >1 for
      atropos) and run fastqc on unmapped data (rather than mapped).  
    * pin atropos version to 1.1.10 and added to requirements.txt
    * Fix parsing of atropos report
    * Update FastQC significantly to use atropos FastqReader instead of pysam.FastxFile
    * documentation for the installation (remove docker, add singularity)
    * rule/module atropos: implement ability to parse json report from atropos
      https://github.com/sequana/sequana/issues/448
    * rule fastqc: the log is now a variable. all pipelines using this rule
      have been updated to save the log in {sample}/logs/ intead of ./logs
    * add polyT in TruSeq adapters

* News:

    * add Singularity container
    * BAM class (bamtools module): add plotting methods (coverage, letters,
      indels)
    * Add Cigar class (cigar module).
    * Sequanix: add option to switch on/off the tooltips
    * rule cutadapt: (1) check whether thread is set to > 1. if not set to 2
      (2) add --report-format to save reports in JSON and TXT

0.5.1
~~~~~~~~~~~~~~~

* BUGS:

   * Set -t thread options correctly in the different rules (e.g. cutadapt)
   * pipeline variant_calling: fix the VCF inputs when snpeff is off .
     See https://github.com/sequana/sequana/issues/471
   * pipeline quality_control. Fix regression bug introduced by the use 
     of sambamba in the bwa_mem_dynamic rule (see 
     ihttps://github.com/sequana/sequana/issues/472)
   * Fix wrong total bases values in summary report of the quality_control
     pipeline computed in FastQC class (see 
     https://github.com/sequana/sequana/issues/470)
   * pipeline pacbio_qc: hard-coded the number of threads to 4 otherwise may
         fail on clusters. Does not change the pipeline or analysis itself
   * sequana_coverage: fix chromosome option.
   * Fix genbank_parser when the genbank contains several concatenated genbank
     entries. This fixes the coverage reports CSV file that had missing
     annotations.
   * Fix regression bug introduced in rule bwa_mem_dynamic that messed 
     up R1 and R2 order as compared to samtools by using sambamba. Fixed by
     using -N parameter.
   * Fix the -p option to be before the input whenever pigz is used in a rules. 
     Indeed -p may be ignored otherwise e.g. on clusters.

* Updates:

   * add pacbio option in the mapping code
   * pacbio_qc: fix pattern to filter input BAM files
   * Speed up fastq_count (https://github.com/sequana/sequana/issues/465)
   * bamtools module: speed up initialisation. add is_sorted method.
   * bedtools: limit number of points to 1,000,000 in plot_coverage and set
     ylimits manually to 6 mean coverage. add __eq__ function. See #464 issue
   * Repeats can handle FastA properly (not limited to first sequence anymore)
   * sequana_mapping: add thread in samtools call



0.5.0 august 2017
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tag a stable release



0.4.2 August 2017
~~~~~~~~~~~~~~~~~~~~~~

* Updates:

  * pipeline: variant calling cleanup and finalised
  * pipeline: denovo updated (busco) and cleanup and finalised
  * pipeline: pacbio_qc finalised 
  * pipeline: rnaseq: finalised
  * module pacbio:  speed up initialisation; add a random_selection method; add a summary method;

* NEWS:

  * Sequanix: can now load cluster config
  * new rules: busco, busco_analysis, canu
  * new pipeline: pacbio_denovo
  * multiqc modules integrated in sequana. See :ref:`developers` for details.
  * module snaketools: new function get_pipeline_stats
  * new gallery example with statistics about the pipelines

* CHANGES:

  * remove random() function from FastQ (useless and will be put in new module
    simulation)  


0.4.1 July 2017
~~~~~~~~~~~~~~~~~~

* Update of Variant calling and denovo pipelines with HTML report creation
* Fix #421 (check for dot command in sequanix)
* Fix #420 (sequanix browser on Mac)
* sequana_coverage #417 division by 0 fixed
* snpeff bugs for special genbank cases fixed


0.4 July 2017
~~~~~~~~~~~~~~~~~~~~~

* Master release for sequanix


0.3 April-June 2017
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* BUG FIXES:

    - sequanix:
        - rulegraph issue on SLURM system. Avoid the os.chdir
    - fastq_samples/ fastq module: fix histogram_gc_content maximum range
    - rulegraph rule: fix issue #405 (spaces in path to snakefile)
    - genome coverage was buggy for multi chromosome and circular option on. Fixed
    - adapters/expdesign modules: fixe the case of design files with same sample
      name and same index but different lanes.
    - sequana_coverage. Fix Issue #416 (float division by zero)

* CHANGES:

    - sequanix:
        - snakemake output is now cleared when pressing RUN
    - quality_control pipeline: default to atropos instead of cutadapt for
          adapter trimming. Kraken: remove classified reads and keep
          unclassified. Unclassified reads are now compressed.
          unclassified reads that are also compressed now.

* NEW:

    - pacbio module: cleanup and add funcion to convert input BAM into Fasta
    - sequence module: Repeats class added
    - new Snakemake pipeline called qc_pacbio to perform quick QC and taxonomy analysis
          for pacbio
    - add ORD, CDS, GC SKEW in sequence module.


0.2. - March - April 2017
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* NEWS:
    - RNA-seq pipeline added (single-end only, paired-end upcoming)
      including all indexes for RNA-seq
    - Hierarchical kraken available
    - add new standalone called **sequana_fox** to expose the pyqt5 browser.
    - Sequanix first release
    - final version of the variant calling, denovo, quality_control and rna-seq
      pipelines.

* CHANGES:

    - Sequanix/Sequana:
      - config file can have the yml extension (in addition to yaml)
      - dropdown widgets in the form based on the docstrings in the config file
      - can import config to override default sequana config file
      - subprocesses killed when the main pipeline is stopped

0.1.21 - Feb 2017
~~~~~~~~~~~~~~~~~~~~~~~~

* NEWS:

    - add sequana_debug_level function at top level to switch verbosity of
      informative messages (default is WARNING).
    - add pacbio module  #351
    - quality control pipeline: atropos can be used in place of cutadapt #346

* CHANGES:

   - Running Median is 10 times faster #345
   - sequana_coverage:  (1) --file1 alone was not working (2) automatically copy
     cluster-config in working directory and update runme.sh accordingly #342
   - sequana standalone:
       - handles cluster_config Snakemake option
       - add error message when adapter name is incorrect
   - sequanix: the help dialog is now created inside designer and has a proper
     scrollable browser dialog. cluster_config Snakemake option is also handle.
   - Remove galleria JS lib and related files (htmltools)
   - sequana_coverage: add --logging-level option

* BUG:

    - Fix #352 : allow gc window size to be even (warning is shown and +1 to
      window size)
    - Fix # 354: cutadapt report that was mixing up R1/R2 trimming in the images.
    - --output-directory in sequana_coverage was failing 
    - in coverage, centralness was buggy (regression) and use number of ROIs
      instead of the total base length #347
    - Fix multi_report summary for single end case #349

0.1.20 - Feb 2017
~~~~~~~~~~~~~~~~~~~~~~~~

* CHANGES: 

    - remove pyquickhelper dependencies and add a simple rest2html function in
      misc module.

0.1.19 - Feb 2017
~~~~~~~~~~~~~~~~~~~~~~~~


* CHANGES:

    - misc module: factorise on_cluster() function used in compressor scripts to
        be used in other tools such as sequanix
    - compressor: limits max number of jobs to 20 (can be bypass manually),
      prevent run on TARS if snakemake-cluster not provided. 
    - rules:
        - dag: now the snakemake is called inside a temporary directory to avoid
          clash with the current snakemake process. This avoid error message. 
          Fixes https://github.com/sequana/sequana/issues/331
    - __init__ was optimized as well as many modules to make use of the lazy
      import mechanism. The reporting package is not part of the exposed module. 
      So::

         from sequana import BAMReport

      is now::

         from sequana.reporting.report_bam import BAMReport

* NEWS:

    - Sequanix stable version
    - add TrueSeq adaptors
    - add lazy import mechanism to speed up the time to import sequana, which 
      speeds up the --help in the standalone 


0.1.17/0.1.18 - Jan 2017
~~~~~~~~~~~~~~~~~~~~~~~~

:Main NEWS: The GUI was completed and the current pipelines stabilised (RNA-seq,
    quality control, variant calling). The test suite was switched from nosetests to
    pytest, in particular to perform tests more eaasily on the Qt GUI. 


* BUG Fixes:

    - experimental design and adapters API simplified fixing a few bugs in the
      process. Doc and tested finalised.
    - Fix cutadapt rules, which was not filling the fwd and rev properly anymore
      when using the design file.
    - in sequana main script, --reference was used by quality_pipeline only.
      Now, available for all.
    - Fix the main script for the reference in variant calling pipeline.


* CHANGES:

    - sequana_compressor: for conversion from e.g gz to bz2, use a pipe instead
        of double IO. Updated docs and tests ready for production.
    - sequana standalone: 
      - --pattern changed to --input-pattern
      - --output-directory changed to --working-directory
    - remove pipetools module (obsolet)
    - GUI revisited with qt designer + can now also read any snakefile/config
      file combo (not just sequana pipelines)
    - RULES: adapters can now use adapter_type without a design (fwd and rev
      gets filled automatically)

* NEWS:

    - add rubicon adapters
    - add ability to read JSON in SequanaConfig

2016
----------
0.1.16
~~~~~~~~~~~

* BUG Fixes:

    - Fix sequana_taxonomy (https://github.com/sequana/sequana/issues/308)
    - Fix typo in sequana_coverage for multiple chromosome (https://github.com/sequana/sequana/issues/307)

* NEWs:

    - SequanaConfig can read back a SequanaConfig instance
    - Added a DummyManager for minimalist manager to create reports


0.1.15
~~~~~~~~~~~

* CHANGES:

    - coverage: https://github.com/sequana/sequana/issues/302
      add histogram, better stats table. add --output-directory
    - Update docker (add bowtie, subread, firefox)
    - snaketools:
          - empty strings are kept as empty strings (not None)
          - remove check() method in SequanaConfig
          - cleanup (removing of templates) ca be switch off

0.1.14
~~~~~~~~~~~

* CHANGES:

    - fastqc.histogram_sequence_lengths (log2 scale to log10)
    - multi_summary fixed and available for the quality_control pipeline
    - sequana_compressor: add --keep-going option by default so that if a file
      fails, other independent files are processed.
    - snaketools:
          - remove SnakeMakeProfile (not used)
          - remove sequana_check_config (not used)
          - remove deprecated __get_tagname
          - remove ExpandedSnakefile since not required anymore
          - Fix sample_file2 option that was not encoded properly
          - PipelineManager and SequanaConfig use new yaml parser
    - sequana_coverage: -- add back the sample name as prefix of the HTML report
      name -- a BED with two coverage columns is now accepted --
      --download-genbank option added
    - sequana_summary works for the quality_control pipeline
    - Simplify combos of input_directory, input_patter, input_samples, the new
      possible mutually exclusive input parameters of sequana standalone and all
      pipelines.

* BUGS:

    - Kraken: if no reads classified at all, errors were raised and
      quality_control summary report would fail. This is fixed now with a "nodata"
      image being shown.

* NEWS

    - GUI (draft version)
    - fq.gz are now allowed in the pipelines and should be supported in the
      future
    - More tests in particular a ./test/pipelines/ new directory


0.1.13
~~~~~~~~~~~

* CHANGES:

    - revisited all pipelines so that they can work of multi samples.
    - quality_phix, quqlity and quality_taxon pipelines merged in
      quality_control pipeline
    - running meadian won't fail anymore with odd window size (we add +1)
    - rulegraph is used as well as dag to create figures of the pipelines

* NEWS:

    - compressor: includes dsrc format in addition to bz2 and gz
    - snakemake rule extension for sphinx
    - add a pipeline manager in snaketools to handle all pipelines
    - a designexp module to handle adapter design files


0.1.12
~~~~~~~~~~~

* BUGS:

   - Fix bug in cutadapt pipeline when there is no adapters. Force a dummy
     adapters (XXX) otherwise trimming is performed on read1 only

* NEWS:

    - compressor rule and script available.
    - coverage annotation
    - multiple_summary draft

0.1.11
~~~~~~~~~~~

* NEWS:

   - add a docker
   - sequana_summary standalone
   - sequana_mapping standalone
   - Module has an overview field

* BUG FIXES:

   - cutadapt report handles single-end tables. Fix the reverse complement
     adapter files for the paired-end case

* CHANGES:

    - sequana_standalone: final version with stats



0.1.10 - July 2016
~~~~~~~~~~~~~~~~~~~~~~~~

* NEWS:

    - sequana_coverage standalone
    - de-novo pipeline

* CHANGES:

    - Remove AdapterDB, a draft version that uses Kraken to detect adapters. Not
      relevant anymore
    - config.yaml is now in each pipeline to have a simplified version
    - sequana can known use single_indexed or multiple_indexed  adapters, which
      are also provided within sequana (Nextera and PCR free cases)
    - Release for production (quality_taxon pipeline)


0.1.7 to 0.1.9 - July 2016
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* NEWS:

  - rule data added and used in phix_removal (fastq_sampling + raw data switch)
  - kmer module
  - sequana_taxonomy standalone

* CHANGES:

  - reports are now in ./sequana/reporting
  - MAJOR refactoring of report/ directories in all pipelines to make them
    independent from the temporary analysis, which can then be removed.

* BUGS:

  - Fix running median issue in bedtools (window size larger than contig size)



0.1.6 - June 2016
~~~~~~~~~~~~~~~~~~~~~~


* NEWS:

  - KrakenDownlad class: download kraken_toydv from sequana/data repository or
    minikraken into a local directotry
  - New method in FastQC to show ACGT content
  - Genomecov renamed into GenomeCov
  - Update main script significantly to create multiruns and handle adapters
  - GC content and plot GC vs coverage added in GenomeCov

* CHANGES:

  - sequana_data by default looks into resources/testing directory
  - in fastq module: FastQC a bit faster andFastQRandom class removed
  - add a moving_average function in misc module

* BUGS:

  - sequana_data was showing __init__ and __pycache__ as possible data sets
  - databases: filelist as a list was not implemented
  - in fastq.FastQ extra_head in gzip mode was missing the last row



prior 0.1.5 June 2016
~~~~~~~~~~~~~~~~~~~~~~

* NEWS

  - sequana_taxonomy standalone available (kraken + krona)
  - sequana standalone available
  - quality_taxon pipeline available
  - module coverage for theoretical computations
  - add gallery in the documentation

* CHANGES:

  - module vcf_to_snpeff renamed as snpeff

* BUG:

  - Fix bug in running median (shift)

