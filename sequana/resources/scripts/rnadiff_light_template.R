library(DESeq2)
library(BiocParallel)
library(limma)

### DEFAULTS
## Normalized counts using counts function in DESEQ (could use rlog and vst as well)

### Notes
## Report generation is taken out of R.
## Enforcing comparisons so condRef not necessary anymore.
## Annotation (gene names, description of each genes) is taken out of R. 
## Enforcing to put the design so batch effect should be included here

### TODO
## Add types of transformations for clustering in python (vst, rlog ...)
## IF NEEDED: Add locfunc to normalize by a list of reference genes (see rnadiff)
## This would need to decompose DEseq to have hand on estimateSizeFactors

## Used only for DESeq wrapper ()
register(MulticoreParam({{threads}}))

## fitType: ["parametric", "local", "mean"]
## BetaPrior: ["FALSE", "TRUE"]


check_counts_meta_tables = function(counts, meta){
    # Verify rownames in meta are the same as colnames in counts

    if (!isTRUE(all.equal(rownames(meta),colnames(counts)))){
        stop("row names of input design do not match count columns. Check both inputs")
    }
    else {
        message("OK: Count and meta tables seems to correspond")
    }
}

make_DR = function(counts, meta, design, fitType, betaPrior){

    check_counts_meta_tables(counts, meta)

    dds = DESeqDataSetFromMatrix(countData=counts, colData=meta, design=design)
    dds = DESeq(dds, parallel=TRUE, fitType=fitType, betaPrior=betaPrior)

    return(dds)
}

pairwise_comparison = function(dds, comps, meta_col_name,
                               independentFiltering, cooksCutoff){
    ## Perform all comparison specified in comps (a list of size 2
    ## lists). Using the groups specified in metadata table column 'meta_col_name'

    compa_res = list()

    for (comp in comps){
        compa_name = paste(c(comp[[1]], comp[[2]]), collapse="_vs_")
        message(compa_name)
        res_i = results(dds, contrast=c(meta_col_name, comp[[1]], comp[[2]]),
                        alpha=0.05, parallel=TRUE,
                        independentFiltering=independentFiltering,
                        cooksCutoff=cooksCutoff)
        summary(res_i)
	shrinked_lfc = lfcShrink(dds, res=res_i, contrast=c(meta_col_name, comp[[1]], comp[[2]]), type="ashr", quiet=TRUE)
        res_i$log2FoldChangeNotShrinked = res_i$log2FoldChange
	res_i$lfcSENotShrinked = res_i$lfcSE
	res_i$log2FoldChange = shrinked_lfc$log2FoldChange
	res_i$lfcSE = shrinked_lfc$lfcSE
        compa_res[[compa_name]] = res_i
    }
    return(compa_res)
}

export_results = function(res, prefix, orderCol){
    ## Export the DE results table

    # Sort by the column specified
    if (orderCol != ''){
        res = res[order(res[,orderCol]), ]
    }


    write.csv(res, file=paste(prefix, 'degs_DESeq2.csv', sep="_"), row.names=TRUE)
    return(res)
}

export_pairwise = function(res, outdir, orderCol='padj'){
    ## For pairwise comparisons
    for (compa in names(res)){
        export_results(res[[compa]], prefix=paste(outdir, compa, sep="/"), orderCol=orderCol)
    }
}

export_counts = function(dds, outdir){

    # Export both raw and VST transformed counts

    ## rlog is accounting for lib size and apparently vst too (see deseq2 manual)
    ## rlog = rlog(dds, blind=FALSE)

    counts = counts(dds)
    norm_counts = counts(dds, normalized=TRUE)
    vst_counts = getVarianceStabilizedData(dds)

    write.csv(counts, paste(outdir, 'counts_raw.csv', sep="/"), row.names=TRUE)
    write.csv(norm_counts, paste(outdir, 'counts_normed.csv', sep="/"), row.names=TRUE)
    write.csv(vst_counts, paste(outdir, 'counts_vst_norm.csv', sep="/"), row.names=TRUE)

    # INFO: vst_counts uses getVarianceStabilizedData that calculates a variance 
    # stabilizing transformation (VST) from the fitted dispersion-mean relation(s) 
    # and then transforms the count data (normalized by division by the size factor), 
    # The batch effect has an impact on the fitted dispersion-mean relation and therefore 
    # on the vst_counts;
    # the counts_vst_batch is derived from limma that does not seem to use the batch effet.
    # to verify that, you can use a batch effect where effect are encoded as numbers or
    # as strings and see that the counts_vst_norm are different indeed. 



    if (grepl("batch", "{{model}}" )) {
        # vst is used here for simplicity, difference with getVarianceStabilizedData ?
        vsd <- vst(dds)
        assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)

        # remove X first character
        vsd_df <- as.data.frame(assay(vsd))
        colnames(vsd_df) <- sub("^X", "", colnames(vsd_df))

        write.table(vsd_df, paste(outdir, "counts_vst_batch.csv", sep="/"), sep=",", row.names = TRUE)

        # also for simplicity, let us save the PDF/PNG of the corrected PCA for comparisons with sequana
        #pdf(paste("{{images_dir}}", "pca_batch.pdf", sep="/"))
        #plotPCA(vsd, intgroup=c("condition"))
        #dev.off()
        png(paste("{{images_dir}}",  "pca_batch.png", sep="/"))
        plotPCA(vsd, intgroup=c("condition"))
        #dev.off()
    } else {
        # if no batch included, let us delete the files. Indeed, if present, we 
        # can include it in the HTML report
        unlink(paste("{{images_dir}}", "pca_batch.pdf", sep="/"), force=TRUE)
        unlink(paste("{{images_dir}}", "pca_batch.png", sep="/"), force=TRUE)
        unlink(paste(outdir, "counts_vst_batch.csv", sep="/"), force=TRUE)
    }
}

export_dds = function(dds, outdir){
    ## Export full dds table
    write.table(mcols(dds), paste(outdir, "overall_dds.csv", sep="/"), sep=",")
}

versions_report = function(outdir){
    write.csv(as.data.frame(installed.packages()[,c(1,3)]), file=paste(outdir, "versions.csv", sep="/"), row.names=FALSE)
}

export_size_factors = function(dds, outdir){
    size_factors = sizeFactors(dds)
    size_factors_df = as.data.frame(size_factors)
    write.csv(size_factors_df, file=paste(outdir, "size_factors.csv", sep="/"), row.names=TRUE)
}


####################
## MAIN
####################

counts = read.csv("{{counts_filename}}", header=T, check.names=F, row.names=1)
target = read.csv("{{design_filename}}", header=T, row.names=1)

dds = make_DR(counts, target, design={{model}}, fitType="{{fit_type}}", betaPrior={{beta_prior}})

res = pairwise_comparison(dds, {{comparisons_str}}, "{{condition}}",
                          independentFiltering={{independent_filtering}},
                          cooksCutoff={{cooks_cutoff}})


export_dds(dds, "{{code_dir}}")
export_pairwise(res, "{{outdir}}")
export_counts(dds, "{{counts_dir}}")
export_size_factors(dds, "{{code_dir}}")

versions_report("{{code_dir}}")





