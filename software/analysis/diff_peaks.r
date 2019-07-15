library(edgeR)
library(stringr)

###function for doing differential peak analysis:
  #reads is a matrix of reads
  #control is the name of the control 
  #toCompare is the names of things to compare to the control
  #padj.Threshold is the adjusted p value threshold for calling binaries
  #log2FC.Threshold is the log2FC threshold for calling binaries

diff_peaks = function (reads, control, toCompare, padj.Threshold=.1, log2FC.Threshold=1) {
 
  experimental.groups = str_replace(colnames(reads), '_.','')
  dge <- DGEList(counts = reads, group = experimental.groups)

  dge.norm <- calcNormFactors(dge)

  # #As of right now I'm deciding not to filter, as it removes almost all the features, 
  # #Which seems sort of silly for this type of data
  # isexpr <- rowSums(cpm(dge.norm) > 20) >= 2
  dge.norm.fltd = dge.norm
  # dge.norm.fltd <- dge.norm[isexpr, ]

  #factor of the treatment groups
  f = factor(experimental.groups)
  #make a design matrix for limma 
  design <- model.matrix(~0 + f)
  #removes the f from the design matrix names
  colnames(design) <- sub("f", "", colnames(design), fixed=T)

  # limma voom to transform count data to log2-counts per million (log2-cpm) with associated weights
  # estimate the mean-variance relationship and use this to compute appropriate observational-level weights.
  # The data are then ready for linear modelling.
  y <- voom(dge.norm.fltd, design)
  #fit linear model
  fit <- lmFit(y, design)

  #make contrast matrix to specify the comparisons to make (everything vs control)
  contrast.matrix <- makeContrasts(contrasts=paste0(toCompare, '-', control), levels=design)
  #compute fit for contrasts
  fit2 <- contrasts.fit(fit, contrast.matrix)

  # Given a linear model fit, compute moderated t-statistics, moderated F-statistic, 
  # and log-odds of differential expression by empirical Bayes moderation of the standard errors towards a common value
  fit.treat <- treat(fit2, log2FC.Threshold)
  
  log2FC.all.binary = decideTests(fit.treat, p.value=padj.Threshold)
  
  log2FC.all = as.data.frame(map_dfc(1:ncol(log2FC.all.binary),
                                     function(coef) topTreat(fit.treat, coef=coef, number=Inf,
                                                             sort.by='none')$logFC))
  rownames(log2FC.all) = rownames(log2FC.all.binary)
  colnames(log2FC.all) = colnames(log2FC.all.binary) = str_replace(colnames(log2FC.all.binary), str_c('-', control), '')
  DE_log2FC.all = log2FC.all[rowSums(log2FC.all.binary != 0) > 0,]
  log2FC.all.binary = log2FC.all.binary[rowSums(log2FC.all.binary != 0) > 0,]
  
  
  return(list(log2FC = DE_log2FC.all, binary = data.frame(log2FC.all.binary),
      fit=fit, design=design, all_log2FC = log2FC.all, dge=dge.norm.fltd))
}


