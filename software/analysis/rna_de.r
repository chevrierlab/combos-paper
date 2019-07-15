#######
# runs the base differential expression analysis on rna data
# data = a dataframe of counts
# experimental.groups = the group for each sample
# others: scalar parameters
#######

rna_de_analysis = function(data, experimental.groups, cpm_threshold = 50, cpm_needed = 2, 
                           padj.Threshold = .01, log2FC.Threshold = 1.5) { 
    
    #edgeR to make DGElist with raw counts and groups
    dge <- DGEList(counts = data, group = experimental.groups)
    
    # edgeR normalization
    dge.norm <- calcNormFactors(dge)
    # filter out genes that do not have >50 reads in at least in 2 samples
    isexpr <- rowSums(cpm(dge.norm) > cpm_threshold) >= cpm_needed
    dge.norm.fltd <- dge.norm[isexpr, ]
    
    ## design matrix
    # create factor with all conditions in experiment
    f <- factor(experimental.groups, levels=unique(experimental.groups))
    # build a design matrix for modeling with limma voom and indicating condition for each sample
    design <- model.matrix(~0 + f)
    # renames column names in the design (remove extra "f" at beginning of each group name)
    colnames(design) <- sub("f", "", colnames(design))
    
    # limma voom transforms count data to log2-counts per million (log2-cpm) with associated weights
    # estimate the mean-variance relationship and use this to compute appropriate observational-level weights.
    # The data are then ready for linear modelling.
    y <- voom(dge.norm.fltd, design, plot=F)
    # fit linear model for each gene given a series of arrays
    fit <- lmFit(y, design)
    
    # list of all non-control conditions
    all.conditions <- levels(experimental.groups)[
        order(str_count(levels(experimental.groups), '_'))]
    all.conditions = all.conditions[all.conditions != 'Ctrl']
    
    # contrast matrix which specifies comparisons to make between RNA samples
    contrast.matrix = makeContrasts(contrasts=str_c(all.conditions, '-Ctrl'), levels=design)
    
    # compute estimated coefficients and standard errors for a given set of contrasts.
    fit2 <- contrasts.fit(fit, contrast.matrix)
    # compute moderated t-statistics, moderated F-statistic, and log-odds of differential 
    # expression by empirical Bayes moderation of the standard errors towards a common value
    fit.treat <- treat(fit2, log2FC.Threshold)
    
    log2FC.all.binary = decideTests(fit.treat, p.value=padj.Threshold)
    
    log2FC.all = as.data.frame(map_dfc(1:length(all.conditions),
                                       function(coef) topTreat(fit.treat, coef=coef, number=Inf,
                                                               sort.by='none')$logFC))
    CI.win = as.data.frame(map_dfc(1:length(all.conditions),
                                   function(coef) {x = topTreat(fit.treat, coef=coef, number=Inf,
                                                                sort.by='none', confint=T);
                                   return((x$logFC - x$CI.L)*2)}))
    rownames(log2FC.all) = rownames(CI.win) = rownames(log2FC.all.binary)
    colnames(log2FC.all) = colnames(CI.win) = colnames(log2FC.all.binary) = all.conditions
    DE_log2FC.all = log2FC.all[rowSums(log2FC.all.binary != 0) > 0,]
    log2FC.all.binary = log2FC.all.binary[rowSums(log2FC.all.binary != 0) > 0,]
    CI.win = CI.win[rowSums(log2FC.all.binary != 0) > 0,]

    
    return(list(log2FC = DE_log2FC.all, binary=data.frame(log2FC.all.binary), 
                fit=fit, design=design, dge=dge.norm.fltd, CIwin = CI.win, contrast_fit = fit2, all_log2FC = log2FC.all))
}



alt_rna_de_analysis = function(data, experimental.groups, cpm_threshold = 50, cpm_needed = 2, 
                           padj.Threshold = .01, log2FC.Threshold = 1.5) { 
    
    #edgeR to make DGElist with raw counts and groups
    dge <- DGEList(counts = data, group = experimental.groups)
    
    # edgeR normalization
    dge.norm <- calcNormFactors(dge)

    dge.norm.fltd <- dge.norm
    
    ## design matrix
    # create factor with all conditions in experiment
    f <- factor(experimental.groups, levels=unique(experimental.groups))
    # build a design matrix for modeling with limma voom and indicating condition for each sample
    design <- model.matrix(~0 + f)
    # renames column names in the design (remove extra "f" at beginning of each group name)
    colnames(design) <- sub("f", "", colnames(design))
    
    # limma voom transforms count data to log2-counts per million (log2-cpm) with associated weights
    # estimate the mean-variance relationship and use this to compute appropriate observational-level weights.
    # The data are then ready for linear modelling.
    y <- voom(dge.norm.fltd, design, plot=F)
    # fit linear model for each gene given a series of arrays
    fit <- lmFit(y, design)
    
    # list of all non-control conditions
    all.conditions <- levels(experimental.groups)[
        order(str_count(levels(experimental.groups), '_'))]
    all.conditions = all.conditions[all.conditions != 'Ctrl']
    
    # contrast matrix which specifies comparisons to make between RNA samples
    contrast.matrix = makeContrasts(contrasts=str_c(all.conditions, '-Ctrl'), levels=design)
    
    # compute estimated coefficients and standard errors for a given set of contrasts.
    fit2 <- contrasts.fit(fit, contrast.matrix)
    # compute moderated t-statistics, moderated F-statistic, and log-odds of differential 
    # expression by empirical Bayes moderation of the standard errors towards a common value
    fit.treat <- treat(fit2, log2FC.Threshold)
    
    log2FC.all.binary = decideTests(fit.treat, p.value=padj.Threshold)
    
    log2FC.all = as.data.frame(map_dfc(1:length(all.conditions),
                                       function(coef) topTreat(fit.treat, coef=coef, number=Inf,
                                                               sort.by='none')$logFC))
    CI.win = as.data.frame(map_dfc(1:length(all.conditions),
                                   function(coef) {x = topTreat(fit.treat, coef=coef, number=Inf,
                                                                sort.by='none', confint=T);
                                   return((x$logFC - x$CI.L)*2)}))
    rownames(log2FC.all) = rownames(CI.win) = rownames(log2FC.all.binary)
    colnames(log2FC.all) = colnames(CI.win) = colnames(log2FC.all.binary) = all.conditions
    log2FC.all.binary = log2FC.all.binary[rowSums(log2FC.all.binary != 0) > 0,]
    CI.win = CI.win[rowSums(log2FC.all.binary != 0) > 0,]
    
    
    return(list(log2FC = log2FC.all, binary=data.frame(log2FC.all.binary), 
                fit=fit, design=design, dge=dge.norm.fltd, CIwin = CI.win, contrast_fit = fit2))
}
