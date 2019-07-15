library(tidyverse)
library(data.table)

###############################
# compare stims
# takes the names of two variables
# runs comparison of them and outputs the binary 
###############################
compare_stims = function(a, b, fit, design, binary, p_thresh) {
    contrast = paste0(a,'-',b)
    contrast.matrix <- makeContrasts(contrasts=contrast, levels = design)
    fit2 = contrasts.fit(fit, contrast.matrix)
    fit.DE = eBayes(fit2)
    log2FC.per.condition <- topTable(fit.DE, number = Inf, 
                                     adjust.method = correction.method, sort.by = "B", 
                                     p.value = p_thresh)
    log2FC.per.condition$binary <- ifelse(log2FC.per.condition$logFC > 0, 1, -1)
    # populate all condition table with binary values (1 = up; -1 = down and 0 = not regulated)
    return(as.numeric(row.names(binary) %in% row.names(log2FC.per.condition)))
}


### conversions ###
##this arguably doesn't belong here, but it's gotta live somewhere
initials_to_names = function (group) {
    initials = c("S", "P", "L", "H", "Z", "cG", 'Cp')
    full_names = c('SeV', 'PAM', 'LPS', 'HMW', 'Zymosan', 'cGAMP', 'CpGB')
    for(i in 1:7) 
        group = str_replace(group, initials[i], str_c(full_names[i], '_'))
    return(stringi::stri_replace_last_fixed(group, '_', ''))
}

names_to_initials = function (group) {
    initials = c("S", "P", "L", "H", "Z", "cG", 'Cp')
    full_names = c('SeV', 'PAM', 'LPS', 'HMW', 'Zymosan', 'cGAMP', 'CpGB')
    for(i in 1:7) 
        group = str_replace(group, full_names[i], initials[i])
    return(str_replace_all(group, '_', ''))
}

#produce_combos() takes a vector of strings (stims) and produces a list with all their possible combinations
#The output is a list where element i is a matrix of the combinations of i stims 
#$strings contains all the combos concatenated as strings, and name contains the largest combo
#size limit allows for a max combination size (say 3)

produce_combos = function(stims, size_lim = Inf, str_break='') {
    
    #make a list where element i is a matrix of the combinations of i stims
    combos = list()
    for (i in 1:min(length(stims), size_lim)) {
        combos[[i]] = combn(stims, i)
    }
    
    counts = unlist(lapply(combos, function(x) dim(x)[2]))
    
    #the strings element of the list contains a list of the strings of *all* the combinations
    combos$strings = unlist(lapply(combos, function(x) apply(x, 2, paste, collapse=str_break)))
    combos$name = paste(stims, collapse='')
    combos$counts = counts
    return(combos)
}

#multiple_combos() takes a list of string vectors (stim_groups) and produces a list that concatenates their
#combinations in the vein of produce_combos
#it's useful for situations in which you don't want all possible combinations of your stims
multiple_combos = function(stim_groups, str_break='') {
    all_combos = lapply(stim_groups, produce_combos, str_break='')
    
    combos= list()
    for(i in 1:length(all_combos[[1]]$counts)){
        combos[[i]] = unique(matrix(unlist(lapply(all_combos, function(x) x[[i]])), nrow=i), MARGIN=2)
    }
    
    counts = unlist(lapply(combos, function(x) dim(x)[2])) 
    
    combos$strings = unlist(lapply(combos, function(x) apply(x, 2, paste, collapse=str_break)))
    combos$name = paste(unlist(lapply(all_combos, function(x) x$name)), collapse='_')
    combos$counts = counts
    return(combos)
}

#compare combos performs a differential peak analysis and compares different combinations
#of stimulations.
#inputs:
#reads: peak read count table
#control: the name of the control treatment
#stims: a character vector of the uncombined treatments
#output:
#a data.frame of whether or not each treatment has new regulation at each feature
compare_combos = function(reads, control, stims, padj.Threshold=.01, log2FC.Threshold=1.5) {
    
    #get combinations
    combos = if (is.list(stims)) multiple_combos(stims) else produce_combos(stims)
    
    #get differential peak binaries, and convert them to proper binaries
    source('software/diff_peaks.r')
    diff = diff_peaks(reads, paste0('combo_', combos$name), control, 
                          combos$strings, padj.Threshold = padj.Threshold, log2FC.Threshold = log2FC.Threshold)
    diff$binary = abs(diff$binary)
    
    #for each of the types of stimulation with more than one stim we need to 
    #un-count each of the genes that appeared in their components
    for (i in 2:length(combos$counts)) {
        #these are all the combos at this level
        tuples = apply(combos[[i]], 2, paste, collapse='')
        
        #for each combo at this level. 
        for (j in 1:length(tuples)) {
            print(tuples[j])
            #subs is the list of all the sub-components for this combo
            subs = produce_combos(combos[[i]][,j])$strings
            subs = subs[-length(subs)]
            
            #How this line works:
            #pulls the binaries for all the subs
            #for each gene in those subs, compute the overall OR (which genes are owned by any of the subs)
            #the new binary for this combo becomes itself AND NOT the overall OR of the subs
            diff$binary[tuples[j]] = as.numeric(diff$binary[tuples[j]] &
                                                    !apply(diff$binary[subs], 1, function(x) Reduce(`||`, x, F)))
        }
    }
    diff$counts = combos$counts
    return(diff)
}

# takes binaries, fit, combos and a design matrix, and calculates the new genes
# based on the components of each stim and comparisons with those components
make_comparison_binaries = function(combos, binary, fit, design, p_thresh, str_break='') {
    #for each of the types of stimulation with more than one stim we need to 
    #un-count each of the genes that appeared in their components
    for (i in 2:length(combos$counts)) {
        #these are all the combos at this level
        tuples = apply(combos[[i]], 2, paste, collapse=str_break)
        
        #for each combo at this level. 
        for (j in 1:length(tuples)) {
            #subs is the list of all the sub-components for this combo
            subs = produce_combos(combos[[i]][,j], str_break=str_break)$strings
            subs = subs[-length(subs)]
            
            new_name = str_c(tuples[j], "_new")
            sub_test = sapply(subs, compare_stims, b = tuples[j], fit = fit, design = design,
                              binary = binary, p_thresh = p_thresh)
            binary[new_name] = 
                as.numeric(binary[tuples[j]] &
                               !apply(binary[subs], 1, function(x) Reduce(`||`, x, F)) &
                               apply(sub_test, 1, function(x) Reduce(`&&`, x, T)))
        }
    }
    return(binary)
}

# compare_combos_contrasts() does the same thing as compare_combos(), but uses the method
# of contrasting triplets with their pairs, etc.
compare_combos_contrasts = 
    function(reads, control, stims, padj.Threshold=.01, 
             log2FC.Threshold=1.5, p.compare.Threshold = .1, lfc.compare.Threshold = 0) {
        #get combinations
        combos = if (is.list(stims)) multiple_combos(stims) else produce_combos(stims)
        
        #get differential peak binaries, and convert them to proper binaries
        source('software/diff_peaks.r')
        diff = diff_peaks(reads, control, 
                              combos$strings, padj.Threshold = padj.Threshold, log2FC.Threshold = log2FC.Threshold)
        diff$binary = abs(diff$binary)
        
        diff$binary = make_comparison_binaries(combos, diff$binary, diff$fit, diff$design, p.compare.Threshold)
        diff$counts = combos$counts
        return(diff)
    }


plot_combos = function(log2FC.all.binary, count_index) {
    # make logical
    true_binary = sapply(log2FC.all.binary, as.logical)
 
    melted_bin = melt(true_binary) %>%
        filter(as.logical(value)) %>% #remove not DE
        transmute(treatment = str_replace(Var2, '_new', ''), #make group uniting new and old
                  treatment = factor(treatment, levels=unique(treatment)),
                  new = str_detect(Var2, 'new'), #new if new
                  number = as.character(count_index[treatment]),
                  temp = value) %>%
        group_by(treatment, new, number) %>%
        summarize(count = sum(temp)) %>%
        mutate(count = log10(count + .1))

    #remove new genes from DE genes for correct plotting 
    for(i in melted_bin$treatment[melted_bin$new]){
        melted_bin$count[melted_bin$treatment == i & !melted_bin$new] =
            melted_bin$count[melted_bin$treatment == i & !melted_bin$new] -
            melted_bin$count[melted_bin$treatment == i & melted_bin$new]
    }

    p = ggplot(melted_bin, aes(treatment, count)) +
        geom_col(aes(alpha=new, fill=number)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_fill_manual(values=c('#a6cee3', '#1f78b4', '#b2df8a')) +
        scale_alpha_manual(values = c(.4, 1)) +
        theme(panel.background = element_rect(fill = "white", colour = "black")) +
        ylab('log10(count)')

    return(p)
}

plot_combo_ratios = function(log2FC.all.binary, count_index) {
    # make logical
    true_binary = sapply(log2FC.all.binary, as.logical)
    
    melted_bin = melt(true_binary) %>%
        filter(as.logical(value)) %>% #remove not DE
        transmute(treatment = str_replace(Var2, '_new', ''), #make group uniting new and old
                  new = str_detect(Var2, 'new'), #new if new
                  number = as.character(count_index[treatment]),
                  temp = value) %>%
        group_by(treatment, new, number) %>%
        summarize(count = sum(temp)) %>%
        filter(number != 1) %>%
        spread(key=new, val=count) %>%
        mutate(ratio = `TRUE` / `FALSE`) %>%
        arrange(number, desc(ratio)) %>%
        ungroup() %>%
        mutate(treatment = factor(treatment, levels = unique(treatment)))
    
    p = ggplot(melted_bin, aes(treatment, ratio)) +
        geom_col(aes(fill=number)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_fill_manual(values=c('#1f78b4', '#b2df8a')) +
        theme(panel.background = element_rect(fill = "white", colour = "black")) +
        ylab('New genes / DE genes')
    
    return(list(p=p, order=levels(melted_bin$treatment), ratios = melted_bin))
}

