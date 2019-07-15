### libraries ###
library(edgeR)
library(limma)
library(tidyverse)
library(data.table)
library(cowplot)

### sources ###
source('software/combo_analysis.r')
source('software/rna_de.r')
source('software/index_functions.r')

### ggtheme ###
theme_set(theme_classic())
theme_update(text = element_text(size = rel(2.25)),
             legend.key.size = unit(.1, 'inches'),
             legend.text = element_text(size = rel(1.5)),
             panel.background = element_rect(fill = "white", colour = "black"),
             axis.line = element_blank())

singles = c('P', 'L', 'H', 'S', 'cG', 'Z', 'Cp')
single_screen = str_c(singles, collapse='|')

#will require at least cpm_threshold in at least cpm_needed samples
cpm_threshold = 50; cpm_needed = 2
padj.Threshold = .01; log2FC.Threshold = log2(1.2); correction.method = 'BH'
p.compare.Threshold = .1; lfc.compare.Threshold = 0

data.filename = 'data/rnaseq/Human/umi/dge/Human.unq.refseq.umi.dat'

#read and prep data
human_data = data.frame(fread(data.filename))
rownames(human_data) = make.names(human_data$V1, unique=T)
human_data$V1 = NULL

#interpret sample names
splits = str_split(colnames(human_data), '_')
treatment = factor(initials_to_names(sapply(splits, function(s) s[[1]])))
indiv = sapply(splits, function(s) s[[2]])

#do DE analysis
human_separate_diffs = map(unique(indiv), function (h) 
    rna_de_analysis(human_data[, indiv==h], treatment[indiv==h], cpm_threshold = cpm_threshold, 
                    cpm_needed = cpm_needed, padj.Threshold = padj.Threshold, 
                    log2FC.Threshold = log2FC.Threshold))
names(human_separate_diffs) = unique(indiv)

#do comparison new-gene filtering
singles = c('LPS', 'HMW', 'Zymosan', 'SeV', 'cGAMP')
combos = produce_combos(singles, size_lim = 3, str_break='_')
for(i in 1:length(human_separate_diffs)){
    human_separate_diffs[[i]]$binary = 
        make_comparison_binaries(combos, human_separate_diffs[[i]]$binary, 
                        human_separate_diffs[[i]]$fit, human_separate_diffs[[i]]$design, 
                        p.compare.Threshold, '_')
}

# make a matrix of lfc combining all the logFC matrices
# has NA for genes that aren't DE for a given individual
i_rep = unlist(map(names(human_separate_diffs), rep, ncol(human_separate_diffs[[1]]$log2FC)))
human_de_genes = unique(c(unlist(map(human_separate_diffs, function(h) rownames(h$log2FC)))))
unified_lfc = matrix(character(), nrow=length(human_de_genes) + 2, ncol=length(i_rep))
rownames(unified_lfc) = c('Treatment', 'individual', human_de_genes)
unified_lfc['Treatment',] = rep(colnames(human_separate_diffs[[1]]$log2FC), 3)
unified_lfc['individual',] = i_rep
for(i in 3:nrow(unified_lfc)) { 
    for(j in 1:ncol(unified_lfc)) {
        if(rownames(unified_lfc)[i] %in% 
           rownames(human_separate_diffs[[unified_lfc['individual', j]]]$log2FC)) {
            unified_lfc[i, j] = 
                human_separate_diffs[[unified_lfc['individual', j]]]$log2FC[rownames(unified_lfc)[i], unified_lfc['Treatment', j]]
        }
    }    
}

write.csv(unified_lfc, 'data/rnaseq/Human/lfc.csv', quote=F)
write.csv(human_separate_diffs[[1]]$log2FC, 'data/rnaseq/Human/27_lfc.csv', quote = F)
write.csv(human_separate_diffs[[2]]$log2FC, 'data/rnaseq/Human/39_lfc.csv', quote = F)
write.csv(human_separate_diffs[[3]]$log2FC, 'data/rnaseq/Human/45_lfc.csv', quote = F)

# This is all a slightly odd way to make a matrix that has lfc all the (non-filtered) genes for the humans
genes = map(unique(indiv), function (h) {
        dge <- DGEList(counts = human_data[, indiv==h], group = treatment[indiv==h])
        # edgeR normalization
        dge.norm <- calcNormFactors(dge)
        # filter out genes that do not have >50 reads in at least in 2 samples
        isexpr <- rowSums(cpm(dge.norm) > cpm_threshold) >= cpm_needed
        return(rownames(dge.norm)[isexpr])
    }
)
genes = unique(unlist(genes))

alt_human_separate_diffs = map(unique(indiv), function (h) 
    alt_rna_de_analysis(human_data[genes, indiv==h], treatment[indiv==h], cpm_threshold = cpm_threshold, 
                    cpm_needed = cpm_needed, padj.Threshold = padj.Threshold, 
                    log2FC.Threshold = log2FC.Threshold))
names(alt_human_separate_diffs) = unique(indiv)
colnames(alt_human_separate_diffs[[1]]$log2FC) = str_c(colnames(alt_human_separate_diffs[[1]]$log2FC), '_27')
colnames(alt_human_separate_diffs[[2]]$log2FC) = str_c(colnames(alt_human_separate_diffs[[2]]$log2FC), '_39')
colnames(alt_human_separate_diffs[[3]]$log2FC) = str_c(colnames(alt_human_separate_diffs[[3]]$log2FC), '_45')
whole_lfc = cbind(alt_human_separate_diffs[[1]]$log2FC,
                      alt_human_separate_diffs[[2]]$log2FC,
                      alt_human_separate_diffs[[3]]$log2FC)
whole_lfc = rbind(str_sub(colnames(whole_lfc), -2, -1), str_sub(colnames(whole_lfc), 1, -4), whole_lfc)
rownames(whole_lfc)[1:2] = c('id', 'treatment')
write.csv(whole_lfc, 'data/rnaseq/Human/whole_lfc.csv', quote = F)

# row annotations for human genes
human_annotation = tibble(gene = human_de_genes) 
for(treat in unique(str_subset(unified_lfc['Treatment',], '_.*_'))) {
    splits = str_split(treat, '_')[[1]]
    pairs = c()
    if (length(splits) == 3)
        pairs = c(str_c(c(splits[1], splits[1], splits[2]), '_', c(splits[2], splits[3], splits[3])))
    
    components = c(splits, pairs, treat)
    for(i in names(human_separate_diffs)) {
        slice = unified_lfc[-1:-2,unified_lfc['Treatment',] %in% components & unified_lfc['individual',] == i]
        slice = array(as.numeric(slice), dim = dim(slice), dimnames = dimnames(slice))
        bin_slice = slice > log2FC.Threshold | slice < -log2FC.Threshold
        name = str_c(i, '_', treat)
        
        human_annotation = mutate(human_annotation, 
                            !!str_c(i, '_', treat) := as.numeric(rowSums(bin_slice) > 0))
    }
    human_annotation = mutate(human_annotation,
                              !!str_c('any_', treat) := as.numeric(!!as.name(str_c('27_', treat)) | 
                                    !!as.name(str_c('39_', treat)) | 
                                    !!as.name(str_c('45_', treat))))
}
write.table(human_annotation, 'data/rnaseq/Human/row_annot.txt', quote=F, row.names = F, sep='\t')

# makes a matrix denoting which genes are new for which treatments/individuals
whole_bin = map_dfr(human_separate_diffs, function(diff){
    bin = diff$binary
    bin = bin[,str_detect(colnames(bin), 'new')]
    colnames(bin) = str_replace(colnames(bin), '_new', '')
    
    bin = rownames_to_column(bin, 'id') %>%
        gather(-id, key='Treatment', value='new')
}, .id='indiv') %>%
    mutate(sample = str_c(Treatment, '_', indiv)) %>%
    select(-Treatment, -indiv) %>%
    spread(sample, new, fill=0)
write_tsv(whole_bin, 'data/rnaseq/Human/whole_bin.tsv')

# new gene bar plot
count_index = str_count(str_replace(colnames(human_separate_diffs$`27`$binary), '_new',''), '_') + 1
names(count_index) = str_replace(colnames(human_separate_diffs$`27`$binary), '_new','')
new_genes = map(human_separate_diffs, function(diff) plot_combo_ratios(diff$binary, count_index))

all_ratios = bind_rows(map(new_genes, ~.$ratios), .id = 'indiv')
p = ggplot(all_ratios, aes(treatment, ratio, fill = interaction(number, indiv))) + 
    geom_col(col = 'black', position = 'dodge', width=.75) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
    scale_fill_manual(values = c("#073D00", "#440000",
                      "#3D7F35", "#8E4B4B", 
                      "#B0E5A9", "#D89797"))
save_plot('plots/paper/3e.pdf', p, base_aspect_ratio = 2)


# calculate amplitude plot data for each human
orders = map(new_genes, ~.$order)
human_amplitudes = map2(
    human_separate_diffs, orders, 
    function(diff, order) {
        treatment_splits = str_split(colnames(diff$log2FC), '_')
        # fill in missing singles with NA (so they all have three)
        treatment_splits = lapply(treatment_splits, function(s) {if (length(s) == 1) c(s, NA, NA)
            else if (length(s) == 2) c(s, NA)
            else s })
        # convert to matrix and name
        treatment_splits = matrix(unlist(treatment_splits), nrow=3)
        colnames(treatment_splits) = colnames(diff$log2FC)
        
        # make a key mapping treatment names to the names of their component combinations
        component_key = tibble(Treatment = colnames(diff$log2FC)) %>%
            mutate(C_1 = treatment_splits[1,],
                   C_2 = treatment_splits[2,],
                   C_3 = treatment_splits[3,],
                   C_12 = str_c(treatment_splits[1,], treatment_splits[2,], sep = '_'),
                   C_13 = str_c(treatment_splits[1,], treatment_splits[3,], sep = '_'),
                   C_23 = str_c(treatment_splits[2,], treatment_splits[3,], sep = '_'),
                   C_123 = str_c(treatment_splits[1,], treatment_splits[2,], treatment_splits[3,], sep = '_'))
        component_key = mutate_all(component_key, function(x) ifelse(is.na(x), component_key$Treatment, x))
        
        # gather components
        pulled_lfc = data.frame()
        for(treat in str_subset(colnames(diff$log2FC), '_')){
            res = diff$log2FC[,unlist(filter(component_key, Treatment == treat)[,2:8])]
            colnames(res) = str_c('C_', c(1,2,3,12,13,23,123))
            res = rownames_to_column(res, 'gene')
            res$Treatment = treat
            
            pulled_lfc = bind_rows(pulled_lfc, res)
        }                       
        
        # calculate data to closest
        pulled_lfc = pulled_lfc %>% 
            mutate(type = str_count(Treatment, '_') + 1,
                   closest = ifelse(type==3,
                                    apply(cbind(C_123 - C_1, C_123 - C_2, C_123 - C_3,
                                                C_123 - C_12, C_123 - C_13, C_123 - C_23), 1, abs_min),
                                    apply(cbind(C_12 - C_1, C_12 - C_2), 1, abs_min)),
                   DE = diff$binary[cbind(gene, Treatment)],
                   type = c('single', 'pair', 'triplet')[type],
                   Treatment = factor(Treatment, levels = order)) %>% 
            filter(DE == 1)
        
        return(pulled_lfc)
    }
)

# plot closest
human_amplitudes = bind_rows(human_amplitudes, .id = 'indiv')
p = ggplot(human_amplitudes, aes(Treatment, closest, fill = interaction(type, indiv))) + 
    geom_boxplot(outlier.shape = NA) + 
    ylim(c(-2.5, 2.5)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
    xlab('Treatment') + 
    ylab('Smallest Divergence') + 
    scale_fill_manual(values = c("#073D00", "#440000",
                                 "#3D7F35", "#8E4B4B", 
                                 "#B0E5A9", "#D89797")) + 
    geom_hline(yintercept = 0, linetype='dashed')
save_plot('plots/paper/3f.pdf', p, base_aspect_ratio = 2)
    

#pca separately by group (run and make matrix of results)
pca_proc = prcomp(t(data.matrix(whole_lfc[-1:-2,])), scale=TRUE,center=TRUE)
pca_res = data.frame(treatment = unlist(whole_lfc['treatment',]), 
                     indiv = unlist(whole_lfc['id',]), pca_proc$x) %>% 
    mutate(number = factor(str_count(treatment, '_') + 1))

#plot pca separately 
p = ggplot(pca_res, aes(PC1, PC2)) + 
    geom_point(aes(fill=number), shape=21, color='black') + 
    geom_text(aes(label=treatment), hjust=-.1, vjust=-.1, size=1.5) +
    scale_fill_manual(values = c('#a6cee3', '#1f78b4', '#b2df8a'),
                      guide = guide_legend(title=NULL),
                      labels = c('singles (5)', 'pair (10)', 'triplet (10)')) +
    xlab(str_c('PC1: ', summary(pca_proc)$importance[2,1] * 100, '% variance')) +
    ylab(str_c('PC2: ', summary(pca_proc)$importance[2,2] * 100, '% variance')) +
    facet_grid(~ indiv)
save_plot('plots/paper/supplemental/12abc.pdf', p, useDingbats = F, base_height=4, base_width=12)

#### correlations between individuals ####

# sort of overkill for just getting cpm but here we are
human_together_diff = rna_de_analysis(human_data, treatment, cpm_threshold = cpm_threshold, cpm_needed = cpm_needed, 
                                      padj.Threshold = padj.Threshold, log2FC.Threshold = log2FC.Threshold)

# make matrix of unified cpm
whole_cpm = as_tibble(t(cpm(human_together_diff$dge))) %>%
    mutate(sample = colnames(human_together_diff$dge),
           treat_indiv = map_chr(str_split(sample, '_'), ~str_c(.[1], '_', .[2]))) %>%
    group_by(treat_indiv) %>%
    select(-sample) %>%
    summarize_all(mean) 
rownames(whole_cpm) = whole_cpm$treat_indiv
whole_cpm$treat_indiv = NULL
whole_cpm = t(as.matrix(whole_cpm))

# calculate correlation matrix of cpm
p = cor(whole_cpm) %>%
    as.data.frame() %>%
    # gather into melted format
    rownames_to_column('first') %>%
    gather(-first, key='second', value = 'cor') %>%
    extract(first, c('treatment1', 'indiv1'), '(.*)_(.*)') %>%
    extract(second, c('treatment2', 'indiv2'), '(.*)_(.*)') %>%
    # remove duplicate points and self-correlations
    filter(treatment1 == treatment2, indiv1 < indiv2) %>%
    select(treatment = treatment1, indiv1, indiv2, cor) %>%
    # square correlations
    mutate(number = str_count(treatment, single_screen),
           r2 = cor ^ 2,
           color = interaction(number, indiv1, indiv2)) %>%
    arrange(number) %>%
    mutate(treatment = factor(treatment, levels = unique(treatment))) %>%
    # plot
    ggplot(aes(treatment, r2, fill = color)) + 
        geom_point(color = 'black', shape = 21) + 
        scale_fill_manual(values = c('darkgrey', "#003182", "#073D00", "#440000",
                                     'grey', "#4D79C0","#3D7F35", "#8E4B4B", 
                                     'white', "#9BC1FF", "#B0E5A9", "#D89797"),
                      guide = guide_legend(title=NULL)) + 
    ylim(0, 1)

save_plot('plots/paper/supplemental/12d.pdf', p, base_aspect_ratio = 1.4)
