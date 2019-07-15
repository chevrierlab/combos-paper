### libraries ###
library(edgeR)
library(limma)
library(tidyverse)
library(gridExtra)
library(data.table)
library(cowplot)

## sources ###
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

#will require at least cpm_threshold in at least cpm_needed samples
cpm_threshold = 50; cpm_needed = 2
padj.Threshold = .01; log2FC.Threshold = log2(1.5); correction.method = 'BH'
p.compare.Threshold = .1; lfc.compare.Threshold = 0

#load data:
SP7_filename = 'data/rnaseq/SP7/umi/Run_2_3/dge/SP7.unq.refseq.umi.dat'
SP7_data = as.data.frame(fread(SP7_filename))
genes = SP7_data$V1
SP7_data$V1 = NULL
SP7_data = as.matrix(SP7_data)
rownames(SP7_data) = genes

SP7_experimental.groups = factor(str_sub(colnames(SP7_data), 1, -7))

#run base DE analysis
SP7_diff = rna_de_analysis(SP7_data, SP7_experimental.groups, cpm_threshold = cpm_threshold, cpm_needed = cpm_needed, 
                           padj.Threshold = padj.Threshold, log2FC.Threshold = log2FC.Threshold)

#make comparisons
singles = c('PAM', 'LPS', 'HMW', 'Zymosan', 'SeV', 'cGAMP', 'CpGB')
single_screen = str_c(singles, collapse = '|')
combos = produce_combos(singles, size_lim = 3, str_break='_')
SP7_diff$binary = make_comparison_binaries(combos, SP7_diff$binary, SP7_diff$fit, SP7_diff$design, 
                                           p.compare.Threshold, '_')
write.csv(SP7_diff$log2FC, 'data/rnaseq/SP7/log2FC.csv', quote=F)
write.csv(cpm(SP7_diff$dge)[rownames(SP7_diff$dge) %in% rownames(SP7_diff$log2FC),], 'data/rnaseq/SP7/cpm.csv', quote=F)


#make count_index and plot
count_index = str_count(str_replace(colnames(SP7_diff$binary), '_new',''), '_') + 1
names(count_index) = str_replace(colnames(SP7_diff$binary), '_new','')

ratio_plot_order = plot_combo_ratios(SP7_diff$binary, count_index)
save_plot('plots/paper/3c.pdf', ratio_plot_order$p,
          base_width = 3, base_height = 2, useDingbats = F)


##### PCA #####

#get log2_cpm data, and include sample and groups
log2_cpm_data = as_tibble(t(cpm(SP7_diff$dge, log = TRUE, prior.count = 3))) %>%
    mutate(sample = colnames(SP7_diff$dge),
           group = SP7_experimental.groups) %>%
    group_by(group)
#calculate mean of each treatment group
SP7_mean_log2_cpm = summarize_at(log2_cpm_data, vars(-sample), mean)

#pca (run and make matrix of results)
pca_proc = prcomp(as.matrix(t(SP7_diff$all_log2FC)), scale=TRUE,center=TRUE)
pca_res = data.frame(group = colnames(SP7_diff$all_log2FC), pca_proc$x) %>% 
    mutate(number = as.character(str_count(group, '_') + 1))

#plot pca
p = ggplot(pca_res, aes(PC1, PC2)) + 
    geom_point(aes(fill=number), shape=21, color='black') + 
    scale_fill_manual(values = c('#a6cee3', '#1f78b4', '#b2df8a'),
                      guide = guide_legend(title=NULL),
                      labels = c('singles (6)', 'pair (15)', 'triplet (20)')) +
    xlab(str_c('PC1: ', summary(pca_proc)$importance[2,1] * 100, '% variance')) +
    ylab(str_c('PC2: ', summary(pca_proc)$importance[2,2] * 100, '% variance'))
save_plot('plots/paper/3b.pdf', p,  useDingbats = F, base_width = 3, base_height = 2.5)


#### histogram materials ####
# pull out 'new' portion of binary matrix
just_news = as.matrix(SP7_diff$binary[,str_detect(colnames(SP7_diff$binary), 'new')])
just_news = array(as.logical(just_news), dim(just_news), dimnames(just_news))

# make dataframe indicating for which rows things are new,
SP7_annotation = tibble(gene = rownames(SP7_diff$binary), new = rowSums(just_news) > 0) %>%
    mutate(new_from = map_chr(gene, function(g) 
        # which treatments they're new in,
        if (sum(just_news[g,]) > 0) str_c(colnames(just_news)[just_news[g,]], collapse=',')
        else ''),
        # if they're new in a pair,
        new_from_pair = map_lgl(str_split(new_from, ','), function (s) any(str_count(s,'_') == 2)),
        # if they're new in a triplet
        new_from_triple = map_lgl(str_split(new_from, ','), function (s) any(str_count(s,'_') == 3)))

# isolate non-new portions of binary matrix (and make properly binary)
proper_bin = abs(SP7_diff$binary)[,!str_detect(colnames(SP7_diff$binary), 'new')]
# make note of genes which are DE in the ATAC-seq treatments (PLS, PSG, ZSG)
for(treatment in c('PAM_LPS_HMW', 'PAM_SeV_cGAMP', 'Zymosan_SeV_cGAMP')) {
    splits = str_split(treatment, '_')[[1]]
    pairs = c()
    if (length(splits) == 3)
        pairs = c(str_c(c(splits[1], splits[1], splits[2]), '_', c(splits[2], splits[3], splits[3])))
    
    components = c(splits, pairs, treatment)
    SP7_annotation = mutate(SP7_annotation, 
                            !!treatment := as.numeric(rowSums(proper_bin[, components]) > 0))
}

atac_trips = "PAM_LPS_HMW|PAM_SeV_cGAMP|Zymosan_SeV_cGAMP"
atac_pairs = "PAM_LPS|LPS_HMW|PAM_HMW|PAM_SeV|SeV_cGAMP|PAM_cGAMP|Zymosan_SeV|Zymosan_cGAMP"


SP7_annotation = SP7_annotation %>% 
    # note which genes are DE in any atac treatment
    mutate(atac = PAM_LPS_HMW | PAM_SeV_cGAMP | Zymosan_SeV_cGAMP,
           # which are new in an atac pair
           atac_new_pair = str_detect(new_from, atac_pairs),
           # which are new in an atac triplet
           atac_new_trip = str_detect(new_from, atac_trips),
           # which are new in an atac treatment
           atac_new = atac_new_pair | atac_new_trip)

write.table(SP7_annotation, 'data/rnaseq/SP7/row_annotation.txt', quote=F, row.names = F, sep='\t')


### dives into amplitudes and newness ###

### amplitudes ###
# split treatment names into singles via regular expressions
treatment_splits = str_split(colnames(SP7_diff$log2FC), '_')
# fill in missing singles with NA (so they all have three)
treatment_splits = lapply(treatment_splits, function(s) {if (length(s) == 1) c(s, NA, NA)
    else if (length(s) == 2) c(s, NA)
    else s })
# convert to matrix and name
treatment_splits = matrix(unlist(treatment_splits), nrow=3)
colnames(treatment_splits) = colnames(SP7_diff$log2FC)

# make a key mapping treatment names to the names of their component combinations
component_key = tibble(Treatment = colnames(SP7_diff$log2FC)) %>%
    mutate(C_1 = treatment_splits[1,],
           C_2 = treatment_splits[2,],
           C_3 = treatment_splits[3,],
           C_12 = str_c(treatment_splits[1,], treatment_splits[2,], sep = '_'),
           C_13 = str_c(treatment_splits[1,], treatment_splits[3,], sep = '_'),
           C_23 = str_c(treatment_splits[2,], treatment_splits[3,], sep = '_'),
           C_123 = str_c(treatment_splits[1,], treatment_splits[2,], treatment_splits[3,], sep = '_'))
component_key = mutate_all(component_key, function(x) ifelse(is.na(x), component_key$Treatment, x))

# gather LFC components
pulled_lfc = data.frame()
for(treat in str_subset(colnames(SP7_diff$log2FC), '_')){
    res = SP7_diff$log2FC[,unlist(filter(component_key, Treatment == treat)[,2:8])]
    colnames(res) = str_c('C_', c(1,2,3,12,13,23,123))
    res = rownames_to_column(res, 'gene')
    res$Treatment = treat
    
    pulled_lfc = bind_rows(pulled_lfc, res)
}                       


pulled_lfc = pulled_lfc %>% 
    mutate(type = str_count(Treatment, '_') + 1,
           # calculate which component is the closest to the treatment, and record the difference
           closest = ifelse(type==3,
                            apply(cbind(C_123 - C_1, C_123 - C_2, C_123 - C_3,
                                        C_123 - C_12, C_123 - C_13, C_123 - C_23), 1, abs_min),
                            apply(cbind(C_12 - C_1, C_12 - C_2), 1, abs_min)),
           # is this gene DE in the treatment of interest
           DE = SP7_diff$binary[cbind(gene, Treatment)],
           type = c('single', 'pair', 'triplet')[type],
           Treatment = factor(Treatment, levels = ratio_plot_order$order)) %>% 
    # only interested in DE genes
    filter(DE == 1)

p = ggplot(pulled_lfc, aes(Treatment, closest, fill=type)) + 
    geom_boxplot(outlier.shape = NA) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
    xlab('Treatment') + 
    ylab('Smallest Divergence') + 
    ylim(c(-1, 2)) + #magic numbers, woo!
    geom_hline(yintercept = 0, linetype='dashed')
save_plot('plots/paper/3d.pdf', p)

#### MA plots ####

# note which genes are new in which triplets
just_new_trips = SP7_diff$binary[, str_count(colnames(SP7_diff$binary), '_') == 3]
colnames(just_new_trips) = str_replace(colnames(just_new_trips), '_new', '')
just_new_trips = rownames_to_column(just_new_trips, 'gene') %>% 
    gather(-gene, key='triplet', value = 'new') %>%
    filter(new == 1) 

# note which genes are new in which pairs
just_new_pairs = SP7_diff$binary[, str_detect(colnames(SP7_diff$binary), str_c('^(', single_screen, ')_(', single_screen, ')_new'))]
colnames(just_new_pairs) = str_replace(colnames(just_new_pairs), '_new', '')
just_new_pairs = rownames_to_column(just_new_pairs, 'gene') %>% 
    gather(-gene, key='pair', value = 'new') %>%
    filter(new == 1) 

# for each triplet with a new gene, make an MA plot
MA_plots = list()
for(trip in unique(just_new_trips$triplet)){
    DE = as_tibble(SP7_diff$binary, rownames = NA)[rownames(SP7_diff$contrast_fit), trip] != 0
    DE[is.na(DE)] = F
    
    plot_mat = 
        tibble(mean_exp = SP7_diff$contrast_fit$Amean, 
               lfc = SP7_diff$contrast_fit$coefficients[,str_c(trip, '-Ctrl')],
               status=factor(ifelse(rownames(SP7_diff$contrast_fit) %in% filter(just_new_trips, triplet == trip)$gene,
                                    'new',
                                    ifelse(DE, 'DE', 'Not DE')), c('Not DE', 'DE', 'new')))
    
    MA_plots[[trip]] = ggplot(plot_mat, aes(mean_exp, lfc, col = status, size = status)) + 
        ggrastr::geom_point_rast() + 
        geom_hline(yintercept = 0) + 
        scale_color_manual(values = c('grey', 'black', 'red')) + 
        scale_size_manual(values = c(.3, .3, 3)) + 
        ggtitle(trip)

}

# for each pair with a new gene, make an MA plot
for(p in unique(just_new_pairs$pair)){
    DE = as_tibble(SP7_diff$binary, rownames = NA)[rownames(SP7_diff$contrast_fit), p] != 0
    DE[is.na(DE)] = F
    
    plot_mat = 
        tibble(mean_exp = SP7_diff$contrast_fit$Amean, 
               lfc = SP7_diff$contrast_fit$coefficients[,str_c(p, '-Ctrl')],
               status=factor(ifelse(rownames(SP7_diff$contrast_fit) %in% filter(just_new_pairs, pair == p)$gene, 
                                           'new',
                                           ifelse(DE, 'DE', 'Not DE')), c('Not DE', 'DE', 'new')))
    
    MA_plots[[p]] = ggplot(plot_mat, aes(mean_exp, lfc, col = status, size = status)) + 
        ggrastr::geom_point_rast() + 
        geom_hline(yintercept = 0) + 
        scale_color_manual(values = c('grey', 'black', 'red')) + 
        scale_size_manual(values = c(.3, .3, 3)) + 
        ggtitle(p)
}

MA_plots = MA_plots[ratio_plot_order$order[ratio_plot_order$order %in% names(MA_plots)]]
save_plot('plots/paper/supplemental/10.pdf', plot_grid(plotlist = MA_plots, ncol = 4), nrow=6, ncol=4)

# makes a table of stats for how many new and proportion of new to DE in the MA plots
MA_proportions = bind_rows(
    map_dfr(unique(just_new_trips$triplet), 
            function(treat) {
                DE = as_tibble(SP7_diff$binary, rownames = NA)[rownames(SP7_diff$contrast_fit), treat] != 0
                DE[is.na(DE)] = F
                status=factor(ifelse(rownames(SP7_diff$contrast_fit) %in% filter(just_new_trips, triplet == treat)$gene,
                                     'new',
                                     ifelse(DE, 'DE', 'Not DE')), c('Not DE', 'DE', 'new'))
                tibble(Treatment = treat, total_new = sum(status == 'new'), 
                       new_to_DE = sum(status == 'new') / (sum(status == 'DE') + sum(status == 'new')))
            }
    ),
    map_dfr(unique(just_new_pairs$pair), 
            function(treat) {
                DE = as_tibble(SP7_diff$binary, rownames = NA)[rownames(SP7_diff$contrast_fit), treat] != 0
                DE[is.na(DE)] = F
                status=factor(ifelse(rownames(SP7_diff$contrast_fit) %in% filter(just_new_pairs, pair == treat)$gene,
                                     'new',
                                     ifelse(DE, 'DE', 'Not DE')), c('Not DE', 'DE', 'new'))
                tibble(Treatment = treat, total_new = sum(status == 'new'), 
                       new_to_DE = sum(status == 'new') / (sum(status == 'DE') + sum(status == 'new')))
            }
    )
)

#### New gene dotplots ####

# prep cpm data for plotting
SP7_log2_cpm_sum = as.data.frame(summarize_at(log2_cpm_data, vars(-sample), list(~mean, ~se)))
rownames(SP7_log2_cpm_sum) = SP7_log2_cpm_sum$group
SP7_log2_cpm_sum = select(SP7_log2_cpm_sum, -group)

splits = str_split(just_new_trips$triplet, '_')

# collect all the data
cpm_for_lines = just_new_trips %>%
    mutate(gene_mean = str_c(gene, '_mean'),
           gene_se = str_c(gene, '_se'),
           ctrl_cpm = SP7_log2_cpm_sum[cbind("Ctrl", gene_mean)],
           ctrl_win = SP7_log2_cpm_sum[cbind("Ctrl", gene_se)],
           s1_cpm = SP7_log2_cpm_sum[cbind(map_chr(splits, function (s) s[[1]]), gene_mean)],
           s1_win = SP7_log2_cpm_sum[cbind(map_chr(splits, function (s) s[[1]]), gene_se)],
           s2_cpm = SP7_log2_cpm_sum[cbind(map_chr(splits, function (s) s[[2]]), gene_mean)],
           s2_win = SP7_log2_cpm_sum[cbind(map_chr(splits, function (s) s[[2]]), gene_se)],
           s3_cpm = SP7_log2_cpm_sum[cbind(map_chr(splits, function (s) s[[3]]), gene_mean)],
           s3_win = SP7_log2_cpm_sum[cbind(map_chr(splits, function (s) s[[3]]), gene_se)],
           d12_cpm = SP7_log2_cpm_sum[cbind(map_chr(splits, function (s) str_c(s[[1]], '_', s[[2]])), gene_mean)],
           d12_win = SP7_log2_cpm_sum[cbind(map_chr(splits, function (s) str_c(s[[1]], '_', s[[2]])), gene_se)],
           d13_cpm = SP7_log2_cpm_sum[cbind(map_chr(splits, function (s) str_c(s[[1]], '_', s[[3]])), gene_mean)],
           d13_win = SP7_log2_cpm_sum[cbind(map_chr(splits, function (s) str_c(s[[1]], '_', s[[3]])), gene_se)],
           d23_cpm = SP7_log2_cpm_sum[cbind(map_chr(splits, function (s) str_c(s[[2]], '_', s[[3]])), gene_mean)],
           d23_win = SP7_log2_cpm_sum[cbind(map_chr(splits, function (s) str_c(s[[2]], '_', s[[3]])), gene_se)],
           t_cpm = SP7_log2_cpm_sum[cbind(map_chr(splits, function (s) str_c(s[[1]], '_', s[[2]], '_', s[[3]])), gene_mean)],
           t_win = SP7_log2_cpm_sum[cbind(map_chr(splits, function (s) str_c(s[[1]], '_', s[[2]], '_', s[[3]])), gene_se)]) %>%
    select(-new, -gene_mean, -gene_se) %>% 
    gather(-gene, -triplet, key='pair_type', value = 'val') %>%
    extract(pair_type, c('pair', 'type'), '(.*)_(.*)') %>%
    spread(type, val) %>%
    mutate(err_top = cpm + win, err_bot = cpm - win,
           pair = factor(pair, levels = c('ctrl', 's1', 's2', 's3', 'd12', 'd13', 'd23', 't')),
           triplet = factor(triplet, levels = unique(triplet)),
           color = str_sub(pair, 1, 1)) 

p = ggplot(cpm_for_lines, aes(pair, cpm, col = color)) + 
    geom_errorbar(aes(ymax = err_top, ymin = err_bot), width = .15) + 
    geom_point(size = .5) +  
    facet_wrap(~interaction(triplet, gene, sep = ', ', lex.order = T), scales = 'free', ncol=5) + 
    theme(text = element_text(size = 6),
          strip.text.x = element_text(margin = margin(.1, 0, .1, 0, 'cm'))) + 
    guides(color = NULL) + 
    ylab('log2 cpm')

save_plot('plots/paper/supplemental/11.pdf', p, 
          nrow=12, ncol=4, base_aspect_ratio = 2, limitsize=F, base_height = NULL, base_width = 2)
