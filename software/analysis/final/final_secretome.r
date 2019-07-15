### libraries ###
library(edgeR)
library(limma)
library(tidyverse)
library(data.table)
library(cowplot)
library(sva)

### sources ###
source('software/combo_analysis.r')
source('software/rna_de.r')
source('software/index_functions.r')

### dependencies ###
source('software/final/final_atac.r')


# single screen is used for pulling apart treatment names using regular expressions
singles = c('Z','S','Cp','P','cG','L', 'H')
single_screen = str_c(singles, collapse='|')

# read in mass-spec data and remove non-needed statistics
all_raw = readxl::read_excel('data/cytokines/One-sample_mod_T_Shiny.xlsx') %>% 
    gather(-id, key = 'key', value = 'val') %>%
    filter(str_detect(key, '^[^\\.]*$')) %>%
    extract(key, c('treatment', 'replicate', 'experiment'), '(.*)_vs_Ctr(.)_(.*)', remove=F) %>%
    rename(sample = key) %>%
    mutate(val = as.numeric(val)) %>%
    group_by(id) %>%
    filter(!is.na(max(val)))

# make matrix
raw_mat = all_raw %>%
    select(id, sample, val) %>%
    spread(sample, val)
rownames(raw_mat) = raw_mat$id
raw_mat$id = NULL
raw_mat = as.matrix(raw_mat) 
raw_mat = raw_mat[rowSums(is.na(raw_mat)) == 0,]

# use combat to correct for batch effect between experiments
batch = str_extract(colnames(raw_mat), '_[^_]*$')
corrected_raw = ComBat(raw_mat, batch)

# average across replicates
full_avg_corrected = corrected_raw %>% 
    as.data.frame() %>%
    rownames_to_column('gene') %>%
    gather(-gene, key='sample', value='fc') %>%
    extract(sample, c('treatment', 'replicate', 'experiment'), '(.*)_vs_Ctr(.)_(.*)', remove=F) %>%
    group_by(treatment, gene) %>%
    summarize(fc = mean(fc)) %>%
    spread(treatment, fc)
rownames(full_avg_corrected) = full_avg_corrected$gene
full_avg_corrected$gene = NULL
full_avg_corrected = as.matrix(full_avg_corrected)
write.csv(full_avg_corrected, 'data/cytokines/secretome_corrected_lfc.csv', quote = F)

# pca on corrected, averaged data
pca_proc = prcomp(t(full_avg_corrected), scale=T, center=T)
pca_res = data.frame(sample = rownames(pca_proc$x), pca_proc$x) %>%
    mutate(number = factor(str_count(sample, single_screen)))

# plot PCA 
p = ggplot(pca_res, aes(PC1, PC2)) + 
    geom_point(aes(fill=number), shape=21, color='black') + 
    geom_text(aes(label=sample), hjust=-.1, vjust=-.1, size=1.5) +
    scale_fill_manual(values = c('#a6cee3', '#1f78b4', '#b2df8a'),
                      guide = 'none',
                      labels = c('singles', 'pair', 'triplet')) +
    xlab(str_c('PC1: ', summary(pca_proc)$importance[2,1] * 100, '% variance')) +
    ylab(str_c('PC2: ', summary(pca_proc)$importance[2,2] * 100, '% variance')) +
    ggtitle('Secretome Mass-Spec')
save_plot('plots/paper/supplemental/13abc.pdf', plot_grid(p, atac_plot, rna_plot, nrow=1, rel_widths = c(2.4, 2.4, 3)),
          base_width = 7, base_height = 2.5, useDingbats = F)

#### comparing to cytokine assay ####

# corresponding gene names to cytokine names
sec_cyto_names = c('Cxcl1', 'Tnf', 'Ccl2',  
                   'Ccl5','Cxcl10', 'Il6')
cytokine_names = c("CXCL1", "TNFa", "CCL2",  
                   "CCL5", "CXCL10", "IL6")

# cytokine concentrations:
cytokines = fread('data/cytokines/6_ligand_antiviral_legendplex_Rtable.txt')
cytokines = mutate(cytokines, Treatment = str_replace(Treatment, '_*.$', '')) %>%
    arrange(str_count(Treatment, single_screen), Treatment)

fc = function(x) x / x[1]

# calculate fold changes
cyto_lfc = cytokines %>%
    group_by(Type, Treatment) %>%
    summarize_all(mean, na.rm=T) %>%
    select(Type, Treatment, cytokine_names) %>%
    ungroup()



# combine secretome and cytokine data
both = full_avg_corrected[sec_cyto_names,] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('Treatment') %>%
    inner_join(cyto_lfc, by='Treatment')

# out_coef * sd of residuals will be used to label outlier dots
out_coef=1

dot_plots = list()
for(i in 1:length(cytokine_names)){
    # ranges for each datatype
    c_min = min(both[cytokine_names[i]])
    c_max = max(both[cytokine_names[i]])
    s_max = max(both[sec_cyto_names[i]])
    s_min = min(both[sec_cyto_names[i]])
    
    # linear model
    mod = summary(lm(as.formula(str_c(sec_cyto_names[i], ' ~ ', cytokine_names[i])), both))
    r2 = mod$r.squared
    intercept = mod$coefficients[1,1]
    slope = mod$coefficients[2,1]
    
    # plot
    c = quo(!!as.name(cytokine_names[i]))
    s = quo(!!as.name(sec_cyto_names[i]))
    p = mutate(both, expected = !!c*coef(mod)[2,1] + coef(mod)[1,1],
               outlier = !!s > expected + mod$sigma*out_coef | 
                   !!s < expected - mod$sigma*out_coef,
               # blank group for non-outliers
               group = ifelse(outlier, Treatment, '')) %>%
        ggplot(aes_(as.name(cytokine_names[i]), as.name(sec_cyto_names[i]))) + 
            # line of fit
            geom_smooth(method='lm', se=F, col='black') +
            # top outlier line
            geom_abline(slope = slope, intercept = intercept + mod$sigma, linetype='dotted') +
            # bottom outlier line
            geom_abline(slope = slope, intercept = intercept - mod$sigma, linetype='dotted') +
            # points
            geom_point(aes(color=Type)) +
            # labels on outliers
            geom_text(aes(label=group), hjust=0, vjust=0, size=2) +
            scale_fill_manual(values = c('#a6cee3', '#1f78b4', '#b2df8a')) +
            labs(title=cytokine_names[i], x='ELISA', y='Mass Spec', color='Type') +
            # label stats
            annotate('text', x=c_min*1.2, y=s_max*.8, label = str_c("italic(R) ^ 2 == ", r2), parse=T, size=2) +
            annotate('text', x=c_max*.9, y=s_min*1.2, label = str_c("y = ", round(slope, 2), "x + ", round(intercept, 2)), size=2)
    dot_plots[[i]] = p
}

save_plot('plots/paper/supplemental/14.pdf', plot_grid(plotlist=dot_plots, ncol=3), 
          ncol = 3, nrow = 2, useDingbats=F)
