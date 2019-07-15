
## Dependencies:
source('software/final/final_SP7.r')

singles = c('P', 'L', 'H', 'S', 'cG', 'Z')
single_screen = str_c(singles, collapse='|')

#run ATAC-seq DA analysis
atac_read_file = 'data/atacseq/020218_SP3/out/combo_PLH_PScG_ZScG/read_table.csv'
atac_reads = read.csv(atac_read_file, row.names = 1)

stims = list(c('P', 'L', 'H'), c('P', 'S', 'cG'), c('Z', 'S', 'cG'))
atac_diff = compare_combos_contrasts(atac_reads, 'Ctrl', stims, 
                                     padj.Threshold, log2FC.Threshold, p.compare.Threshold, lfc.compare.Threshold)
write.csv(rbind(group = initials_to_names(colnames(atac_diff$log2FC)), atac_diff$log2FC),
          'data/atacseq/020218_SP3/log2FC.csv', quote = F)

#prep count_index and plot
count_index = rep(c(1,2,3), atac_diff$counts)
names(count_index) = colnames(atac_diff$log2FC)
atac_plot_order = plot_combo_ratios(atac_diff$binary, count_index)
save_plot('plots/paper/supplemental/13g.pdf', atac_plot_order$p, 
           base_width = 3, base_height = 2, useDingbats = F)

#plot subsetted SP7
SP7_by_atac_bin = SP7_diff$binary[,initials_to_names(colnames(atac_diff$binary))]
names(count_index) = initials_to_names(names(count_index))
SP7_atac_plot_order = plot_combo_ratios(SP7_by_atac_bin, count_index)
save_plot('plots/paper/supplemental/13h.pdf', SP7_atac_plot_order$p, 
           base_width = 3, base_height = 2, useDingbats = F)

#pull out actual reads data from ATACseq analysis
combos = if (is.list(stims)) multiple_combos(stims) else produce_combos(stims)
reads = merge_and_count(atac_data_dir, paste0('combo_', combos$name), 'Ctrl', combos$strings)$reads

#get log CPM counts
experimental.groups = names(colnames(reads))
dge <- DGEList(counts = reads, group = experimental.groups)
dge.norm <- calcNormFactors(dge)
log2_cpm_atac = t(cpm(dge.norm, log = TRUE, prior.count = 3))

#calculate mean for each treatment group
if(length(list.files(atac_data_dir, pattern='^mean_cpm.RData$')) && !overwrite) {
    load(str_c(atac_data_dir, '/mean_cpm.RData'))
} else { 
    atac_mean_log2_cpm = map_dfr(unique(names(rownames(log2_cpm_atac))), function (n)
        mean = as.data.frame(t(colMeans(as.matrix(log2_cpm_atac[names(rownames(log2_cpm_atac)) == n,])))),
        .id = 'group')
    atac_mean_log2_cpm$group = unique(names(rownames(log2_cpm_atac)))
    save(atac_mean_log2_cpm, file=str_c(atac_data_dir, '/mean_cpm.RData'))
    write.csv(t(atac_mean_log2_cpm), str_c(atac_data_dir, '/mean_cpm.csv'), quote=F)
}


#pca on atacseq
atac_pca_proc = prcomp(as.matrix(t(atac_diff$all_log2FC)), scale=TRUE,center=TRUE)
atac_pca_res = data.frame(group = colnames(atac_diff$all_log2FC), atac_pca_proc$x) %>% 
    mutate(number = as.character(ifelse(group=='Ctrl', 0, str_count(group, single_screen))))

#reduce rna to just the treatments included in the atacseq data
rna_stims = list(c('PAM', 'LPS', 'HMW'), c('PAM', 'SeV', 'cGAMP'), c('Zymosan', 'SeV', 'cGAMP'))
rna_combos = multiple_combos(rna_stims, str_break='_') 
rna_by_atac_lfc = SP7_diff$all_log2FC[,rna_combos$strings]

#run pca on subsetted rna
rna_by_atac_pca = prcomp(as.matrix(t(rna_by_atac_lfc)), scale=TRUE,center=TRUE)
rna_by_atac_pca_res = data.frame(group = names_to_initials(colnames(rna_by_atac_lfc)), rna_by_atac_pca$x) %>% 
    mutate(number = as.character(str_count(group, single_screen)))

#atac_plot
atac_plot = ggplot(atac_pca_res, aes(PC1, PC2)) + 
    geom_point(aes(fill=number), shape=21, color='black') + 
    geom_text(aes(label=group), hjust=-.1, vjust=-.1, size=1.5) +
    scale_fill_manual(values = c('#a6cee3', '#1f78b4', '#b2df8a'),
                      guide = 'none') +
    xlab(str_c('PC1: ', summary(atac_pca_proc)$importance[2,1] * 100, '% variance')) +
    ylab(str_c('PC2: ', summary(atac_pca_proc)$importance[2,2] * 100, '% variance')) + 
    ggtitle('ATAC-seq')
#rna_plot
rna_plot = ggplot(rna_by_atac_pca_res, aes(PC1, PC2)) + 
    geom_point(aes(fill=number), shape=21, color='black') + 
    geom_text(aes(label=group), hjust=-.1, vjust=-.1, size=1.5) +
    scale_fill_manual(values = c('#a6cee3', '#1f78b4', '#b2df8a'),
                      guide = guide_legend(title=NULL),
                      labels = c('singles (6)', 'pair (8)', 'triplet (3)')) +
    xlab(str_c('PC1: ', summary(rna_by_atac_pca)$importance[2,1] * 100, '% variance')) +
    ylab(str_c('PC2: ', summary(rna_by_atac_pca)$importance[2,2] * 100, '% variance')) + 
    ggtitle('RNA-seq')


single_screen = "S|cG|Z|P|L|H"

#heatmap setup
write.csv(t(atac_mean_log2_cpm), 'data/atacseq/020218_SP3/mean_cpm.csv', quote=F)
just_news = as.matrix(atac_diff$binary[,str_detect(colnames(atac_diff$binary), 'new')])
just_news = array(as.logical(just_news), dim(just_news), dimnames(just_news))
atac_annotation = tibble(gene = rownames(atac_diff$binary), new = rowSums(just_news) > 0) %>%
    mutate(new_from = map_chr(gene, function(g) 
            if (sum(just_news[g,]) > 0) str_c(colnames(just_news)[just_news[g,]], collapse=',')
            else ''),
        new_from_pair = map_lgl(str_split(new_from, ','), function (s) any(str_count(s,single_screen) == 2)),
        new_from_triple = map_lgl(str_split(new_from, ','), function (s) any(str_count(s,single_screen) == 3)))
proper_bin = abs(atac_diff$bin)[,!str_detect(colnames(atac_diff$binary), 'new')]
for(splits in list(c('P', 'S', 'cG'), c('Z', 'S', 'cG'), c('P', 'L', 'H'))) {
    treatment = str_c(splits[1], splits[2], splits[3])
    pairs = c(str_c(c(splits[1], splits[1], splits[2]), c(splits[2], splits[3], splits[3])))
    components = c(splits, pairs, treatment)
    atac_annotation = mutate(atac_annotation, 
                            !!treatment := as.numeric(rowSums(proper_bin[, components]) > 0))
}
write.table(atac_annotation, 'data/atacseq/020218_SP3/row_annotation.txt', quote=F, row.names = F, sep='\t')
