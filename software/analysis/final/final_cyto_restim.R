### libraries ###
library(tidyverse)
library(data.table)
library(cowplot)
library(parallel)

### sources ###
source('software/index_functions.r')

##### PCA minus 1 #####

# read data, select only variables of interest, calculate mean by treatment
restim_file = 'data/cytokines/121118_restim_legendplex_E28.txt'
restim_data = fread(restim_file) %>% 
    select(Treatment, IFNg_EL:IL13_LP) %>%
    group_by(Treatment) %>%
    summarise_at(vars(IFNg_EL:IL13_LP), mean)

# standard errors if needed
restim_se = fread(restim_file) %>% 
    select(Treatment, IFNg_EL:IL13_LP) %>%
    group_by(Treatment) %>%
    summarise_at(vars(IFNg_EL:IL13_LP), se)

# make matrix and do pca
as_mat = as.matrix(restim_data[,-1])
rownames(as_mat) = restim_data$Treatment
pca_proc = prcomp(as_mat, scale = T, center = T)
to_plot = as.tibble(pca_proc$x) %>% mutate(Treatment = rownames(pca_proc$x))

# plot pca
p=ggplot(to_plot, aes(PC1, PC2)) + 
    geom_point() + 
    geom_text(aes(label=Treatment), hjust=0, vjust=0) + 
    ggtitle('whole PCA')
save_plot('plots/paper/supplemental/16a.pdf', p, useDingbats = F, base_width = 2, base_height = 2)


##### IL17A isserlis #####

singles = c('cG', 'H', 'L', 'P', 'S', 'Z', 'Ova')
single_screen = str_c(singles, collapse = '|')

elisa_file = 'data/cytokines/ELISA_restim_E14-22-24-26_121718.txt'
elisa_data = fread(elisa_file) %>% 
    select(Treatment, Experiment, OVA_conc, LN, Mouse_ID, IL17A)  %>% 
    mutate(type = str_count(Treatment, single_screen),
           type = ifelse(Treatment == 'Ova', 0, type)) %>% 
    group_by(Treatment, OVA_conc, Experiment) 


# split treatment names into singles via regular expressions
treatment_splits = str_extract_all(unique(elisa_data$Treatment), single_screen)
# fill in missing singles with NA (so they all have three)
treatment_splits = lapply(treatment_splits, function(s) { if (length(s) == 1) c(s, NA, NA)
    else if (length(s) == 2) c(s, NA)
    else s })

# setup component key
treatment_splits = matrix(unlist(treatment_splits), nrow=3)
colnames(treatment_splits) = unique(elisa_data$Treatment)
component_key = tibble(Treatment = unique(elisa_data$Treatment)) %>%
    mutate(C_1 = treatment_splits[1,],
           C_2 = treatment_splits[2,],
           C_3 = treatment_splits[3,],
           C_12 = str_c(treatment_splits[1,], treatment_splits[2,]),
           C_13 = str_c(treatment_splits[1,], treatment_splits[3,]),
           C_23 = str_c(treatment_splits[2,], treatment_splits[3,]),
           C_123 = str_c(treatment_splits[1,], treatment_splits[2,], treatment_splits[3,]))

# calculate means and se within Treatment/OVA groups
raw_summaries = summarize(elisa_data, mean = mean(IL17A, na.rm=T), se = se(IL17A)) %>%
    mutate(Screen_rep_ID = interaction(OVA_conc, Experiment))

# make component matrix for mean and se
pulled_summaries = bind_cols(raw_summaries, 
                         as.tibble(values_from_key(raw_summaries, 'mean', component_key, T)),
                         as.tibble(values_from_key(raw_summaries, 'se', component_key, T))) 

pulled_triplets = pulled_summaries %>%
    # keep just triplets
    filter(str_count(Treatment, single_screen) == 3) %>%
    # gather means and se, such that type contains the metric (mean/se) and the component
    gather(mean_1:se_123, key='type', value='val') %>% 
    # extract type into a metric and component column
    extract(type, c('metric', 'component'), '(.*)_(.*)') %>% 
    # spread metric and value so that there is a mean and se column
    spread(metric, val) %>%
    # calculate error bars, type (single, pair, triplet) and the name from the component
    mutate(err_top = mean + se,
           err_bot = mean - se,
           type = as.character(str_length(component)),
           name = map2_chr(Treatment, component, function(treat, comp)
               str_c(str_extract_all(treat, single_screen)[[1]]
                     [as.numeric(str_split(comp, '')[[1]])], collapse = ''))) %>%
    # order by type (single, pair, triplet)
    arrange(type)

# factor the name (so it's in order by type)
pulled_triplets$name = factor(pulled_triplets$name, levels = unique(pulled_triplets$name))

# make a bar plot for each OVA and Treatment combination
plot_types = expand.grid(unique(pulled_triplets$OVA_conc), unique(pulled_triplets$Treatment))
raw_plots = map2(plot_types$Var1, plot_types$Var2,
     function(OVA, treat) {
         data = filter(pulled_triplets, Treatment == treat, OVA_conc == OVA)
         ymax = max(data$err_top)
         rounded_max = 200*(ymax%/%200 + as.logical(ymax%%200))
         ggplot(data, aes(name, mean)) + 
             geom_col(aes(fill = type), width = .4) +
             geom_errorbar(aes(ymin = err_bot, ymax=err_top), width = .2) +
             ggtitle(str_c(treat, '_', OVA)) + 
             scale_y_continuous(breaks = c(0, rounded_max/2, rounded_max), limits = c(0, rounded_max))
     }
)

save_plot('plots/paper/5c_sup17ab.pdf', 
          plot_grid(plotlist = raw_plots[plot_types$Var2 != 'LHcG'], nrow = 4, ncol = 3),
          nrow = 4, ncol = 4, base_aspect_ratio = 1.5, base_height = 1)


### combined boostrap ###
experiments = c(12, 14, 16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 29)

# read in data
raw_sum = read_csv('data/cytokines/restim_ELISA_all_052819.txt') %>%
    # remove outliers and bad experiments
    filter(IL17A_outliers == 0, Experiment %in% experiments, OVA_conc != 50, !is.na(IL17A)) %>%
    # avg within experiment/treatments (across replicates/sides)
    group_by(Experiment, Treatment, OVA_conc) %>%
    summarize(mean = mean(IL17A), se = se(IL17A)) %>%
    # calculate error bars
    mutate(err_top = mean + se, err_bot = mean - se,
           number = factor(str_count(Treatment, single_screen)))

treatments = unique(raw_sum$Treatment)
treatments = treatments[treatments != 'Ova']

# split treatment names into singles via regular expressions
treatment_splits = str_extract_all(treatments, single_screen)
# fill in missing singles with NA (so they all have three)
treatment_splits = lapply(treatment_splits, function(s) { if (length(s) == 1) c(s, NA, NA)
    else if (length(s) == 2) c(s, NA)
    else s })

# setup component key
treatment_splits = matrix(unlist(treatment_splits), nrow=3)
colnames(treatment_splits) = treatments
component_key = tibble(Treatment = treatments) %>%
    mutate(C_1 = treatment_splits[1,],
           C_2 = treatment_splits[2,],
           C_3 = treatment_splits[3,],
           C_12 = str_c(treatment_splits[1,], treatment_splits[2,]),
           C_13 = str_c(treatment_splits[1,], treatment_splits[3,]),
           C_23 = str_c(treatment_splits[2,], treatment_splits[3,]),
           C_123 = str_c(treatment_splits[1,], treatment_splits[2,], treatment_splits[3,]))

### normalization ###
pulled_raw = raw_sum %>% 
    # only interested in mean
    select(Experiment:OVA_conc, IL17A = mean) %>%
    # make screen_rep_ID for pulling components (just Exp, OVA)
    mutate(Screen_rep_ID = str_c(Experiment, OVA_conc, sep = '_')) %>%
    # don't want Ova
    filter(Treatment != 'Ova')

# pull components
pulled_raw = bind_cols(pulled_raw, as.tibble(values_from_key(pulled_raw, 'IL17A', component_key, T))) %>%
    # only keep triplets
    filter(str_count(Treatment, single_screen) == 3)


norm = pulled_raw %>%
    # original IL17A column is the same as IL17A_123
    select(-IL17A) %>%
    gather(IL17A_1:IL17A_123, key = 'component', value = 'IL17A') %>%
    group_by(Treatment, Experiment) %>%
    # divide by max in Triplet/Experiment group
    mutate(IL17A = IL17A/max(IL17A))

norm_factor = 2
for_bootstrap = norm %>%
    # keep only triplets we want for bootstrapping
    filter(Treatment %in% c('HZcG', 'ZScG', 'LScG', 'PScG')) %>%
    ungroup() %>%
    # calculate one_minus_norm
    mutate(one_minus_norm = 1 - IL17A / norm_factor,
           component = str_replace(component, 'IL17A_', ''),
           triplet = Treatment, 
           Treatment = map2_chr(str_extract_all(Treatment, single_screen), str_split(component, ''), 
                                function(names, needed) 
                                    str_c(names[as.numeric(needed)], collapse='')),
           Screen_rep_ID = str_c(Screen_rep_ID, '_', triplet)) 


# calculate appropriate isserlis differences distribution 
filtered_diffs = for_bootstrap %>%
    group_by(Treatment, OVA_conc, component, triplet) %>%
    summarize(one_minus_norm = mean(one_minus_norm)) %>%
    ungroup() %>%
    select(-Treatment) %>%
    mutate(component = str_c('one_minus_norm_', component)) %>%
    spread(component, one_minus_norm)
filtered_diffs$iss = unlist(isserlis(filtered_diffs, 'one_minus_norm'))
filtered_diffs = filtered_diffs %>%
    mutate(diff = abs(one_minus_norm_123 - iss)) %>%
    select(triplet, OVA_conc, one_minus_norm_123, iss, diff)

### bootstrap ###
raw_triplets = filter(for_bootstrap, str_length(component) == 3)
not_triplets = filter(for_bootstrap, str_length(component) != 3)

n_it = 1000
n_cores = 10 #for parallelization
replace=F


bootstrap_it = function(dummy) {
    shuffled = not_triplets
    order = sample.int(nrow(not_triplets), replace=replace)
    shuffled$one_minus_norm = shuffled$one_minus_norm[order]
    shuffled = bind_rows(shuffled, raw_triplets)
    
    shuffled = cbind(raw_triplets, split_values_from_key(raw_triplets, shuffled, 'one_minus_norm', component_key, T))
    shuffled = shuffled %>%
        filter(str_count(Treatment, single_screen) == 3) %>%
        group_by(Treatment, OVA_conc) %>%
        summarize_at(vars(one_minus_norm_1:one_minus_norm_123), mean)
    shuffled$iss = unlist(isserlis(shuffled, 'one_minus_norm'))

    diffs = abs(shuffled$one_minus_norm_123 - shuffled$iss)
    
    return(diffs)
}

diff_dist = unlist(mclapply(1:n_it, bootstrap_it, mc.cores = n_cores))
p = bind_rows(bootstrap = tibble(diff = diff_dist), original = tibble(diff = filtered_diffs$diff), .id = 'source') %>%
    ggplot(aes(diff, fill = source)) + 
    geom_density(alpha = .5, bw=.05) 

save_plot('plots/paper/5d.pdf', p, base_aspect_ratio = 2)