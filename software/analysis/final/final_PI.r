### libraries ###
library(data.table)
library(tidyverse)
library(parallel)
library(cowplot)

### sources ###
source('software/index_functions.r')

#plot theming 
theme_set(theme_classic())
theme_update(text = element_text(size = rel(2.25)),
             legend.key.size = unit(.1, 'inches'),
             legend.text = element_text(size = rel(1.5)),
             panel.background = element_rect(fill = "white", colour = "black"),
             axis.line = element_blank())


# read in data, remove excess information
PI_data_file = 'data/PI/dat_all_isserlis_select_120718.txt'
PI_data = as.tibble(fread(PI_data_file)) %>% 
    select(Sample, Treatment, Screen_rep_ID, Treatment_type, Peak_0:PI_norm)

# single screen is used for pulling apart treatment names using regular expressions
singles = c('Z','S','Cp','P','cG','L', 'H')
single_screen = str_c(singles, collapse='|')

##### CSFE barplots #####
peak_summaries = PI_data %>%
    # summarize peaks by treatment
    group_by(Treatment, Treatment_type) %>%
    summarize_at(vars(Peak_0:Peak_6), funs(mean(., na.rm=T), se)) %>%
    ungroup() %>%
    # set order of treatment_type
    mutate(Treatment_type = factor(Treatment_type, levels = c('Single', 'Pair', 'Triplet')),
           Treatment = str_replace(Treatment, 'cG', 'G'),
           Treatment = str_replace(Treatment, 'Cp', 'C')) %>%
    # gather mean and se's
    gather(-Treatment, -Treatment_type, key='key', value='val') %>%
    # separate key into peak and mean/se
    extract('key', c('Peak', 'metric'), '(.*_.*)_(.*)') %>%
    # spread mean/se into separate columns
    spread(metric, val) %>%
    # order by treatment type
    arrange(Treatment_type) %>%
    # make errorbars
    mutate(err_top = mean + se, err_bot = mean - se)
colors = c(Single = 'darkgreen', Pair = 'blue', Triplet = 'red')

plots = list()
for(treat in unique(peak_summaries$Treatment)) {
    treat_data = filter(peak_summaries, Treatment == treat)
    yMax = max(treat_data$err_top)
    plots[[treat]] = ggplot(treat_data, aes(desc(Peak), mean)) + 
        geom_errorbar(aes(ymin=err_bot, ymax=err_top, width=.4)) +
        geom_col(fill = colors[treat_data$Treatment_type[1]]) +
        ggtitle(treat) +
        scale_y_continuous(breaks = c(0, yMax/2, yMax), limits = c(0, yMax)) +
        theme(plot.title = element_text(size = 10, hjust = 0.5),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 6),   
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = 'black', fill = NA))
}
save_plot('plots/paper/supplemental/4abc.pdf', plot_grid(plotlist=plots, ncol = 7), ncol=7, nrow=9, base_height = 2)

# calculate proportion cells divided
prop_div = peak_summaries %>% 
    mutate(p0 = Peak == 'Peak_0') %>%
    group_by(Treatment, p0) %>% 
    summarize(sum = sum(mean)) %>%
    spread(p0, sum) %>%
    transmute(prop_divided = `FALSE` / (`FALSE` + `TRUE`))
write_csv(prop_div, 'data/PI/prop_div.csv')

# histogram of proportion divided
p = prop_div %>%
    mutate(type = c('Single', 'Pair', 'Triplet')[str_length(Treatment)],
           type = factor(type, levels = c('Single', 'Pair', 'Triplet'))) %>%
    ggplot(aes(prop_divided)) + 
        geom_histogram(bins = 10) +
        facet_grid(. ~ type )
save_plot('plots/paper/1d.pdf', p, ncol = 3)

# slightly different take on calculating proportion divided
raw_prop_div = PI_data %>%
    select(Sample, Treatment, Treatment_type, Peak_0:Peak_6) %>%
    gather(Peak_0:Peak_6, key='Peak', value = 'count') %>%
    mutate(p0 = Peak == 'Peak_0',
           Treatment_type = factor(Treatment_type, levels = c('Single', 'Pair', 'Triplet')),
           Treatment = str_replace(Treatment, 'cG', 'G'),
           Treatment = str_replace(Treatment, 'Cp', 'C')) %>%
    group_by(Sample, Treatment, Treatment_type, p0) %>% 
    summarize(sum = sum(count)) %>%
    spread(p0, sum) %>%
    transmute(prop_divided = `FALSE` / (`FALSE` + `TRUE`)) %>%
    group_by(Treatment, Treatment_type) %>%
    summarize(SEM = se(prop_divided), AVG = mean(prop_divided, na.rm=T)) %>%
    filter(!is.na(AVG)) %>%
    mutate(Datatype = "Prop_div") %>%
    arrange(Treatment_type, desc(AVG))

# compare IFNg and Viability to prop divided
big_elisa = read_tsv('data/cytokines/110918_ELISA_IFNg_IL17_Viability_norm_formatted.txt') %>%
    mutate(Treatment = str_replace(Treatment, 'cG', 'G'),
           Treatment = str_replace(Treatment, 'Cp', 'C')) %>%
    select(Treatment, Treatment_type = Type, Datatype = Cytokine, AVG, SEM) %>%
    filter(Datatype != 'IL17')

p = bind_rows(big_elisa, raw_prop_div) %>%
    mutate(Treatment = factor(Treatment, levels = unique(raw_prop_div$Treatment)),
           err_top = AVG + SEM, err_bot = AVG - SEM,
           Datatype = factor(Datatype, levels = c('Prop_div', 'IFNg', 'Viability'))) %>%
    ggplot(aes(Treatment, AVG)) + 
    geom_errorbar(aes(ymin = err_bot, ymax = err_top), width = .5) + 
    geom_col(aes(fill= Treatment_type)) + 
    facet_grid(Datatype ~ ., scales='free')

save_plot('plots/paper/supplemental/5abc.pdf', p, nrow=3, base_aspect_ratio = 4)

##### Setup for building out component matrices #####

# split treatment names into singles via regular expressions
treatment_splits = str_extract_all(unique(PI_data$Treatment), single_screen)
# fill in missing singles with NA (so they all have three)
treatment_splits = lapply(treatment_splits, function(s) {if (length(s) == 1) c(s, NA, NA)
                                                               else if (length(s) == 2) c(s, NA)
                                                               else s })
# convert to matrix and name
treatment_splits = matrix(unlist(treatment_splits), nrow=3)
colnames(treatment_splits) = unique(PI_data$Treatment)

# make a key mapping treatment names to the names of their component combinations
component_key = tibble(Treatment = unique(PI_data$Treatment)) %>%
    mutate(C_1 = treatment_splits[1,],
           C_2 = treatment_splits[2,],
           C_3 = treatment_splits[3,],
           C_12 = str_c(treatment_splits[1,], treatment_splits[2,]),
           C_13 = str_c(treatment_splits[1,], treatment_splits[3,]),
           C_23 = str_c(treatment_splits[2,], treatment_splits[3,]),
           C_123 = str_c(treatment_splits[1,], treatment_splits[2,], treatment_splits[3,]))


##### Making isserlis dot plots #####

### Function for making dot plots comparing experimental values to isserlis predicted values ###
### Takes a data_frame, the name of the metric to be compared, and whether or not to scale (center) the data ###

dot_plot = function(data, name, scale) {
    # the x value will be the normalized real, and the y the isserlis prediction
    x = str_c(name, '_norm')
    y = str_c(name, '_iss')
    
    # calculate mean for centering
    x_mean = mean(unlist(data[, str_c(x, '_mean')]))
    y_mean = mean(unlist(data[, str_c(y, '_mean')]))
    
    # calculate centering of the data and error bars
    annot = data %>% mutate(x_centered = if (scale) !!as.name(str_c(x, '_mean')) - x_mean + .5 else !!as.name(str_c(x, '_mean')),
                            y_centered = if (scale) !!as.name(str_c(y, '_mean')) - y_mean + .5 else !!as.name(str_c(y, '_mean')),
                            x_err_top = x_centered + !!as.name(str_c(x, '_se')),
                            x_err_bot = x_centered - !!as.name(str_c(x, '_se')),
                            y_err_top = y_centered + !!as.name(str_c(y, '_se')),
                            y_err_bot = y_centered - !!as.name(str_c(y, '_se')))
    
    # model the relationship between isserlis and experimental results 
    model = summary(lm(y_centered ~ x_centered, annot))
    r_squared = model$r.squared
    slope = model$coefficients[2,1]
    
    # max and mins for labeling plot with r-squared and model
    xmin = min(annot$x_centered)
    xmax = max(annot$x_centered)
    ymax = max(annot$y_centered)
    
    # make plot
    p = ggplot(annot, aes(x_centered, y_centered)) +
        geom_point(size=.75) + 
        annotate('text', x=xmin*1.2, y=ymax*.8, label = str_c("italic(R) ^ 2 == ", round(r_squared, 2)), parse=T, size=2) + 
        geom_errorbar(aes(ymin=y_err_bot, ymax=y_err_top), size=.25) + 
        geom_errorbarh(aes(xmin = x_err_bot, xmax=x_err_top), size=.25) + 
        xlab(x) + 
        ylab(y)
    
    if (scale) {
        # if scaled just plot x=y line
        p = p + 
            annotate('segment', x=0, y=0, xend=1, yend=1) + 
            annotate('text', x = xmax*.9, ymax*.9, label = "y = x", size=3)
    } else {
        # if not scaled plot line of best fit
        p = p +
            geom_abline(slope = slope, intercept = 0) +
            annotate('text', x = xmax*.9, ymax*.9, 
                     label = str_c("y = ", round(slope, 2), "x"), size=2) +
            geom_hline(yintercept = 0) + 
            geom_vline(xintercept = 0)
    }
    
    return(p)
}

# calculate DI, EI and RI, and normalize
index_data = PI_data %>%
    mutate(DI = divisions/cell_start,
           EI = cell_total/cell_start,
           RI = cell_activated/cell_start) %>%
    group_by(Screen_rep_ID) %>%
    mutate(DI_norm = DI/max(DI, na.rm = T),
           EI_norm = EI/max(EI, na.rm = T),
           RI_norm = RI/max(RI, na.rm = T)) %>%
    ungroup()

# build out component matrix for indices
index_data = cbind(index_data, values_from_key(index_data, 'PI_norm', component_key, T),
                   values_from_key(index_data, 'DI_norm', component_key, T),
                   values_from_key(index_data, 'EI_norm', component_key, T),
                   values_from_key(index_data, 'RI_norm', component_key, T)) %>%
    # calculate isserliss values for matrices
    mutate(PI_iss = isserlis(., 'PI_norm'),
           DI_iss = isserlis(., 'DI_norm'),
           EI_iss = isserlis(., 'EI_norm'),
           RI_iss = isserlis(., 'RI_norm')) %>%
    # calculate interaction values for PI
    make_interaction() %>% 
    #group for summarization
    group_by(Treatment, Treatment_type)

# calculate mean and se for each index 
summaries = index_data %>%
    summarize_at(vars(PI_norm, DI_norm, EI_norm, RI_norm, PI_iss:I_iss), funs(mean(., na.rm=T), se))

# filter for just triplets
triplets = filter(summaries, Treatment_type == 'Triplet')

# calculate R^2 value for correlation between each index's experimental and predicted values
r_squareds = numeric(4)
names(r_squareds) = c('PI', 'DI', 'EI', 'RI')
for (name in names(r_squareds)) {
    r_squareds[name] = summary(lm(as.formula(str_c(name, '_norm_mean ~ ', name, '_iss_mean')), triplets))$r.squared
}

# plot PI values and interactions dot plots
PI_plot = dot_plot(triplets, 'PI', T)
interaction_plot = dot_plot(rename(triplets, I_norm_mean = I_123_mean, I_norm_se = I_123_se), 'I', F)
save_plot('plots/paper/2d_e.pdf', plot_grid(PI_plot, interaction_plot, nrow=1), useDingbats = F, base_width = 5, base_height = 3)

# plot dot plots for other indices
for(name in c('DI', 'EI', 'RI')){
    p = dot_plot(triplets, name, T)
   save_plot(str_c('plots/paper/supplemental/7_', name, '_dotplot.pdf'), p, useDingbats=F)
}

#function for taking ligand initials and returning a vector of ligand families
initials_to_families = function(Treatment) {
    initials = c('L', 'Cp','P', 'H', 'Z','S', 'cG')
    families = c( rep('TLR', 4),'CLR', 'RLR', 'CDS')
    
    #we make a matrix of counts for each family because we need uniform ordering
    counts = matrix(0,nrow = length(Treatment), ncol=4, dimnames = list(NULL, unique(families)))
    for(i in 1:length(initials)) {
        counts[,families[i]] = counts[,families[i]] + str_count(Treatment, initials[i])
    }
    
    #turn counts into strings
    treatment_families = apply(counts, 1, function(cs) str_c(rep(unique(families), cs), collapse='_'))
    return(treatment_families)
}

# plot PI experimental vs isserlis color coded by ligand families
exp_mean = mean(triplets$PI_norm_mean)
iss_mean = mean(triplets$PI_iss_mean)
triplets = mutate(triplets, ligand_families = initials_to_families(Treatment),
                  centered_exp = PI_norm_mean - exp_mean + .5,
                  centered_iss = PI_iss_mean - iss_mean + .5)

p = ggplot(triplets, aes(centered_exp, centered_iss)) + 
    geom_point(aes(col=ligand_families)) +
    annotate('segment', x=0, y=0, xend=1, yend=1) + 
    geom_point(aes(col=ligand_families)) + 
    xlab("Normalized PI") +
    ylab("Predicted Normalized PI") + 
    geom_text(aes(label = Treatment), size = 1)
save_plot('plots/paper/supplemental/9b.pdf', 
          ggExtra::ggMarginal(p, type='density', groupFill = T))

# plot PI experimental vs isserlis color coded by number of unique ligand families
triplets = mutate(triplets, family_count = map_chr(str_split(ligand_families, '_'), function(s) length(unique(s))))
p = ggplot(triplets, aes(centered_exp, centered_iss)) + 
    geom_point(aes(col=family_count)) +
    annotate('segment', x=0, y=0, xend=1, yend=1) + 
    geom_point(aes(col=family_count)) + 
    xlab("Normalized PI") +
    ylab("Predicted Normalized PI") + 
    geom_text(aes(label = Treatment), size = 1)
save_plot('plots/paper/supplemental/9a.pdf',
          ggExtra::ggMarginal(p, type='density', groupFill = T))

##### Making peak line plots and interaction barplots #####

# calculate averages and standard errors for each peak
peak_averages = PI_data %>%
    select(Treatment, Treatment_type, Peak_1:Peak_6) %>%
    group_by(Treatment, Treatment_type) %>%
    summarise_at(vars(Peak_1:Peak_6), funs(mean(., na.rm=T), se))

# for each peak, build out the component matrix (pull the values using the key)
pulled = map(names(peak_averages)[c(-1,-2)], values_from_key, 
                data=peak_averages, key=component_key, match=F)
pulled = as.tibble(do.call(cbind, pulled))

# convert to tidy data
tidied_averages = bind_cols(peak_averages, pulled) %>% 
    #keep only the triplets
    filter(Treatment_type == 'Triplet') %>%
    #keep only treatment and peak values
    select(Treatment, Peak_1_mean_1:Peak_6_se_123) %>%
    # gather, such that each entry has a Treatment, type (containing Peak, mean/se, and component), and value
    gather(Peak_1_mean_1:Peak_6_se_123, key='type', value = 'val') %>%
    # separate type into peak, type (mean/se) and component
    extract(type, c('peak', 'type', 'component'), 'Peak_(.)_(.*)_(.*)') %>%
    #only interested in single and pair components
    filter(str_length(component) < 3) %>%
    # spread by type, so that now we have mean and se columns
    spread(type, val) %>%
    # add name by pulling component names from Treatment
    mutate(name = map2_chr(Treatment, component, function(treat, comp)
        str_c(str_extract_all(treat, single_screen)[[1]]
                        [as.numeric(str_split(comp, '')[[1]])], collapse = '')),
        # add error bar tops and bottoms
        err_top = mean + se,
        err_bot = mean - se,
        # correct types
        peak = as.numeric(peak),
        component = factor(component, levels = c('1','2','3','12','13','23')))

# the plot we'd like to make has 3 'panes' per triplet: one for each pair. 
# each single appears in two pairs, so we need to duplicate them
# the bind_rows binds the original to a set of duplicates, and marks the duplicates in the 'copy' column
# the mutate uses copy and component to assign a pane to the singles (pairs are their pane)
tidied_averages = bind_rows('1' = tidied_averages, 
                            '2' = filter(tidied_averages, str_length(component) == 1),
                            .id = 'copy') %>%
    mutate(pane = ifelse(component == '1', ifelse(copy == '1', '12', '13'),
                         ifelse(component == '2', ifelse(copy == '1', '12', '23'),
                                ifelse(component == '3', ifelse(copy == '1', '13', '23'),
                                       as.character(component)))))

# tidying triplets into data for interaction box plots
bar_plot_melt = triplets %>%
    # keep just interactions
    select(Treatment, Treatment_type, I_12_mean:I_iss_mean, I_12_se:I_iss_se) %>%
    # gather columns, type describes the component and the type of metric (mean/se)
    gather(I_12_mean:I_iss_se, key='type', value='val') %>% 
    # separate type into component and type (mean/se)
    extract(type, c('component', 'type'), 'I_(.*)_(.*)') %>% 
    # spread val by type into separate columns for mean and se
    spread(type, val) %>% 
    # calculate error bars, factor component, color based on component
    mutate(err_top = mean + se,
           err_bot = mean - se,
           component = factor(component, levels = c('12', '13', '23', '123', 'iss')),
           color = ifelse(component %in% c('12', '13', '23'), 'double', component))

# keep just pairs and their interaction values
pairs = summaries %>%
    filter(Treatment_type == 'Pair') %>%
    select(Treatment, I_12_mean, I_12_se)

# thresholds for synergistic and antagonistic are the average of standard errors
pair_threshold = mean(pairs$I_12_se)
triplet_threshold = mean(triplets$I_123_se)

### function for making peak line plots and interaction bar plots ###
# inputs:
    # trips: a character vector of the triplets to include
    # bp_trip: the triplets as tidied and prepared for bar plots (always bar_plot_melt)
    # tidied_averages: the peaks as tidied and averaged for line plots (always tidied_averages)
    # pair_threshold: the Synergistic/Antagonisitc threshold for pairs
    # triplet_threshold: the Synergisitc/Antagonistic threshold for triplets

plot_lines = function(trips, bp_trip, tidied_averages, pair_threshold, triplet_threshold) {
    plots = list()
    # make a plot for each triplet
    for(treat in trips) {
        # keep just this triplet for lineplot
        this_treat = filter(tidied_averages, Treatment == treat)
        
        # calculate maximum and minimums for lineplots
        ymin = min(this_treat$err_bot)
        ymax = ceiling(max(this_treat$err_top))
        rounded_max = 100*(ymax%/%100 + as.logical(ymax%%100))
        
        # method for getting name from component (for better labeling)
        comp_to_name = unique(this_treat$name)
        names(comp_to_name) = unique(this_treat$component)
        # names of components for each pair
        comp_names = map(c('12', '13', '23'), function(p) 
            comp_to_name[c(str_sub(p,1,1), str_sub(p,2,2), p)])
        
        # plot lines for each pair
        lines = map2(c('12', '13', '23'), comp_names, function(p, cn)
            ggplot(filter(this_treat, pane==p), aes(peak, mean, col = component)) + 
                geom_errorbar(aes(ymin = err_bot, ymax = err_top), width=.2, alpha=.5) + 
                geom_line() +
                theme(legend.position = c(.9,.8), legend.title=element_blank(),
                      legend.background = element_blank()) + 
                scale_x_continuous(breaks = 1:6) +
                scale_y_continuous(breaks = c(0, rounded_max/2, rounded_max), limits = c(0, ymax)) + 
                coord_cartesian(ylim = c(0, rounded_max)) +
                scale_color_hue(labels = cn))
        
        # make bar plots for just the pair interactions
        bar_pairs = bp_trip %>% 
            filter(Treatment == treat, component %in% c('12', '13', '23')) %>%
            ggplot(aes(component, mean)) + 
                geom_col(aes(fill = color), width = .4) +
                geom_errorbar(aes(ymin = err_bot, ymax=err_top), width = .2) +
                geom_hline(yintercept = pair_threshold, linetype='dotted') + 
                geom_hline(yintercept = -pair_threshold, linetype='dotted') + 
                geom_hline(yintercept = 0) +
                ylim(c(-.25, .25)) +
                guides(fill='none') +
                scale_x_discrete(labels = comp_to_name[c('12','13','23')])
        
        # make bar plots for the experimental and isserlis triplet interactions
        bar_trips = bp_trip %>% 
            filter(Treatment == treat, !component %in% c('12', '13', '23')) %>%
            ggplot(aes(component, mean)) + 
                geom_col(aes(fill = color), width = .4) +
                geom_errorbar(aes(ymin = err_bot, ymax=err_top), width = .2) +
                geom_hline(yintercept = triplet_threshold, linetype='dotted') + 
                geom_hline(yintercept = -triplet_threshold, linetype='dotted') + 
                geom_hline(yintercept = 0) + 
                ylim(c(-.25, .25)) + 
                guides(fill='none')
        
        #combine plot components into one row
        plots[[treat]] = plot_grid(plotlist = c(lines, list(bar_pairs, bar_trips)), nrow = 1, rel_widths = c(2,2,2,1,4/5))
    }
    return(plot_grid(plotlist=plots, ncol=1))
}

# plot lines for just three triplets
example_triplets = c('ZSCp', 'PScG', 'PLZ')
p=plot_lines(example_triplets, bar_plot_melt, tidied_averages, pair_threshold, triplet_threshold)
save_plot('plots/paper/2abc.pdf', p, useDingbats = F, base_width = 7, base_height = 4)


##### Making bar plots for triplet and pair interactions #####

# calculate errorbars and if synergistic/antagonistic
pairs = pairs %>% 
    mutate(err_top = I_12_mean + I_12_se,
           err_bot = I_12_mean - I_12_se,
           syn_ant_add = ifelse(err_bot > pair_threshold, 'syn',
                                ifelse(err_top < -pair_threshold, 'ant', 'add')))
pairs$rank = row_number(desc(pairs$I_12_mean))

# calculate proportions of synergistic/antagonisitc
syn = mean(pairs$syn_ant_add == 'syn')
ant = mean(pairs$syn_ant_add == 'ant')
add = 1 - syn - ant

# plot pair bars
p = ggplot(pairs, aes(rank, I_12_mean)) + 
    geom_col(aes(fill = syn_ant_add), col='black', size=.25) + 
    geom_errorbar(aes(ymin = err_bot, ymax = err_top)) +
    annotate('text', x = nrow(pairs)*.8, y = max(pairs$I_12_mean)*.6, 
             label=str_c('SYN: ', round(syn, 2), ', ANTAG: ', 
                         round(ant, 2), ', ADD: ', round(add, 2)), size=1) + 
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = pair_threshold, linetype='dotted') + 
    geom_hline(yintercept = -pair_threshold, linetype='dotted') + 
    scale_x_continuous(breaks = sort(pairs$rank), labels = arrange(pairs, rank)$Treatment) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
save_plot('plots/paper/supplemental/6a.pdf', p, base_width=3, base_height=3, useDingbats=F)

# calculate errorbars and if synergistic/antagonistic
triplets = triplets %>%
    mutate(err_top = I_123_mean + I_123_se,
           err_bot = I_123_mean - I_123_se,
           syn_ant_add = ifelse(err_bot > triplet_threshold, 'syn',
                                ifelse(err_top < -triplet_threshold, 'ant', 'add')))
triplets$rank = row_number(desc(triplets$I_123_mean))

# calculate proportions of synergistic/antagonisitc
syn = mean(triplets$syn_ant_add == 'syn')
ant = mean(triplets$syn_ant_add == 'ant')
add = 1 - syn - ant

# plot triplet bars
p = ggplot(triplets, aes(rank, I_123_mean)) + 
    geom_col(aes(fill = syn_ant_add)) + 
    geom_errorbar(aes(ymin = err_bot, ymax = err_top)) +
    annotate('text', x = nrow(triplets)*.8, y = max(triplets$I_12_mean)*.6, 
         label=str_c('SYN: ', round(syn, 2), ', ANTAG: ', 
                     round(ant, 2), ', ADD: ', round(add, 2)), size = 1) + 
    geom_hline(yintercept = triplet_threshold, linetype='dotted') + 
    geom_hline(yintercept = -triplet_threshold, linetype='dotted') + 
    scale_x_continuous(breaks = sort(triplets$rank), labels = arrange(triplets, rank)$Treatment) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

save_plot('plots/paper/supplemental/6b.pdf', p, base_width=3, base_height=3, useDingbats=F)

##### Bootstraping #####
## Goal: produce a distribution of R-squared values for PI, DI, EI, and RI which can be plotted
    ## and used to produce p-values for the original R-squared values
## Method: Shuffle single and pair values, recalculate isserlis predictions, and calculate an R-squared

# select only necessary data
for_bootstrap = index_data %>%
    select(Treatment, Screen_rep_ID, Treatment_type, PI_norm, DI_norm, EI_norm, RI_norm) %>%
    ungroup()

# Separate triplets and not triplets so we can shuffle not triplets
not_triplets = filter(for_bootstrap, Treatment_type != 'Triplet')
raw_triplets = filter(for_bootstrap, Treatment_type == 'Triplet')


n_it = 1000 # how many iterations to use
n_cores = 10 # how many threads for parallelization
replace=F # whether or not to sample with replacement

### function that runs one iteration of bootstrapping ###
### takes a dummy variable (needed for mclapply), returns vector of R-squared for PI, RI, DI and EI ###
### This version is depricated due to slowness ###
### split_bootstrap_it() is used instead ###

split_bootstrap_it = function(dummy) {
    
    # shuffle the indices of all the pairs and singles, and rebind them to the triplets
    shuffled = not_triplets
    order = sample.int(nrow(not_triplets), replace=replace)
    shuffled$PI_norm = shuffled$PI_norm[order]
    shuffled$DI_norm = shuffled$DI_norm[order]
    shuffled$EI_norm = shuffled$EI_norm[order]
    shuffled$RI_norm = shuffled$RI_norm[order]
    shuffled = bind_rows(shuffled, raw_triplets)
    
    # build out the component matrix
    shuffled = cbind(raw_triplets, split_values_from_key(raw_triplets, shuffled, 
                                                         c('PI_norm', 'DI_norm', 'EI_norm', 'RI_norm'), component_key, T)) %>%
        # calculate isserlis
        mutate(PI_iss = isserlis(., 'PI_norm'),
               DI_iss = isserlis(., 'DI_norm'),
               EI_iss = isserlis(., 'EI_norm'),
               RI_iss = isserlis(., 'RI_norm')) %>%
        # keep just triplets (I think this is redundant but I'm not positive)
        filter(Treatment_type == 'Triplet') %>%
        # group for summarizing 
        group_by(Treatment)
    
    # summarize (calculating standard error I think is unnecessary)
    shuffled_summaries = shuffled %>%
        summarize_at(vars(PI_norm, DI_norm, EI_norm, RI_norm, PI_iss:RI_iss), funs(mean(., na.rm=T), se))
    
    # calculate the correlations
    r_squareds = numeric(4)
    names(r_squareds) = c('PI', 'DI', 'EI', 'RI')
    for (name in names(r_squareds)) {
        r_squareds[name] = cor(shuffled_summaries[,str_c(name, '_norm_mean')], shuffled_summaries[,str_c(name, '_iss_mean')])^2
    }
    return(r_squareds)
}

# Actual bootstrap code has been commented out so that this code doesn't take an eternity to run:

############################## Bootstrap run ################################

# r2_dist = data.frame()

# run 1000 iterations at a time in parallel
# collect at every 1000 and write to a file
# for(i in 1:(n_it/1000)) {
#     these_r2s = mclapply(1:1000, split_bootstrap_it, mc.cores=n_cores)
#     these_r2s = as.data.frame(do.call(rbind, these_r2s))
#     r2_dist = bind_rows(r2_dist, these_r2s)
#     write.table(r2_dist, 'data/PI/r_squareds.txt')
# }

# plot distributions of each index, labelled with p-values
# plots = list()
# for(stat in colnames(r2_dist)) {
#     pval = mean(r2_dist[,stat] > r_squareds[stat])
#     plots[[stat]] = ggplot(r2_dist, aes(!!as.name(stat))) +
#         stat_density(geom='line') +
#         geom_vline(xintercept = r_squareds[stat], col = 'red') +
#         ggtitle(str_c(stat, ' p = ', pval)) +
#         xlim(0, 1)
# }
# save_plot(str_c('plots/paper/supplemental/7_bootstrap_', n_it, '.pdf'),
#                 plot_grid(plotlist=plots, ncol=1), ncol=1, nrow=4, useDingbats=F)

############################## End bootstrap run ################################

#### Second Round ####
smarta_coculture = read_tsv('data/PI/coculture_CFSE_smarta_all_042319.txt')

smarta_indices = smarta_coculture %>%
    mutate(cell_total = (Peak_1 + Peak_2 + Peak_3 + Peak_4 + Peak_5 + Peak_6 + Peak_7),
          divisions = (Peak_1 * 1 / (2 ^ 1)) + 
                (Peak_2 * 2 / (2 ^ 2)) + 
                (Peak_3 * 3 / (2 ^ 3)) + 
                (Peak_4 * 4 / (2 ^ 4)) +
                (Peak_5 * 5 / (2 ^ 5)) +
                (Peak_6 * 6 / (2 ^ 6)) + 
                (Peak_7 * 7 / (2 ^ 7)),
           cell_start = (Peak_0 / 2 ^ 0) +
                (Peak_1 / (2 ^ 1)) + 
                (Peak_2 / (2 ^ 2)) + 
                (Peak_3 / (2 ^ 3)) + 
                (Peak_4 / (2 ^ 4)) +
                (Peak_5 / (2 ^ 5)) +
                (Peak_6 / (2 ^ 6)) + 
                (Peak_7 / (2 ^ 7)),
           cell_activated = cell_start - Peak_0,
           PI = divisions / cell_activated,
           DI = divisions/cell_start,
           EI = cell_total/cell_start,
           RI = cell_activated/cell_start) %>%
    filter(Treatment != 'Ova') %>%
    group_by(Screen_rep_ID) %>%
    mutate(PI_norm = PI/max(PI, na.rm = T),
           DI_norm = DI/max(DI, na.rm = T),
           EI_norm = EI/max(EI, na.rm = T),
           RI_norm = RI/max(RI, na.rm = T)) %>%
    ungroup()

smarta_indices = cbind(smarta_indices, values_from_key(smarta_indices, 'PI_norm', component_key, T),
                   values_from_key(smarta_indices, 'DI_norm', component_key, T),
                   values_from_key(smarta_indices, 'EI_norm', component_key, T),
                   values_from_key(smarta_indices, 'RI_norm', component_key, T)) %>%
    # calculate isserliss values for matrices
    mutate(PI_iss = isserlis(., 'PI_norm'),
           DI_iss = isserlis(., 'DI_norm'),
           EI_iss = isserlis(., 'EI_norm'),
           RI_iss = isserlis(., 'RI_norm')) %>%
    # calculate interaction values for PI
    make_interaction() %>% 
    #group for summarization
    group_by(Treatment, Treatment_type)

# calculate mean and se for each index 
smarta_summaries = smarta_indices %>%
    summarize_at(vars(PI_norm, DI_norm, EI_norm, RI_norm, PI_iss:I_iss), funs(mean(., na.rm=T), se))

# filter for just triplets
smarta_triplets = filter(smarta_summaries, Treatment_type == 'T')

PI_plot = dot_plot(smarta_triplets, 'PI', T)
interaction_plot = dot_plot(rename(smarta_triplets, I_norm_mean = I_123_mean, I_norm_se = I_123_se), 'I', F)
save_plot('plots/paper/supplemental/8ab.pdf', 
          plot_grid(PI_plot + geom_text(aes(label=Treatment), size=1, hjust=-.2, vjust=-.5), 
                    interaction_plot + geom_text(aes(label=Treatment), size=1, hjust=-.2, vjust=-.5), nrow=1), 
          useDingbats = F, base_width = 5, base_height = 3)
