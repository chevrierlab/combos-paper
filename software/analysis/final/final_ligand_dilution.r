### libraries ###
library(tidyverse)
library(cowplot)

### sources ###
source('software/index_functions.r')

# read in dilution information
dilutions = read_tsv('data/Ligand_dilution_all.txt') %>%
    filter(!is.na(Treatment), Treatment != 'Ctrl', Exp_ID == 3) %>%
    mutate(Viability = Live_cells_number / All_cells_number * 100,
           Ligand_concentration = as.numeric(str_replace(Ligand_concentration, '_', '.'))) %>%
    select(Treatment, Ligand_concentration, Proliferation_percent, Viability) %>%
    group_by(Treatment, Ligand_concentration)

# values selected for final use
reds = c(PAM = 250, LPS = 100, Zymosan = 30, SeV = 10, cGAMP = 20, HMW = 20, CpGB = 10)

# calculate means and standard errors
dilution_sum = dilutions %>%
    summarise(prolif_mean = mean(Proliferation_percent), prolif_se = se(Proliferation_percent), 
              viability_mean = mean(Viability), viability_se = se(Viability)) %>%
    mutate(prolif_top = prolif_mean + prolif_se, prolif_bot = prolif_mean - prolif_se,
           viability_top = viability_mean + viability_se, viability_bot = viability_mean - viability_se,
           red = Ligand_concentration == reds[Treatment])

# dot plots
# mostly uses logorithmic fitting, a few we decided on linear or cubic
dot_plots = map(unique(dilution_sum$Treatment), 
                function(treat) filter(dilution_sum, Treatment == treat) %>%
                    ggplot(aes(Ligand_concentration, prolif_mean)) +
                    geom_smooth(formula = (if (treat == 'Zymosan') y ~ poly(x, 3) else if (treat == 'HMW') y ~ x else y ~ log(x)), 
                                           method='lm', se=F, col='cadetblue4') + 
                    geom_pointrange(aes(ymin=prolif_bot, ymax=prolif_top, col = red), size=1) + 
                    scale_color_manual(values = c('black', 'red')) +
                    guides(col = 'none'))

# bar plots
bar_plots = map(unique(dilution_sum$Treatment),
                function(treat) filter(dilution_sum, Treatment == treat) %>%
                    ggplot(aes(factor(Ligand_concentration), viability_mean)) +
                    geom_col(aes(fill = red), width = .75) + 
                    geom_errorbar(aes(ymin = viability_bot, max = viability_top), width = .5) + 
                    scale_fill_manual(values = c('black', 'red')) +
                    guides(fill='none'))

zipped_plots = map2(dot_plots, bar_plots, function(d, b) plot_grid(d, b, nrow=1))

save_plot('plots/paper/supplemental/2b.pdf',
          plot_grid(plotlist=zipped_plots, ncol=1, labels = unique(dilution_sum$Treatment)), nrow=7, ncol=2)
