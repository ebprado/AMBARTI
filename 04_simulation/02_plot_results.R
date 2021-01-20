library(ggplot2)

save_file = "/Users/estevaoprado/Documents/GitHub/AMBARTI/04_simulation/results/"
load(paste(save_file, '00_results_consolidated.RData',    sep=''))

tab = save_results
tab$id = factor(tab$id, levels = c('classical AMMI', 'Bayesian AMMI (postproc)', 'Bayesian AMMI (NO postproc)', 'AMBARTI'),
                        labels = c('AMMI', ' B-AMMI (postproc)', 'B-AMMI (no postproc)', 'AMBARTI'))

tab$I = as.factor(tab$I)
tab$J = as.factor(tab$J)
tab$lambda = as.factor(tab$lambda)
tab$sa = as.factor(tab$sa)
tab$sb = as.factor(tab$sb)
tab$sy = as.factor(tab$sy)

tab %>%
  ggplot(aes(x = id, y=rmse_test, colour=algorithm2)) +
  geom_boxplot(outlier.shape = 1) +
  labs(#title = '(a)',
    # title = paste(varA, sep=''),
    # subtitle = 'Descriptive statistics',
    colour = '',
    y = 'RMSE') +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        legend.position = 'bottom',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # facet_grid(. ~ num_trees_class, scales='free')
  # facet_grid(.~sample_size) +
  facet_wrap(~ ncov, scales='free', nrow=1) +
  labs(colour='') +
  guides(col = guide_legend(nrow=2))