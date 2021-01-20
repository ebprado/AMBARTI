library(ggplot2)
library(tidyverse)

save_file = "/Users/estevaoprado/Documents/GitHub/AMBARTI/04_simulation/results/"
load(paste(save_file, '00_results_consolidated.RData',    sep=''))

tab = save_results
tab$id = factor(tab$id, levels = c('classical AMMI', 'Bayesian AMMI (postproc)', 'Bayesian AMMI (NO postproc)', 'AMBARTI'),
                        labels = c('AMMI', ' B-AMMI (postproc)', 'B-AMMI (no postproc)', 'AMBARTI'))

tab$Q = factor(sapply(strsplit(as.character(tab$lambda), split = ' '), function(x) length(x)), levels = c('1','2','3'), labels = c('Q = 1', 'Q = 2', 'Q = 3'))
tab$I = as.factor(tab$I)
tab$J = as.factor(tab$J)
tab$lambda = as.factor(tab$lambda)
tab$sa = as.factor(tab$sa)
tab$sb = as.factor(tab$sb)
tab$sy = as.factor(tab$sy)

tab2$lambda <- factor(tab2$lambda, levels = c('50'),
                      labels = c(expression(paste(lambda, '= 50'))))

myplot <- function(varA, varB){

  # pdf(paste(SaveFigures,'P_by_lambda.pdf', sep = ''), width = 8, height = 6)

  tab %>%
    ggplot(aes_string(x = varA , y = varB, colour='id')) +
    geom_boxplot(outlier.shape = 1) +
    labs(
      title = paste(varB, sep=''),
      # subtitle = 'Descriptive statistics',
      colour = '',
      y = 'RMSE') +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 15, hjust = 0.5),
          legend.position = 'bottom',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    # facet_grid(. ~ num_trees_class, scales='free')
    # facet_grid(.~Q + sa) +
    facet_wrap(~ Q, scales='free', nrow=1) +
    labs(colour='') +
    guides(col = guide_legend(nrow=2))

  # dev.off()

}
myplot('sb','y_test_rmse')
myplot('sa','y_test_rmse')
myplot('sa','rrmse_alpha')
myplot('sb','rrmse_beta')
myplot('sb','lambda_rrmse')
myplot('sa','lambda_rrmse')
myplot('sb','gamma_rrmse')
myplot('sa','gamma_rrmse')
myplot('sb','delta_rrmse')
myplot('sa','delta_rrmse')

mean(apply(ambarti$beta_hat, 2, mean))
