library(ggplot2)
library(tidyverse)

# Where the file containing the consolidated results is
save_file = "/Users/estevaoprado/Documents/GitHub/AMBARTI/04_simulation/01_simulation_AMMI/results/"

# Load the consolidated results
load(paste(save_file, '00_results_consolidated.RData', sep=''))

# Some preprocessing
tab = save_results
tab$id = factor(tab$id, levels = c('AMBARTI', 'classical AMMI', 'Bayesian AMMI (NO postproc)', 'Bayesian AMMI (postproc)'),
                        labels = c('AMBARTI', 'AMMI', 'B-AMMI (no PP)', 'B-AMMI (PP)'))
tab$Q = factor(tab$Q, levels = c('1','2','3'), labels = c('Q = 1', 'Q = 2', 'Q = 3'))
tab$I = factor(tab$I, levels=c('10'), labels=c('I = 10'))
tab$J = factor(tab$J, levels=c('10'), labels=c('J = 10'))
tab$lambda = as.factor(tab$lambda)
tab$sa = factor(tab$sa, levels=c('1','5'), labels=c(expression(paste(sigma[alpha],' = 1')), expression(paste(sigma[alpha],' = 5'))))
tab$sb = factor(tab$sb, levels=c('1','5'), labels=c(expression(paste(sigma[beta],' = 1')), expression(paste(sigma[beta],' = 5'))))
tab$sy = as.factor(tab$sy)

# Generate plots
save_plots = "/Users/estevaoprado/Documents/GitHub/AMBARTI/04_simulation/01_simulation_AMMI/"
myplot <- function(varA, varB, varC){

  pdf(paste(save_plots, varB, varA, '.pdf', sep=''), width = 8, height = 6)
  if (varA == 'sa'){aux = expression(sigma[alpha])}
  if (varA == 'sb'){aux = expression(sigma[beta])}

  if(varB=='rrmse_alpha'){varC = expression("RRMSE - "~alpha[i])}
  if(varB=='rrmse_beta'){varC = expression("RRMSE - "~beta[j])}
  if(varB=='lambda_rrmse'){varC = expression("RRMSE - "~lambda[q])}
  if(varB=='gamma_rrmse'){varC = expression("RRMSE - "~gamma[iq])}
  if(varB=='delta_rrmse'){varC = expression("RRMSE - "~delta[jq])}

  if (varB %in% c('rrmse_alpha','rrmse_beta','lambda_rrmse','gamma_rrmse','delta_rrmse')){aux_y = 'RRMSE'}
  if (varB %in% c('rmse_blinear','y_test_rmse','y_train_rmse')){aux_y = 'RMSE'}

  xxx = tab %>%
    group_by(Q) %>%
    filter((!!sym(varB)) < quantile(!!sym(varB), 0.9, na.rm = TRUE)) %>%
    ggplot(aes_string(x = varA , y = varB, colour='id')) +
    geom_boxplot(outlier.shape = 1) +
    labs(x = aux,
      title = varC,
      colour = '',
      y = aux_y) +
    scale_x_discrete(labels = parse(text = levels(tab$sa))) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 15, hjust = 0.5),
          strip.text.y = element_text(angle = 0),
          legend.position = 'bottom',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    # facet_wrap(~ Q, scales='free_x', nrow=1) +
    # facet_wrap(Q ~ sb, scales='free_x', nrow=2) +
    facet_grid(sb ~ Q, labeller = labeller(sb = label_parsed, Q=label_value)) +
    labs(colour='') +
    guides(col = guide_legend(nrow=2))

  print(xxx)

  dev.off()
}

myplot('sa','rrmse_alpha', 'RRMSE - alpha')
myplot('sa','rrmse_beta',  'RRMSE - beta')
myplot('sa','lambda_rrmse', 'RRMSE - lambda')
# myplot('sa','lambda_rrmse', 'RRMSE - lambda')
myplot('sa','gamma_rrmse', 'RRMSE - gamma')
# myplot('sa','gamma_rrmse', 'RRMSE - gamma')
myplot('sa','delta_rrmse', 'RRMSE - delta')
# myplot('sa','delta_rrmse', 'RRMSE - delta')
myplot('sa','rmse_blinear', 'RMSE - Interaction')
# myplot('sb','rmse_blinear', 'RMSE - Bilinear part')
myplot('sa','y_test_rmse', 'RMSE - y test')
# myplot('sa','y_test_rmse', 'RMSE - y test')

# Check -----------
# they're quite similar, but still different
myplot('sa','lambda_rrmse', 'RRMSE - lambda')

tab %>%
  filter(id %in% c('AMMI', 'B-AMMI (postproc)'), Q == 'Q = 1', sa == '1') %>%
  select(id,
         lambda_rrmse,
         sb)

myplot('sa','rmse_blinear', 'RMSE - Bilinear part')
tab %>%
  filter(id %in% c('AMMI', 'B-AMMI (postproc)'), Q == 'Q = 1', sa == '1') %>%
  select(id,
         rmse_blinear,
         sb)
