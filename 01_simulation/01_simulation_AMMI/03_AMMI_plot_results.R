library(ggplot2)
library(tidyverse)

# Where the file containing the consolidated results is
save_file = "~/R/AMBARTI/01_simulation_AMMI/results/"

# Load the consolidated results
load(paste(save_file, '00_results_consolidated.RData', sep=''))

# Some preprocessing
tab = save_results
tab$id = factor(tab$id, levels = c('AMBARTI', 'classical AMMI', 'Bayesian AMMI (NO postproc)', 'Bayesian AMMI (postproc)'),
                        labels = c('AMBARTI', 'AMMI', 'B-AMMI (no PP)', 'B-AMMI (PP)'))
tab$Q = factor(tab$Q, levels = c('1','2','3'), labels = c('Q = 1', 'Q = 2', 'Q = 3'))
tab$I = factor(tab$I, levels=c('10', '25', '50'), labels=c('I = 10', 'I = 25', 'I = 50'))
tab$J = factor(tab$J, levels=c('10', '25', '50'), labels=c('J = 10', 'J = 25', 'I = 50'))
tab$lambda = as.factor(tab$lambda)
tab$sa = factor(tab$sa, levels=c('1','5'), labels=c(expression(paste(sigma[g],' = 1')), expression(paste(sigma[g],' = 5'))))
tab$sb = factor(tab$sb, levels=c('1','5'), labels=c(expression(paste(sigma[e],' = 1')), expression(paste(sigma[e],' = 5'))))
tab$sy = as.factor(tab$sy)

tab = tab %>% filter(id %in% c('AMBARTI', 'AMMI'))

# Generate plots
save_plots = "~/R/AMBARTI/01_simulation_AMMI/"
myplot <- function(varA, varB, varC, varD){

  pdf(paste(save_plots, 'AMMI_', gsub(' = ', '', varD), '_' , varB, varA, '.pdf', sep=''), width = 8, height = 6)
  if (varA == 'sa'){aux = expression(sigma[g])}
  if (varA == 'sb'){aux = expression(sigma[e])}

  if(varB=='rrmse_g'){varC = expression("RRMSE - "~g[i])}
  if(varB=='rrmse_e'){varC = expression("RRMSE - "~e[j])}
  if(varB=='lambda_rrmse'){varC = expression("RRMSE - "~lambda[q])}
  if(varB=='gamma_rrmse'){varC = expression("RRMSE - "~gamma[iq])}
  if(varB=='delta_rrmse'){varC = expression("RRMSE - "~delta[jq])}

  if (varB %in% c('rrmse_g','rrmse_e','lambda_rrmse','gamma_rrmse','delta_rrmse')){aux_y = 'RRMSE'}
  if (varB %in% c('rmse_blinear','y_test_rmse','y_train_rmse')){aux_y = 'RMSE'}

  xxx = tab %>%
    group_by(Q) %>%
    filter((!!sym(varB)) < quantile(!!sym(varB), 0.9, na.rm = TRUE)) %>%
    filter(I == varD) %>% 
    ggplot(aes_string(x = varA , y = varB, colour='id')) +
    geom_boxplot(outlier.shape = 1) +
    labs(x = aux,
      # title = varC,
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
    guides(col = guide_legend(nrow=1))

  print(xxx)

  dev.off()
}
num_gen_env = 'I = 10'
myplot('sa','rrmse_g', 'RRMSE - g', num_gen_env)
myplot('sa','rrmse_e',  'RRMSE - e', num_gen_env)
myplot('sa','lambda_rrmse', 'RRMSE - lambda', num_gen_env)
# myplot('sa','lambda_rrmse', 'RRMSE - lambda')
myplot('sa','gamma_rrmse', 'RRMSE - gamma', num_gen_env)
# myplot('sa','gamma_rrmse', 'RRMSE - gamma')
myplot('sa','delta_rrmse', 'RRMSE - delta', num_gen_env)
# myplot('sa','delta_rrmse', 'RRMSE - delta')
myplot('sa','rmse_blinear', 'RMSE - Interaction', num_gen_env)
# myplot('sb','rmse_blinear', 'RMSE - Bilinear part')
myplot('sa','y_test_rmse', 'RMSE - y test', num_gen_env)
myplot('sa','y_train_rmse', 'RMSE - y train', num_gen_env)

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
