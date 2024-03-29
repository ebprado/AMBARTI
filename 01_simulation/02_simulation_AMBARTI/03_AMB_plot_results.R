library(ggplot2)
library(tidyverse)

# Where the file containing the consolidated results is
save_file = "~/R/AMBARTI/02_simulation_AMBARTI/results/"

# Load the consolidated results
load(paste(save_file, '101000_AMB_results_consolidated.RData', sep='')); save_results10 = save_results
load(paste(save_file, '252500_AMB_results_consolidated.RData', sep='')); save_results25 = save_results

# Some preprocessing
# tab = save_results10
tab = rbind(save_results10, save_results25)
tab$id = factor(tab$id, levels = c('classical AMMI', 'Bayesian AMMI (postproc)', 'Bayesian AMMI (NO postproc)', 'AMBARTI', 'AMBARTI (PP)'),
                        labels = c('AMMI', 'B-AMMI', 'B-AMMI (No PP)', 'AMBARTI', 'AMBARTI (PP)'))

tab$Q = factor(tab$Q, levels=c('1', '2', '3'), labels=c('Q=1', 'Q=2', 'Q=3'))
tab$id = paste(tab$id,' (',tab$Q, ')', sep='')
tab$id = gsub(' (NA)','', tab$id, fixed = TRUE)
tab$id = factor(tab$id, levels=sort(unique(tab$id)), labels = sort(unique(tab$id)))
tab$I = factor(tab$I, levels=c('10', '25', '50'), labels=c('I = 10', 'I = 25', 'I = 50'))
tab$J = factor(tab$J, levels=c('10', '25', '50'), labels=c('J = 10', 'I = 25', 'I = 50'))
tab$lambda = as.factor(tab$lambda)
tab$sa = factor(tab$sa, levels=c('1','5'), labels=c(expression(paste(sigma[g],' = 1')), expression(paste(sigma[g],' = 5'))))
                tab$sb = factor(tab$sb, levels=c('1','5'), labels=c(expression(paste(sigma[e],' = 1')), expression(paste(sigma[e],' = 5'))))
tab$sy = as.factor(tab$sy)

tab = tab %>% filter(id %in% c('AMMI (Q=1)', 'AMMI (Q=2)', 'AMMI (Q=3)','B-AMMI (Q=1)', 'B-AMMI (Q=2)', 'B-AMMI (Q=3)', 'AMBARTI'), I != 'I = 50')

# Generate plots
save_plots = "~/R/AMBARTI/02_simulation_AMBARTI/"
myplot <- function(varA, varB, varC){

  pdf(paste(save_plots, 'new_AMBARTI_', varB, varA, '.pdf', sep=''), width = 9, height = 6)
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
    filter((!!sym(varB)) < quantile(!!sym(varB), 0.95, na.rm = TRUE)) %>%
    ggplot(aes_string(x = varA , y = varB, colour='id')) +
    geom_boxplot(outlier.shape = 1) +
    labs(
      x = '',
      # x = aux,
      # title = varC,
      colour = '',
      y = aux_y) +
    scale_x_discrete(labels = parse(text = levels(tab$sa))) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 15, hjust = 0.5),
          legend.position = 'bottom',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    # facet_grid(~ sa, scales='free_x', nrow=1,labeller = label_parsed) +
    facet_grid(I ~ sb, scales='free_x', labeller = labeller(sb=label_parsed)) +
    labs(colour='') +
    guides(col = guide_legend(nrow=1))

  print(xxx)

  dev.off()
}

myplot('sa','rrmse_g', 'RRMSE - g')
myplot('sa','rrmse_e',  'RRMSE - e')
myplot('sa','rmse_blinear', 'RMSE - Interaction')
#myplot('sa','y_train_rmse', 'RMSE - y train')
myplot('sa','y_test_rmse', 'RMSE - y test')

