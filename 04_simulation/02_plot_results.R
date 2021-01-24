library(ggplot2)
library(tidyverse)

# Where the file containing the consolidated results is
save_file = "/Users/estevaoprado/Documents/GitHub/AMBARTI/04_simulation/results/"

# Load the consolidated results
load(paste(save_file, '00_results_consolidated.RData', sep=''))

# Some preprocessing
tab = save_results
tab$id = factor(tab$id, levels = c('classical AMMI', 'Bayesian AMMI (postproc)', 'Bayesian AMMI (NO postproc)', 'AMBARTI'),
                        labels = c('AMMI', 'B-AMMI (postproc)', 'B-AMMI (no postproc)', 'AMBARTI'))

tab$Q = factor(sapply(strsplit(as.character(tab$lambda), split = ' '), function(x) length(x)), levels = c('1','2','3'), labels = c('Q = 1', 'Q = 2', 'Q = 3'))
tab$I = factor(tab$I, levels=c('10'), labels=c('I = 10'))
tab$J = factor(tab$J, levels=c('10'), labels=c('J = 10'))
tab$lambda = as.factor(tab$lambda)
tab$sa = as.factor(tab$sa)
tab$sb = as.factor(tab$sb)
tab$sy = as.factor(tab$sy)

# Generate plots
myplot <- function(varA, varB, varC){

  if (varA == 'sa'){aux = expression(sigma[alpha])}
  if (varA == 'sb'){aux = expression(sigma[beta])}

  if(varB=='rrmse_alpha'){varC = expression("RRMSE - "~alpha[i])}
  if(varB=='rrmse_beta'){varC = expression("RRMSE - "~beta[j])}
  if(varB=='lambda_rrmse'){varC = expression("RRMSE - "~lambda[q])}
  if(varB=='gamma_rrmse'){varC = expression("RRMSE - "~gamma[iq])}
  if(varB=='delta_rrmse'){varC = expression("RRMSE - "~delta[jq])}

  tab %>%
    group_by(Q) %>%
    filter((!!sym(varB)) < quantile(!!sym(varB), 0.9, na.rm = TRUE)) %>%
    ggplot(aes_string(x = varA , y = varB, colour='id')) +
    geom_boxplot(outlier.shape = 1) +
    labs(x = aux,
      title = varC,
      colour = '',
      y = 'RMSE') +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 15, hjust = 0.5),
          legend.position = 'bottom',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    facet_wrap(~ Q, scales='free', nrow=1) +
    # facet_grid(I ~ Q) +
    labs(colour='') +
    guides(col = guide_legend(nrow=2))

}

myplot('sa','rrmse_alpha', 'RRMSE - alpha')
myplot('sb','rrmse_beta',  'RRMSE - beta')
myplot('sb','lambda_rrmse', 'RRMSE - lambda')
# myplot('sa','lambda_rrmse', 'RRMSE - lambda')
myplot('sb','gamma_rrmse', 'RRMSE - gamma')
# myplot('sa','gamma_rrmse', 'RRMSE - gamma')
myplot('sb','delta_rrmse', 'RRMSE - delta')
# myplot('sa','delta_rrmse', 'RRMSE - delta')
myplot('sa','rmse_blinear', 'RMSE - Bilinear part')
# myplot('sb','rmse_blinear', 'RMSE - Bilinear part')
myplot('sb','y_test_rmse', 'RMSE - y test')
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

#-------------------------------------------------
# Individual plots
#-------------------------------------------------
library(reshape2)
setwd(save_file)
filename = 'I10J10sa1sb1sy1L1012r10'

data_filename    = paste(filename, '_data.RData',           sep='')
ammi_filename    = paste(filename, '_classical_AMMI.RData', sep='')
bammi_filename   = paste(filename, '_bayesian_AMMI.RData',  sep='')
ambarti_filename = paste(filename, '_AMBARTI.RData',        sep='')

load(paste(save_file, data_filename,      sep=''))
load(paste(save_file, ammi_filename,      sep=''))
load(paste(save_file, bammi_filename,     sep=''))
load(paste(save_file, ambarti_filename,   sep=''))

# Get parameter estimates from the classical AMMI
res_AMMI = organise_classical_AMMI(classical_AMMI, data)



bAMMI_help_plot <- function(object, data){
  # Get training info
  x_train = data$x
  y_train = data$y

  # Get test info
  x_test = data$x
  y_test = data$y_test

  # Get the number of PCs
  PC = length(data$lambda)

  # Get some MCMC info
  nburn      = object$BUGSoutput$n.burnin
  niter      = object$BUGSoutput$n.iter
  npost      = niter - nburn
  seq_burn   = seq(1, nburn, by=1)

  # Get estimates info
  estimate   = as.data.frame(object$BUGSoutput$sims.matrix)
  estimate$id= 'Bayesian AMMI'
  #mu_hat     = estimate[-seq_burn,c('id',colnames(estimate)[grepl('mu_all', colnames(estimate))])]
  alpha_hat  = estimate[-seq_burn,c('id',colnames(estimate)[grepl('alpha', colnames(estimate))])]
  beta_hat   = estimate[-seq_burn,c('id',colnames(estimate)[grepl('beta', colnames(estimate))])]
  delta_hat  = estimate[-seq_burn,c('id',colnames(estimate)[grepl('delta', colnames(estimate))])]
  gamma_hat  = estimate[-seq_burn,c('id',colnames(estimate)[grepl('gamma', colnames(estimate))])]
  lambda_hat = estimate[-seq_burn,c('id',colnames(estimate)[grepl('lambda', colnames(estimate))])]

  alpha_hat  = melt(alpha_hat, measure.vars = colnames(alpha_hat)[-1])
  beta_hat   = melt(beta_hat, measure.vars = colnames(beta_hat)[-1])
  delta_hat  = melt(delta_hat, measure.vars = colnames(delta_hat)[-1])
  gamma_hat  = melt(gamma_hat, measure.vars = colnames(gamma_hat)[-1])
  lambda_hat = melt(lambda_hat, measure.vars = colnames(lambda_hat)[-1])

  alpha_hat$true  = rep(as.numeric(data[['alpha']]),  each=npost)
  beta_hat$true   = rep(as.numeric(data[['beta']]),   each=npost)
  delta_hat$true  = rep(as.numeric(data[['delta']]),  each=npost)
  gamma_hat$true  = rep(as.numeric(data[['gamma']]),  each=npost)
  lambda_hat$true = rep(as.numeric(data[['lambda']]), each=npost)

  return(list(alpha_hat   = alpha_hat,
              beta_hat    = beta_hat,
              delta_hat   = delta_hat,
              gamma_hat   = gamma_hat,
              lambda_hat  = lambda_hat))
}

bAMMI_save_info = bAMMI_help_plot(bayesian_AMMI, data)

alpha_hat = bAMMI_save_info$alpha_hat
beta_hat = bAMMI_save_info$beta_hat
lambda_hat = bAMMI_save_info$lambda_hat
gamma_hat = bAMMI_save_info$gamma_hat
delta_hat = bAMMI_save_info$delta_hat

AMBARTI_help_plot <- function(object, data){

  # Get training info
  x_train = data$x
  y_train = data$y

  # Get test info
  x_test = data$x
  x_test$g = as.factor(x_test$g)
  x_test$e = as.factor(x_test$e)
  y_test = data$y_test

  # Get estimates info
  estimate = as.data.frame(object$beta_hat)
  alpha_hat = estimate[,grepl('g', names(estimate))]
  beta_hat  = estimate[,grepl('e', names(estimate))]
  names(alpha_hat) <- paste(gsub('g','alpha[', names(alpha_hat)), ']', sep='')
  names(beta_hat) <- paste(gsub('e','beta[', names(beta_hat)), ']', sep='')
  alpha_hat$id = 'AMBARTI'
  beta_hat$id = 'AMBARTI'

  alpha_hat = melt(alpha_hat, measure.vars = colnames(alpha_hat)[grepl('alpha', colnames(alpha_hat))])
  beta_hat = melt(beta_hat, measure.vars = colnames(beta_hat)[grepl('beta', colnames(beta_hat))])

  alpha_hat$true = rep(as.numeric(data[['alpha']]), each=ambarti$npost)
  beta_hat$true  = rep(as.numeric(data[['beta']]), each=ambarti$npost)

  return(list(alpha_hat   = alpha_hat,
              beta_hat    = beta_hat))
}

AMBARTI_save_info = AMBARTI_help_plot(ambarti, data)

alpha_hat = rbind(alpha_hat, AMBARTI_save_info$alpha_hat)
beta_hat = rbind(beta_hat, AMBARTI_save_info$beta_hat)

new_parse_format <- function(text) {
  out <- vector("expression", length(text))
  for (i in seq_along(text)) {
    expr <- parse(text = text[[i]])
    out[[i]] <- if (length(expr) == 0)
      NA
    else expr[[1]]
  }
  out
}

plot_individual <- function(object, data){

  db = object
  names(db) = c('Method', 'Parameter', 'value', 'true')
  aux_name = strsplit(deparse(substitute(object)), split = '_')[[1]][1]
  orig_labels = as.character(unique(object$variable))
  fixed_labels = new_parse_format(gsub(',','', orig_labels))

  db %>%
    ggplot(aes(x=Parameter, y=value)) +
    geom_boxplot(aes(colour=Method))+
    geom_point(data=db, aes(y=true), colour='darkgray', shape=4, size=4) +
    theme_bw() +
    labs(title = bquote(Comparison~of~.(sym(aux_name)))) +
    theme(axis.text.x = element_text(size=12),
          axis.title = element_blank(),
          plot.title = element_text(size = 18, hjust = 0.5),
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.position = "bottom") +
    scale_x_discrete(limit = orig_labels,
                     labels = fixed_labels) +
    scale_color_discrete(labels=unique(db$Method))
}
plot_individual(alpha_hat, data)
plot_individual(beta_hat, data)
plot_individual(lambda_hat, data)
plot_individual(gamma_hat, data)
plot_individual(delta_hat, data)
