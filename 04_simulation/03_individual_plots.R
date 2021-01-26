save_file = "/Users/estevaoprado/Documents/GitHub/AMBARTI/04_simulation/results/"
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

AMMI_help_plot <- function(object, data){

  # Get training info

  x_train = data$x
  y_train = data$y

  # Get test info
  x_test = data$x
  y_test = data$y_test

  g = as.factor(x_train[,"g"])
  e = as.factor(x_train[,"e"])

  # Fit the linear model
  linear_mod = aov(y_train ~ g + e + g:e)

  aux_mod = lm(y_train ~ g + e)
  aux_mod_sd = summary(aux_mod)$coefficients[,'Std. Error'][-1]
  alpha_sd = max(aux_mod_sd[grepl('g', names(aux_mod_sd))])# if there's no repetition, they'll be the same
  beta_sd = max(aux_mod_sd[grepl('e', names(aux_mod_sd))])

  # Get the residuals for the interaction gen:env
  #interaction_tab = matrix(residuals(linear_mod), ncol = length(unique(gen)), length(unique(env)))
  interaction_tab = model.tables(linear_mod, type='effects', cterms = 'g:e')
  interaction_tab = interaction_tab$tables$`g:e`

  # Get the number of PCs
  PC = length(data$lambda)

  # Run the Singular Value Decomposition (SVD) to compute lambda, gamma, and delta
  sv_dec <- svd(interaction_tab, nu = PC, nv = PC)

  # Get parameter estimates
  # mu_hat     = linear_mod$coefficients[1] # slightly biased compared to mean(y_train)
  mu_hat     = mean(y_train)
  alpha_hat  = aggregate(x = y_train - mu_hat, by = list(g), FUN = "mean")[,2]
  beta_hat   = aggregate(x = y_train - mu_hat, by = list(e), FUN = "mean")[,2]
  lambda_hat = sv_dec$d[1:PC]
  gamma_hat  = -1*sv_dec$u
  delta_hat  = -1*sv_dec$v

  # Set up the bilinear term
  blin_train = rep(0, length(y_train))
  blin_test  = rep(0, length(y_test))

  for (k in 1:PC) {
    blin_train = blin_train + lambda_hat[k]*gamma_hat[x_train[,'g'],k]*delta_hat[x_train[,'e'],k]
    blin_test  = blin_test + lambda_hat[k]*gamma_hat[x_test[,'g'],k]*delta_hat[x_test[,'e'],k]
  }

  # Compute the predicted response for the TRAINING data
  y_hat_train = mu_hat + alpha_hat[x_train[,'g']] + beta_hat[x_train[,'e']] + blin_train

  # Compute the predicted response for the TEST data
  y_hat_test = mu_hat + alpha_hat[x_test[,'g']] + beta_hat[x_test[,'e']] + blin_test

  #
  g_names = rownames(interaction_tab)
  e_names = colnames(interaction_tab)

  alpha_hat  = data.frame(id = 'AMMI', variable = paste('alpha[',g_names ,']', sep=''), value=alpha_hat, true=data$alpha)
  beta_hat   = data.frame(id = 'AMMI', variable = paste('beta[',e_names ,']', sep=''), value=beta_hat, true=data$beta)
  lambda_hat = data.frame(id = 'AMMI', variable = paste('lambda[',1:PC,']',sep=''), value=lambda_hat, true=data$lambda)
  gamma_hat  = data.frame(id = 'AMMI', variable = paste('gamma[',g_names,',',rep(1:PC, each=length(g_names)),']',sep=''), value=as.numeric(gamma_hat), true=as.numeric(data$gamma))
  delta_hat  = data.frame(id = 'AMMI', variable = paste('delta[',g_names,',',rep(1:PC, each=length(e_names)),']',sep=''), value=as.numeric(delta_hat), true=as.numeric(data$delta))

  return(list(alpha_hat   = alpha_hat,
              beta_hat    = beta_hat,
              delta_hat   = delta_hat,
              gamma_hat   = gamma_hat,
              lambda_hat  = lambda_hat,
              y_hat_train = y_hat_train,
              y_hat_test  = y_hat_test,
              blinear_hat = blin_train,
              id          = 'classical AMMI'))
}

res_AMMI = AMMI_help_plot(classical_AMMI, data)

alpha_hat  = res_AMMI$alpha_hat
beta_hat   = res_AMMI$beta_hat
lambda_hat = res_AMMI$lambda_hat
gamma_hat  = res_AMMI$gamma_hat
delta_hat  = res_AMMI$delta_hat

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

alpha_hat  = rbind(alpha_hat,  bAMMI_save_info$alpha_hat)
beta_hat   = rbind(beta_hat,   bAMMI_save_info$beta_hat)
lambda_hat = rbind(lambda_hat, bAMMI_save_info$lambda_hat)
gamma_hat  = rbind(gamma_hat,  bAMMI_save_info$gamma_hat)
delta_hat  = rbind(delta_hat,  bAMMI_save_info$delta_hat)

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
