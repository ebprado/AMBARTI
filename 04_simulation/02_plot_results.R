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


object = bayesian_AMMI

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
estimate   = object$BUGSoutput$sims.matrix
mu_hat     = estimate[-seq_burn,colnames(estimate)[grepl('mu_all', colnames(estimate))]]
alpha_hat  = estimate[-seq_burn,colnames(estimate)[grepl('alpha',  colnames(estimate))]]
beta_hat   = estimate[-seq_burn,colnames(estimate)[grepl('beta',   colnames(estimate))]]
delta_hat  = estimate[-seq_burn,colnames(estimate)[grepl('delta',  colnames(estimate))]]
gamma_hat  = estimate[-seq_burn,colnames(estimate)[grepl('gamma',  colnames(estimate))]]
lambda_hat = estimate[-seq_burn,colnames(estimate)[grepl('lambda', colnames(estimate))]]
lambda_hat = as.matrix(lambda_hat)

# Alpha ------------
alpha_hat2 = melt(alpha_hat)
names(alpha_hat2) <- c('id', 'Parameter', 'value')
alpha_hat2$true = rep(data$alpha, each=npost)
alpha_hat2 %>%
  ggplot(aes(x=Parameter, y=value, group = Parameter, colour=Parameter)) +
  geom_boxplot() +
  geom_point(data=alpha_hat2, aes(y=true), shape=4) +
  labs(title = expression(~'Comparison of '~alpha[i])) +
  theme(axis.title = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5)) +
  scale_x_discrete(labels = c('alpha[1]' = expression(alpha[1]),
                              'alpha[2]' = expression(alpha[2]),
                              'alpha[3]' = expression(alpha[3]),
                              'alpha[4]' = expression(alpha[4]),
                              'alpha[5]' = expression(alpha[5]),
                              'alpha[6]' = expression(alpha[6]),
                              'alpha[7]' = expression(alpha[7]),
                              'alpha[8]' = expression(alpha[8]),
                              'alpha[9]' = expression(alpha[9]),
                              'alpha[10]' = expression(alpha[10]))) +
  scale_color_discrete(labels = c('alpha[1]' = expression(alpha[1]),
                                'alpha[2]' = expression(alpha[2]),
                                'alpha[3]' = expression(alpha[3]),
                                'alpha[4]' = expression(alpha[4]),
                                'alpha[5]' = expression(alpha[5]),
                                'alpha[6]' = expression(alpha[6]),
                                'alpha[7]' = expression(alpha[7]),
                                'alpha[8]' = expression(alpha[8]),
                                'alpha[9]' = expression(alpha[9]),
                                'alpha[10]' = expression(alpha[10])))







mu_pos <- MCMCchains(bayesian_AMMI, params = 'mu_all') %>% data.frame()
alpha_pos <- MCMCchains(bayesian_AMMI, params = 'alpha') %>% data.frame()
beta_pos <- MCMCchains(bayesian_AMMI, params = 'beta') %>% data.frame()
lambda_pos <- MCMCchains(bayesian_AMMI, params = 'lambda') %>% data.frame()
gamma_pos <- MCMCchains(bayesian_AMMI, params = 'gamma') %>% data.frame()
delta_pos <- MCMCchains(bayesian_AMMI, params = 'delta') %>% data.frame()
sigma_E_pos <- MCMCchains(bayesian_AMMI, params = 'sigma_E') %>% data.frame()

p_mu_all <- ggplot(data = melt(mu_pos), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable),fill = "gray93", color = "black") +
  geom_point(data = data.frame(x = names(mu_pos), y = 100), aes(x=x, y = y), color = 'red') +
  theme(axis.title = element_blank())

p_alpha <- ggplot(data = melt(alpha_pos), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable),fill = "gray93", color = "black") +
  geom_point(data = data.frame(x = names(alpha_pos), y = alpha), aes(x=x, y = y), color = 'red') +
  theme(axis.title = element_blank())

p_beta <- ggplot(data = melt(beta_pos), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable),fill = "gray93", color = "black") +
  geom_point(data = data.frame(x = names(beta_pos), y = beta), aes(x=x, y = y), color = 'red') +
  theme(axis.title = element_blank())

p_lambda <- ggplot(data = melt(lambda_pos), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable),fill = "gray93", color = "black") +
  geom_point(data = data.frame(x = names(lambda_pos), y = lambda), aes(x=x, y = y), color = 'red') +
  theme(axis.title = element_blank())

p_gamma <- ggplot(data = melt(gamma_pos), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable),fill = "gray93", color = "black") +
  geom_point(data = data.frame(x = names(gamma_pos), y = gamma), aes(x=x, y = y), color = 'red') +
  theme(axis.title = element_blank())

p_delta <- ggplot(data = melt(delta_pos), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable),fill = "gray93", color = "black") +
  geom_point(data = data.frame(x = names(delta_pos), y = delta), aes(x=x, y = y), color = 'red') +
  theme(axis.title = element_blank())

p_sigma_E <- ggplot(data = melt(sigma_E_pos), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable),fill = "gray93", color = "black") +
  geom_point(data = data.frame(x = names(sigma_E_pos), y = sigma_E), aes(x=x, y = y), color = 'red') +
  theme(axis.title = element_blank())

all_par <- list(
  p_alpha,
  p_beta,
  p_lambda,
  p_gamma,
  p_delta,
  p_sigma_E
)
