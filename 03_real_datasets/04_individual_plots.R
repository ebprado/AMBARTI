library(devtools)
install_github("ebprado/AMBARTI/R package",
               ref = 'main',
               auth_token = '363d4ad84eaa25d3eb26752732cc208f7e698086')

library(AMBARTI)
save_file = "~/Documents/GitHub/AMBARTI/03_real_datasets/results/"

#-------------------------------------------------
# Individual plots
#-------------------------------------------------
library(reshape2)
setwd(save_file)
year = 2019
filename = paste('Ireland_VCU_',
                 year,
                 sep='')

data_filename    = paste(filename, '_data.RData',              sep='')
ammi_filenameQ1  = paste(filename, '_classical_AMMI_Q1.RData', sep='')
ammi_filenameQ2  = paste(filename, '_classical_AMMI_Q2.RData', sep='')
ammi_filenameQ3  = paste(filename, '_classical_AMMI_Q3.RData', sep='')
ambarti_filename = paste(filename, '_AMBARTI.RData',           sep='')

load(paste(save_file, data_filename,    sep=''))
load(paste(save_file, ammi_filenameQ1,  sep=''))
load(paste(save_file, ammi_filenameQ2,  sep=''))
load(paste(save_file, ammi_filenameQ3,  sep=''))
load(paste(save_file, ambarti_filename, sep=''))

alpha_hat   = NULL
beta_hat    = NULL
lambda_hat  = NULL
gamma_hat   = NULL
delta_hat   = NULL
blinear_hat = NULL
y_train_hat = NULL
y_test_hat  = NULL

# Get parameter estimates from classical AMMI
for (Q in 1:3){
  aux_AMMI    = get(paste('classical_AMMI_Q',Q, sep=''))
  res_AMMI    = AMMI_help_plot(aux_AMMI, data, Q = Q)

  alpha_hat   = rbind(alpha_hat,    res_AMMI$alpha_hat)
  beta_hat    = rbind(beta_hat,     res_AMMI$beta_hat)
  lambda_hat  = rbind(lambda_hat,   res_AMMI$lambda_hat)
  gamma_hat   = rbind(gamma_hat,    res_AMMI$gamma_hat)
  delta_hat   = rbind(delta_hat,    res_AMMI$delta_hat)
  blinear_hat = rbind(blinear_hat,  res_AMMI$blinear_hat)
  y_train_hat = rbind(y_train_hat,  res_AMMI$y_hat_train)
  y_test_hat  = rbind(y_test_hat,   res_AMMI$y_hat_test)

}

# Get parameter estimates from AMBARTI
AMBARTI_save_info = AMBARTI_help_plot(ambarti, data)

alpha_hat   = rbind(alpha_hat,   AMBARTI_save_info$alpha_hat)
beta_hat    = rbind(beta_hat,    AMBARTI_save_info$beta_hat)
blinear_hat = rbind(blinear_hat, AMBARTI_save_info$blinear_hat)
y_train_hat = rbind(y_train_hat, AMBARTI_save_info$y_hat_train)
y_test_hat  = rbind(y_test_hat,  AMBARTI_save_info$y_hat_test)

# Some postprocessing

alpha_hat$id   = factor(alpha_hat$id,   levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI (PP)', 'B-AMMI (No PP)', 'AMBARTI'))
beta_hat$id    = factor(beta_hat$id,    levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI (PP)', 'B-AMMI (No PP)', 'AMBARTI'))
lambda_hat$id  = factor(lambda_hat$id,  levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI (PP)', 'B-AMMI (No PP)', 'AMBARTI'))
gamma_hat$id   = factor(gamma_hat$id,   levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI (PP)', 'B-AMMI (No PP)', 'AMBARTI'))
delta_hat$id   = factor(delta_hat$id,   levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI (PP)', 'B-AMMI (No PP)', 'AMBARTI'))
blinear_hat$id = factor(blinear_hat$id, levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI (PP)', 'B-AMMI (No PP)', 'AMBARTI'))
y_train_hat$id = factor(y_train_hat$id, levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI (PP)', 'B-AMMI (No PP)', 'AMBARTI'))
y_test_hat$id  = factor(y_test_hat$id,  levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI (PP)', 'B-AMMI (No PP)', 'AMBARTI'))

alpha_hat$id   = ifelse(alpha_hat$id   != 'AMBARTI', paste(alpha_hat$id,   ' (Q = ', alpha_hat$Q, ')',    sep=''),as.character(alpha_hat$id))
beta_hat$id    = ifelse(beta_hat$id    != 'AMBARTI', paste(beta_hat$id,    ' (Q = ', beta_hat$Q, ')',     sep=''),as.character(beta_hat$id))
lambda_hat$id  = ifelse(lambda_hat$id  != 'AMBARTI', paste(lambda_hat$id,   ' (Q = ', lambda_hat$Q, ')',    sep=''),as.character(lambda_hat$id))
gamma_hat$id   = ifelse(gamma_hat$id   != 'AMBARTI', paste(gamma_hat$id,    ' (Q = ', gamma_hat$Q, ')',     sep=''),as.character(gamma_hat$id))
delta_hat$id   = ifelse(delta_hat$id   != 'AMBARTI', paste(delta_hat$id,    ' (Q = ', delta_hat$Q, ')',     sep=''),as.character(delta_hat$id))
blinear_hat$id = ifelse(blinear_hat$id != 'AMBARTI', paste(blinear_hat$id, ' (Q = ', blinear_hat$Q, ')',  sep=''),as.character(blinear_hat$id))
y_train_hat$id = ifelse(y_train_hat$id != 'AMBARTI', paste(y_train_hat$id, ' (Q = ', y_train_hat$Q, ')',  sep=''),as.character(y_train_hat$id))
y_test_hat$id  = ifelse(y_test_hat$id  != 'AMBARTI', paste(y_test_hat$id,  ' (Q = ', y_test_hat$Q, ')',   sep=''),as.character(y_test_hat$id))

alpha_hat$id   = factor(alpha_hat$id,   levels = sort(as.character(unique(alpha_hat$id))),   labels = sort(as.character(unique(alpha_hat$id))))
beta_hat$id    = factor(beta_hat$id,    levels = sort(as.character(unique(beta_hat$id))),    labels = sort(as.character(unique(beta_hat$id))))
blinear_hat$id = factor(blinear_hat$id, levels = sort(as.character(unique(blinear_hat$id))), labels = sort(as.character(unique(blinear_hat$id))))

lambda_hat$id = factor(lambda_hat$id,   levels = sort(as.character(unique(lambda_hat$id))),   labels = sort(as.character(unique(lambda_hat$id))))
gamma_hat$id  = factor(gamma_hat$id,    levels = sort(as.character(unique(gamma_hat$id))),    labels = sort(as.character(unique(gamma_hat$id))))
delta_hat$id  = factor(delta_hat$id,    levels = sort(as.character(unique(delta_hat$id))), labels = sort(as.character(unique(delta_hat$id))))

y_train_hat$id = factor(y_train_hat$id, levels = sort(as.character(unique(y_train_hat$id))), labels = sort(as.character(unique(y_train_hat$id))))
y_test_hat$id  = factor(y_test_hat$id,  levels = sort(as.character(unique(y_test_hat$id))),  labels = sort(as.character(unique(y_test_hat$id))))

# Box plots for the parameter estimates ----------
# alpha_hat = alpha_hat %>% group_by(id, Q, variable) %>% summarise(value=mean(value), true=NA)
plot_individual_boxplots <- function(object){

  db = object
  names(db) = c('Method', 'Q', 'Parameter', 'value', 'true')
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
                     labels = fixed_labels)
}
plot_individual_boxplots(alpha_hat)
plot_individual_boxplots(beta_hat)
# plot_individual_boxplots(lambda_hat)
# plot_individual_boxplots(gamma_hat)
# plot_individual_boxplots(delta_hat)

alpha_hat %>% group_by(id) %>% summarise(mean=mean(value)) # Bayes AMMI (no postproc) has mean > 0
beta_hat %>% group_by(id) %>% summarise(mean=mean(value)) # Bayes AMMI (no postproc) has mean > 0

## Density plots for the DIFFERENCE (true - predicted)

plot_individual_density <- function(object, data){

  if (deparse(substitute(object)) == 'blinear_hat') {aux_title = expression(Sigma[q]~lambda[q]~gamma[iq]~delta[jq]~-~Sigma[q]~hat(lambda)[q]~hat(gamma)[iq]~hat(delta)[jq])}
  if (deparse(substitute(object)) == 'y_train_hat') {aux_title = expression('Training data:'~y - hat(y))}
  if (deparse(substitute(object)) == 'y_test_hat') {aux_title = expression('Test data:'~y - hat(y))}

  db = object
  names(db) = c('Method', 'Q', 'Parameter', 'value', 'true')

  db %>%
    ggplot(aes(x=value, colour=Method)) +
    geom_density(alpha=0.4)+
    theme_bw() +
    labs(title = aux_title) +
    theme(axis.text.x = element_text(size=12),
          axis.title = element_blank(),
          plot.title = element_text(size = 18, hjust = 0.5),
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.position = "bottom")
}
# plot_individual_density(blinear_hat, data)
plot_individual_density(y_train_hat, data)
# plot_individual_density(y_test_hat, data)

# Calculate correlation ------
yhat_ammi_q1 = y_train_hat %>% filter(Q==1) %>% mutate(y_hat = value + true)
yhat_ammi_q2 = y_train_hat %>% filter(Q==2) %>% mutate(y_hat = value + true)
yhat_ammi_q3 = y_train_hat %>% filter(Q==3) %>% mutate(y_hat = value + true)
cor(data$y, yhat_ammi_q1$y_hat)
cor(data$y, yhat_ammi_q2$y_hat)
cor(data$y, yhat_ammi_q3$y_hat)
cor(data$y, apply(ambarti$y_hat,2,mean))
