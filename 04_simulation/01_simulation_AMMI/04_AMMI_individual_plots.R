library(devtools)
install_github("ebprado/AMBARTI/R package",
               ref = 'main',
               auth_token = '363d4ad84eaa25d3eb26752732cc208f7e698086')

library(AMBARTI)
save_file = "/Users/estevaoprado/Documents/GitHub/AMBARTI/04_simulation/01_simulation_AMMI/results/"

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

# Get parameter estimates from classical AMMI
res_AMMI = AMMI_help_plot(classical_AMMI, data)

g_hat       = res_AMMI$g_hat
beta_hat    = res_AMMI$beta_hat
lambda_hat  = res_AMMI$lambda_hat
gamma_hat   = res_AMMI$gamma_hat
delta_hat   = res_AMMI$delta_hat
blinear_hat = res_AMMI$blinear_hat
y_train_hat = res_AMMI$y_hat_train
y_test_hat  = res_AMMI$y_hat_test

# Get parameter estimates from Bayesian AMMI (WITHOUT postprocessing)
bAMMI_save_info = bAMMI_help_plot(bayesian_AMMI, data)

g_hat       = rbind(g_hat,    bAMMI_save_info$g_hat)
beta_hat    = rbind(beta_hat,     bAMMI_save_info$beta_hat)
lambda_hat  = rbind(lambda_hat,   bAMMI_save_info$lambda_hat)
gamma_hat   = rbind(gamma_hat,    bAMMI_save_info$gamma_hat)
delta_hat   = rbind(delta_hat,    bAMMI_save_info$delta_hat)
blinear_hat = rbind(blinear_hat,  bAMMI_save_info$blinear_hat)
y_train_hat = rbind(y_train_hat,  bAMMI_save_info$y_hat_train)
y_test_hat  = rbind(y_test_hat,   bAMMI_save_info$y_hat_test)

# Get parameter estimates from Bayesian AMMI (WITH postprocessing)
bAMMI_save_info_WITHPOS = bAMMI_help_plot_WITHPOS(bayesian_AMMI, data)

g_hat       = rbind(g_hat,   bAMMI_save_info_WITHPOS$g_hat)
beta_hat    = rbind(beta_hat,    bAMMI_save_info_WITHPOS$beta_hat)
lambda_hat  = rbind(lambda_hat,  bAMMI_save_info_WITHPOS$lambda_hat)
gamma_hat   = rbind(gamma_hat,   bAMMI_save_info_WITHPOS$gamma_hat)
delta_hat   = rbind(delta_hat,   bAMMI_save_info_WITHPOS$delta_hat)
blinear_hat = rbind(blinear_hat, bAMMI_save_info_WITHPOS$blinear_hat)
y_train_hat = rbind(y_train_hat, bAMMI_save_info_WITHPOS$y_hat_train)
y_test_hat  = rbind(y_test_hat,  bAMMI_save_info_WITHPOS$y_hat_test)

# Get parameter estimates from AMBARTI
AMBARTI_save_info = AMBARTI_help_plot(ambarti, data)

g_hat       = rbind(g_hat, AMBARTI_save_info$g_hat)
beta_hat    = rbind(beta_hat, AMBARTI_save_info$beta_hat)
blinear_hat = rbind(blinear_hat, AMBARTI_save_info$blinear_hat)
y_train_hat = rbind(y_train_hat, AMBARTI_save_info$y_hat_train)
y_test_hat  = rbind(y_test_hat,  AMBARTI_save_info$y_hat_test)

# Some postprocessing

g_hat$id       = factor(g_hat$id,   levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI (PP)', 'B-AMMI (No PP)', 'AMBARTI'))
beta_hat$id    = factor(beta_hat$id,    levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI (PP)', 'B-AMMI (No PP)', 'AMBARTI'))
blinear_hat$id = factor(blinear_hat$id, levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI (PP)', 'B-AMMI (No PP)', 'AMBARTI'))
y_train_hat$id = factor(y_train_hat$id, levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI (PP)', 'B-AMMI (No PP)', 'AMBARTI'))
y_test_hat$id  = factor(y_test_hat$id,  levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI (PP)', 'B-AMMI (No PP)', 'AMBARTI'))


g_hat$id       = factor(g_hat$id,   levels = sort(as.character(unique(g_hat$id))),   labels = sort(as.character(unique(g_hat$id))))
beta_hat$id    = factor(beta_hat$id,    levels = sort(as.character(unique(beta_hat$id))),    labels = sort(as.character(unique(beta_hat$id))))
blinear_hat$id = factor(blinear_hat$id, levels = sort(as.character(unique(blinear_hat$id))), labels = sort(as.character(unique(blinear_hat$id))))
y_train_hat$id = factor(y_train_hat$id, levels = sort(as.character(unique(y_train_hat$id))), labels = sort(as.character(unique(y_train_hat$id))))
y_test_hat$id  = factor(y_test_hat$id,  levels = sort(as.character(unique(y_test_hat$id))),  labels = sort(as.character(unique(y_test_hat$id))))

# Box plots for the parameter estimates ----------

plot_individual_boxplots <- function(object, data){
  
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
                     labels = fixed_labels) # +
  # scale_color_discrete(labels=unique(db$Method))
}
plot_individual_boxplots(g_hat, data)
plot_individual_boxplots(beta_hat, data)
plot_individual_boxplots(lambda_hat, data)
plot_individual_boxplots(gamma_hat, data)
plot_individual_boxplots(delta_hat, data)

g_hat %>% group_by(id) %>% summarise(mean=mean(value)) # Bayes AMMI (no postproc) has mean > 0
beta_hat %>% group_by(id) %>% summarise(mean=mean(value)) # Bayes AMMI (no postproc) has mean > 0

## Density plots for the DIFFERENCE (true - predicted)

plot_individual_density <- function(object){
  
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
plot_individual_density(blinear_hat)
plot_individual_density(y_train_hat)
plot_individual_density(y_test_hat)
