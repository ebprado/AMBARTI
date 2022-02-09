library(devtools)
library(ggplot2)
install_github("ebprado/AMBARTI/R package",
               ref = 'main')

library(AMBARTI)
save_file = "~/R/AMBARTI/04_real_data_sets/results/"

#-------------------------------------------------
# Individual plots
#-------------------------------------------------
library(reshape2)
setwd(save_file)
year = 2015
filename = paste('Ireland_VCU_',
                 year,
                 sep='')

data_filename    = paste(filename, '_data.RData',              sep='')
ammi_filenameQ1  = paste(filename, '_classical_AMMI_Q1.RData', sep='')
ammi_filenameQ2  = paste(filename, '_classical_AMMI_Q2.RData', sep='')
ammi_filenameQ3  = paste(filename, '_classical_AMMI_Q3.RData', sep='')
bammi_filenameQ1  = paste(filename, '_bayesian_AMMI_Q1.RData', sep='')
bammi_filenameQ2  = paste(filename, '_bayesian_AMMI_Q2.RData', sep='')
bammi_filenameQ3  = paste(filename, '_bayesian_AMMI_Q3.RData', sep='')
ambarti_filename = paste(filename, '_AMBARTI.RData',           sep='')

load(paste(save_file, data_filename,    sep=''))
load(paste(save_file, ammi_filenameQ1,  sep=''))
load(paste(save_file, ammi_filenameQ2,  sep=''))
load(paste(save_file, ammi_filenameQ3,  sep=''))
load(paste(save_file, bammi_filenameQ1,  sep=''))
load(paste(save_file, bammi_filenameQ2,  sep=''))
load(paste(save_file, bammi_filenameQ3,  sep=''))
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
  
  # x_train <- data$x
  # y_train <- data$y
  # 
  # g = as.factor(x_train[,"g"])
  # e = as.factor(x_train[,"e"])
  # 
  # # Fit the linear model
  # linear_mod = aov(y_train ~ g + e + g:e)
  # linear_mo2 = lm(y_train ~ g + e, contrasts = list(
  #   g = contr.treatment(n = 18, contrasts = FALSE),
  #   e = contr.treatment(n = 9, contrasts = FALSE)))
  # summary(linear_mo2)

  alpha_hat   = rbind(alpha_hat,    res_AMMI$g_hat)
  beta_hat    = rbind(beta_hat,     res_AMMI$e_hat)
  lambda_hat  = rbind(lambda_hat,   res_AMMI$lambda_hat)
  gamma_hat   = rbind(gamma_hat,    res_AMMI$gamma_hat)
  delta_hat   = rbind(delta_hat,    res_AMMI$delta_hat)
  blinear_hat = rbind(blinear_hat,  res_AMMI$blinear_hat)
  y_train_hat = rbind(y_train_hat,  res_AMMI$y_hat_train)
  y_test_hat  = rbind(y_test_hat,   res_AMMI$y_hat_test)
  
}

# Helps create the boxplots for AMMI based on the assumption that the coefficients in a linear model are normally distributed

for(qq in 1:3){
  for (k in 1:length(unique(res_AMMI$g_hat[,3]))){
    alpha_hat = rbind(alpha_hat, data.frame(id = 'AMMI', Q=qq, variable=paste('g[', k, ']', sep=''), value=rnorm(1000, alpha_hat[k,4], sd=0.18373), true = NA))
  }
  for (k in 1:length(unique(res_AMMI$e_hat[,3]))){
    beta_hat = rbind(beta_hat, data.frame(id = 'AMMI', Q=qq, variable=paste('e[', k, ']', sep=''), value=rnorm(1000, beta_hat[k,4], sd=0.12993), true = NA))
  }
}

for (Q in 1:3){
  aux_bAMMI    = get(paste('bayesian_AMMI_Q',Q, sep=''))
  res_bAMMI    = bAMMI_help_plot_WITHPOS(aux_bAMMI, data, Q = Q)
  # res_bAMMI    = bAMMI_help_plot(aux_bAMMI, data, Q = Q)
  
  alpha_hat   = rbind(alpha_hat,    res_bAMMI$g_hat)
  beta_hat    = rbind(beta_hat,     res_bAMMI$e_hat)
  lambda_hat  = rbind(lambda_hat,   res_bAMMI$lambda_hat)
  gamma_hat   = rbind(gamma_hat,    res_bAMMI$gamma_hat)
  delta_hat   = rbind(delta_hat,    res_bAMMI$delta_hat)
  blinear_hat = rbind(blinear_hat,  res_bAMMI$blinear_hat)
  y_train_hat = rbind(y_train_hat,  res_bAMMI$y_hat_train)
  y_test_hat  = rbind(y_test_hat,   res_bAMMI$y_hat_test)
  
}

for (Q in 1:3){
  aux_bAMMI    = get(paste('bayesian_AMMI_Q',Q, sep=''))
  res_bAMMI    = bAMMI_help_plot(aux_bAMMI, data, Q = Q) # NO POST-PROCESSING
  
  alpha_hat   = rbind(alpha_hat,    res_bAMMI$g_hat)
  beta_hat    = rbind(beta_hat,     res_bAMMI$e_hat)
  lambda_hat  = rbind(lambda_hat,   res_bAMMI$lambda_hat)
  gamma_hat   = rbind(gamma_hat,    res_bAMMI$gamma_hat)
  delta_hat   = rbind(delta_hat,    res_bAMMI$delta_hat)
  blinear_hat = rbind(blinear_hat,  res_bAMMI$blinear_hat)
  y_train_hat = rbind(y_train_hat,  res_bAMMI$y_hat_train)
  y_test_hat  = rbind(y_test_hat,   res_bAMMI$y_hat_test)
  
}

# Get parameter estimates from AMBARTI
AMBARTI_save_info = AMBARTI_help_plot(ambarti, data)

alpha_hat   = rbind(alpha_hat,   AMBARTI_save_info$g_hat)
beta_hat    = rbind(beta_hat,    AMBARTI_save_info$e_hat)
blinear_hat = rbind(blinear_hat, AMBARTI_save_info$blinear_hat)
y_train_hat = rbind(y_train_hat, AMBARTI_save_info$y_hat_train)
y_test_hat  = rbind(y_test_hat,  AMBARTI_save_info$y_hat_test)

# Some postprocessing

alpha_hat$id   = factor(alpha_hat$id,   levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI', 'B-AMMI (No PP)', 'AMBARTI'))
beta_hat$id    = factor(beta_hat$id,    levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI', 'B-AMMI (No PP)', 'AMBARTI'))
lambda_hat$id  = factor(lambda_hat$id,  levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI', 'B-AMMI (No PP)', 'AMBARTI'))
gamma_hat$id   = factor(gamma_hat$id,   levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI', 'B-AMMI (No PP)', 'AMBARTI'))
delta_hat$id   = factor(delta_hat$id,   levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI', 'B-AMMI (No PP)', 'AMBARTI'))
blinear_hat$id = factor(blinear_hat$id, levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI', 'B-AMMI (No PP)', 'AMBARTI'))
y_train_hat$id = factor(y_train_hat$id, levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI', 'B-AMMI (No PP)', 'AMBARTI'))
y_test_hat$id  = factor(y_test_hat$id,  levels = c('AMMI', 'Bayesian AMMI (posproc)', 'Bayesian AMMI (no postproc)', 'AMBARTI'), labels = c('AMMI', 'B-AMMI', 'B-AMMI (No PP)', 'AMBARTI'))

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
    theme_bw(base_size=18) +
    # labs(title = bquote(Comparison~of~.(sym(aux_name)))) +
    theme(axis.text.x = element_text(size=15),
          axis.title = element_blank(),
          plot.title = element_text(size = 20, hjust = 0.5),
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_x_discrete(limit = orig_labels,
                     labels = fixed_labels) + 
    guides(color = guide_legend(nrow = 1))
}

alpha_hat$id = gsub('[(Q = 1)]', '', alpha_hat$id) # only for the paper
beta_hat$id = gsub('[(Q = 1)]', '', beta_hat$id) # only for the paper

plot_individual_boxplots(alpha_hat %>% filter(id %in% c('AMBARTI','AMMI')))
plot_individual_boxplots(beta_hat %>% filter(id %in% c('AMBARTI','AMMI')))
# plot_individual_boxplots(lambda_hat)
# plot_individual_boxplots(gamma_hat)
# plot_individual_boxplots(delta_hat)