
library(devtools)
install_github("ebprado/AMBARTI/R package",
               ref = 'main',
               auth_token = '363d4ad84eaa25d3eb26752732cc208f7e698086')

library(AMBARTI)
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

# Get parameter estimates from classical AMMI
res_AMMI = AMMI_help_plot(classical_AMMI, data)

alpha_hat  = res_AMMI$alpha_hat
beta_hat   = res_AMMI$beta_hat
lambda_hat = res_AMMI$lambda_hat
gamma_hat  = res_AMMI$gamma_hat
delta_hat  = res_AMMI$delta_hat
blinear_hat = res_AMMI$blinear_hat
y_train_hat = res_AMMI$y_hat_train
y_test_hat = res_AMMI$y_hat_test

# Get parameter estimates from Bayesian AMMI (WITHOUT postprocessing)
bAMMI_save_info = bAMMI_help_plot(bayesian_AMMI, data)

alpha_hat   = rbind(alpha_hat,    bAMMI_save_info$alpha_hat)
beta_hat    = rbind(beta_hat,     bAMMI_save_info$beta_hat)
lambda_hat  = rbind(lambda_hat,   bAMMI_save_info$lambda_hat)
gamma_hat   = rbind(gamma_hat,    bAMMI_save_info$gamma_hat)
delta_hat   = rbind(delta_hat,    bAMMI_save_info$delta_hat)
blinear_hat = rbind(blinear_hat,  bAMMI_save_info$blinear_hat)
y_train_hat = rbind(y_train_hat,  bAMMI_save_info$y_hat_train)
y_test_hat  = rbind(y_test_hat,   bAMMI_save_info$y_test_hat)

# Get parameter estimates from Bayesian AMMI (WITH postprocessing)
bAMMI_save_info_WITHPOS = bAMMI_help_plot_WITHPOS(bayesian_AMMI, data)

alpha_hat   = rbind(alpha_hat,   bAMMI_save_info_WITHPOS$alpha_hat)
beta_hat    = rbind(beta_hat,    bAMMI_save_info_WITHPOS$beta_hat)
lambda_hat  = rbind(lambda_hat,  bAMMI_save_info_WITHPOS$lambda_hat)
gamma_hat   = rbind(gamma_hat,   bAMMI_save_info_WITHPOS$gamma_hat)
delta_hat   = rbind(delta_hat,   bAMMI_save_info_WITHPOS$delta_hat)
blinear_hat = rbind(blinear_hat, bAMMI_save_info_WITHPOS$blinear_hat)
y_train_hat = rbind(y_train_hat, bAMMI_save_info_WITHPOS$y_hat_train)
y_test_hat  = rbind(y_test_hat,  bAMMI_save_info_WITHPOS$y_test_hat)

# Get parameter estimates from AMBARTI
AMBARTI_save_info = AMBARTI_help_plot(ambarti, data)

alpha_hat   = rbind(alpha_hat, AMBARTI_save_info$alpha_hat)
beta_hat    = rbind(beta_hat, AMBARTI_save_info$beta_hat)
blinear_hat = rbind(blinear_hat, AMBARTI_save_info$blinear_hat)
y_train_hat = rbind(y_train_hat, AMBARTI_save_info$y_hat_train)
y_test_hat  = rbind(y_test_hat,  AMBARTI_save_info$y_test_hat)

# Plot ----------

plot_individual_boxplots <- function(object, data){

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
plot_individual_boxplots(alpha_hat, data)
plot_individual_boxplots(beta_hat, data)
plot_individual_boxplots(lambda_hat, data)
plot_individual_boxplots(gamma_hat, data)
plot_individual_boxplots(delta_hat, data)

## ------------

plot_individual_density <- function(object, data){

  db = object
  names(db) = c('Method', 'Parameter', 'value', 'true')

  db %>%
    ggplot(aes(x=value, colour=Method)) +
    geom_density(alpha=0.4)+
    theme_bw() +
    labs(title = 'bla bla') +
    theme(axis.text.x = element_text(size=12),
          axis.title = element_blank(),
          plot.title = element_text(size = 18, hjust = 0.5),
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.position = "bottom")  +
    # scale_x_discrete(limit = orig_labels,
                     # labels = fixed_labels) +
    scale_color_discrete(labels=unique(db$Method))
}
plot_individual_density(blinear_hat, data)




