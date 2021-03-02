library(AMBARTI)
library(dplyr)
library(ggplot2)
save_file = "~/Documents/GitHub/AMBARTI/03_real_datasets/results/"

# Relative Root Mean Squared Error (RRMSE)
RRMSE <- function(true, predicted){
  sqrt(mean(((true - predicted)/true)^2))
}

RMSE <- function(true, predicted){
  sqrt(mean(((true - predicted))^2))
}

library(reshape2)
setwd(save_file)
years = 2010:2019
save=NULL

set.seed(001)

for (year in years){
  filename = paste('Ireland_VCU_',
                   year,
                   sep='')

  data_filename    = paste(filename, '_data.RData',              sep='')
  ammi_filenameQ3  = paste(filename, '_classical_AMMI_Q3.RData', sep='')
  ambarti_filename = paste(filename, '_AMBARTI.RData',           sep='')

  load(paste(save_file, data_filename,    sep=''))
  load(paste(save_file, ammi_filenameQ3,  sep=''))
  load(paste(save_file, ambarti_filename, sep=''))

  file_path = "~/Documents/GitHub/AMBARTI/03_real_datasets/01_data_sets/"
  file_name = 'Historical_data_Ireland_VCU.csv'

  ireland = read.csv(paste(file_path, file_name, sep=''), sep = ';', dec = ',')
  ireland$Year = as.factor(ireland$Year)

  ireland_year = ireland %>%
    filter(Year==year) %>%
    mutate(g = factor(Genotype, levels = unique(Genotype), labels=paste('g', 1:length(unique(Genotype)), sep='')),
           e = factor(Location, levels = unique(Location), labels=paste('e', 1:length(unique(Location)), sep='')),
           y = yld_ton_ha) %>%
    select(g, Genotype, e, Location, y)

  ireland_year = ireland_year %>%
    select(g, e, y) %>%
    group_by(g, e) %>%
    sample_n(2, ) %>%
    data.frame()

  y_hat_ambarti = predict_ambarti(ambarti, ireland_year, 'mean')

  x_test = ireland_year

  predict_ammi = function(object, data) {
    x_test$g = as.numeric(gsub('g', '', x_test$g))
    x_test$e = as.numeric(gsub('e', '', x_test$e))
    object = classical_AMMI_Q3
    mu_hat = object$mu_hat
    # g_hat      = object$g_hat
    # e_hat      = object$e_hat
    g_hat      = object$alpha_hat
    e_hat      = object$beta_hat
    lambda_hat = object$lambda_hat
    gamma_hat  = object$gamma_hat
    delta_hat  = object$delta_hat
    Q = ncol(object$gamma_hat)
    blin_test = rep(0, nrow(x_test))
    for (k in 1:Q) {
      blin_test = blin_test + lambda_hat[k]*gamma_hat[x_test[,'g'],k]*delta_hat[x_test[,'e'],k]
    }
    y_hat_test = mu_hat + g_hat[x_test[,'g']] + e_hat[x_test[,'e']] + blin_test
    return(y_hat_test)
  }

  y_hat_test = predict_ammi(classical_AMMI_Q3, x_test)
  save = rbind(save,
               cbind(year = year,
                     RMSE_AMMI = RMSE(ireland_year$y, y_hat_test),
                     RMSE_ambarti  = RMSE(ireland_year$y, y_hat_ambarti)))
}
