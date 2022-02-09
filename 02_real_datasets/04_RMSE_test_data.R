library(AMBARTI)
library(dplyr)
library(ggplot2)
library(BASS)
library(gss)
library(xtable)
library(earth)
library(xtable)
save_file = "/users/research/eprado/R/AMBARTI/04_real_data_sets/results/"

# Relative Root Mean Squared Error (RRMSE)
RRMSE <- function(true, predicted){
  sqrt(mean(((true - predicted)/true)^2))
}

RMSE <- function(true, predicted){
  sqrt(mean(((true - predicted))^2))
}

predict_ammi = function(object, data) {
  data$g = as.numeric(gsub('g', '', data$g))
  data$e = as.numeric(gsub('e', '', data$e))
  mu_hat = object$mu_hat
  g_hat      = object$g_hat
  e_hat      = object$e_hat
  # g_hat      = object$alpha_hat
  # e_hat      = object$beta_hat
  lambda_hat = object$lambda_hat
  gamma_hat  = object$gamma_hat
  delta_hat  = object$delta_hat
  Q = ncol(object$gamma_hat)
  blin_test = rep(0, nrow(data))
  for (k in 1:Q) {
    blin_test = blin_test + lambda_hat[k]*gamma_hat[data[,'g'],k]*delta_hat[data[,'e'],k]
  }
  y_hat_test = mu_hat + g_hat[data[,'g']] + e_hat[data[,'e']] + blin_test
  return(y_hat_test)
}
predict_bammi <- function(object, traindata, newdata, Q = NULL){
  # order the test set
  # newdata = newdata[order(newdata$g,newdata$e),]
  
  # Get some MCMC info
  nburn      = object$BUGSoutput$n.burnin
  niter      = object$BUGSoutput$n.iter
  npost      = niter - nburn
  seq_burn   = seq(1, nburn, by=1)
  
  # Get estimates info
  estimate   = object$BUGSoutput$sims.matrix
  mu_hat     = estimate[-seq_burn,colnames(estimate)[grepl('mu_all', colnames(estimate))]]
  g_hat      = estimate[-seq_burn,colnames(estimate)[grepl('^g\\[',  colnames(estimate))]]
  e_hat      = estimate[-seq_burn,colnames(estimate)[grepl('^e\\[',   colnames(estimate))]]
  delta_hat  = estimate[-seq_burn,colnames(estimate)[grepl('delta',  colnames(estimate))]]
  gamma_hat  = estimate[-seq_burn,colnames(estimate)[grepl('gamma',  colnames(estimate))]]
  lambda_hat = estimate[-seq_burn,colnames(estimate)[grepl('lambda', colnames(estimate))]]
  lambda_hat = as.matrix(lambda_hat)
  
  # Compute the bilinear term
  # blin_train = matrix(0, nrow = npost, ncol = nrow(x_train))
  blin_train  = matrix(0, nrow = npost, ncol = nrow(traindata))
  
  for (k in 1:Q) {
    # blin_train = blin_train + lambda_hat[,k] * gamma_hat[,x_train[,'g']]* delta_hat[,x_train[,'e']]
    blin_train  = blin_train +  lambda_hat[,k] * gamma_hat[,traindata[,'g']] * delta_hat[, traindata[,'e']]
  }
  
  # Compute the predicted response for the TRAINING traindata
  mu_ij = mu_hat + g_hat[,traindata[,'g']] + e_hat[,traindata[,'e']] + blin_train
  colnames(mu_ij) = NULL
  
  # -------------------------------------------------------------
  # Start postprocessing
  # -------------------------------------------------------------
  
  n_gen = length(unique(traindata[,'g']))
  n_env = length(unique(traindata[,'e']))
  
  # Create matrices/lists to store the postprocessing results
  snew_mu_hat     = matrix(NA, nrow=npost, ncol=1)
  snew_g_hat      = matrix(NA, nrow=npost, ncol=n_gen)
  snew_e_hat      = matrix(NA, nrow=npost, ncol=n_env)
  snew_lambda_hat = matrix(NA, nrow=npost, ncol=Q)
  snew_gamma_hat  = list()
  snew_delta_hat  = list()
  
  for (i in 1:nrow(mu_ij)){
    
    # Get the matrix with mu_ij for each genotype and environment
    matrix_mu_ij  = matrix(mu_ij[i,], nrow=n_gen, ncol=n_env)
    new_mu_hat    = mean(matrix_mu_ij)
    new_g_hat     = rowMeans(matrix_mu_ij) - new_mu_hat # their sum is 0
    new_e_hat     = colMeans(matrix_mu_ij) - new_mu_hat # their sum is 0
    
    # Center matrix by row and column
    resA = sweep(matrix_mu_ij*0, 2, -colMeans(matrix_mu_ij))
    resB = matrix_mu_ij*0 + rowMeans(matrix_mu_ij)
    res_double_centered = matrix_mu_ij - resA - resB + new_mu_hat # sum zero by row and column
    
    # Run the Singular Value Decomposition (SVD) to compute lambda, gamma, and delta
    sv_dec <- svd(res_double_centered, nu = Q, nv = Q)
    
    # Get new parameter estimates
    new_lambda_hat = sv_dec$d[1:Q]
    new_gamma_hat  = -1*sv_dec$u
    new_delta_hat  = -1*sv_dec$v
    
    # Store the new results
    snew_mu_hat[i,]       = new_mu_hat
    snew_g_hat[i,]        = new_g_hat
    snew_e_hat[i,]        = new_e_hat
    snew_lambda_hat[i,]   = new_lambda_hat
    snew_gamma_hat[[i]]   = new_gamma_hat
    snew_delta_hat[[i]]   = new_delta_hat
  }
  # Clean column names
  colnames(snew_mu_hat) = NULL
  colnames(snew_g_hat) = NULL
  colnames(snew_e_hat) = NULL
  colnames(snew_lambda_hat) = NULL
  
  # Summarise the new posterior results
  new_mu_hat     = apply(snew_mu_hat,2,mean)
  new_g_hat      = apply(snew_g_hat,2,mean)
  new_e_hat      = apply(snew_e_hat,2,mean)
  new_lambda_hat = apply(snew_lambda_hat ,2,mean)
  new_gamma_hat  = apply(simplify2array(snew_gamma_hat), 1:2, mean)
  new_delta_hat  = apply(simplify2array(snew_delta_hat), 1:2, mean)
  
  # Compute the bilinear term
  new_blin_train = rep(0, nrow(traindata))
  new_blin_test  = rep(0, nrow(newdata))

  for (k in 1:Q) {
    new_blin_train = new_blin_train + new_lambda_hat[k]*new_gamma_hat[traindata[,'g'],k]*new_delta_hat[traindata[,'e'],k]
    new_blin_test  = new_blin_test + new_lambda_hat[k]*new_gamma_hat[newdata[,'g'],k]*new_delta_hat[newdata[,'e'],k]
  }
  
  # Compute the predicted values for the TRAINING and TEST traindata sets
  new_mu_ij_train = new_mu_hat + new_g_hat[traindata[,'g']] + new_e_hat[traindata[,'e']] + new_blin_train
  new_mu_ij_test = new_mu_hat + new_g_hat[newdata[,'g']] + new_e_hat[newdata[,'e']]  + new_blin_test
  
  return(new_mu_ij_test)
}
predict_gss = function(object, newdata){
  return(predict(object, newdata))
}
predict_bmars = function(object, newdata){
  return(predict(object, newdata))
}
predict_mars = function(object, newdata){
  return(predict(object, newdata))
}
file_path = "~/R/AMBARTI/04_real_data_sets/01_data_sets/"
file_name = 'Historical_data_Ireland_VCU.csv'
ireland = read.csv(paste(file_path, file_name, sep=''), sep = ';', dec = ',')

# Set up where to save the results ---- 
save_path = "~/R/AMBARTI/04_real_data_sets/results/"

# Get real data sets ------ 

library(reshape2)
setwd(save_file)
save = NULL

validation <- function(year){
  
  set.seed(001)
  data_filename    = paste('Ireland_VCU_',  year,  '_data.RData',             sep='') # train data (only for BAMMI)
  ammi_filenameQ3  = paste('Ireland_VCU_',  year, '_classical_AMMI_Q3.RData', sep='')
  bammi_filenameQ3 = paste('Ireland_VCU_',  year, '_bayesian_AMMI_Q3.RData',  sep='')
  gss_filename     = paste('Ireland_VCU_',  year, '_gss.RData',               sep='')
  bmars_filename   = paste('Ireland_VCU_',  year, '_bmars.RData',             sep='')
  mars_filename    = paste('Ireland_VCU_',  year, '_mars.RData',              sep='')
  ambarti_filename = paste('Ireland_VCU_',  year, '_AMBARTI.RData',           sep='')
  
  load(paste(save_file, data_filename,    sep=''))
  load(paste(save_file, ammi_filenameQ3,  sep=''))
  load(paste(save_file, bammi_filenameQ3, sep=''))
  load(paste(save_file, gss_filename,     sep=''))
  load(paste(save_file, bmars_filename,   sep=''))
  load(paste(save_file, mars_filename,    sep=''))
  load(paste(save_file, ambarti_filename, sep=''))
  
  ireland_year = ireland %>% 
    filter(Year==year) %>% 
    mutate(g = factor(Genotype, levels = unique(Genotype), labels=paste('g', 1:length(unique(Genotype)), sep='')),
           e = factor(Location, levels = unique(Location), labels=paste('e', 1:length(unique(Location)), sep='')),
           y = yld_ton_ha) %>% 
    dplyr::select(g, e, y) %>% 
    data.frame()
  
  ireland_year = ireland_year %>% 
    group_by(g,e) %>% 
    mutate(id = 1:n())
    
  ireland_year = ireland_year %>% 
    filter(id %in% c(2,4)) %>% 
    dplyr::select(g,e,y) %>% 
    arrange(g,e) %>% 
    data.frame()
  
  x_test = ireland_year[,c('g', 'e')]
  
  y_hat_ambarti = predict_ambarti(ambarti, x_test, 'mean')
  
  y_hat_test_ammi = predict_ammi(classical_AMMI_Q3, x_test)
  y_hat_test_bammi = predict_bammi(bayesian_AMMI_Q3, data$x, x_test, Q=3)
  y_hat_test_gss = predict_gss(gss, x_test)
  y_hat_test_bmars = predict_bmars(bmars, x_test)
  y_hat_test_mars = predict_mars(mars, x_test)
  
  save = rbind(save,
               cbind(dataset = year,
                     RMSE_bAMMI = RMSE(ireland_year$y, y_hat_test_bammi),
                     RMSE_ambarti  = RMSE(ireland_year$y, y_hat_ambarti$npostproc),
                     RMSE_gss = RMSE(ireland_year$y, y_hat_test_gss),
                     RMSE_bmars = RMSE(ireland_year$y, apply(y_hat_test_bmars,2,mean)),
                     RMSE_mars = RMSE(ireland_year$y, y_hat_test_mars),
                     RMSE_AMMI = RMSE(ireland_year$y, y_hat_test_ammi)
               )
  )
  save(save, file = paste(save_file, 'Prediction_test_data_VCU_ireland_' , year, '.RData',sep=''))
  return(save)
}

validation(2010)
validation(2011)
validation(2012)
validation(2013)
validation(2014)
validation(2015)
validation(2016)
validation(2017)
validation(2018)
validation(2019)

#### Table 1 in the manuscript ------ 
results = NULL
for (year in 2010:2019){
  load(paste(save_file, 'Prediction_test_data_VCU_ireland_', year, '.RData',  sep=''))
  results = rbind(results, save)
}
head(results)
print(xtable(results[,c(1,3,7,2,4,5)]))
