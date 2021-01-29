# Relative Root Mean Squared Error (RRMSE)
RRMSE <- function(true, predicted){
  sqrt(mean(((true - predicted)/true)^2))
}

RMSE <- function(true, predicted){
  sqrt(mean(((true - predicted))^2))
}

#' @export
#'
get_metrics = function(object, data, rep){

  id          = object$id
  y_train     = data$y
  blinear     = data$blinear
  y_test      = data$y_test
  y_hat_train = object$y_hat_train
  y_hat_test  = object$y_hat_test
  blinear_hat = object$blinear_hat

  rmse_y_train = RMSE(y_train, y_hat_train)
  rmse_y_test  = RMSE(y_test, y_hat_test)
  rmse_blinear = RMSE(blinear, blinear_hat)

  alpha  = data$alpha
  beta   = data$beta
  lambda = data$lambda
  gamma  = data$gamma
  delta  = data$delta

  alpha_hat  = object$alpha_hat
  beta_hat   = object$beta_hat

  rrmse_alpha = RRMSE(alpha, alpha_hat)
  rrmse_beta  = RRMSE(beta, beta_hat)

  if (id != 'AMBARTI') {

    lambda_hat = object$lambda_hat
    gamma_hat  = object$gamma_hat
    delta_hat  = object$delta_hat

    rrmse_lambda = RRMSE(lambda, lambda_hat)
    rrmse_gamma  = RRMSE(gamma, gamma_hat)
    rrmse_delta  = RRMSE(delta, delta_hat)

  } else{

    rrmse_lambda = NA
    rrmse_gamma = NA
    rrmse_delta = NA

  }

  I  = data$I
  J  = data$J
  sa = data$s_alpha
  sb = data$s_beta
  sy = data$s_y
  lambda = gsub(', ', ' ', toString(data$lambda))

aux = data.frame(
                 id           = id,
                 rep          = rep,
                 I            = I,
                 J            = J,
                 sa           = sa,
                 sb           = sb,
                 sy           = sy,
                 lambda       = lambda,
                 y_train_rmse = rmse_y_train,
                 y_test_rmse  = rmse_y_test,
                 rmse_blinear = rmse_blinear,
                 rrmse_alpha  = rrmse_alpha,
                 rrmse_beta   = rrmse_beta,
                 lambda_rrmse = rrmse_lambda,
                 gamma_rrmse  = rrmse_gamma,
                 delta_rrmse  = rrmse_delta
                 )

return(aux)

}

#' @export
#' @importFrom stats 'aov' 'model.tables' 'lm'
organise_classical_AMMI <- function(object, data){

  # Get training info
  x_train = data$x
  y_train = data$y

  # Get test info
  x_test = data$x
  y_test = data$y_test

  mu_hat     = object$mu_hat
  alpha_hat  = object$alpha_hat
  beta_hat   = object$beta_hat
  lambda_hat = object$lambda_hat
  gamma_hat  = object$gamma_hat
  delta_hat  = object$delta_hat
  if (is.null(data$Q) == FALSE) {Q = data$Q}

  # Set up the bilinear term
  blin_train = rep(0, length(y_train))
  blin_test  = rep(0, length(y_test))

  for (k in 1:Q) {
    blin_train = blin_train + lambda_hat[k]*gamma_hat[x_train[,'g'],k]*delta_hat[x_train[,'e'],k]
    blin_test  = blin_test + lambda_hat[k]*gamma_hat[x_test[,'g'],k]*delta_hat[x_test[,'e'],k]
  }

  # Compute the predicted response for the TRAINING data
  y_hat_train = mu_hat + alpha_hat[x_train[,'g']] + beta_hat[x_train[,'e']] + blin_train

  # Compute the predicted response for the TEST data
  y_hat_test = mu_hat + alpha_hat[x_test[,'g']] + beta_hat[x_test[,'e']] + blin_test

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

#' @export
#'
organise_bayesian_AMMI_WITH_postprocessing <- function(object, data){

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

  # Compute the bilinear term
  blin_train = matrix(0, nrow = npost, ncol = nrow(x_train))
  blin_test  = matrix(0, nrow = npost, ncol = nrow(x_test))

  for (k in 1:PC) {
    blin_train = blin_train + lambda_hat[,k] * gamma_hat[,x_train[,'g']]* delta_hat[,x_train[,'e']]
    blin_test  = blin_test +  lambda_hat[,k] * gamma_hat[,x_test[,'g']] * delta_hat[, x_test[,'e']]
  }

  # Compute the predicted response for the TRAINING data
  mu_ij = mu_hat + alpha_hat[,x_train[,'g']] + beta_hat[,x_train[,'e']] + blin_train
  colnames(mu_ij) = NULL

  # -------------------------------------------------------------
  # Start postprocessing
  # -------------------------------------------------------------

  n_gen = length(unique(x_train[,'g']))
  n_env = length(unique(x_train[,'e']))

  # Create matrices/lists to store the postprocessing results
  snew_mu_hat     = matrix(NA, nrow=npost, ncol=1)
  snew_alpha_hat  = matrix(NA, nrow=npost, ncol=n_gen)
  snew_beta_hat   = matrix(NA, nrow=npost, ncol=n_env)
  snew_lambda_hat = matrix(NA, nrow=npost, ncol=PC)
  snew_gamma_hat  = list()
  snew_delta_hat  = list()

  for (i in 1:nrow(mu_ij)){

    # Get the matrix with mu_ij for each genotype and environment
    matrix_mu_ij  = matrix(mu_ij[i,], nrow=n_gen, ncol=n_env)
    new_mu_hat    = mean(matrix_mu_ij)
    new_alpha_hat = new_mu_hat - rowMeans(matrix_mu_ij) # their sum is 0
    new_beta_hat  = new_mu_hat - colMeans(matrix_mu_ij) # their sum is 0

    # Center matrix by row and column
    # Thank https://stackoverflow.com/questions/43639063/double-centering-in-r
    resA = t(matrix_mu_ij*0 + colMeans(matrix_mu_ij))
    resB = matrix_mu_ij*0 + rowMeans(matrix_mu_ij)
    res_double_centered = matrix_mu_ij - resA - resB + new_mu_hat # sum zero by row and column

    # Run the Singular Value Decomposition (SVD) to compute lambda, gamma, and delta
    sv_dec <- svd(res_double_centered, nu = PC, nv = PC)

    # Get new parameter estimates
    new_lambda_hat = sv_dec$d[1:PC]
    new_gamma_hat  = -1*sv_dec$u
    new_delta_hat  = -1*sv_dec$v

    # Store the new results
    snew_mu_hat[i,]       = new_mu_hat
    snew_alpha_hat[i,]    = new_alpha_hat
    snew_beta_hat[i,]     = new_beta_hat
    snew_lambda_hat[i,]   = new_lambda_hat
    snew_gamma_hat[[i]]   = new_gamma_hat
    snew_delta_hat[[i]]   = new_delta_hat
  }
  # Clean column names
  colnames(snew_mu_hat) = NULL
  colnames(snew_alpha_hat) = NULL
  colnames(snew_beta_hat) = NULL
  colnames(snew_lambda_hat) = NULL

  # Summarise the new posterior results
  new_mu_hat     = apply(snew_mu_hat     ,2,mean)
  new_alpha_hat  = apply(snew_alpha_hat  ,2,mean)
  new_beta_hat   = apply(snew_beta_hat   ,2,mean)
  new_lambda_hat = apply(snew_lambda_hat ,2,mean)
  new_gamma_hat  = apply(simplify2array(snew_gamma_hat), 1:2, mean)
  new_delta_hat  = apply(simplify2array(snew_delta_hat), 1:2, mean)

  # Compute the bilinear term
  new_blin_train = rep(0, length(y_train))
  new_blin_test  = rep(0, length(y_test))

  for (k in 1:PC) {
    new_blin_train = new_blin_train + new_lambda_hat[k]*new_gamma_hat[x_train[,'g'],k]*new_delta_hat[x_train[,'e'],k]
    new_blin_test  = new_blin_test + new_lambda_hat[k]*new_gamma_hat[x_test[,'g'],k]*new_delta_hat[x_test[,'e'],k]
  }

  # Compute the predicted values for the TRAINING and TEST data sets
  new_mu_ij_train = new_mu_hat + new_alpha_hat[x_train[,'g']] + new_beta_hat[x_train[,'e']] + new_blin_train
  new_mu_ij_test = new_mu_hat + new_alpha_hat[x_test[,'g']] + new_beta_hat[x_test[,'e']] + new_blin_test

  return(list(alpha_hat   = new_alpha_hat,
              beta_hat    = new_beta_hat,
              delta_hat   = new_delta_hat,
              gamma_hat   = new_gamma_hat,
              lambda_hat  = new_lambda_hat,
              y_hat_train = new_mu_ij_train,
              y_hat_test  = new_mu_ij_test,
              blinear_hat = new_blin_train,
              id          = 'Bayesian AMMI (postproc)'))
}

# plot(new_alpha_hat, ylim=c(-8,8), main='alpha - WITH postprocessing')
# points(data$alpha, col=2)
# points(alpha_hat[1,], col=3)
#
# plot(new_beta_hat, ylim=c(-8,8), main='beta - WITH postprocessing')
# points(data$beta, col=2)
# points(beta_hat[1,], col=3)
#
# plot(new_delta_hat, ylim=c(-2,2), main='delta - WITH postprocessing')
# points(data$delta, col=2)
# points(delta_hat[1,], col=3)
#
# plot(new_gamma_hat, ylim=c(-2,2), main='gamma - WITH postprocessing')
# points(data$gamma, col=2)
# points(gamma_hat[1,], col=3)
#' @export
#'
organise_AMBARTI <- function(object, data){

  # Get training info
  x_train = data$x
  y_train = data$y

  # Get test info
  x_test = data$x
  x_test$g = as.factor(x_test$g)
  x_test$e = as.factor(x_test$e)
  y_test = data$y_test

  # Get estimates info
  estimate = apply(object$beta_hat, 2, mean)
  alpha_hat = estimate[grepl('g', names(estimate))]
  beta_hat = estimate[grepl('e', names(estimate))]
  y_hat_train = apply(object$y_hat, 2, mean)
  y_hat_test = as.numeric(predict_ambarti(object, x_test, type = 'mean'))
  blinear_hat = apply(object$y_hat_bart,2,mean)


  return(list(alpha_hat   = alpha_hat,
              beta_hat    = beta_hat,
              y_hat_train = y_hat_train,
              y_hat_test  = y_hat_test,
              blinear_hat = blinear_hat,
              id          = 'AMBARTI'))

}

#' @export
#'
organise_bayesian_AMMI_WITHOUT_postprocessing <- function(object, data){

  # Get training info
  x_train = data$x
  y_train = data$y

  # Get test info
  x_test = data$x
  y_test = data$y_test

  # Get the number of PCs
  PC = length(data$lambda)

  # Get estimates info
  estimate   = object$BUGSoutput$mean
  mu_hat     = as.numeric(estimate$mu_all)
  alpha_hat  = estimate$alpha
  beta_hat   = estimate$beta
  delta_hat  = estimate$delta
  gamma_hat  = estimate$gamma
  lambda_hat = estimate$lambda

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

  return(list(alpha_hat  = alpha_hat,
              beta_hat    = beta_hat,
              delta_hat   = delta_hat,
              gamma_hat   = gamma_hat,
              lambda_hat  = lambda_hat,
              y_hat_train = y_hat_train,
              y_hat_test  = y_hat_test,
              blinear_hat = blin_train,
              id          = 'Bayesian AMMI (NO postproc)'))
}




###### Plot functions ##########
###### Plot functions ##########
###### Plot functions ##########
###### Plot functions ##########

#' @export
#'
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

  id = 'AMMI'
  alpha_hat  = data.frame(id = id, variable = paste('alpha[',g_names ,']', sep=''), value=alpha_hat, true=data$alpha)
  beta_hat   = data.frame(id = id, variable = paste('beta[',e_names ,']', sep=''), value=beta_hat, true=data$beta)
  lambda_hat = data.frame(id = id, variable = paste('lambda[',1:PC,']',sep=''), value=lambda_hat, true=data$lambda)
  gamma_hat  = data.frame(id = id, variable = paste('gamma[',g_names,',',rep(1:PC, each=length(g_names)),']',sep=''), value=as.numeric(gamma_hat), true=as.numeric(data$gamma))
  delta_hat  = data.frame(id = id, variable = paste('delta[',g_names,',',rep(1:PC, each=length(e_names)),']',sep=''), value=as.numeric(delta_hat), true=as.numeric(data$delta))

  aux_blinear_hat   = data.frame(id = id, variable = 'blinear', value = data$blinear - blin_train, true = data$blinear)
  aux_y_hat_train   = data.frame(id = id, variable = 'blinear', value = y_train - y_hat_train, true = data$y)
  aux_y_hat_test    = data.frame(id = id, variable = 'blinear', value = y_test - y_hat_test, true = data$y_test)

  return(list(alpha_hat   = alpha_hat,
              beta_hat    = beta_hat,
              delta_hat   = delta_hat,
              gamma_hat   = gamma_hat,
              lambda_hat  = lambda_hat,
              y_hat_train = aux_y_hat_train,
              y_hat_test  = aux_y_hat_test,
              blinear_hat = aux_blinear_hat,
              id          = id))
}

#' @export
#'
bAMMI_help_plot <- function(object, data){
# Get training info
x_train = data$x
y_train = data$y

# Get test info
x_test = data$x
y_test = data$y_test

# Get the number of PCs
PC = length(data$lambda)

# Get estimates info
estimate   = object$BUGSoutput$mean
mu_hat     = as.numeric(estimate$mu_all)
alpha_hat  = estimate$alpha
beta_hat   = estimate$beta
delta_hat  = estimate$delta
gamma_hat  = estimate$gamma
lambda_hat = estimate$lambda

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

# Get some MCMC info
nburn      = object$BUGSoutput$n.burnin
niter      = object$BUGSoutput$n.iter
npost      = niter - nburn
seq_burn   = seq(1, nburn, by=1)

# Get estimates info
id = 'Bayesian AMMI (no postproc)'
estimate   = as.data.frame(object$BUGSoutput$sims.matrix)
estimate$id= id
#mu_hat     = estimate[-seq_burn,c('id',colnames(estimate)[grepl('mu_all', colnames(estimate))])]
alpha_hat  = estimate[-seq_burn,c('id',colnames(estimate)[grepl('alpha', colnames(estimate))])]
beta_hat   = estimate[-seq_burn,c('id',colnames(estimate)[grepl('beta', colnames(estimate))])]
gamma_hat  = estimate[-seq_burn,c('id',colnames(estimate)[grepl('gamma', colnames(estimate))])]
delta_hat  = estimate[-seq_burn,c('id',colnames(estimate)[grepl('delta', colnames(estimate))])]
lambda_hat = estimate[-seq_burn,c('id',colnames(estimate)[grepl('lambda', colnames(estimate))])]

alpha_hat  = melt(alpha_hat, measure.vars = colnames(alpha_hat)[-1])
beta_hat   = melt(beta_hat, measure.vars = colnames(beta_hat)[-1])
lambda_hat = melt(lambda_hat, measure.vars = colnames(lambda_hat)[-1])
gamma_hat  = melt(gamma_hat, measure.vars = colnames(gamma_hat)[-1])
delta_hat  = melt(delta_hat, measure.vars = colnames(delta_hat)[-1])

alpha_hat$true  = rep(as.numeric(data[['alpha']]),  each=npost)
beta_hat$true   = rep(as.numeric(data[['beta']]),   each=npost)
lambda_hat$true = rep(as.numeric(data[['lambda']]), each=npost)
gamma_hat$true  = rep(as.numeric(data[['gamma']]),  each=npost)
delta_hat$true  = rep(as.numeric(data[['delta']]),  each=npost)

aux_blinear_hat   = data.frame(id = id, variable = 'blinear', value = data$blinear - blin_train, true = data$blinear)
aux_y_hat_train   = data.frame(id = id, variable = 'blinear', value = y_train - y_hat_train, true = data$y)
aux_y_hat_test    = data.frame(id = id, variable = 'blinear', value = y_test - y_hat_test, true = data$y_test)

return(list(alpha_hat   = alpha_hat,
            beta_hat    = beta_hat,
            delta_hat   = delta_hat,
            gamma_hat   = gamma_hat,
            lambda_hat  = lambda_hat,
            blinear_hat = aux_blinear_hat,
            y_hat_train = aux_y_hat_train,
            y_hat_test  = aux_y_hat_test))
}

#' @export
#' @importFrom reshape2 'melt'
bAMMI_help_plot_WITHPOS <- function(object, data){

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

  # Compute the bilinear term
  blin_train = matrix(0, nrow = npost, ncol = nrow(x_train))
  blin_test  = matrix(0, nrow = npost, ncol = nrow(x_test))

  for (k in 1:PC) {
    blin_train = blin_train + lambda_hat[,k] * gamma_hat[,x_train[,'g']]* delta_hat[,x_train[,'e']]
    blin_test  = blin_test +  lambda_hat[,k] * gamma_hat[,x_test[,'g']] * delta_hat[, x_test[,'e']]
  }

  # Compute the predicted response for the TRAINING data
  mu_ij = mu_hat + alpha_hat[,x_train[,'g']] + beta_hat[,x_train[,'e']] + blin_train
  colnames(mu_ij) = NULL

  # -------------------------------------------------------------
  # Start postprocessing
  # -------------------------------------------------------------

  n_gen = length(unique(x_train[,'g']))
  n_env = length(unique(x_train[,'e']))

  # Create matrices/lists to store the postprocessing results
  snew_mu_hat     = matrix(NA, nrow=npost, ncol=1)
  snew_alpha_hat  = matrix(NA, nrow=npost, ncol=n_gen)
  snew_beta_hat   = matrix(NA, nrow=npost, ncol=n_env)
  snew_lambda_hat = matrix(NA, nrow=npost, ncol=PC)
  snew_gamma_hat  = matrix(NA, nrow=npost, ncol=PC*n_gen)
  snew_delta_hat  = matrix(NA, nrow=npost, ncol=PC*n_env)

  for (i in 1:nrow(mu_ij)){

    # Get the matrix with mu_ij for each genotype and environment
    matrix_mu_ij  = matrix(mu_ij[i,], nrow=n_gen, ncol=n_env)
    new_mu_hat    = mean(matrix_mu_ij)
    new_alpha_hat = new_mu_hat - rowMeans(matrix_mu_ij) # their sum is 0
    new_beta_hat  = new_mu_hat - colMeans(matrix_mu_ij) # their sum is 0

    # Center matrix by row and column
    # Thank https://stackoverflow.com/questions/43639063/double-centering-in-r
    resA = t(matrix_mu_ij*0 + colMeans(matrix_mu_ij))
    resB = matrix_mu_ij*0 + rowMeans(matrix_mu_ij)
    res_double_centered = matrix_mu_ij - resA - resB + new_mu_hat # sum zero by row and column

    # Run the Singular Value Decomposition (SVD) to compute lambda, gamma, and delta
    sv_dec <- svd(res_double_centered, nu = PC, nv = PC)

    # Get new parameter estimates
    new_lambda_hat = sv_dec$d[1:PC]
    new_gamma_hat  = -1*sv_dec$u
    new_delta_hat  = -1*sv_dec$v

    # Store the new results
    snew_mu_hat[i,]     = new_mu_hat
    snew_alpha_hat[i,]  = new_alpha_hat
    snew_beta_hat[i,]   = new_beta_hat
    snew_lambda_hat[i,] = new_lambda_hat
    snew_gamma_hat[i,]  = as.numeric(new_gamma_hat)
    snew_delta_hat[i,]  = as.numeric(new_delta_hat)
  }
  # Clean column names
  colnames(snew_mu_hat)     = colnames(mu_hat)
  colnames(snew_alpha_hat)  = colnames(alpha_hat)
  colnames(snew_beta_hat)   = colnames(beta_hat)
  colnames(snew_lambda_hat) = colnames(lambda_hat)
  colnames(snew_gamma_hat)  = colnames(gamma_hat)
  colnames(snew_delta_hat)  = colnames(delta_hat)

  # Summarise the new posterior results
  new_mu_hat     = apply(snew_mu_hat     ,2,mean)
  new_alpha_hat  = apply(snew_alpha_hat  ,2,mean)
  new_beta_hat   = apply(snew_beta_hat   ,2,mean)
  new_lambda_hat = apply(snew_lambda_hat ,2,mean)
  new_gamma_hat  = apply(snew_gamma_hat,  2,mean)
  new_delta_hat  = apply(snew_delta_hat,  2,mean)

  # Compute the bilinear term
  new_blin_train = rep(0, length(y_train))
  new_blin_test  = rep(0, length(y_test))

  aux_new_gamma_hat = matrix(new_gamma_hat,ncol=PC)
  aux_new_delta_hat = matrix(new_delta_hat,ncol=PC)

  for (k in 1:PC) {
    new_blin_train = new_blin_train + new_lambda_hat[k]*aux_new_gamma_hat[x_train[,'g'],k]*aux_new_delta_hat[x_train[,'e'],k]
    new_blin_test  = new_blin_test + new_lambda_hat[k]*aux_new_gamma_hat[x_test[,'g'],k]*aux_new_delta_hat[x_test[,'e'],k]
  }

  # Compute the predicted values for the TRAINING and TEST data sets
  new_mu_ij_train = new_mu_hat + new_alpha_hat[x_train[,'g']] + new_beta_hat[x_train[,'e']] + new_blin_train
  new_mu_ij_test = new_mu_hat + new_alpha_hat[x_test[,'g']] + new_beta_hat[x_test[,'e']] + new_blin_test

  id = 'Bayesian AMMI (posproc)'
  snew_alpha_hat  = as.data.frame(snew_alpha_hat)
  snew_beta_hat   = as.data.frame(snew_beta_hat)
  snew_lambda_hat = as.data.frame(snew_lambda_hat)
  snew_gamma_hat  = as.data.frame(snew_gamma_hat)
  snew_delta_hat  = as.data.frame(snew_delta_hat)

  snew_alpha_hat$id  = id
  snew_beta_hat$id   = id
  snew_lambda_hat$id = id
  snew_gamma_hat$id  = id
  snew_delta_hat$id  = id

  snew_alpha_hat  = melt(snew_alpha_hat,  measure.vars  = colnames(snew_alpha_hat)[-length(colnames(snew_alpha_hat))])
  snew_beta_hat   = melt(snew_beta_hat,   measure.vars  = colnames(snew_beta_hat)[-length(colnames(snew_beta_hat))])
  snew_lambda_hat = melt(snew_lambda_hat, measure.vars  = colnames(snew_lambda_hat)[-length(colnames(snew_lambda_hat))])
  snew_gamma_hat  = melt(snew_gamma_hat,  measure.vars  = colnames(snew_gamma_hat)[-length(colnames(snew_gamma_hat))])
  snew_delta_hat  = melt(snew_delta_hat,  measure.vars  = colnames(snew_delta_hat)[-length(colnames(snew_delta_hat))])

  snew_alpha_hat$true   = rep(as.numeric(data[['alpha']]),  each=npost)
  snew_beta_hat$true    = rep(as.numeric(data[['beta']]),   each=npost)
  snew_lambda_hat$true  = rep(as.numeric(data[['lambda']]),  each=npost)
  snew_gamma_hat$true   = rep(as.numeric(data[['gamma']]),  each=npost)
  snew_delta_hat$true   = rep(as.numeric(data[['delta']]), each=npost)

  aux_blinear_hat   = data.frame(id = id, variable = 'blinear', value = data$blinear - new_blin_train, true = data$blinear)
  aux_y_hat_train   = data.frame(id = id, variable = 'blinear', value = y_train - new_mu_ij_train, true = data$y)
  aux_y_hat_test    = data.frame(id = id, variable = 'blinear', value = y_test - new_mu_ij_test, true = data$y_test)


  return(list(alpha_hat   = snew_alpha_hat,
              beta_hat    = snew_beta_hat,
              delta_hat   = snew_delta_hat,
              gamma_hat   = snew_gamma_hat,
              lambda_hat  = snew_lambda_hat,
              y_hat_train = aux_y_hat_train,
              y_hat_test  = aux_y_hat_test,
              blinear_hat = aux_blinear_hat,
              id          = id))
}
#' @export

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
  id = 'AMBARTI'
  estimate    = as.data.frame(object$beta_hat)
  alpha_hat   = estimate[,grepl('g', names(estimate))]
  beta_hat    = estimate[,grepl('e', names(estimate))]
  y_hat_train = apply(object$y_hat, 2, mean)
  y_hat_test  = as.numeric(predict_ambarti(object, x_test, type = 'mean'))
  blinear_hat = apply(object$y_hat_bart,2,mean)
  names(alpha_hat) = paste(gsub('g','alpha[', names(alpha_hat)), ']', sep='')
  names(beta_hat) = paste(gsub('e','beta[', names(beta_hat)), ']', sep='')
  alpha_hat$id = id
  beta_hat$id = id

  alpha_hat = melt(alpha_hat, measure.vars = colnames(alpha_hat)[grepl('alpha', colnames(alpha_hat))])
  beta_hat  = melt(beta_hat, measure.vars = colnames(beta_hat)[grepl('beta', colnames(beta_hat))])

  alpha_hat$true = rep(as.numeric(data[['alpha']]), each=object$npost)
  beta_hat$true  = rep(as.numeric(data[['beta']]), each=object$npost)

  aux_blinear_hat   = data.frame(id = id, variable = 'blinear', value = data$blinear - blinear_hat, true = data$blinear)
  aux_y_hat_train   = data.frame(id = id, variable = 'blinear', value = y_train - y_hat_train, true = data$y)
  aux_y_hat_test    = data.frame(id = id, variable = 'blinear', value = y_test - y_hat_test, true = data$y_test)


  return(list(alpha_hat   = alpha_hat,
              beta_hat    = beta_hat,
              y_hat_train = aux_y_hat_train,
              y_hat_test  = aux_y_hat_test,
              blinear_hat = aux_blinear_hat,
              id          = 'AMBARTI'))
}

#' @export
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