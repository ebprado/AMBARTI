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
#' @importFrom stats 'aov' 'model.tables'
organise_classical_AMMI <- function(object, data){

  # Get training info

  x_train = data$x
  y_train = data$y

  # Get test info
  x_test = data$x
  y_test = data$y_test

  gen = as.factor(x_train[,"gen"])
  env = as.factor(x_train[,"env"])

  # Fit the linear model
  linear_mod = aov(y_train ~ gen + env + gen:env)
  # linear_mod = lm(y_train ~ gen + env)

  # Get the residuals for the interaction gen:env
  #interaction_tab = matrix(residuals(linear_mod), ncol = length(unique(gen)), length(unique(env)))
  interaction_tab = model.tables(linear_mod, type='effects', cterms = 'gen:env')
  interaction_tab = interaction_tab$tables$`gen:env`

  # Get the number of PCs
  PC = length(data$lambda)

  # Run the Singular Value Decomposition (SVD) to compute lambda, gamma, and delta
  sv_dec <- svd(interaction_tab, nu = PC, nv = PC)

  # Get parameter estimates
  # mu_hat     = linear_mod$coefficients[1] # slightly biased compared to mean(y_train)
  mu_hat     = mean(y_train)
  alpha_hat  = aggregate(x = y_train - mu_hat, by = list(gen), FUN = "mean")[,2]
  beta_hat   = aggregate(x = y_train - mu_hat, by = list(env), FUN = "mean")[,2]
  lambda_hat = sv_dec$d[1:PC]
  gamma_hat  = -1*sv_dec$u
  delta_hat  = -1*sv_dec$v

  # Set up the bilinear term
  blin_train = rep(0, length(y_train))
  blin_test  = rep(0, length(y_test))

  for (k in 1:length(lambda_hat)) {
    blin_train = blin_train + lambda_hat[k]*gamma_hat[x_train[,'gen'],k]*delta_hat[x_train[,'env'],k]
    blin_test  = blin_test + lambda_hat[k]*gamma_hat[x_test[,'gen'],k]*delta_hat[x_test[,'env'],k]
  }

  # Compute the predicted response for the TRAINING data
  y_hat_train = mu_hat + alpha_hat[x_train[,'gen']] + beta_hat[x_train[,'env']] + blin_train

  # Compute the predicted response for the TEST data
  y_hat_test = mu_hat + alpha_hat[x_test[,'gen']] + beta_hat[x_test[,'env']] + blin_test

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

# plot(y_train, bb$y_hat_train);abline(0,1)
# cor(y_train, bb$y_hat_train)
#
# plot(y_test, bb$y_hat_test);abline(0,1)
# cor(y_test, bb$y_hat_test)
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
    blin_train = blin_train + lambda_hat[,k] * gamma_hat[,x_train[,'gen']]* delta_hat[,x_train[,'env']]
    blin_test  = blin_test +  lambda_hat[,k] * gamma_hat[,x_test[,'gen']] * delta_hat[, x_test[,'env']]
  }

  # Compute the predicted response for the TRAINING data
  mu_ij = mu_hat + alpha_hat[,x_train[,'gen']] + beta_hat[,x_train[,'env']] + blin_train
  colnames(mu_ij) = NULL

  # -------------------------------------------------------------
  # Start postprocessing
  # -------------------------------------------------------------

  n_gen = length(unique(x_train[,'gen']))
  n_env = length(unique(x_train[,'env']))

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
    new_blin_train = new_blin_train + new_lambda_hat[k]*new_gamma_hat[x_train[,'gen'],k]*new_delta_hat[x_train[,'env'],k]
    new_blin_test  = new_blin_test + new_lambda_hat[k]*new_gamma_hat[x_test[,'gen'],k]*new_delta_hat[x_test[,'env'],k]
  }

  # Compute the predicted values for the TRAINING and TEST data sets
  new_mu_ij_train = new_mu_hat + new_alpha_hat[x_train[,'gen']] + new_beta_hat[x_train[,'env']] + new_blin_train
  new_mu_ij_test = new_mu_hat + new_alpha_hat[x_test[,'gen']] + new_beta_hat[x_test[,'env']] + new_blin_test

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
  names(x_test) = c('g', 'e')
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
    blin_train = blin_train + lambda_hat[k]*gamma_hat[x_train[,'gen'],k]*delta_hat[x_train[,'env'],k]
    blin_test  = blin_test + lambda_hat[k]*gamma_hat[x_test[,'gen'],k]*delta_hat[x_test[,'env'],k]
  }

  # Compute the predicted response for the TRAINING data
  y_hat_train = mu_hat + alpha_hat[x_train[,'gen']] + beta_hat[x_train[,'env']] + blin_train
  # Compute the predicted response for the TEST data
  y_hat_test = mu_hat + alpha_hat[x_test[,'gen']] + beta_hat[x_test[,'env']] + blin_test

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
