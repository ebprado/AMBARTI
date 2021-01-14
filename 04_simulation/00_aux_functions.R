organise_classical_AMMI <- function(object){

  classical_AMMI$ANOVA
  classical_AMMI$genXenv
  classical_AMMI$analysis
  classical_AMMI$means
  classical_AMMI$biplot
  train_data$y
}

organise_bayesian_AMMI <- function(object, train_data){

  x_train = train_data$x
  y_train = train_data$y

  estimate = bayesian_AMMI$BUGSoutput$mean

  mu_hat = estimate$mu_all
  alpha_hat = estimate$alpha
  beta_hat = estimate$beta
  delta_hat = estimate$delta
  gamma_hat = estimate$gamma
  lambda_hat = estimate$lambda

  y_hat =

}