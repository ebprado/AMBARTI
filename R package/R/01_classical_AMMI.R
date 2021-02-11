#' @export
run_classical_AMMI <- function(data, Q){

  x_train <- data$x
  y_train <- data$y

  g = as.factor(x_train[,"g"])
  e = as.factor(x_train[,"e"])

  # Fit the linear model
  linear_mod = aov(y_train ~ g + e + g:e)
  # linear_mod = lm(y_train ~ g + e)

  # Get the residuals for the interaction g:e
  #interaction_tab = matrix(residuals(linear_mod), ncol = length(unique(g)), length(unique(env)))
  interaction_tab = model.tables(linear_mod, type='effects', cterms = 'g:e')
  interaction_tab = interaction_tab$tables$`g:e`

  # Get the number of PCs
  if (is.null(data$Q) == FALSE) {Q = data$Q}

  # Run the Singular Value Decomposition (SVD) to compute lambda, gamma, and delta
  sv_dec <- svd(interaction_tab, nu = Q, nv = Q)

  # Get parameter estimates
  # mu_hat     = linear_mod$coefficients[1] # slightly biased compared to mean(y_train)
  mu_hat     = mean(y_train)
  alpha_hat  = aggregate(x = y_train - mu_hat, by = list(g), FUN = "mean")[,2]
  beta_hat   = aggregate(x = y_train - mu_hat, by = list(e), FUN = "mean")[,2]
  lambda_hat = sv_dec$d[1:Q]
  gamma_hat  = -1*sv_dec$u
  delta_hat  = -1*sv_dec$v

  return(list(mu_hat     = mu_hat,
              alpha_hat  = alpha_hat,
              beta_hat   = beta_hat,
              lambda_hat = lambda_hat,
              gamma_hat  = gamma_hat,
              delta_hat  = delta_hat))
}
