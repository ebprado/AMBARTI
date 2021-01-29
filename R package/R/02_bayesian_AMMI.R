#' @export
#' @importFrom R2jags 'jags'
#'

run_bayesian_AMMI <- function(data,
                              Q,
                              s_mu = 20,
                              s_lambda = 10,
                              s_alpha_hyperpar = 10,
                              s_beta_hyperpar = 10,
                              S_ME = 2,
                              n.thin = 1,
                              n.burnin = 2000,
                              n.iter = 3000){

  # Specify the Bayesian AMMI model similar to Josse et al (JABES, 2014)
  model_code = '
  model
  {
  # Likelihood
  for (k in 1:N) {
  Y[k] ~ dnorm(mu[k], sigma_E^-2)
  mu[k] = mu_all + alpha[genotype[k]] + beta[environment[k]] + sum(lambda * gamma[genotype[k],1:Q] * delta[environment[k],1:Q])
  }
  # Priors
  mu_all ~ dnorm(0, s_mu^-2) # Prior on grand mean
  for(i in 1:I) {
  alpha[i] ~ dnorm(0, s_alpha^-2) # Prior on genotype effect
  }
  for(j in 1:J) {
  beta[j] ~ dnorm(0, s_beta^-2) # Prior on environment effect
  }
  # Priors on gamma
  for(q in 1:Q) {
  gamma[1, q] ~ dnorm(0, 1)T(0,) # First one is restriced to be positive
  for(i in 2:I) {
  gamma[i, q] ~ dnorm(0, 1) # Prior on genotype interactions
  }
  }
  # Priors on delta
  for(q in 1:Q) {
  for(j in 1:J) {
  delta[j, q] ~ dnorm(0, 1) # Prior on environment interactions
  }
  }
  # Prior on eigenvalues
  for(q in 1:Q) {
  lambda_raw[q] ~ dnorm(0, s_lambda^-2)T(0,)
  }
  lambda = sort(lambda_raw)
  # Prior on residual standard deviation
  sigma_E ~ dunif(0, S_ME)
  }
  '

  # Get the quantities needed in the JAGS model list
  N = data$I * data$J
  Y = data$y
  I = data$I
  J = data$J
  if (is.null(data$Q) == FALSE) {Q = data$Q}
  genotype = data$x[,'g']
  environment = data$x[,'e']

  # Set up the data
  model_data = list(N = N,
                    Y = Y,
                    I = I,
                    J = J,
                    Q = Q,
                    genotype = genotype,
                    environment = environment,
                    s_mu = s_mu,
                    s_alpha = s_alpha_hyperpar,
                    s_beta = s_beta_hyperpar,
                    s_lambda = s_lambda,
                    S_ME = S_ME)

  # Choose the parameters to watch
  model_parameters =  c("alpha", "beta", "lambda", "gamma", "delta",
                        'sigma_E', 'mu_all')

  # Run the model
  model_run = jags(data = model_data,
                   parameters.to.save = model_parameters,
                   model.file=textConnection(model_code),
                   progress.bar = 'none',
                   n.thin = n.thin,
                   n.burnin = n.burnin,
                   n.iter = n.iter)

  return(model_run)

}
