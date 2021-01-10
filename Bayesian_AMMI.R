
rm(list = ls())
library(R2jags)
library(ggplot2)
library(tidyverse)
# Simulate data -----------------------------------------------------------

# We will follow the simulation strategy detailed in Section 3.1 of the Josse et al paper

# Specify fixed values
Q = 1 # Number of components
I = 10 # Number of genotypes
J = 10# Number of environments
N = I*J # Total number of obs
m = 90
s_mu = 20
s_alpha = 10
s_beta = 10
s_lambda = 10
S_ME = 10

# Some further fixed values
mu = 10
sigma_E = 3/2 # Not sure why S_ME was specified if they're also giving sigma_E
sigma_alpha = 2
sigma_beta = 2
alpha = rnorm(I, 0, sigma_alpha)
beta = rnorm(J, 0, sigma_beta)
lambda_1 = 12
gamma = seq(2, -2,length.out = I)/sqrt(10)
delta = seq(-0.5, 0.5,length.out = J)

# Now simulate the values
set.seed(123)
G_by_E = expand.grid(1:I, 1:J) ## setting the interaction matrix
mu_ij = mu + alpha[G_by_E[,1]] + beta[G_by_E[,2]] + lambda_1 * gamma[G_by_E[,1]] * delta[G_by_E[,2]] ## maybe insert lambda2
Y = rnorm(N, mu_ij, sigma_E) ## response variable

# Can create some plots
qplot(x = G_by_E[,1], y = Y, geom = 'boxplot', group = G_by_E[,1], xlab = 'Genotype')
qplot(x = G_by_E[,2], y = Y, geom = 'boxplot', group = G_by_E[,2], xlab = 'Environment')

# Second model - general Q ------------------------------------------------

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

# Set up the data
model_data = list(N = N,
                  Y = Y,
                  I = I,
                  J = J,
                  Q = 2, # Set Q to be 2 even though the simulation was for Q = 1
                  genotype = G_by_E[,1],
                  environment = G_by_E[,2],
                  s_mu = s_mu,
                  s_alpha = s_alpha,
                  s_beta = s_beta,
                  s_lambda = s_lambda,
                  S_ME = S_ME)

# Choose the parameters to watch
model_parameters =  c("alpha", "beta", "lambda", "gamma", "delta",
                      'sigma_E')

# Run the model
model_run = jags(data = model_data,
                 parameters.to.save = model_parameters,
                 model.file=textConnection(model_code))

# Plot the results
print(model_run)

#
post_means = model_run$BUGSoutput$mean

# Plot the main effects estimates and add the true values
alpha_hat = post_means$alpha
plot(1:length(alpha), alpha, col=2, cex=2, main='AMMI - Genotype', ylim= c(-4,4)) # true values
points(alpha_hat, cex=2, pch = 3, col=4) # estimates
legend(8,3,'AMMI', pch = 3, col=4)
legend(8,4,'True', col=2, pch = 1, cex=1)

beta_hat = post_means$beta
plot(1:length(beta), beta, col=2, cex=2, main='AMMI - Environment', ylim= c(-4,4)) # true values
points(beta_hat, cex=2, pch = 3, col=4) # estimates
legend(8,3,'AMMI', pch = 3, col=4)
legend(8,4,'True', col=2, pch = 1, cex=1)


