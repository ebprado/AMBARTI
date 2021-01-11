library(devtools)
check()
document()
build()
install()

install_github("ebprado/semibart") # I just added the BART predictions to the output
library(semibart)

install_github("ebprado/AMBARTI") # This my implementation of the semibart idea
library(AMBARTI)

# Simulate data -----------------------------------------------------------

# We will follow the simulation strategy detailed in Section 3.1 of the
# Josse et al paper
set.seed(001)
# Specify fixed values
Q = 1 # Number of components
I = 10 # Number of genotypes
J = 10# Number of environments
N = I*J # Total number of obs

# Some further fixed values
mu = 10
sigma_E = 1
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
mu_ij = mu + alpha[G_by_E[,1]] + beta[G_by_E[,2]]  +
  lambda_1 * gamma[G_by_E[,1]] * delta[G_by_E[,2]] ## maybe insert lambda2
Y = rnorm(N, mu_ij, sigma_E) ## response variable

# ---------------------------------------
# AMBARTI
# ---------------------------------------

# Some pre-processing
x.ambarti = G_by_E
names(x.ambarti) = c('g', 'e')
x.ambarti$g = as.factor(x.ambarti$g)
x.ambarti$e = as.factor(x.ambarti$e)
y = Y
set.seed(101)

# Run AMBARTI
fit.ambarti = ambarti(x.ambarti, y, ntrees = 50, nburn = 100, npost = 100, sparse= FALSE)

# Get the final prediction (y hat)
yhat_ambarti = apply(fit.ambarti$y_hat, 2, mean)
yhat_ambarti2 = predict_ambarti(fit.ambarti, newdata = x.ambarti, type = 'mean')
cor(y, yhat_ambarti);

# Get the prediction specifically from BART
yhat_bart = apply(fit.ambarti$y_hat_bart, 2, mean);
cor(y, yhat_bart); # correlation btw y and the BART component (from AMBARTI)

# Plot the main effects estimates and add the true values
alpha_hat = apply(fit.ambarti$beta_hat[,1:10], 2, mean)
plot(1:length(alpha), alpha, col=2, cex=2, main='AMBARTI-Genotype', ylim=c(-5,5), ylab=expression(g[i]), xlab = expression(i)) # true values
points(alpha_hat, cex=2, pch = 2) # estimates
legend(8,4,'AMBARTI', col=1, pch = 2, bty = 'n')
legend(8,5,'True', col=2, pch = 1, cex=1, bty = 'n')

# Plot the main effects estimates and add the true values
beta_hat = apply(fit.ambarti$beta_hat[,11:20], 2, mean)
plot(1:length(beta), beta, col=2, cex=2, main='AMBARTI-Environment', ylim=c(-5,5), ylab=expression(e[j]), xlab = expression(j)) # true values
points(beta_hat, cex=2, pch = 2) # estimates
legend(8,4,'AMBARTI', col=1, pch = 2, bty='n')
legend(8,5,'True', col=2, pch = 1, cex=1, bty='n')

fit.ambarti$trees[[100]][[1]] # show the first tree in the last (100) MCMC iteration.

# ---------------------------------------
# AMMI
# ---------------------------------------

model_code = '
model
{
  # Likelihood
  for (k in 1:N) {
  Y[k] ~ dnorm(mu[k], sigma_E^-2)
  mu[k] = mu_all + alpha[genotype[k]] + beta[environment[k]] + lambda_1 * gamma[genotype[k]] * delta[environment[k]]
  }
  # Priors
  mu_all ~ dnorm(0, s_mu^-2) # Prior on grand mean
  for(i in 1:I) {
  alpha[i] ~ dnorm(0, s_alpha^-2) # Prior on genotype effect
  }
  gamma[1] ~ dnorm(0, 1)T(0,) # First one is restriced to be positive
  for(i in 2:I) {
  gamma[i] ~ dnorm(0, 1) # Prior on genotype interactions
  }
  for(j in 1:J) {
  beta[j] ~ dnorm(0, s_beta^-2) # Prior on environment effect
  delta[j] ~ dnorm(0, 1) # Prior on environment interactions
  }
  # Prior on first (and only) eigenvalue
  lambda_1 ~ dnorm(0, s_lambda^-2)T(0,)
  # Prior on residual standard deviation
  sigma_E ~ dunif(0, S_ME)
}
'
s_mu = 20
s_alpha = 10
s_beta = 10
s_lambda = 10
S_ME = 10
# Set up the data
model_data = list(N = N,
                  Y = Y,
                  I = I,
                  J = J,
                  genotype = G_by_E[,1],
                  environment = G_by_E[,2],
                  s_mu = s_mu,
                  s_alpha = s_alpha,
                  s_beta = s_beta,
                  s_lambda = s_lambda,
                  S_ME = S_ME)

# Choose the parameters to watch
model_parameters =  c("alpha", "beta", "lambda_1", "gamma", "delta",
                      'sigma_E')

# Run the model
model_run = jags(data = model_data,n.chains=5,
                 parameters.to.save = model_parameters,
                 model.file=textConnection(model_code))

#
post_means = model_run$BUGSoutput$mean

# Plot the main effects estimates and add the true values
alpha_hat = post_means$alpha
plot(1:length(alpha), alpha, col=2, cex=2, main='AMMI - Genotype', ylim= c(-4,4)) # true values
points(alpha_hat, cex=2, pch = 3, col=4) # estimates
legend(8,3,'AMMI', pch = 3, col=4, bty='n')
legend(8,4,'True', col=2, pch = 1, cex=1)

beta_hat = post_means$beta
plot(1:length(beta), beta, col=2, cex=2, main='AMMI - Environment', ylim= c(-4,4)) # true values
points(beta_hat, cex=2, pch = 3, col=4) # estimates
legend(8,3,'AMMI', pch = 3, col=4, bty='n')
legend(8,4,'True', col=2, pch = 1, cex=1, bty='n')

delta_hat = post_means$delta
gamma_hat = post_means$gamma
lambda_hat = as.numeric(post_means$lambda)

yhat_ammi = mu + alpha_hat[G_by_E[,1]] + beta_hat[G_by_E[,2]]  +
  lambda_hat * gamma_hat[G_by_E[,1]] * delta_hat[G_by_E[,2]]

cor(y, yhat_ammi)

##################################
# BART (just to have a benchmark)
##################################
library(BART)
bart = BART::wbart(x.ambarti, y)
cor(y, bart$yhat.train.mean) # BART and semibart are quite similar. That's fine.

##################################
# Semiparametric BART
##################################

# Some pre-processing
x = G_by_E
names(x) = c('g', 'e')
x$g = as.factor(x$g)
x$e = as.factor(x$e)
y = Y

cov_g = x[,'g']
cov_e = x[,'e']

classes_g = sort(unique(cov_g))
classes_e = sort(unique(cov_e))

ng = tapply(cov_g, cov_g, length)
ne = tapply(cov_e, cov_e, length)

x <- model.matrix(~ -1 + g + e, data=x,
                  contrasts.arg=list(g=contrasts(as.factor(x$g), contrasts=F),
                                     e=contrasts(as.factor(x$e), contrasts=F)))
set.seed(001)

# Run Semiparametric BART
semib = semibart(x.train = x, y.train = y, a.train = x)

# Get the main effects estimates
betahat = apply(semib$beta,2,mean)[1:10] # The first 10 are associated to the covariate g (genotype)

# Plot the main effects estimates and add the true values
plot(betahat, cex=2, ylim = c(-5,5), main='Genotype - semibart') # estimates (black)
points(1:length(alpha), alpha, col=2, cex=2) # true values (red). Looks not too bad.
legend(8,4,'semi BART', col=1, pch = 2)
legend(8,5,'True', col=2, pch = 1, cex=1)

# Plot the main effects estimates and add the true values
alphahat = apply(semib$beta,2,mean)[11:20] # The remaining 10 are associated to the covariate e (environment)
plot(alphahat, cex=2, ylim = c(-5,5), main='Environment - semibart') # estimates
points(1:length(beta), beta, col=2, cex=2) # true values. Looks fine.
legend(8,4,'semi BART', col=1, pch = 2)
legend(8,5,'True', col=2, pch = 1, cex=1)

# Correlation btw y and BART estimate
  cor(y, apply(semib$bartfit, 2, mean)) # ~0.31

# Compute the final prediction (y hat)
yhat = x%*%apply(semib$beta,2,mean) + apply(semib$bartfit, 2, mean)
plot(y, yhat, main = 'semibart - y versus y hat'); abline(0,1) # Looks fine
cor(y, yhat); # ~0.89
