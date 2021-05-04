# AMBARTI models for agricultural experiments

This repository contains R scripts and data sets that can be used to reproduce the simulation and real datasets results in Danilo A. Sarti, Estevão B. Prado, Alan N. Inglis, Antônia A. L. dos Santos, Catherine B. Hurley, Rafael A. Moral, Andrew C. Parnel. _Bayesian Additive Regression Trees for Genotype by Environment Interaction Models_. biorxiv (2021).

In addition, it provides an implementation of AMBARTI in the format of an R package named ```AMBARTI```.

## Installation
``` r
library(devtools)
install_github("ebprado/AMBARTI/R package")
```
Below, we fit AMBARTI to a simulated example from the equation of the Bayesian Additive Main effects and Multiplicative Interaction (AMMI) effects model presented in [Josse et al (JABES, 2014)](https://link.springer.com/content/pdf/10.1007/s13253-014-0168-z.pdf).

## Example
``` r
library(AMBARTI)
rm(list = ls())

# Simulate data -----------------------------------------------------------

# Specify fixed values
Q = 1 # Number of components
I = 10 # Number of genotypes
J = 10 # Number of environments
N = I*J # Total number of obs

# Some further fixed values
mu = 10 # Grand mean
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
plot(1:length(alpha), alpha, col=2, cex=2, main='AMBARTI-Genotype', ylim=c(-5,5)) # true values
points(alpha_hat, cex=2, pch = 2) # estimates
legend(8,4,'AMBARTI', col=1, pch = 2)
legend(8,5,'True', col=2, pch = 1, cex=1)

# Plot the main effects estimates and add the true values
beta_hat = apply(fit.ambarti$beta_hat[,11:20], 2, mean)
plot(1:length(beta), beta, col=2, cex=2, main='AMBARTI-Environment', ylim=c(-5,5)) # true values
points(beta_hat, cex=2) # estimates
legend(8,4,'AMBARTI', col=1, pch = 2)
legend(8,5,'True', col=2, pch = 1, cex=1)

fit.ambarti$trees[[100]][[1]] # show the first tree in the last (100) MCMC iteration.

# ---------------------------------------
# BART (just to have a benchmark)
# ---------------------------------------
library(BART)
bart = BART::wbart(x.ambarti, y)
cor(y, bart$yhat.train.mean)
```
