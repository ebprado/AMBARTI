# AMBARTI models for agricultural experiments

This repository contains R scripts and data sets that can be used to fit the Additive Main Effect Bayesian Additive Regression Tree interaction (AMBARTI) models for agricultural experiments.

In addition, it provides an implementation of AMBARTI in the format of an R package named ```AMBARTI```.

## Installation
``` r
library(devtools)
install_github("ebprado/AMBARTI")
```
Below, we generate a simulated example of the Bayesian version of the Additive Main effects and Multiplicative Interaction (AMMI) effects model presented in the [Josse et al (JABES, 2014)](https://link.springer.com/content/pdf/10.1007/s13253-014-0168-z.pdf). The simulation code was written by [Andrew Parnell](https://github.com/andrewcparnell) and [Danilo Sarti](https://github.com/danilosarti).

## Example
``` r
rm(list = ls())

# Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file

# Likelihood
# Y_{ij} ~ N(mu_{ij}, sigma^2_E)
# with
# mu_{ij} = mu + alpha_i + beta_j + sum_{q=1}^Q lambda_q*gamma_iq*delta_jq
# Our idea is to estimate alpha_i + beta_j parametrically and the component 
# "sum_{q=1}^Q lambda_q*gamma_iq*delta_jq" via BART.

# Notation
# Y_ij = response (e.g. yield) for genotype i and environment j, i = 1, ..., I
# genotypes and j = 1, ..., J environments
# mu is the grand mean
# alpha_i is the genotype effect
# beta_j is the environment effect
# lambda_q is the q-th eigenvalue q = 1,.., Q of the interaction matrix
# Q is the number of components used to model the interaction. Usually Q is fixed 
# at a small number, e.g. 2
# gamma_{iq} is the interaction effect for the q-th eigenvector for genotype i
# delta_{iq} is the interaction effect for the q-th eigenvector for environment j
# E_{ij} is a residual term with E_{ij} ~ N(0, sigma^2_E)
# Usually these models have quite complicated restrictions on the gamma/delta/lambda
# values but Josse et al show that these are not fully necessary

# Priors
# alpha_i ~ N(0, s_alpha^2)
# beta_j ~ N(0, s_beta^2)

# Simulate data -----------------------------------------------------------

# We will follow the simulation strategy detailed in Section 3.1 of the
# Josse et al paper

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
fit.ambarti = ambarti(x.ambarti, y, ntrees = 50, skip_trees = FALSE, nburn = 100, npost = 100, sparse= FALSE)

# Get the final prediction (y hat)
yhat_ambarti = apply(fit.ambarti$y_hat, 2, mean)
yhat_ambarti2 = predict_ambarti(fit.ambarti, newdata = x.ambarti, type = 'mean')
cor(y, yhat_ambarti);

# Get the prediction specifically from BART
yhat_bart = apply(fit.ambarti$y_hat_bart, 2, mean);
cor(y, yhat_bart); # correlation btw y and BART (AMBARTI package)

# Plot the main effects estimates and add the true values
plot(1:length(alpha), alpha, col=2, cex=2, main='AMBARTI - Genotype', ylim= c(-5,5)) # true values
points(apply(fit.ambarti$beta_hat[,1:10], 2, mean), cex=2, pch = 2) # estimates
legend(8,4,'AMBARTI', col=1, pch = 2)
legend(8,5,'True', col=2, pch = 1, cex=1)

# Plot the main effects estimates and add the true values
plot(1:length(beta), beta, col=2, cex=2, main='AMBARTI - Environment', ylim = c(-5,5)) # true values
points(apply(fit.ambarti$beta_hat[,11:20], 2, mean), cex=2) # estimates
legend(8,4,'AMBARTI', col=1, pch = 2)
legend(8,5,'True', col=2, pch = 1, cex=1)

fit.ambarti$trees[[100]][[1]] # shows the tree 1 in the last (100) MCMC iteration.

# ---------------------------------------
# BART (just to have a benchmark)
# ---------------------------------------
library(BART)
bart = BART::wbart(x.ambarti, y)
cor(y, bart$yhat.train.mean) # BART and semibart are quite similar. That's fine.
```
