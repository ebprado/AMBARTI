# AMBARTI models for agricultural experiments

This repository contains R scripts and data sets that can be used to reproduce the simulation and real datasets results in [Sarti*, D.A., Prado*, E.B., Inglis, A.N., dos Santos, A.A.L, Hurley, C., Moral, R.A. \& Parnell, A.C. _Bayesian additive regression trees for genotype by environment interaction models_. Annals of Applied Statistics, 17, 1 (2023)](https://www.biorxiv.org/content/10.1101/2021.05.07.442731v5.full.pdf). 

  * joint first authors.

In addition, it provides an implementation of AMBARTI in the format of an R package named ```AMBARTI```.

## Installation
``` r
library(devtools)
install_github("ebprado/AMBARTI/R package", ref='main')
```
Below, we fit AMBARTI to a simulated example from the equation of the Bayesian Additive Main effects and Multiplicative Interaction (AMMI) effects model presented in [Josse et al (JABES, 2014)](https://link.springer.com/content/pdf/10.1007/s13253-014-0168-z.pdf).

## Example
``` r
library(AMBARTI)
rm(list = ls())

# Simulate data -----------------------------------------------------------
I          = 10 # Number of genotypes
J          = 10 # Number of environments
s_alpha    = 1 # standard deviation of alpha
s_beta     = 1 # standard deviation of alpha
s_y        = 1 # standard deviation of y
lambda     = c(10, 8) # values for lambda

# Set a seed to make it reproducible
set.seed(001)

# Generate data from the AMMI model  
data = generate_data_AMMI(I, J, s_alpha, s_beta, s_y, lambda)
data_test = generate_data_AMMI(I, J, s_alpha, s_beta, s_y, lambda)

# run classical AMMI
classical_AMMI = run_classical_AMMI(data)

# run AMBARTI
fit.ambarti = run_AMBARTI(data, ntrees = 50, nburn = 1000, npost = 1000) # it takes a little while

# Get the final prediction (y hat)
yhat_ambarti = apply(fit.ambarti$y_hat, 2, mean)
yhat_ambarti2 = predict_ambarti(fit.ambarti, newdata = data_test$x , type = 'mean')
cor(data$y, yhat_ambarti);

# Get the prediction specifically from BART
yhat_bart = apply(fit.ambarti$y_hat_bart, 2, mean)
cor(data$y, yhat_bart); # correlation btw y and the BART component (from AMBARTI)

# Plot the main effects estimates and add the true values
g_hat = apply(fit.ambarti$g_hat, 2, mean)
plot(1:length(data$g), data$g, col=2, cex=2, main='AMBARTI-Genotype', ylim=c(-5,5)) # true values
points(g_hat, cex=2, pch = 2) # estimates
legend(7,4,'AMBARTI', col=1, pch = 2, bty='n')
legend(7,5,'True', col=2, pch = 1, cex=1, bty='n')

# Plot the main effects estimates and add the true values
e_hat = apply(fit.ambarti$e_hat, 2, mean)
plot(1:length(data$e), data$e, col=2, cex=2, main='AMBARTI-Environment', ylim=c(-5,5)) # true values
points(e_hat, cex=2, pch = 2) # estimates
legend(7,4,'AMBARTI', col=1, pch = 2, bty='n')
legend(7,5,'True', col=2, pch = 1, cex=1, bty='n')

fit.ambarti$trees[[100]][[1]] # show the first tree in the 100th MCMC iteration.

# ---------------------------------------
# BART (just to have a benchmark)
# ---------------------------------------
library(dbarts)
bart = dbarts::bart2(data$x, data$y)
cor(data$y, bart$yhat.train.mean)
```
