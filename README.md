# AMBARTI models for agricultural experiments

This repository contains R scripts and data sets that can be used to fit the Additive Main Effect Bayesian Additive Regression Tree interaction (AMBARTI) models for agricultural experiments.

In addition, it provides an implementation of AMBARTI in the format of an R package named ```AMBARTI```.

## Installation
``` r
library(devtools)
install_github("ebprado/AMBARTI")
```
## Example
``` r
library(AMBARTI)

# Simulate a Friedman data set
friedman_data = function(n, num_cov, sd_error){
  x = matrix(runif(n*num_cov),n,num_cov)
  y = 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5] + rnorm(n, sd=sd_error)
  return(list(y = y,
              x = x))
}
# Training data
data = friedman_data(200, 10, 1)
y = data$y
x = data$x

# Test data
data_test = friedman_data(100, 10, 1)
y.test = data_test$y
x.test = data_test$x

# Run MOTR-BART
set.seed(99)
fit.motr.bart = motr_bart(x, y, ntrees = 10, nburn = 100, npost = 100)
y.test.hat = predict_motr_bart(fit.motr.bart, x.test, 'mean')
plot(y.test, y.test.hat); abline(0, 1)
cor(y.test, y.test.hat)

# Run MOTR-BART for classification
set.seed(01)
y = ifelse(y > median(y), 1, 0)
y.test = ifelse(y.test > median(y.test), 1, 0)
fit.motr.bart = motr_bart_class(x, y, ntrees = 10, nburn = 100, npost = 100)
y.test.hat = predict_motr_bart_class(fit.motr.bart, x.test, 'mean')
```
