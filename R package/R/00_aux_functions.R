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
  Q           = object$Q

  rmse_y_train = RMSE(y_train, y_hat_train)
  rmse_y_test  = RMSE(y_test, y_hat_test)
  rmse_blinear = RMSE(blinear, blinear_hat)

  g  = data$g
  e  = data$e
  if (is.null(data$lambda) == FALSE) {lambda = data$lambda}
  if (is.null(data$gamma) == FALSE) {gamma   = data$gamma}
  if (is.null(data$delta) == FALSE) {delta   = data$delta}

  g_hat  = object$g_hat
  e_hat  = object$e_hat

  rrmse_g = RRMSE(g, g_hat)
  rrmse_e = RRMSE(e, e_hat)

  if (id != 'AMBARTI' && is.null(data$lambda) == FALSE) {

    lambda_hat = object$lambda_hat
    gamma_hat  = object$gamma_hat
    delta_hat  = object$delta_hat

    rrmse_lambda = RRMSE(lambda, lambda_hat)
    rrmse_gamma  = RRMSE(gamma, gamma_hat)
    rrmse_delta  = RRMSE(delta, delta_hat)
    lambda = gsub(', ', ' ', toString(data$lambda))

  } else{

    rrmse_lambda = NA
    rrmse_gamma = NA
    rrmse_delta = NA
    lambda = NA

  }

  I  = data$I
  J  = data$J
  sa = data$s_g
  sb = data$s_e
  sy = data$s_y

aux = data.frame(
                 id           = id,
                 rep          = rep,
                 I            = I,
                 J            = J,
                 sa           = sa,
                 sb           = sb,
                 sy           = sy,
                 Q            = Q,
                 lambda       = lambda,
                 y_train_rmse = rmse_y_train,
                 y_test_rmse  = rmse_y_test,
                 rmse_blinear = rmse_blinear,
                 rrmse_g      = rrmse_g,
                 rrmse_e      = rrmse_e,
                 lambda_rrmse = rrmse_lambda,
                 gamma_rrmse  = rrmse_gamma,
                 delta_rrmse  = rrmse_delta
                 )

return(aux)

}

#' @export
#' @importFrom stats 'aov' 'model.tables' 'lm'
organise_classical_AMMI <- function(object, data, Q = NULL){

  # Get training info
  x_train = data$x
  y_train = data$y

  # Get test info
  x_test = data$x
  y_test = data$y_test

  mu_hat     = object$mu_hat
  # g_hat      = object$g_hat
  # e_hat      = object$e_hat
  g_hat      = object$alpha_hat
  e_hat      = object$beta_hat
  lambda_hat = object$lambda_hat
  gamma_hat  = object$gamma_hat
  delta_hat  = object$delta_hat

  if (is.null(data$Q) == FALSE) {Q = data$Q}

  # Set up the bilinear term
  blin_train = rep(0, length(y_train))
  blin_test  = rep(0, length(y_test))

  for (k in 1:Q) {
    blin_train = blin_train + lambda_hat[k]*gamma_hat[x_train[,'g'],k]*delta_hat[x_train[,'e'],k]
    blin_test  = blin_test + lambda_hat[k]*gamma_hat[x_test[,'g'],k]*delta_hat[x_test[,'e'],k]
  }

  # Compute the predicted response for the TRAINING data
  y_hat_train = mu_hat + g_hat[x_train[,'g']] + e_hat[x_train[,'e']] + blin_train

  # Compute the predicted response for the TEST data
  y_hat_test = mu_hat + g_hat[x_test[,'g']] + e_hat[x_test[,'e']] + blin_test

  return(list(g_hat       = g_hat,
              e_hat       = e_hat,
              delta_hat   = delta_hat,
              gamma_hat   = gamma_hat,
              lambda_hat  = lambda_hat,
              y_hat_train = y_hat_train,
              y_hat_test  = y_hat_test,
              blinear_hat = blin_train,
              Q           = Q,
              id          = 'classical AMMI'))
}

#' @export
#'
organise_bayesian_AMMI_WITH_postprocessing <- function(object, data, Q = NULL){

  # Get training info
  x_train = data$x
  y_train = data$y

  # Get test info
  x_test = data$x
  y_test = data$y_test

  # Get the number of PCs
  if (is.null(data$Q) == FALSE) {Q = data$Q}

  # Get some MCMC info
  nburn      = object$BUGSoutput$n.burnin
  niter      = object$BUGSoutput$n.iter
  npost      = niter - nburn
  seq_burn   = seq(1, nburn, by=1)

  # Get estimates info
  estimate   = object$BUGSoutput$sims.matrix
  mu_hat     = estimate[-seq_burn,colnames(estimate)[grepl('mu_all', colnames(estimate))]]
  g_hat      = estimate[-seq_burn,colnames(estimate)[grepl('^g*?(\\d+).*',  colnames(estimate))]]
  e_hat      = estimate[-seq_burn,colnames(estimate)[grepl('^e*?(\\d+).',   colnames(estimate))]]
  delta_hat  = estimate[-seq_burn,colnames(estimate)[grepl('delta',  colnames(estimate))]]
  gamma_hat  = estimate[-seq_burn,colnames(estimate)[grepl('gamma',  colnames(estimate))]]
  lambda_hat = estimate[-seq_burn,colnames(estimate)[grepl('lambda', colnames(estimate))]]
  lambda_hat = as.matrix(lambda_hat)

  # Compute the bilinear term
  blin_train = matrix(0, nrow = npost, ncol = nrow(x_train))
  blin_test  = matrix(0, nrow = npost, ncol = nrow(x_test))

  for (k in 1:Q) {
    blin_train = blin_train + lambda_hat[,k] * gamma_hat[,x_train[,'g']]* delta_hat[,x_train[,'e']]
    blin_test  = blin_test +  lambda_hat[,k] * gamma_hat[,x_test[,'g']] * delta_hat[, x_test[,'e']]
  }

  # Compute the predicted response for the TRAINING data
  mu_ij = mu_hat + g_hat[,x_train[,'g']] + e_hat[,x_train[,'e']] + blin_train
  colnames(mu_ij) = NULL

  # -------------------------------------------------------------
  # Start postprocessing
  # -------------------------------------------------------------

  n_gen = length(unique(x_train[,'g']))
  n_env = length(unique(x_train[,'e']))

  # Create matrices/lists to store the postprocessing results
  snew_mu_hat     = matrix(NA, nrow=npost, ncol=1)
  snew_g_hat      = matrix(NA, nrow=npost, ncol=n_gen)
  snew_e_hat      = matrix(NA, nrow=npost, ncol=n_env)
  snew_lambda_hat = matrix(NA, nrow=npost, ncol=Q)
  snew_gamma_hat  = list()
  snew_delta_hat  = list()

  for (i in 1:nrow(mu_ij)){

    # Get the matrix with mu_ij for each genotype and environment
    matrix_mu_ij  = matrix(mu_ij[i,], nrow=n_gen, ncol=n_env)
    new_mu_hat    = mean(matrix_mu_ij)
    new_g_hat     = new_mu_hat - rowMeans(matrix_mu_ij) # their sum is 0
    new_e_hat     = new_mu_hat - colMeans(matrix_mu_ij) # their sum is 0

    # Center matrix by row and column
    # Thank https://stackoverflow.com/questions/43639063/double-centering-in-r
    resA = t(matrix_mu_ij*0 + colMeans(matrix_mu_ij))
    resB = matrix_mu_ij*0 + rowMeans(matrix_mu_ij)
    res_double_centered = matrix_mu_ij - resA - resB + new_mu_hat # sum zero by row and column

    # Run the Singular Value Decomposition (SVD) to compute lambda, gamma, and delta
    sv_dec <- svd(res_double_centered, nu = Q, nv = Q)

    # Get new parameter estimates
    new_lambda_hat = sv_dec$d[1:Q]
    new_gamma_hat  = -1*sv_dec$u
    new_delta_hat  = -1*sv_dec$v

    # Store the new results
    snew_mu_hat[i,]       = new_mu_hat
    snew_g_hat[i,]        = new_g_hat
    snew_e_hat[i,]        = new_e_hat
    snew_lambda_hat[i,]   = new_lambda_hat
    snew_gamma_hat[[i]]   = new_gamma_hat
    snew_delta_hat[[i]]   = new_delta_hat
  }
  # Clean column names
  colnames(snew_mu_hat) = NULL
  colnames(snew_g_hat) = NULL
  colnames(snew_e_hat) = NULL
  colnames(snew_lambda_hat) = NULL

  # Summarise the new posterior results
  new_mu_hat     = apply(snew_mu_hat,2,mean)
  new_g_hat      = apply(snew_g_hat,2,mean)
  new_e_hat      = apply(snew_e_hat,2,mean)
  new_lambda_hat = apply(snew_lambda_hat ,2,mean)
  new_gamma_hat  = apply(simplify2array(snew_gamma_hat), 1:2, mean)
  new_delta_hat  = apply(simplify2array(snew_delta_hat), 1:2, mean)

  # Compute the bilinear term
  new_blin_train = rep(0, length(y_train))
  new_blin_test  = rep(0, length(y_test))

  for (k in 1:Q) {
    new_blin_train = new_blin_train + new_lambda_hat[k]*new_gamma_hat[x_train[,'g'],k]*new_delta_hat[x_train[,'e'],k]
    new_blin_test  = new_blin_test + new_lambda_hat[k]*new_gamma_hat[x_test[,'g'],k]*new_delta_hat[x_test[,'e'],k]
  }

  # Compute the predicted values for the TRAINING and TEST data sets
  new_mu_ij_train = new_mu_hat + new_g_hat[x_train[,'g']] + new_e_hat[x_train[,'e']] + new_blin_train
  new_mu_ij_test = new_mu_hat + new_g_hat[x_test[,'g']] + new_e_hat[x_test[,'e']] + new_blin_test

  return(list(g_hat       = new_g_hat,
              e_hat       = new_e_hat,
              delta_hat   = new_delta_hat,
              gamma_hat   = new_gamma_hat,
              lambda_hat  = new_lambda_hat,
              y_hat_train = new_mu_ij_train,
              y_hat_test  = new_mu_ij_test,
              blinear_hat = new_blin_train,
              Q           = Q,
              id          = 'Bayesian AMMI (postproc)'))
}

# plot(new_g_hat, ylim=c(-8,8), main='g - WITH postprocessing')
# points(data$g, col=2)
# points(g_hat[1,], col=3)
#
# plot(new_e_hat, ylim=c(-8,8), main='e - WITH postprocessing')
# points(data$e, col=2)
# points(e_hat[1,], col=3)
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
  x_test   = data$x
  x_test$g = as.factor(x_test$g)
  x_test$e = as.factor(x_test$e)
  y_test   = data$y_test

  # Get estimates info
  # estimate_g  = apply(object$g_hat, 2, mean)
  # estimate_e  = apply(object$e_hat, 2, mean)
  estimate  = apply(object$beta_hat, 2, mean)
  # g_hat       = estimate_g[grepl('^g', names(estimate_g))]
  # e_hat       = estimate_e[grepl('^e', names(estimate_e))]
  g_hat       = estimate[grepl('^g', names(estimate))]
  e_hat       = estimate[grepl('^e', names(estimate))]
  y_hat_train = apply(object$y_hat, 2, mean)
  y_hat_test  = as.numeric(predict_ambarti(object, x_test, type = 'mean'))
  blinear_hat = apply(object$y_hat_bart,2,mean)

  return(list(g_hat       = g_hat,
              e_hat       = e_hat,
              y_hat_train = y_hat_train,
              y_hat_test  = y_hat_test,
              blinear_hat = blinear_hat,
              Q           = NA,
              id          = 'AMBARTI'))

}

#' @export
#'
organise_bayesian_AMMI_WITHOUT_postprocessing <- function(object, data, Q = NULL){

  # Get training info
  x_train = data$x
  y_train = data$y

  # Get test info
  x_test = data$x
  y_test = data$y_test

  # Get the number of PCs
  if (is.null(data$Q) == FALSE) {Q = data$Q}

  # Get estimates info
  estimate   = object$BUGSoutput$mean
  mu_hat     = as.numeric(estimate$mu_all)
  g_hat      = estimate$g
  e_hat      = estimate$e
  delta_hat  = estimate$delta
  gamma_hat  = estimate$gamma
  lambda_hat = estimate$lambda

  # Set up the bilinear term
  blin_train = rep(0, length(y_train))
  blin_test  = rep(0, length(y_test))

  for (k in 1:Q) {
    blin_train = blin_train + lambda_hat[k]*gamma_hat[x_train[,'g'],k]*delta_hat[x_train[,'e'],k]
    blin_test  = blin_test + lambda_hat[k]*gamma_hat[x_test[,'g'],k]*delta_hat[x_test[,'e'],k]
  }

  # Compute the predicted response for the TRAINING data
  y_hat_train = mu_hat + g_hat[x_train[,'g']] + e_hat[x_train[,'e']] + blin_train
  # Compute the predicted response for the TEST data
  y_hat_test = mu_hat + g_hat[x_test[,'g']] + e_hat[x_test[,'e']] + blin_test

  return(list(g_hat       = g_hat,
              e_hat       = e_hat,
              delta_hat   = delta_hat,
              gamma_hat   = gamma_hat,
              lambda_hat  = lambda_hat,
              y_hat_train = y_hat_train,
              y_hat_test  = y_hat_test,
              blinear_hat = blin_train,
              Q           = Q,
              id          = 'Bayesian AMMI (NO postproc)'))
}




###### Plot functions ##########
###### Plot functions ##########
###### Plot functions ##########
###### Plot functions ##########

#' @export
#'
AMMI_help_plot <- function(object, data, Q = NULL){

  # Get training info
  x_train      = data$x
  y_train      = data$y

  # Get test info
  x_test = data$x

  if(is.null(data$y_test) == FALSE) {y_test      = data$y_test} else {y_test=NA}
  if(is.null(data$g) == FALSE)      {g_true      = data$g} else {g_true=NA}
  if(is.null(data$e) == FALSE)      {e_true      = data$e} else {e_true=NA}
  if(is.null(data$lambda) == FALSE) {lambda_true = data$lambda} else {lambda_true=NA}
  if(is.null(data$gamma) == FALSE)  {gamma_true  = data$gamma}  else {gamma_true=NA}
  if(is.null(data$delta) == FALSE)  {delta_true  = data$delta}  else {delta_true=NA}
  if(is.null(data$blinear) == FALSE)  {blinear_true  = data$blinear}  else {blinear_true=NA}

  # Get parameter estimates
  mu_hat     = object$mu_hat
  # g_hat      = object$g_hat
  # e_hat      = object$e_hat
  g_hat      = object$alpha_hat
  e_hat      = object$beta_hat
  lambda_hat = object$lambda_hat
  gamma_hat  = object$gamma_hat
  delta_hat  = object$delta_hat

  # Get the number of PCs
  if(is.null(data$Q) == FALSE){Q = data$Q}

  # Set up the bilinear term
  blin_train = rep(0, length(y_train))
  blin_test  = rep(0, length(y_test))

  for (k in 1:Q) {
    blin_train = blin_train + lambda_hat[k]*gamma_hat[x_train[,'g'],k]*delta_hat[x_train[,'e'],k]
    blin_test  = blin_test + lambda_hat[k]*gamma_hat[x_test[,'g'],k]*delta_hat[x_test[,'e'],k]
  }

  # Compute the predicted response for the TRAINING data
  y_hat_train = mu_hat + g_hat[x_train[,'g']] + e_hat[x_train[,'e']] + blin_train

  # Compute the predicted response for the TEST data
  y_hat_test = mu_hat + g_hat[x_test[,'g']] + e_hat[x_test[,'e']] + blin_test

  #
  g_names = unique(x_train$g)
  e_names = unique(x_train$e)

  if (is.null(data$y_test)==TRUE) {
    g_names = gsub('g', '', unique(x_train$g))
    e_names = gsub('e', '', unique(x_train$e))
  }

  id = 'AMMI'
  g_hat      = data.frame(id = id, Q = Q, variable = paste('g[',g_names ,']', sep=''), value=g_hat, true=g_true)
  e_hat      = data.frame(id = id, Q = Q, variable = paste('e[',e_names ,']', sep=''), value=e_hat, true=e_true)
  lambda_hat = data.frame(id = id, Q = Q, variable = paste('lambda[',1:Q,']',sep=''), value=lambda_hat, true=lambda_true)
  gamma_hat  = data.frame(id = id, Q = Q, variable = paste('gamma[',g_names,',',rep(1:Q, each=length(g_names)),']',sep=''), value=as.numeric(gamma_hat), true=as.numeric(gamma_true))
  delta_hat  = data.frame(id = id, Q = Q, variable = paste('delta[',e_names,',',rep(1:Q, each=length(e_names)),']',sep=''), value=as.numeric(delta_hat), true=as.numeric(delta_true))

  aux_blinear_hat = data.frame(id = id, Q = Q, variable = 'blinear', value = blinear_true - blin_train, true = blinear_true)
  aux_y_hat_train = data.frame(id = id, Q = Q, variable = 'yhattrain', value = y_train - y_hat_train, true = y_train)
  aux_y_hat_test  = data.frame(id = id, Q = Q, variable = 'yhattest', value = y_test - y_hat_test, true = y_test)

  return(list(g_hat       = g_hat,
              e_hat       = e_hat,
              delta_hat   = delta_hat,
              gamma_hat   = gamma_hat,
              lambda_hat  = lambda_hat,
              y_hat_train = aux_y_hat_train,
              y_hat_test  = aux_y_hat_test,
              blinear_hat = aux_blinear_hat,
              Q           = Q,
              id          = id))
}

#' @export
#'
bAMMI_help_plot <- function(object, data, Q = NULL){

# Get training info
x_train = data$x
y_train = data$y

g_true = data$g
e_true = data$e
if(is.null(data$lambda) == FALSE) {lambda_true = data$lambda} else {lambda_true=NA}
if(is.null(data$gamma) == FALSE)  {gamma_true  = data$gamma}  else {gamma_true=NA}
if(is.null(data$delta) == FALSE)  {delta_true  = data$delta}  else {delta_true=NA}
blinear_true = data$blinear

# Get test info
x_test = data$x
y_test = data$y_test

# Get the number of PCs
if(is.null(data$Q) == FALSE){Q = data$Q}

# Get estimates info
estimate   = object$BUGSoutput$mean
mu_hat     = as.numeric(estimate$mu_all)
g_hat      = estimate$g
e_hat      = estimate$e
delta_hat  = estimate$delta
gamma_hat  = estimate$gamma
lambda_hat = estimate$lambda

# Set up the bilinear term
blin_train = rep(0, length(y_train))
blin_test  = rep(0, length(y_test))

for (k in 1:Q) {
  blin_train = blin_train + lambda_hat[k]*gamma_hat[x_train[,'g'],k]*delta_hat[x_train[,'e'],k]
  blin_test  = blin_test + lambda_hat[k]*gamma_hat[x_test[,'g'],k]*delta_hat[x_test[,'e'],k]
}

# Compute the predicted response for the TRAINING data
y_hat_train = mu_hat + g_hat[x_train[,'g']] + e_hat[x_train[,'e']] + blin_train
# Compute the predicted response for the TEST data
y_hat_test = mu_hat + g_hat[x_test[,'g']] + e_hat[x_test[,'e']] + blin_test

# Get some MCMC info
nburn      = object$BUGSoutput$n.burnin
niter      = object$BUGSoutput$n.iter
npost      = niter - nburn
seq_burn   = seq(1, nburn, by=1)

# Get estimates info
id = 'Bayesian AMMI (no postproc)'
estimate    = as.data.frame(object$BUGSoutput$sims.matrix)
estimate$id = id
estimate$Q  = Q
#mu_hat     = estimate[-seq_burn,c('id',colnames(estimate)[grepl('mu_all', colnames(estimate))])]
g_hat      = estimate[-seq_burn,c('id','Q',colnames(estimate)[grepl('^g*?(\\d+).*', colnames(estimate))])]
e_hat      = estimate[-seq_burn,c('id','Q',colnames(estimate)[grepl('^e*?(\\d+).*', colnames(estimate))])]
gamma_hat  = estimate[-seq_burn,c('id','Q',colnames(estimate)[grepl('gamma', colnames(estimate))])]
delta_hat  = estimate[-seq_burn,c('id','Q',colnames(estimate)[grepl('delta', colnames(estimate))])]
lambda_hat = estimate[-seq_burn,c('id','Q',colnames(estimate)[grepl('lambda', colnames(estimate))])]
if(Q==1) {colnames(lambda_hat)[3] = 'lambda[1]'}

# Remove columns id and Q
m_vars_g      = colnames(g_hat)[-which(colnames(g_hat) %in% c('id', 'Q'))]
m_vars_e      = colnames(e_hat)[-which(colnames(e_hat) %in% c('id', 'Q'))]
m_vars_lambda = colnames(lambda_hat)[-which(colnames(lambda_hat) %in% c('id', 'Q'))]
m_vars_gamma  = colnames(gamma_hat)[-which(colnames(gamma_hat) %in% c('id', 'Q'))]
m_vars_delta  = colnames(delta_hat)[-which(colnames(delta_hat) %in% c('id', 'Q'))]

g_hat      = melt(g_hat,  measure.vars = m_vars_g)
e_hat      = melt(e_hat,   measure.vars = m_vars_e)
lambda_hat = melt(lambda_hat, measure.vars = m_vars_lambda)
gamma_hat  = melt(gamma_hat,  measure.vars = m_vars_gamma)
delta_hat  = melt(delta_hat,  measure.vars = m_vars_delta)

g_hat$true      = rep(g_true,  each=npost)
e_hat$true      = rep(e_true,   each=npost)
lambda_hat$true = rep(lambda_true, each=npost)
gamma_hat$true  = rep(gamma_true,  each=npost)
delta_hat$true  = rep(delta_true,  each=npost)

aux_blinear_hat = data.frame(id = id, Q = Q, variable = 'blinear', value = data$blinear - blin_train, true = blinear_true)
aux_y_hat_train = data.frame(id = id, Q = Q, variable = 'yhattrain', value = y_train - y_hat_train, true = y_train)
aux_y_hat_test  = data.frame(id = id, Q = Q, variable = 'yhattest', value = y_test - y_hat_test, true = y_test)

return(list(g_hat       = g_hat,
            e_hat       = e_hat,
            delta_hat   = delta_hat,
            gamma_hat   = gamma_hat,
            lambda_hat  = lambda_hat,
            blinear_hat = aux_blinear_hat,
            y_hat_train = aux_y_hat_train,
            y_hat_test  = aux_y_hat_test,
            Q           = Q,
            id          = id))
}

#' @export
#' @importFrom reshape2 'melt'
bAMMI_help_plot_WITHPOS <- function(object, data, Q = NULL){

  # Get training info
  x_train = data$x
  y_train = data$y
  g_true  = data$g
  e_true  = data$e
  if(is.null(data$lambda) == FALSE) {lambda_true = data$lambda} else {lambda_true=NA}
  if(is.null(data$gamma) == FALSE)  {gamma_true  = data$gamma}  else {gamma_true=NA}
  if(is.null(data$delta) == FALSE)  {delta_true  = data$delta}  else {delta_true=NA}
  blinear_true = data$blinear

  # Get test info
  x_test = data$x
  y_test = data$y_test

  # Get the number of PCs
  if(is.null(data$Q) == FALSE){Q = data$Q}

  # Get some MCMC info
  nburn      = object$BUGSoutput$n.burnin
  niter      = object$BUGSoutput$n.iter
  npost      = niter - nburn
  seq_burn   = seq(1, nburn, by=1)

  # Get estimates info
  estimate   = object$BUGSoutput$sims.matrix
  mu_hat     = estimate[-seq_burn,colnames(estimate)[grepl('mu_all', colnames(estimate))]]
  g_hat      = estimate[-seq_burn,colnames(estimate)[grepl('^g*?(\\d+).*',  colnames(estimate))]]
  e_hat      = estimate[-seq_burn,colnames(estimate)[grepl('^e*?(\\d+).',   colnames(estimate))]]
  delta_hat  = estimate[-seq_burn,colnames(estimate)[grepl('delta',  colnames(estimate))]]
  gamma_hat  = estimate[-seq_burn,colnames(estimate)[grepl('gamma',  colnames(estimate))]]
  lambda_hat = estimate[-seq_burn,colnames(estimate)[grepl('lambda', colnames(estimate))]]
  lambda_hat = as.matrix(lambda_hat)

  # Compute the bilinear term
  blin_train = matrix(0, nrow = npost, ncol = nrow(x_train))
  blin_test  = matrix(0, nrow = npost, ncol = nrow(x_test))

  for (k in 1:Q) {
    blin_train = blin_train + lambda_hat[,k] * gamma_hat[,x_train[,'g']]* delta_hat[,x_train[,'e']]
    blin_test  = blin_test +  lambda_hat[,k] * gamma_hat[,x_test[,'g']] * delta_hat[, x_test[,'e']]
  }

  # Compute the predicted response for the TRAINING data
  mu_ij = mu_hat + g_hat[,x_train[,'g']] + e_hat[,x_train[,'e']] + blin_train
  colnames(mu_ij) = NULL

  # -------------------------------------------------------------
  # Start postprocessing
  # -------------------------------------------------------------

  n_gen = length(unique(x_train[,'g']))
  n_env = length(unique(x_train[,'e']))

  # Create matrices/lists to store the postprocessing results
  snew_mu_hat     = matrix(NA, nrow=npost, ncol=1)
  snew_g_hat      = matrix(NA, nrow=npost, ncol=n_gen)
  snew_e_hat      = matrix(NA, nrow=npost, ncol=n_env)
  snew_lambda_hat = matrix(NA, nrow=npost, ncol=Q)
  snew_gamma_hat  = matrix(NA, nrow=npost, ncol=Q*n_gen)
  snew_delta_hat  = matrix(NA, nrow=npost, ncol=Q*n_env)

  for (i in 1:nrow(mu_ij)){

    # Get the matrix with mu_ij for each genotype and environment
    matrix_mu_ij = matrix(mu_ij[i,], nrow=n_gen, ncol=n_env)
    new_mu_hat   = mean(matrix_mu_ij)
    new_g_hat    = new_mu_hat - rowMeans(matrix_mu_ij) # their sum is 0
    new_e_hat    = new_mu_hat - colMeans(matrix_mu_ij) # their sum is 0

    # Center matrix by row and column
    # Thank https://stackoverflow.com/questions/43639063/double-centering-in-r
    resA = t(matrix_mu_ij*0 + colMeans(matrix_mu_ij))
    resB = matrix_mu_ij*0 + rowMeans(matrix_mu_ij)
    res_double_centered = matrix_mu_ij - resA - resB + new_mu_hat # sum zero by row and column

    # Run the Singular Value Decomposition (SVD) to compute lambda, gamma, and delta
    sv_dec <- svd(res_double_centered, nu = Q, nv = Q)

    # Get new parameter estimates
    new_lambda_hat = sv_dec$d[1:Q]
    new_gamma_hat  = -1*sv_dec$u
    new_delta_hat  = -1*sv_dec$v

    # Store the new results
    snew_mu_hat[i,]     = new_mu_hat
    snew_g_hat[i,]      = new_g_hat
    snew_e_hat[i,]      = new_e_hat
    snew_lambda_hat[i,] = new_lambda_hat
    snew_gamma_hat[i,]  = as.numeric(new_gamma_hat)
    snew_delta_hat[i,]  = as.numeric(new_delta_hat)
  }
  # Clean column names
  colnames(snew_mu_hat) = colnames(mu_hat)
  colnames(snew_g_hat)  = colnames(g_hat)
  colnames(snew_e_hat)  = colnames(e_hat)
  if(Q==1) {colnames(snew_lambda_hat) = 'lambda[1]'} else{colnames(snew_lambda_hat) = colnames(lambda_hat)}
  colnames(snew_gamma_hat)  = colnames(gamma_hat)
  colnames(snew_delta_hat)  = colnames(delta_hat)

  # Summarise the new posterior results
  new_mu_hat     = apply(snew_mu_hat     ,2,mean)
  new_g_hat      = apply(snew_g_hat  ,2,mean)
  new_e_hat      = apply(snew_e_hat   ,2,mean)
  new_lambda_hat = apply(snew_lambda_hat ,2,mean)
  new_gamma_hat  = apply(snew_gamma_hat,  2,mean)
  new_delta_hat  = apply(snew_delta_hat,  2,mean)

  # Compute the bilinear term
  new_blin_train = rep(0, length(y_train))
  new_blin_test  = rep(0, length(y_test))

  aux_new_gamma_hat = matrix(new_gamma_hat,ncol=Q)
  aux_new_delta_hat = matrix(new_delta_hat,ncol=Q)

  for (k in 1:Q) {
    new_blin_train = new_blin_train + new_lambda_hat[k]*aux_new_gamma_hat[x_train[,'g'],k]*aux_new_delta_hat[x_train[,'e'],k]
    new_blin_test  = new_blin_test + new_lambda_hat[k]*aux_new_gamma_hat[x_test[,'g'],k]*aux_new_delta_hat[x_test[,'e'],k]
  }

  # Compute the predicted values for the TRAINING and TEST data sets
  new_mu_ij_train = new_mu_hat + new_g_hat[x_train[,'g']] + new_e_hat[x_train[,'e']] + new_blin_train
  new_mu_ij_test = new_mu_hat + new_g_hat[x_test[,'g']] + new_e_hat[x_test[,'e']] + new_blin_test

  id = 'Bayesian AMMI (posproc)'
  snew_g_hat      = as.data.frame(snew_g_hat)
  snew_e_hat      = as.data.frame(snew_e_hat)
  snew_lambda_hat = as.data.frame(snew_lambda_hat)
  snew_gamma_hat  = as.data.frame(snew_gamma_hat)
  snew_delta_hat  = as.data.frame(snew_delta_hat)

  snew_g_hat$id      = id
  snew_e_hat$id      = id
  snew_lambda_hat$id = id
  snew_gamma_hat$id  = id
  snew_delta_hat$id  = id

  snew_g_hat$Q      = Q
  snew_e_hat$Q      = Q
  snew_lambda_hat$Q = Q
  snew_gamma_hat$Q  = Q
  snew_delta_hat$Q  = Q

  # Remove columns id and Q
  m_vars_g      = colnames(snew_g_hat)[-which(colnames(snew_g_hat) %in% c('id', 'Q'))]
  m_vars_e      = colnames(snew_e_hat)[-which(colnames(snew_e_hat) %in% c('id', 'Q'))]
  m_vars_lambda = colnames(snew_lambda_hat)[-which(colnames(snew_lambda_hat) %in% c('id', 'Q'))]
  m_vars_gamma  = colnames(snew_gamma_hat)[-which(colnames(snew_gamma_hat) %in% c('id', 'Q'))]
  m_vars_delta  = colnames(snew_delta_hat)[-which(colnames(snew_delta_hat) %in% c('id', 'Q'))]

  # Transpose the data
  snew_g_hat      = melt(snew_g_hat,  measure.vars  = m_vars_g)
  snew_e_hat      = melt(snew_e_hat,   measure.vars  = m_vars_e)
  snew_lambda_hat = melt(snew_lambda_hat, measure.vars  = m_vars_lambda)
  snew_gamma_hat  = melt(snew_gamma_hat,  measure.vars  = m_vars_gamma)
  snew_delta_hat  = melt(snew_delta_hat,  measure.vars  = m_vars_delta)

  # Add true values
  snew_g_hat$true       = rep(g_true,  each=npost)
  snew_e_hat$true       = rep(e_true,   each=npost)
  snew_lambda_hat$true  = rep(lambda_true, each=npost)
  snew_gamma_hat$true   = rep(gamma_true,  each=npost)
  snew_delta_hat$true   = rep(delta_true,  each=npost)

  aux_blinear_hat   = data.frame(id = id, Q = Q, variable = 'blinear', value = data$blinear - new_blin_train, true = blinear_true)
  aux_y_hat_train   = data.frame(id = id, Q = Q, variable = 'yhattrain', value = y_train - new_mu_ij_train, true = y_train)
  aux_y_hat_test    = data.frame(id = id, Q = Q, variable = 'yhattest', value = y_test - new_mu_ij_test, true = y_test)

  return(list(g_hat       = snew_g_hat,
              e_hat       = snew_e_hat,
              delta_hat   = snew_delta_hat,
              gamma_hat   = snew_gamma_hat,
              lambda_hat  = snew_lambda_hat,
              y_hat_train = aux_y_hat_train,
              y_hat_test  = aux_y_hat_test,
              blinear_hat = aux_blinear_hat,
              Q           = Q,
              id          = id))
}
#' @export

AMBARTI_help_plot <- function(object, data){

  # Get training info
  x_train = data$x
  y_train = data$y
  if(is.null(data$g) == FALSE)  {g_true  = data$g} else {g_true=NA}
  if(is.null(data$e) == FALSE)   {e_true = data$e} else {e_true=NA}
  if(is.null(data$blinear) == FALSE)  {blinear_true  = data$blinear}  else {blinear_true=NA}

  # Get test info
  x_test = data$x
  x_test$g = as.factor(x_test$g)
  x_test$e = as.factor(x_test$e)
  if(is.null(data$y_test) == FALSE)  {y_test  = data$y_test}  else {y_test=NA}

  # Get estimates info
  id = 'AMBARTI'
  # estimate_g  = as.data.frame(object$g_hat)
  # estimate_e  = as.data.frame(object$e_hat)
  # g_hat       = estimate_g[,grepl('^g', names(estimate_g))]
  # e_hat       = estimate_e[,grepl('^e', names(estimate_e))]
  estimate  = as.data.frame(object$beta_hat)
  g_hat       = estimate[,grepl('^g', names(estimate))]
  e_hat       = estimate[,grepl('^e', names(estimate))]
  y_hat_train = apply(object$y_hat, 2, mean)

  if(is.null(data$y_test) == FALSE) {
    y_test       = data$y_test
    y_hat_test   = as.numeric(predict_ambarti(object, x_test, type = 'mean'))
    names(g_hat) = paste(gsub('g','g[', names(g_hat)), ']', sep='')
    names(e_hat) = paste(gsub('e','e[', names(e_hat)), ']', sep='')
  } else {
    y_test=NA
    y_hat_test=NA
    names(g_hat) = paste(gsub('gg','g[', names(g_hat)), ']', sep='')
    names(e_hat) = paste(gsub('ee','e[', names(e_hat)), ']', sep='')
  }

  blinear_hat = apply(object$y_hat_bart,2,mean)
  g_hat$id = id
  g_hat$Q = NA
  e_hat$id = id
  e_hat$Q = NA

  g_hat = melt(g_hat, measure.vars = colnames(g_hat)[grepl('g', colnames(g_hat))])
  e_hat = melt(e_hat, measure.vars = colnames(e_hat)[grepl('e', colnames(e_hat))])

  g_hat$true = rep(as.numeric(g_true), each=object$npost)
  e_hat$true = rep(as.numeric(e_true), each=object$npost)

  aux_blinear_hat = data.frame(id = id, Q = NA, variable = 'blinear', value = blinear_true - blinear_hat, true = blinear_true)
  aux_y_hat_train = data.frame(id = id, Q = NA, variable = 'yhattrain', value = y_train - y_hat_train, true = y_train)
  aux_y_hat_test  = data.frame(id = id, Q = NA, variable = 'yhattest', value = y_test - y_hat_test, true = y_test)


  return(list(g_hat       = g_hat,
              e_hat       = e_hat,
              y_hat_train = aux_y_hat_train,
              y_hat_test  = aux_y_hat_test,
              blinear_hat = aux_blinear_hat,
              Q           = NA,
              id          = 'AMBARTI'))
}


#' @export
new_parse_format <- function(text) {
  out <- vector("expression", length(text))
  for (i in seq_along(text)) {
    expr <- parse(text = text[[i]])
    out[[i]] <- if (length(expr) == 0)
      NA
    else expr[[1]]
  }
  out
}