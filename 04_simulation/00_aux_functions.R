organise_classical_AMMI <- function(object, train_data){

  # Get training info
  x_train = train_data$x
  y_train = train_data$y

  # Get test info
  x_test = test_data$x
  y_test = test_data$y

  gen = as.factor(x_train[,"gen"])
  env = as.factor(x_train[,"env"])

  # Fit the linear model
  linear_mod = aov(y_train ~ gen + env + gen:env)

  # Get the residuals for the interaction gen:env
  interaction_tab = model.tables(linear_mod, type='effects', cterms = 'gen:env')
  interaction_tab = interaction_tab$tables$`gen:env`

  # Get the number of PCs
  PC = length(train_data$lambda)

  # Run the Singular Value Decomposition (SVD) to compute lambda, gamma, and delta
  sv_dec <- svd(interaction_tab, nu = PC, nv = PC)

  # Get parameter estimates
  # mu_hat     = linear_mod$coefficients[1] # slightly biased compared to mean(y_train)
  mu_hat = mean(y_train)
  alpha_hat  = aggregate(x = y_train - mu_hat, by = list(gen), FUN = "mean")[,2]
  beta_hat   = aggregate(x = y_train - mu_hat, by = list(env), FUN = "mean")[,2]
  lambda_hat = sv_dec$d[1:PC]
  gamma_hat  = -1*sv_dec$u
  delta_hat  = -1*sv_dec$v

  # Set up the bilinear term
  blin_train = rep(0, length(y_train))
  blin_test  = rep(0, length(y_test))

  for (k in 1:length(lambda_hat)) {
    blin_train = blin_train + lambda_hat[k]*gamma_hat[x_train[,'gen'],k]*delta_hat[x_train[,'env'],k]
    blin_test  = blin_test + lambda_hat[k]*gamma_hat[x_test[,'gen'],k]*delta_hat[x_test[,'env'],k]
  }

  # Compute the predicted response for the TRAINING data
  y_hat_train = mu_hat + alpha_hat[x_train[,'gen']] + beta_hat[x_train[,'env']] + blin_train
  plot(y_train, y_hat_train);abline(0,1)
  # Compute the predicted response for the TEST data
  y_hat_test = mu_hat + alpha_hat[x_test[,'gen']] + beta_hat[x_test[,'env']] + blin_test
  plot(y_test, y_hat_test);abline(0,1)

  aa = list(alpha_hat   = alpha_hat,
            beta_hat    = beta_hat,
            delta_hat   = delta_hat,
            gamma_hat   = gamma_hat,
            lambda_hat  = lambda_hat,
            y_hat_train = y_hat_train,
            y_hat_test  = y_hat_test)

  return(list(alpha_hat   = alpha_hat,
              beta_hat    = beta_hat,
              delta_hat   = delta_hat,
              gamma_hat   = gamma_hat,
              lambda_hat  = lambda_hat,
              y_hat_train = y_hat_train,
              y_hat_test  = y_hat_test))
}

plot(y_train, bb$y_hat_train);abline(0,1)
cor(y_train, bb$y_hat_train)

plot(y_test, bb$y_hat_test);abline(0,1)
cor(y_test, bb$y_hat_test)

organise_bayesian_AMMI <- function(object, data){

  # Get training info
  x_train = data$x
  y_train = data$y_train

  # Get test info
  x_test = data$x
  y_test = data$y_test

  # Get estimates info
  estimate   = bayesian_AMMI$BUGSoutput$mean
  mu_hat     = as.numeric(estimate$mu_all)
  alpha_hat  = estimate$alpha
  beta_hat   = estimate$beta
  delta_hat  = estimate$delta
  gamma_hat  = estimate$gamma
  lambda_hat = estimate$lambda

  # Set up the bilinear term
  blin_train = rep(0, length(y_train))
  blin_test  = rep(0, length(y_test))

  for (k in 1:length(lambda_hat)) {
    blin_train = blin_train + lambda_hat[k]*gamma_hat[x_train[,'gen'],k]*delta_hat[x_train[,'env'],k]
    blin_test  = blin_test + lambda_hat[k]*gamma_hat[x_test[,'gen'],k]*delta_hat[x_test[,'env'],k]
  }

  # Compute the predicted response for the TRAINING data
  y_hat_train = mu_hat + alpha_hat[x_train[,'gen']] + beta_hat[x_train[,'env']] + blin_train

  # Compute the predicted response for the TEST data
  y_hat_test = mu_hat + alpha_hat[x_test[,'gen']] + beta_hat[x_test[,'env']] + blin_test

  bb = list(alpha_hat  = alpha_hat,
            beta_hat    = beta_hat,
            delta_hat   = delta_hat,
            gamma_hat   = gamma_hat,
            lambda_hat  = lambda_hat,
            y_hat_train = y_hat_train,
            y_hat_test  = y_hat_test)

  return(list(alpha_hat  = alpha_hat,
              beta_hat    = beta_hat,
              delta_hat   = delta_hat,
              gamma_hat   = gamma_hat,
              lambda_hat  = lambda_hat,
              y_hat_train = y_hat_train,
              y_hat_test  = y_hat_test))
}

plot(train_data$alpha, ylim=c(-4,4))
points(bb$alpha_hat, col=2)
points(aa$alpha_hat, col=3)
