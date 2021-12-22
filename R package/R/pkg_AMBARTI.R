#' @export
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm' 'aggregate' 'contrasts' 'model.matrix' 'as.formula'
#' @importFrom MCMCpack 'rdirichlet'
#' @importFrom truncnorm 'rtruncnorm'

# x
# y
# ntrees = 10
# node_min_size = 5
# alpha = 0.95
# beta = 2
# nu = 3
# lambda = 0.1
# mu_mu = 0
# mu_g = 0
# mu_e = 0
# sigma2 = 1
# sigma2_mu = 1
# sigma2_psi = 1
# a_g = 1
# b_g = 1
# a_e = 1
# b_e = 1
# nburn = 1000
# npost = 1000
# nthin = 1
# data = generate_data_AMMI(20, 20, 1,1,1,c(8, 10, 15))
# x = data$x
# names(x) <- c('g','e')
# x$g <- as.factor(x$g)
# x$e <- as.factor(x$e)
# y = data$y

# ambarti = ambarti(x,y,ntrees=200, nburn=100, npost=100)
# newdata = data$x
# object = ambarti
# bb = predict_ambarti(object, newdata, type = 'mean')
# aa = apply(ambarti$y_hat,2,mean)
# plot(aa,bb)
# plot(y, aa); abline(0,1)
# cor(y, aa)
# cc = run_classical_AMMI(data)
# dd = organise_classical_AMMI(cc, data)
# dd
# points(y, dd$y_hat_train, col=2)
# cor(y, dd$y_hat_train)
ambarti = function(x,
                   y,
                   ntrees = 10,
                   node_min_size = 5,
                   alpha = 0.95,
                   beta = 2,
                   nu = 3,
                   lambda = 0.1,
                   mu_mu = 0,
                   mu_g = 0,
                   mu_e = 0,
                   sigma2 = 1,
                   sigma2_mu = 1,
                   sigma2_psi = 1,
                   a_g = 1,
                   b_g = 1,
                   a_e = 1,
                   b_e = 1,
                   nburn = 1000,
                   npost = 1000,
                   nthin = 1) {

  # Extract the categories for genotype and environment

  cov_g = x[,'g']
  cov_e = x[,'e']

  classes_g = sort(unique(cov_g))
  classes_e = sort(unique(cov_e))

  ng = as.numeric(tapply(cov_g, cov_g, length)) # the number of obs within each g_i
  ne = as.numeric(tapply(cov_e, cov_e, length)) # the number of obs within each e_j
  ncov = length(ng) + length(ne)

  x <- model.matrix(~ -1 + g + e, data=x,
                    contrasts.arg=list(g=contrasts(x$g, contrasts=F),
                                       e=contrasts(x$e, contrasts=F)))
  x_e = x[,grepl('^e', colnames(x))]
  x_g = x[,grepl('^g', colnames(x))]

  ### aux

  n_class_e   = length(classes_e) # number of distinct levels
  aux_comb_e  = 2:floor(n_class_e/2) # levels needed to create the combinations without redundances
  num_comb_e  = choose(n_class_e, aux_comb_e) # total of combinations 2x2, 3x3, up to floor(ne/2)xfloor(ne/2).
  prob_comb_e = num_comb_e/sum(num_comb_e) # probability of observating a combination 2x2, 3x3, up to floor(ne/2)xfloor(ne/2).

  n_class_g   = length(classes_g) # number of distinct levels
  aux_comb_g  = 2:floor(n_class_g/2) # levels needed to create the combinations without redundances
  num_comb_g  = choose(n_class_g, aux_comb_g) # total of combinations 2x2, 3x3, up to floor(ne/2)xfloor(ne/2).
  prob_comb_g = num_comb_g/sum(num_comb_g) # probability of observating a combination 2x2, 3x3, up to floor(ne/2)xfloor(ne/2).

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  sigma2_g_store = rep(NA, store_size)
  sigma2_e_store = rep(NA, store_size)
  bart_store = matrix(NA, ncol = length(y), nrow = store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  var_count = rep(0, ncol(x))
  var_count_store = matrix(0, ncol = ncol(x), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(x), nrow = store_size)
  g_e_hat_store = matrix(0, ncol = ncol(x), nrow = store_size)
  colnames(g_e_hat_store) <- colnames(x)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  p = ncol(x)
  s = rep(1/p, p)
  p_g = ncol(x_g)
  p_e = ncol(x_e)
  s_g = rep(1/p_g, p_g)
  s_e = rep(1/p_e, p_e)
  sigma2_psi = diag(p)*sigma2_psi
  sigma2_psi_inv = solve(sigma2_psi)
  linear_effects = rep(0, length(y))
  yhat_bart = rep(0, length(y))

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = y_scale)
  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  yhat_bart = get_predictions(curr_trees, x, single_tree = ntrees == 1, internal=TRUE)

  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = TotIter,
                             style = 3, width = 60,
                             title = 'Running rBART...')

  # Start the MCMC iterations loop
  for (i in 1:TotIter) {

    utils::setTxtProgressBar(pb, i)

    # If at the right place, store everything
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      curr                 = (i - nburn)/nthin
      tree_store[[curr]]   = curr_trees
      sigma2_store[curr]   = sigma2
      sigma2_g_store[curr] = sigma2_g
      sigma2_e_store[curr] = sigma2_e
      bart_store[curr,]    = yhat_bart
      y_hat_store[curr,]   = y_hat
      g_e_hat_store[curr,] = g_e_hat
    }

    # Update the random effects alpha_i and beta_j
    g_e_hat = update_linear_component(y_scale, 0, x, sigma2, sigma2_psi_inv)
    linear_effects = x%*%g_e_hat

    # Start looping through trees
    for (j in 1:ntrees) {

      # partial residuals for the trees
      current_partial_residuals = y_scale - (yhat_bart + linear_effects) + tree_fits_store[,j]

      # Propose a new tree via grow/change/prune/swap
      type = sample(c('grow', 'prune'), 1)
      if(i < max(floor(0.1*nburn), 10)) type = 'grow' # Grow for the first few iterations

      # Generate a new tree based on the current

      if (type == 'grow' || (type=='prune' && nrow(curr_trees[[j]]$tree_matrix) == 1)){

        x_e_inter   = create_interaction(x_e, n_class_e, aux_comb_e, prob_comb_e)
        x_g_inter   = create_interaction(x_g, n_class_g, aux_comb_g, prob_comb_g)
        x_e_g_inter = cbind(x_e_inter, x_g_inter)

        # Below, there are two calls because we need to add an interaction of genotype and then
        # add to the same tree an interaction of environment, otherwise we run the risk of allowing
        # confunding.

        new_trees[[j]] = update_tree(X             = x_e_inter,
                                     type          = type,
                                     curr_tree     = curr_trees[[j]],
                                     node_min_size = node_min_size,
                                     s             = 1,
                                     index         = colnames(x_e_inter))

        new_trees[[j]] = update_tree(X             = x_g_inter,
                                     type          = type,
                                     curr_tree     = new_trees[[j]],
                                     node_min_size = node_min_size,
                                     s             = 1,
                                     index         = colnames(x_g_inter))
      } else {

        # We can't prune a terminal node because we run the risk of removing either
        # a genotype or an environment. If we remove one of them, the predicted values
        # from the tree will be confounded with the main effect associated to the enviroment/genotype
        # removed.

        new_trees[[j]] = create_stump(num_trees = 1,
                                      y = y_scale)[[1]]
      }

      # NEW TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_new = tree_full_conditional(new_trees[[j]],
                                    current_partial_residuals,
                                    sigma2,
                                    sigma2_mu) +
        get_tree_prior(new_trees[[j]], alpha, beta)

      # CURRENT TREE: compute the log of the marginalised likelihood + log of the tree prior
      l_old = tree_full_conditional(curr_trees[[j]],
                                    current_partial_residuals,
                                    sigma2,
                                    sigma2_mu) +
        get_tree_prior(curr_trees[[j]], alpha, beta)

      # Exponentiate the results above
      a = exp(l_new - l_old)

      # The current tree "becomes" the new tree, if the latter is better
      if(a > runif(1)) {
        curr_trees[[j]] = new_trees[[j]]
      }

      # Update mu whether tree accepted or not
      curr_trees[[j]] = simulate_mu(curr_trees[[j]],
                                    current_partial_residuals,
                                    sigma2,
                                    sigma2_mu)

      current_fit = get_predictions(curr_trees[[j]], x_e_g_inter, single_tree = TRUE, internal=TRUE)
      yhat_bart = yhat_bart - tree_fits_store[,j] # subtract the old fit
      yhat_bart = yhat_bart + current_fit # add the new fit
      tree_fits_store[,j] = current_fit # update the new fit

    } # End loop through trees

    # Updating the final predictions
    y_hat = linear_effects + yhat_bart

    gi = g_e_hat[grepl('^g', colnames(x))]
    ej = g_e_hat[grepl('^e', colnames(x))]

    sum_of_squares   = sum((y_scale - y_hat)^2)
    sum_of_squares_g = sum((gi - mu_g)^2)
    sum_of_squares_e = sum((ej - mu_e)^2)

    # Update sigma2 (variance of the residuals)
    sigma2   = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)
    sigma2_g = update_sigma2_g(sum_of_squares_g, length(gi), a_g, b_g)
    sigma2_e = update_sigma2_e(sum_of_squares_e, length(ej), a_e, b_e)
    aux_s2eg = c(rep(1/sigma2_g, length(gi)), rep(1/sigma2_e, length(ej)))
    sigma2_psi_inv = diag(aux_s2eg)

  } # End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  return(list(trees      = tree_store,
              sigma2     = sigma2_store*y_sd^2,
              sigma2_g   = sigma2_g_store,
              sigma2_e   = sigma2_e_store,
              y_hat      = y_hat_store*y_sd + y_mean,
              y_hat_bart = bart_store*y_sd,
              npost      = npost,
              nburn      = nburn,
              nthin      = nthin,
              ntrees     = ntrees,
              y_mean     = y_mean,
              y_sd       = y_sd,
              g_hat      = g_e_hat_store[, grepl('^g', colnames(x))]*y_sd,
              e_hat      = g_e_hat_store[, grepl('^e', colnames(x))]*y_sd,
              x          = x))

} # End main function
