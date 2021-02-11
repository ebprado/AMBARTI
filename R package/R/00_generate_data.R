#' @export
#' @importFrom truncnorm 'rtruncnorm'
#'
generate_data_AMMI <- function(I, # Number of genotypes
                               J, # Number of environments
                               s_alpha, # standard deviation of alpha
                               s_beta, # standard deviation of alpha
                               s_y, # standard deviation of y
                               lambda # values for lambda (number of Q)
){

  # Total number of observations
  N = I*J

  # Number of components in the bilinear part
  Q = length(lambda)

  # Generate alpha (genotypes)
  alpha = rnorm(I, 0, s_alpha)

  # Generate beta (environments)
  beta = rnorm(J, 0, s_beta)

  # Set the grand mean
  mu = 100

  # Generate gamma
  # gamma <- matrix(NA, nrow = I ,ncol = Q)
  # gamma[1,] <- truncnorm::rtruncnorm(Q, a=0)
  # gamma[-1,] <- rnorm((I-1)*Q)
  gamma = generate_gamma_delta(I, Q)

  # Generate delta
  # delta <- matrix(rnorm(J*Q), nrow = J ,ncol = Q)
  delta = generate_gamma_delta(J, Q)
  # Generate the "design matrix"
  x = expand.grid(1:I, 1:J)
  names(x) <- c('g', 'e') # g = genotype and e = envorinment
  x$g = as.factor(x$g)
  x$e = as.factor(x$e)

  # Generate the interaction/bilinear part
  blin = rep(0, I*J)
  for (k in 1:length(lambda)) {
    blin <- blin + lambda[k]*gamma[x[,'g'],k]*delta[x[,'e'],k]
  }

  # Now simulate the response
  mu_ij = mu + alpha[x[,'g']] + beta[x[,'e']] + blin

  # Compute the response for the TRAINING set
  y = rnorm(N, mu_ij, s_y)

  # Compute the response for the TEST set
  y_test = rnorm(N, mu_ij, s_y)
  # I1 <- rep(1,I)
  # J1 <- rep(1,J)

  # This is a matrix representation from Alessandra (it works fine)
  # mu_ij <- mu*I1%*%t(J1) + kronecker(alpha,t(J1)) + kronecker(t(beta), (I1)) + gamma%*%diag(lambda)%*%t(delta)
  # y <- rnorm(N, c(mu.Y), s_y)

  return(list(y       = y,
              y_test  = y_test,
              x       = x,
              I       = I,
              J       = J,
              Q       = Q,
              s_alpha = s_alpha,
              s_beta  = s_beta,
              s_y     = s_y,
              lambda  = lambda,
              alpha   = alpha,
              beta    = beta,
              gamma   = gamma,
              delta   = delta,
              blinear = blin))
}

square_root_matrix <- function(x){

  # When Q = 1, x will be a scalar
  if (nrow(x) == 1) {return(sqrt(x))}

  # When Q > 1, then x will be a matrix
  if (nrow(x) > 1) {
    # Jordan normal form
    X = eigen(x)
    P = X$vectors
    A = diag(X$values)

    A_sqrt = diag(sqrt(X$values))
    P_inv = solve(P)
    x_sqrt = P %*% A_sqrt %*%  P_inv
    return(x_sqrt)
  }
}

generate_gamma_delta <- function(INDEX, Q) {

  first_row = TRUE

  while(first_row) {
    raw_par = matrix(rnorm(INDEX*Q), ncol=Q)
    par_mean  = matrix(rep(apply(raw_par,2,mean), each = nrow(raw_par)), ncol=Q)
    par_aux  = raw_par - par_mean

    # Constraints ----
    # apply(par_aux,2,sum)
    parTpar = solve(t(par_aux)%*%(par_aux))
    A = square_root_matrix(parTpar)
    samples = par_aux%*%A

    # Force the first to be positive
    for (i in 1:nrow(samples)){
      row1 = samples[1, ]
      if (all(samples[i, ] > 0)) {
        aux = samples[i, ]
        samples[1,] = aux
        samples[i,] = row1
        return(samples)
      }
    }
    # t(samples)%*%samples == 0
    # apply(samples,2,sum) == diag(Q)
  }
}


#' @export
#' @importFrom stats 'rnorm' 'aggregate' 'contrasts' 'model.matrix' 'as.formula'

generate_data_AMBARTI = function(I,
                                 J,
                                 s_alpha, # standard deviation of alpha
                                 s_beta, # standard deviation of alpha
                                 s_y, # standard deviation of y
                                 ntrees = 200,
                                 node_min_size = 5,
                                 mu_mu = 0,
                                 sigma2_mu = 3) {

  # Generate the "design matrix"
  N = I*J
  x = expand.grid(1:I, 1:J)
  names(x) <- c('g', 'e') # g = genotype and e = envorinment
  x$g = as.factor(x$g)
  x$e = as.factor(x$e)
  x_orig = x

  # Extract the categories for genotype and environment
  cov_g = x[,'g']
  cov_e = x[,'e']

  classes_g = unique(cov_g)
  classes_e = unique(cov_e)

  ng = as.numeric(tapply(cov_g, cov_g, length)) # the number of obs within each g_i
  ne = as.numeric(tapply(cov_e, cov_e, length)) # the number of obs within each e_j
  ncov = length(ng) + length(ne)

  x <- model.matrix(~ -1 + g + e, data=x,
                    contrasts.arg=list(g=contrasts(x$g, contrasts=F),
                                       e=contrasts(x$e, contrasts=F)))

  x_e = x[,grepl('e', colnames(x))]
  x_g = x[,grepl('g', colnames(x))]

  ### aux

  n_class_e   = length(classes_e) # number of distinct levels
  aux_comb_e  = 2:floor(n_class_e/2) # levels needed to create the combinations without redundances
  num_comb_e  = choose(n_class_e, aux_comb_e) # total of combinations 2x2, 3x3, up to floor(ne/2)xfloor(ne/2).
  prob_comb_e = num_comb_e/sum(num_comb_e) # probability of observating a combination 2x2, 3x3, up to floor(ne/2)xfloor(ne/2).

  n_class_g   = length(classes_g) # number of distinct levels
  aux_comb_g  = 2:floor(n_class_g/2) # levels needed to create the combinations without redundances
  num_comb_g  = choose(n_class_g, aux_comb_g) # total of combinations 2x2, 3x3, up to floor(ne/2)xfloor(ne/2).
  prob_comb_g = num_comb_g/sum(num_comb_g) # probability of observating a combination 2x2, 3x3, up to floor(ne/2)xfloor(ne/2).

  # Extract control parameters
  tree_store = vector('list', 1)
  p_g = ncol(x_g)
  p_e = ncol(x_e)
  s_g = rep(1/p_g, p_g)
  s_e = rep(1/p_e, p_e)

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees, y = 1:nrow(x))

  # Initialise the new trees as current one
  new_trees = curr_trees

  # Start looping through trees
  for (j in 1:ntrees) {

    # Propose a new tree
    type = 'grow'

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

    curr_trees[[j]] = new_trees[[j]]

    # Generate the predicted values
    curr_trees[[j]] = simulate_mu_bart_prior(curr_trees[[j]], mu_mu, sigma2_mu)

  } # End loop through trees

  # We do that because we assume that the sum of BART predictions within g_i and e_j is zero.
  bart_part = get_predictions(curr_trees, x_e_g_inter, single_tree = ntrees == 1, internal=TRUE)
  bart_part = bart_part - mean(bart_part)
  mean_by_g = aggregate(bart_part, by=list(x_orig[,'g']), mean)[,2]
  mean_by_e = aggregate(bart_part, by=list(x_orig[,'e']), mean)[,2]
  bart_part = bart_part - mean_by_g[x_orig[,'g']]
  bart_part = bart_part - mean_by_e[x_orig[,'e']]

  tree_store[[1]] = curr_trees

  # Generate alpha_i and beta_j
  alpha = rnorm(I, 0, sd = s_alpha)
  beta  = rnorm(J, 0, sd = s_beta)

  # Generate y's
  mu = 100
  y_train = rnorm(N, mu + alpha[cov_g] + beta[cov_e] + bart_part, sd=s_y)
  y_test  = rnorm(N, mu + alpha[cov_g] + beta[cov_e] + bart_part, sd=s_y)

  return(list(y       = y_train,
              y_test  = y_test,
              alpha   = alpha,
              beta    = beta,
              blinear = bart_part,
              x       = x_orig,
              I       = I,
              J       = J,
              s_alpha = s_alpha,
              s_beta  = s_beta,
              s_y     = s_y,
              ntrees  = ntrees,
              trees   = tree_store))

} # End main function



#' @export
#' @importFrom stats 'rnorm'
generate_data_full_model <- function(I, # Number of genotypes
                                     J, # Number of environments
                                     s_alpha, # standard deviation of alpha
                                     s_beta, # standard deviation of alpha
                                     s_y # standard deviation of y
){

  # Total number of observations
  N = I*J

  # Generate alpha (genotypes)
  alpha = rnorm(I, 0, s_alpha)

  # Generate beta (environments)
  beta = rnorm(J, 0, s_beta)

  # Set the grand mean
  mu = 100

  # Generate the "design matrix"
  x = expand.grid(1:I, 1:J)
  names(x) <- c('g', 'e') # g = genotype and e = envorinment
  x$g = as.factor(x$g)
  x$e = as.factor(x$e)

  # Now simulate the response
  inter = alpha[x[,'g']] * beta[x[,'e']]
  mu_ij = mu + alpha[x[,'g']] + beta[x[,'e']] + inter

  # Compute the response for the TRAINING set
  y = rnorm(N, mu_ij, s_y)

  # Compute the response for the TEST set
  y_test = rnorm(N, mu_ij, s_y)

  return(list(y       = y,
              y_test  = y_test,
              x       = x,
              I       = I,
              J       = J,
              s_alpha = s_alpha,
              s_beta  = s_beta,
              s_y     = s_y,
              alpha   = alpha,
              beta    = beta,
              blinear = inter))
}
