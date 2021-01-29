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
  gamma <- matrix(NA, nrow = I ,ncol = Q)
  gamma[1,] <- truncnorm::rtruncnorm(Q, a=0)
  gamma[-1,] <- rnorm((I-1)*Q)

  # Generate delta
  delta <- matrix(rnorm(J*Q), nrow = J ,ncol = Q)

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

  ###########################################
  #### Genotype
  ###########################################
  number_geno = length(ng)
  num_comb_g = max(2, floor(number_geno/2))
  formula_g = as.formula(paste('y', "~", '.^', num_comb_g))
  x_all_iter_g <- model.matrix( formula_g, data = data.frame(y = rep(0, nrow(x)), x[, grepl("g", colnames(x))]))
  individual_g = (1:(number_geno + 1))
  name_all_comb_g = colnames(x_all_iter_g)
  name_all_comb_g = name_all_comb_g[-c(individual_g)] # remove the individual effects

  if (number_geno%%2 == 0){ #even
    repeated_comb_g = choose(number_geno, num_comb_g)/2
    name_all_comb_g = name_all_comb_g[-c(repeated_comb_g)] # remove some equivalent columns
  }

  x_g = matrix(NA, ncol=length(name_all_comb_g), nrow=nrow(x))
  colnames(x_g) = name_all_comb_g

  for (k in 1:ncol(x_g)){
    name_col_g = unlist(strsplit(name_all_comb_g[k],':'))
    x_g[,k] = apply(x[,name_col_g],1,sum)
  }

  ###########################################
  #### Environment
  ###########################################

  number_env = length(ne)
  num_comb_e = max(2, floor(number_env/2))
  formula_e = as.formula(paste('y', "~", '.^', num_comb_e))
  x_all_iter_e <- model.matrix(formula_e, data = data.frame(y = rep(0, nrow(x)), x[, grepl("e", colnames(x))]))
  individual_e = (1:(number_env + 1))
  repeated_comb_e = choose(number_env, num_comb_e)/2
  name_all_comb_e = colnames(x_all_iter_e)
  name_all_comb_e = name_all_comb_e[-c(individual_e)] # remove individual effects

  if (number_env%%2 == 0){ #even
    repeated_comb_e = choose(number_env, num_comb_e)/2
    name_all_comb_e = name_all_comb_e[-c(repeated_comb_e)] # remove some equivalent columns
  }

  x_e = matrix(NA, ncol=length(name_all_comb_e), nrow=nrow(x))
  colnames(x_e) = name_all_comb_e

  for (k in 1:ncol(x_e)){
    name_col_e = unlist(strsplit(name_all_comb_e[k],':'))
    x_e[,k] = apply(x[,name_col_e],1,sum)
  }
  # Put x_g and x_e into a data frame and get the column indices
  if (number_env <= 2 || number_geno <= 2) {stop('Number of genotypes and environments needs to be >= 3.')}
  x_g_e = as.data.frame(cbind(x_g, x_e))
  ind_x_g = 1:ncol(x_g)
  ind_x_e = (ncol(x_g) + 1):ncol(x_g_e)

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

    # Below, there are two calls because we need to add an interaction of genotype and then
    # add to the same tree an interaction of environment, otherwise we run the risk of allowing
    # confunding.

    new_trees[[j]] = update_tree(X = x_g_e,
                                 type = type,
                                 curr_tree = curr_trees[[j]],
                                 node_min_size = node_min_size,
                                 s = s_g,
                                 index = ind_x_g)

    new_trees[[j]] = update_tree(X = x_g_e,
                                 type = type,
                                 curr_tree = new_trees[[j]],
                                 node_min_size = node_min_size,
                                 s = s_e,
                                 index = ind_x_e)

    curr_trees[[j]] = new_trees[[j]]

    # Generate the predicted values
    curr_trees[[j]] = simulate_mu_bart_prior(curr_trees[[j]], mu_mu, sigma2_mu)

  } # End loop through trees

  bart_part = get_predictions(curr_trees, x_g_e, single_tree = ntrees == 1)
  # hist(bart_part - mean(bart_part))
  # We take the mean because we assume that the sum of BART predictions within g_i and e_j is zero.
  bart_part = bart_part - mean(bart_part)

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
