#' @export
#' @importFrom truncnorm 'rtruncnorm'
#'
generate_data <- function(I, # Number of genotypes
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
  for(k in 2:nrow(gamma)){
    gamma[k,] <- rnorm(Q)
  }

  # Generate delta
  delta <- matrix(rnorm(J*Q), nrow = J ,ncol = Q)

  # Generate the "design matrix"
  x = expand.grid(1:I, 1:J)
  names(x) <- c('gen', 'env') # gen = genotype and env = envorinment

  # Generate the interaction/bilinear part
  blin = rep(0, I*J)
  for (k in 1:length(lambda)) {
    blin <- blin + lambda[k]*gamma[x[,1],k]*delta[x[,2],k]
  }

  # Now simulate the response
  mu_ij = mu + alpha[x[,1]] + beta[x[,2]] + blin

  # Compute the response for the TRAINING set
  y_train = rnorm(N, mu_ij, s_y)

  # Compute the response for the TEST set
  y_test = rnorm(N, mu_ij, s_y)
  # I1 <- rep(1,I)
  # J1 <- rep(1,J)

  # This is a matrix representation from Alessandra (it works fine)
  # mu_ij <- mu*I1%*%t(J1) + kronecker(alpha,t(J1)) + kronecker(t(beta), (I1)) + gamma%*%diag(lambda)%*%t(delta)
  # y <- rnorm(N, c(mu.Y), s_y)

  return(list(y_train = y_train,
              y_test  = y_test,
              x       = x,
              I       = I,
              J       = J,
              s_alpha = s_alpha,
              s_beta  = s_beta,
              s_y     = s_y,
              lambda  = lambda,
              alpha   = alpha,
              beta    = beta,
              gamma   = gamma,
              delta   = delta))
}
