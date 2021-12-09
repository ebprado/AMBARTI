#' @export
predict_ambarti = function(object, newdata,
                           type = c('all', 'median', 'mean')) {

  x_train_names = colnames(object$x)
  x_e_train_names = x_train_names[grepl('^e', colnames(x_train_names))]
  x_g_train_names = x_train_names[grepl('^g', colnames(x_train_names))]

  new_x <- model.matrix(~ -1 + g + e, data=newdata,
                        contrasts.arg=list(g=contrasts(newdata$g, contrasts=F),
                                           e=contrasts(newdata$e, contrasts=F)))

  x_e_test = new_x[,grepl('^e', colnames(new_x))]
  x_g_test = new_x[,grepl('^g', colnames(new_x))]

  g_estimates_mean = apply(object$g_hat,2,'mean')
  e_estimates_mean = apply(object$e_hat,2,'mean')
  g_estimates_median = apply(object$g_hat,2,'median')
  e_estimates_median = apply(object$e_hat,2,'median')
  n_its = object$npost
  y_hat_mat = matrix(NA, nrow = n_its, ncol = nrow(new_x))
  new_y_bart_mat = matrix(NA, nrow=n_its, ncol=nrow(new_x))
  I = length(unique(newdata$g))
  J = length(unique(newdata$e))

  # Now loop through iterations and get predictions
  for (i in 1:n_its) {

    # Get current set of trees
    curr_trees = object$trees[[i]]

    # Create all covariates that were used in the training
    get_cov = unlist(lapply(curr_trees, function(x) x$tree_matrix[,'split_variable']))
    get_cov = get_cov[!is.na(get_cov)]
    new_x_g_e = create_covariates_prediction(get_cov, new_x)

    # Use get_predictions function to get predictions
    y_hat_mat[i,] = get_predictions(curr_trees,
                                    new_x_g_e,
                                    single_tree = length(curr_trees) == 2,
                                    internal=FALSE)

    # imposing the sum-to-zero constraint by row and column (only on the BART component)
    matBART = matrix(y_hat_mat[i,], nrow=I, ncol=J) # turns it into a g by e matrix
    row.mean = sweep(matBART*0,2, -colMeans(matBART)) # take column means
    col.mean = matBART*0 + rowMeans(matBART) # take row means
    matBARTdc = mean(matBART) - row.mean - col.mean + matBART # double center
    y_hat_bart_dc = as.numeric(matBARTdc) # turns it back into a (now constrained) vector
    new_y_bart_mat[i,] = y_hat_bart_dc

  }

  new_g_hat = colMeans(object$g_hat) - mean(colMeans(object$g_hat)) # constraint: g's sum to zero
  new_e_hat = colMeans(object$e_hat) - mean(colMeans(object$e_hat)) # constraint: e's sum to zero
  new_beta_hat = c(new_g_hat,new_e_hat)

  out = switch(type,
               mean =
                 list(wpostproc = object$y_mean + new_x%*%new_beta_hat + object$y_sd*apply(new_y_bart_mat,2,'mean'),
                      npostproc = object$y_mean + new_x%*%c(g_estimates_mean, e_estimates_mean) + object$y_sd * apply(y_hat_mat,2,'mean')),
               median = list(wpostproc = object$y_mean + new_x%*%new_beta_hat + object$y_sd* apply(new_y_bart_mat,2,'median'),
                             npostproc = object$y_mean + new_x%*%c(g_estimates_median, e_estimates_median) + object$y_sd * apply(y_hat_mat,2,'median'))
  )

  return(out)

} # end of predict function
