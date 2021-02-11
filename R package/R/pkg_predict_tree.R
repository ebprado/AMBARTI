#' @export
predict_ambarti = function(object, newdata,
                           type = c('all', 'median', 'mean')) {

  x_train_names = colnames(object$x)
  x_e_train_names = x_train_names[grepl('e', colnames(x_train_names))]
  x_g_train_names = x_train_names[grepl('g', colnames(x_train_names))]

  new_x <- model.matrix(~ -1 + g + e, data=newdata,
                    contrasts.arg=list(g=contrasts(newdata$g, contrasts=F),
                                       e=contrasts(newdata$e, contrasts=F)))

  x_e_test = new_x[,grepl('e', colnames(new_x))]
  x_g_test = new_x[,grepl('g', colnames(new_x))]

  n_its = object$npost
  y_hat_mat = matrix(NA, nrow = n_its,
                     ncol = nrow(new_x))

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
  }

  # Sort out what to return
  out = switch(type,
               all = object$y_mean + object$y_sd * y_hat_mat,
               mean = object$y_mean + new_x%*%apply(object$beta_hat,2,'mean') + object$y_sd * apply(y_hat_mat,2,'mean'),
               median = object$y_mean + new_x%*%apply(object$beta_hat,2,'median') + object$y_sd * apply(y_hat_mat,2,'median'))

  return(out)

} # end of predict function
