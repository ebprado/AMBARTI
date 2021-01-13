
#' @export
predict_ambarti = function(object, newdata,
                           type = c('all', 'median', 'mean')) {

  # Extract the categories for genotype and environment

  cov_g = newdata[,'g']
  cov_e = newdata[,'e']

  classes_g = sort(unique(cov_g))
  classes_e = sort(unique(cov_e))

  ng = as.numeric(tapply(cov_g, cov_g, length)) # the number of obs within each g_i
  ne = as.numeric(tapply(cov_e, cov_e, length)) # the number of obs within each e_j
  ncov = length(ng) + length(ne)

  new_x <- model.matrix(~ -1 + g + e, data=newdata,
                        contrasts.arg=list(g=contrasts(as.factor(newdata$g), contrasts=F),
                                           e=contrasts(as.factor(newdata$e), contrasts=F)))

  y_aux = rep(0, nrow(newdata))

  ###########################################
  #### Genotype
  ###########################################
  number_geno = length(ng)
  num_comb_g = floor(number_geno/2)
  formula_g = as.formula(paste('y', "~", '.^', num_comb_g))
  x_all_iter_g <- model.matrix( formula_g, data = data.frame(y = y_aux, new_x[, grepl("g", colnames(new_x))]))
  individual_g = (1:(number_geno + 1))
  name_all_comb_g = colnames(x_all_iter_g)
  name_all_comb_g = name_all_comb_g[-c(individual_g)] # remove the individual effects

  if (number_geno%%2 == 0){ #even
    repeated_comb_g = choose(number_geno, num_comb_g)/2
    name_all_comb_g = name_all_comb_g[-c(repeated_comb_g)] # remove some equivalent columns
  }

  new_x_g = matrix(NA, ncol=length(name_all_comb_g), nrow=nrow(new_x))
  colnames(new_x_g) = name_all_comb_g

  for (k in 1:ncol(new_x_g)){
    name_col_g = unlist(strsplit(name_all_comb_g[k],':'))
    new_x_g[,k] = apply(new_x[,name_col_g],1,sum)
  }

  ###########################################
  #### Environment
  ###########################################

  number_env = length(ne)
  num_comb_e = floor(number_env/2)
  formula_e = as.formula(paste('y', "~", '.^', num_comb_e))
  x_all_iter_e <- model.matrix(formula_e, data = data.frame(y = y_aux, new_x[, grepl("e", colnames(new_x))]))
  individual_e = (1:(number_env + 1))
  repeated_comb_e = choose(number_env, num_comb_e)/2
  name_all_comb_e = colnames(x_all_iter_e)
  name_all_comb_e = name_all_comb_e[-c(individual_e)] # remove individual effects

  if (length(ne)%%2 == 0){ #even
    repeated_comb_e = choose(number_env, num_comb_e)/2
    name_all_comb_e = name_all_comb_e[-c(repeated_comb_e)] # remove some equivalent columns
  }

  new_x_e = matrix(NA, ncol=length(name_all_comb_e), nrow=nrow(new_x))
  colnames(new_x_e) = name_all_comb_e

  for (k in 1:ncol(new_x_e)){
    name_col_e = unlist(strsplit(name_all_comb_e[k],':'))
    new_x_e[,k] = apply(new_x[,name_col_e],1,sum)
  }
  # Put x_g and x_e into a data frame and get the column indices
  new_x_g_e = as.data.frame(cbind(new_x_g, new_x_e))

  # Make sure the columns in the train and test data are in the same order due to the main effects
  if (ncol(new_x_g_e) == ncol(object$x_g_e)){
    new_x_g_e = new_x_g_e[,colnames(object$x_g_e)]
  } else {
    stop('The train and test data sets have either different structures or different factor levels.')
  }

  n_its = object$npost
  y_hat_mat = matrix(NA, nrow = n_its,
                     ncol = nrow(new_x))

  # Now loop through iterations and get predictions
  for (i in 1:n_its) {
    # Get current set of trees
    curr_trees = object$trees[[i]]

    # Use get_predictions function to get predictions
    y_hat_mat[i,] = get_predictions(curr_trees,
                                    new_x_g_e,
                                    single_tree = length(curr_trees) == 2)
  }

  # Sort out what to return
  out = switch(type,
               all = object$y_mean + object$y_sd * y_hat_mat,
               mean = object$y_mean + new_x%*%apply(object$beta_hat,2,'mean') + object$y_sd * apply(y_hat_mat,2,'mean'),
               median = object$y_mean + new_x%*%apply(object$beta_hat,2,'median') + object$y_sd * apply(y_hat_mat,2,'median'))

  return(out)

} # end of predict function
