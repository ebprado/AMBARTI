# -------------------------------------------------------------------------#
# Description: this script contains auxiliar functions needed to update    #
# the trees with details and to map the predicted values to each obs       #
# -------------------------------------------------------------------------#

# 1. fill_tree_details: takes a tree matrix and returns the number of obs in each node in it and the indices of each observation in each terminal node
# 2. get_predictions: gets the predicted values from a current set of trees
# 3. get_children: it's a function that takes a node and, if the node is terminal, returns the node. If not, returns the children and calls the function again on the children
# 4. resample: an auxiliar function
# 5. create_interactions: create a new column that represents an interaction between two columns of a given X.
# 5. create_covariates_prediction: for each MCMC iteration, it creates the covariates needed for each tree.
# Fill_tree_details -------------------------------------------------------

fill_tree_details = function(curr_tree, X, node_to_split = NULL) {

  # Collect right bits of tree
  tree_matrix = curr_tree$tree_matrix

  # tree_matrix = curr_tree$tree_matrix

  # Create a new tree matrix to overwrite
  new_tree_matrix = tree_matrix

  # Start with dummy node indices

  if (is.null(node_to_split)== FALSE) {
    loop_indices = which(tree_matrix[,'parent'] == node_to_split) # only nodes that were just created
    node_indices = curr_tree$node_indices
  } else {
    loop_indices = 2:nrow(tree_matrix)
    node_indices = rep(1, nrow(X))
  }

  # For all but the top row, find the number of observations falling into each one
  for(i in loop_indices) {

    # Get the parent
    curr_parent = as.numeric(tree_matrix[i,'parent'])

    # Find the split variable and value of the parent
    split_var = tree_matrix[curr_parent,'split_variable']
    split_val = as.numeric(tree_matrix[curr_parent, 'split_value'])

    # Find whether it's a left or right terminal node
    left_or_right = ifelse(tree_matrix[curr_parent,'child_left'] == i,
                           'left', 'right')
    if(left_or_right == 'left') {
      # If left use less than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] < split_val)
      node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] < split_val] = i
    } else {
      # If right use greater than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] >= split_val)
      node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] >= split_val] = i
    }
  } # End of loop through table

  return(list(tree_matrix = new_tree_matrix,
              node_indices = node_indices))

} # End of function


# Get predictions ---------------------------------------------------------

get_predictions = function(trees, X, single_tree = FALSE, internal) {

  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]

  # Normally trees will be a list of lists but just in case
  if(single_tree) {
    # Deal with just a single tree
    if(nrow(trees$tree_matrix) == 1) {
      predictions = rep(trees$tree_matrix[1, 'mu'], nrow(X))
    } else {
      # Loop through the node indices to get predictions
      predictions = rep(NA, nrow(X))
      unique_node_indices = unique(trees$node_indices)
      # Get the node indices for the current X matrix
      if (internal==TRUE) {
        curr_X_node_indices = trees$node_indices
      } else{
        curr_X_node_indices = fill_tree_details(trees, X)$node_indices
        }

      # Now loop through all node indices to fill in details
      for(i in 1:length(unique_node_indices)) {
        predictions[curr_X_node_indices == unique_node_indices[i]] =
          as.numeric(trees$tree_matrix[unique_node_indices[i], 'mu'])
      }
    }
    # More here to deal with more complicated trees - i.e. multiple trees
  } else {
    # Do a recursive call to the function
    partial_trees = trees
    partial_trees[[1]] = NULL # Blank out that element of the list
    predictions = get_predictions(trees[[1]], X, single_tree = TRUE, internal)  +
      get_predictions(partial_trees, X,
                      single_tree = length(partial_trees) == 1, internal)
    #single_tree = !is.null(names(partial_trees)))
    # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)
  }

  return(predictions)
}

# get_children ------------------------------------------------------------

get_children = function(tree_mat, parent) {
  # Create a holder for the children
  all_children = NULL
  if(as.numeric(tree_mat[parent,'terminal']) == 1) {
    # If the node is terminal return the list so far
    return(c(all_children, parent))
  } else {
    # If not get the current children
    curr_child_left = as.numeric(tree_mat[parent, 'child_left'])
    curr_child_right = as.numeric(tree_mat[parent, 'child_right'])
    # Return the children and also the children of the children recursively
    return(c(all_children,
             get_children(tree_mat,curr_child_left),
             get_children(tree_mat,curr_child_right)))
  }
}

# Sample function ----------------------------------------------------------

resample <- function(x, ...) x[sample.int(length(x), size=1), ...]

# Create interactions -------

create_interaction = function(X, nclass, aux_comb, prob){
  num_cov = sample(aux_comb, 1, prob = prob) # number of covariates to be sampled
  s_covs  = sample(x = 1:nclass, size = num_cov, replace = FALSE) # sampled covariates
  # new_cov = matrix(apply(X[,s_covs], 1,sum), ncol=1) # create the interaction from the sampled covariates
  new_cov = matrix(rowSums(X[,s_covs]), ncol=1) # create the interaction from the sampled covariates
  colnames(new_cov) = paste(sort(colnames(X)[s_covs]), collapse = ':')
  return(new_cov)
}

create_covariates_prediction = function(variables, data){
  if (length(variables)) {
    mat_covs = matrix(NA, nrow=nrow(data), ncol=length(variables))
    colnames(mat_covs) = variables
    for (k in 1:length(variables)){
      columns = unlist(strsplit(variables[k], ':'))
      mat_covs[,k] = apply(data[,columns],1,sum)
    }
    return(mat_covs)
  }
  mat_covs = matrix(NA, nrow=nrow(data), ncol=1)
  return(mat_covs)

}
