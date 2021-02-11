# -------------------------------------------------------------------------#
# Description: this script contains functions that are used to perform the #
#              growing, pruning, changing, and swapping moves. It also has #
#              a function to initialise the trees to stumps                #
# -------------------------------------------------------------------------#

# 01. create_stump: initialises the trees to a stump
# 02. update_tree: calls the corresponding function associated to the move grow, prune, change, or swap.
# 03. grow_tree: grows a tree
# 04. prune_tree: prunes a tree
# 05. change_tree: changes the splitting rule that defines a pair of terminal nodes
# 06. swap_tree: exchanges the splitting rules that define two pair of terminal nodes

# Function to create stump ------------------------------------------------

create_stump = function(num_trees,
                        y) {

  # Each tree is a list of 2 elements
  # The 2 elements are the tree matrix (8 columns), and the node indices
  # The columns of the tree matrix are:
  # Terminal (0 = no, 1 = yes)
  # Child left
  # Child right
  # Node parents
  # Split variable
  # Split value
  # mu
  # Node size

  # Create holder for trees
  all_trees = vector('list', length = num_trees)
  # Loop through trees
  for (j in 1:num_trees) {
    # Set up each tree to have two elements in the list as described above
    all_trees[[j]] = vector('list', length = 2)
    # Give the elements names
    names(all_trees[[j]]) = c('tree_matrix',
                              'node_indices')
    # Create the two elements: first is a matrix
    all_trees[[j]][[1]] = matrix(NA, ncol = 8, nrow = 1)

    # Second is the assignment to node indices
    all_trees[[j]][[2]] = rep(1, length(y))

    # Create column names
    colnames(all_trees[[j]][[1]]) = c('terminal',
                                      'child_left',
                                      'child_right',
                                      'parent',
                                      'split_variable',
                                      'split_value',
                                      'mu',
                                      'node_size')

    # Set values for stump
    all_trees[[j]][[1]][1,] = c(1, NA, NA, NA, NA, NA, 0 , length(y))

  } # End of loop through trees

  return(all_trees)

} # End of function

# Function to update trees ------------------------------------------------

update_tree = function(X, # Feature matrix
                       type = c('grow',   # Grow existing tree
                                'prune',  # Prune existing tree
                                'change', # Change existing tree - change split variable and value for an internal node
                                'swap'),  # Swap existing tree - swap splitting rules for two pairs of terminal nodes
                       curr_tree,         # The current set of trees (not required if type is stump)
                       node_min_size,     # The minimum size of a node to grow
                       s,                 # probability vector to be used during the growing process
                       index)
  {

  # Call the appropriate function to get the new tree
  new_tree = switch(type,
                    grow = grow_tree(X, curr_tree, node_min_size, s, index),
                    prune = prune_tree(X, curr_tree))

  # Return the new tree
  return(new_tree)

} # End of update_tree function

# Grow_tree function ------------------------------------------------------

grow_tree = function(X, curr_tree, node_min_size, s, index) {

  # Set up holder for new tree
  new_tree = curr_tree

  # Get the list of terminal nodes
  terminal_nodes = as.numeric(which(new_tree$tree_matrix[,'terminal'] == 1))

  # Find terminal node sizes
  terminal_node_size = as.numeric(new_tree$tree_matrix[terminal_nodes,'node_size'])

  available_values = NULL
  max_bad_trees = 10
  count_bad_trees = 0
  bad_trees = TRUE

  while (bad_trees ){

    # Set up holder for new tree
    new_tree = curr_tree

    # Add two extra rows to the tree in question
    new_tree$tree_matrix = rbind(new_tree$tree_matrix,
                                 c(1, NA, NA, NA, NA, NA, NA, NA), # Make sure they're both terminal
                                 c(1, NA, NA, NA, NA, NA, NA, NA))

    # Choose a random terminal node to split
    node_to_split = sample(terminal_nodes, 1,
                           prob = as.integer(terminal_node_size > node_min_size)) # Choose which node to split, set prob to zero for any nodes that are too small

    # Choose a split variable uniformly from all columns (the first one is the intercept)
    split_variable = sample(index, 1, prob = s)

    # Alternatively follow BARTMachine and choose a split value using sample on the internal values of the available
    available_values = sort(unique(X[new_tree$node_indices == node_to_split,
                                     split_variable]))

    if(length(available_values) == 1){
      split_value = available_values[1]
    } else if (length(available_values) == 2){
      split_value = available_values[2]
    }  else {
      # split_value = sample(available_values[-c(1,length(available_values))], 1)
      split_value = resample(available_values[-c(1,length(available_values))])
    }

    curr_parent = new_tree$tree_matrix[node_to_split, 'parent'] # Make sure to keep the current parent in there. Will be NA if at the root node
    new_tree$tree_matrix[node_to_split,1:6] = c(0, # Now not temrinal
                                                nrow(new_tree$tree_matrix) - 1, # child_left is penultimate row
                                                nrow(new_tree$tree_matrix),  # child_right is penultimate row
                                                curr_parent,
                                                split_variable,
                                                split_value)

    #  Fill in the parents of these two nodes
    new_tree$tree_matrix[nrow(new_tree$tree_matrix),'parent'] = node_to_split
    new_tree$tree_matrix[nrow(new_tree$tree_matrix)-1,'parent'] = node_to_split

    # Now call the fill function on this tree
    #new_tree = fill_tree_details(new_tree, X)
    new_tree = fill_tree_details(new_tree, X, node_to_split)

    # Store the covariate name to use it to update the Dirichlet prior of Linero (2016).
    new_tree$var = split_variable

    # Check for bad tree
    if(any(as.numeric(new_tree$tree_matrix[,'node_size']) <= node_min_size)) {
      count_bad_trees = count_bad_trees + 1
    } else {
      bad_trees = FALSE
    }

    if(count_bad_trees == max_bad_trees) {
      curr_tree$var = 0
      return(curr_tree)
      }
  }
  # Return new_tree
  return(new_tree)

} # End of grow_tree function

# Prune_tree function -----------------------------------------------------

prune_tree = function(X, curr_tree) {

  # Create placeholder for new tree
  new_tree = curr_tree

  if(nrow(new_tree$tree_matrix) == 1) { # No point in pruning a stump!
    new_tree$var = 0
    return(new_tree)
  }

  # Get the list of terminal nodes
  terminal_nodes = which(as.numeric(new_tree$tree_matrix[,'terminal']) == 1)

  # Pick a random terminal node to prune
  # ONLY PICK NODES WHERE BOTH LEFT AND RIGHT CHILD ARE TERMINAL
  bad_node_to_prune = TRUE # Assume a bad node pick
  while(bad_node_to_prune) {

    # Choose a random terminal node
    node_to_prune = sample(terminal_nodes, 1)

    # Find the parent of this terminal node
    parent_pick = as.numeric(new_tree$tree_matrix[node_to_prune, 'parent'])
    var_pruned_nodes = as.numeric(new_tree$tree_matrix[parent_pick, 'split_variable'])

    # Get the two children of this parent
    child_left = as.numeric(new_tree$tree_matrix[parent_pick, 'child_left'])
    child_right = as.numeric(new_tree$tree_matrix[parent_pick, 'child_right'])

    # See whether either are terminal
    child_left_terminal = as.numeric(new_tree$tree_matrix[child_left, 'terminal'])
    child_right_terminal = as.numeric(new_tree$tree_matrix[child_right, 'terminal'])

    # If both are terminal then great
    if( (child_left_terminal == 1) & (child_right_terminal == 1) ) {
      bad_node_to_prune = FALSE # Have chosen a pair of terminal nodes so exist while loop
    }

  }# End of bad node to prune while loop

  # Delete these two rows from the tree matrix
  new_tree$tree_matrix = new_tree$tree_matrix[-c(child_left,child_right),,
                                              drop = FALSE]
  # Make this node terminal again with no children or split values
  new_tree$tree_matrix[parent_pick,c('terminal',
                                     'child_left',
                                     'child_right',
                                     'split_variable',
                                     'split_value')] = c(1, NA, NA, NA, NA)

  # If we're back to a stump no need to call fill_tree_details
  if(nrow(new_tree$tree_matrix) == 1) {
    new_tree$var = var_pruned_nodes
    new_tree$node_indices = rep(1, nrow(X))
  } else {
    # If we've removed some nodes from the middle we need to re-number all the child_left and child_right values - the parent values will still be correct
    if(node_to_prune <= nrow(new_tree$tree_matrix)) { # Only need do this if we've removed some observations from the middle of the tree matrix
      # If you're pruning any nodes which affect parent indices further down the tree then make sure to shift the parent values
      bad_parents = which(as.numeric(new_tree$tree_matrix[,'parent'])>=node_to_prune)
      # Shift them back because you have removed two rows
      new_tree$tree_matrix[bad_parents,'parent'] = as.numeric(new_tree$tree_matrix[bad_parents,'parent']) - 2

      for(j in node_to_prune:nrow(new_tree$tree_matrix)) {
        # Find the current parent
        curr_parent = as.numeric(new_tree$tree_matrix[j,'parent'])
        # Find both the children of this node
        curr_children = which(as.numeric(new_tree$tree_matrix[,'parent']) == curr_parent)
        # Input these children back into the parent
        new_tree$tree_matrix[curr_parent,c('child_left','child_right')] = sort(curr_children)
      } # End for loop of correcting parents and children
    } # End if statement to fill in tree details

    # Call the fill function on this tree
    # new_tree = fill_tree_details(new_tree, X, node_to_prune)
    # Get the old indices
    get_old_indices = which(new_tree$node_indices == c(child_left, child_right))
    # Replace the old node indices by the indice of their parent node
    new_tree$node_indices[get_old_indices] = parent_pick

    # Store the covariate name that was used in the splitting rule of the terminal nodes that were just pruned
    new_tree$var = var_pruned_nodes

  }

  # Return new_tree
  return(new_tree)

} # End of prune_tree function