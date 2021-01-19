library(devtools)
install_github("ebprado/AMBARTI/R package",
               ref = 'main',
               auth_token = '363d4ad84eaa25d3eb26752732cc208f7e698086')
library(AMBARTI)

save_file = "/Users/estevaoprado/Documents/GitHub/AMBARTI/04_simulation/results/"

I = 10 # c(5, 15, 30) # Number of genotypes
J = 10 # c(5, 15, 30) # Number of environments
s_alpha = 1 # c(1, 5) # standard deviation of alpha
s_beta = 1 # c(1, 5) # standard deviation of beta
s_y = c(1, 5) # standard deviation of y
lambda = c(8, 12) # c('8', '12', '8, 12', '10, 12','8, 10, 12')
n_rep = 10 # Number of Monte Carlo repetition

# Get all combinations of the quantities above
all_comb = expand.grid(I = I,
                       J = J,
                       s_alpha = s_alpha,
                       s_beta = s_beta,
                       s_y = s_y,
                       lambda = lambda)

# Get the number of combinations
n_comb = nrow(all_comb)
save_results = NULL

for (i in 1:n_comb){

  comb = all_comb[i,] # Get the row of the combination i
  I = comb$I # Number of genotypes
  J = comb$J # Number of environments
  s_alpha = comb$s_alpha # standard deviation of alpha
  s_beta = comb$s_beta # standard deviation of alpha
  s_y = comb$s_y # standard deviation of y
  aux_lambda = as.character(comb$lambda)
  lambda = as.numeric(unlist(strsplit(aux_lambda,','))) # values for lambda

  for (j in 1:n_rep){
    # load files
    filename = paste('I', I,
                     'J', J,
                     'sa', s_alpha,
                     'sb', s_beta,
                     'sy', s_y,
                     'L', gsub(', ', '', toString(lambda)),
                     'r', j,
                     sep='')

    data_filename    = paste(filename, '_data.RData',           sep='')
    ammi_filename    = paste(filename, '_classical_AMMI.RData', sep='')
    bammi_filename   = paste(filename, '_bayesian_AMMI.RData',  sep='')
    ambarti_filename = paste(filename, '_AMBARTI.RData',        sep='')

    load(paste(save_file, data_filename,      sep=''))
    load(paste(save_file, ammi_filename,      sep=''))
    load(paste(save_file, bammi_filename,     sep=''))
    load(paste(save_file, ambarti_filename,   sep=''))

    # Get parameter estimates from the classical AMMI
    res_AMMI = organise_classical_AMMI(classical_AMMI, data)

    # Get parameter estimates from the Bayesian AMMI
    res_bAMMI_with_post = organise_bayesian_AMMI_WITH_postprocessing(bayesian_AMMI, data)
    res_bAMMI_without_post = organise_bayesian_AMMI_WITHOUT_postprocessing(bayesian_AMMI, data)

    # Get parameter estimates from AMBARTI
    res_ambarti = organise_AMBARTI(ambarti, data)

    # Compute the Relative Root Mean Square Error (RRMSE)
    metrics_AMMI       = get_metrics(res_AMMI, data, j)
    metrics_bAMMI_post = get_metrics(res_bAMMI_with_post, data, j)
    metrics_bAMMI      = get_metrics(res_bAMMI_without_post, data, j)
    metrics_ambarti    = get_metrics(res_ambarti, data, j)

    save_results = rbind(save_results,
                         metrics_AMMI,
                         metrics_bAMMI_post,
                         metrics_bAMMI,
                         metrics_ambarti)
  }
}


