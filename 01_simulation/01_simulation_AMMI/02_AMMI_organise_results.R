library(devtools)
install_github("ebprado/AMBARTI/R package")
library(AMBARTI)

save_file = "~/R/AMBARTI/01_simulation_AMMI/results/"

I = c(10) # c(5, 25, 50) # Number of genotypes
J = c(10) # c(5, 25, 50) # Number of environments
s_alpha = c(1, 5) # c(1, 5) # standard deviation of alpha
s_beta = c(1,5) # c(1, 5) # standard deviation of beta
s_y = 1 # c(1, 5) # standard deviation of y
lambda = c('8', '12', '12, 8', '12, 10','12, 10, 8')
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

for (j in 1:n_rep){

  for (i in 1:n_comb){

    comb       = all_comb[i,] # Get the row of the combination i
    I          = comb$I # Number of genotypes
    J          = comb$J # Number of environments
    s_alpha    = comb$s_alpha # standard deviation of alpha
    s_beta     = comb$s_beta # standard deviation of alpha
    s_y        = comb$s_y # standard deviation of y
    aux_lambda = as.character(comb$lambda)
    lambda     = as.numeric(unlist(strsplit(aux_lambda,','))) # values for lambda
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
    # bammi_filename   = paste(filename, '_bayesian_AMMI.RData',  sep='')
    ambarti_filename = paste(filename, '_AMBARTI.RData',        sep='')

    load(paste(save_file, data_filename,      sep=''))
    load(paste(save_file, ammi_filename,      sep=''))
    # load(paste(save_file, bammi_filename,     sep=''))
    load(paste(save_file, ambarti_filename,   sep=''))

    # Get parameter estimates from the classical AMMI
    res_AMMI = organise_classical_AMMI(classical_AMMI, data)

    # Get parameter estimates from the Bayesian AMMI
    # res_bAMMI_with_post = organise_bayesian_AMMI_WITH_postprocessing(bayesian_AMMI, data)
    # res_bAMMI_without_post = organise_bayesian_AMMI_WITHOUT_postprocessing(bayesian_AMMI, data)

    # Get parameter estimates from AMBARTI
    res_ambarti = organise_AMBARTI(ambarti, data)

    # Compute the RMSE and RRMSE
    metrics_AMMI       = get_metrics(res_AMMI, data, j)
    metrics_bAMMI_post = get_metrics(res_bAMMI_with_post, data, j)
    metrics_bAMMI      = get_metrics(res_bAMMI_without_post, data, j)
    metrics_ambarti    = get_metrics(res_ambarti, data, j)

    save_results = rbind(save_results,
                         metrics_AMMI,
                         metrics_bAMMI_post,
                         metrics_bAMMI,
                         metrics_ambarti
                         )

    print(paste('comb = ', i, ' out of ', n_comb , ' (rep = ', j, ')', sep=''))
  }
}

save(save_results, file = paste(save_file, '00_results_consolidated.RData',    sep=''))
