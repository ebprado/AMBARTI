library(devtools)
install_github("ebprado/AMBARTI/R package",
               ref = 'main',
               auth_token = '363d4ad84eaa25d3eb26752732cc208f7e698086')
library(AMBARTI)

save_file = "/Users/estevaoprado/Documents/GitHub/AMBARTI/04_simulation/02_simulation_AMBARTI/results/"

I = c(10) # c(5, 15, 30) # Number of genotypes
J = c(10) # c(5, 15, 30) # Number of environments
s_alpha = c(1, 5) # c(1, 5) # standard deviation of alpha
s_beta = c(1,5) # c(1, 5) # standard deviation of beta
s_y = 1 # c(1, 5) # standard deviation of y
n_rep = 10 # Number of Monte Carlo repetition

# Get all combinations of the quantities above
all_comb = expand.grid(I = I,
                       J = J,
                       s_alpha = s_alpha,
                       s_beta = s_beta,
                       s_y = s_y)

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

  for (j in 1:n_rep){
    # load files
    filename = paste('AMBARTI',
                     'I', I,
                     'J', J,
                     'sa', s_alpha,
                     'sb', s_beta,
                     'sy', s_y,
                     'r', j,
                     sep='')

    data_filename    = paste(filename, '_data.RData',              sep='')
    ammi_filenameQ1  = paste(filename, '_classical_AMMI_Q1.RData', sep='')
    ammi_filenameQ2  = paste(filename, '_classical_AMMI_Q2.RData', sep='')
    ammi_filenameQ3  = paste(filename, '_classical_AMMI_Q3.RData', sep='')
    bammi_filenameQ1 = paste(filename, '_bayesian_AMMI_Q1.RData',  sep='')
    bammi_filenameQ2 = paste(filename, '_bayesian_AMMI_Q2.RData',  sep='')
    bammi_filenameQ3 = paste(filename, '_bayesian_AMMI_Q3.RData',  sep='')
    ambarti_filename = paste(filename, '_AMBARTI.RData',           sep='')

    load(paste(save_file, data_filename,    sep=''))
    load(paste(save_file, ammi_filenameQ1,  sep=''))
    load(paste(save_file, ammi_filenameQ2,  sep=''))
    load(paste(save_file, ammi_filenameQ3,  sep=''))
    load(paste(save_file, bammi_filenameQ1, sep=''))
    load(paste(save_file, bammi_filenameQ2, sep=''))
    load(paste(save_file, bammi_filenameQ3, sep=''))
    load(paste(save_file, ambarti_filename, sep=''))

    # Get parameter estimates from the classical AMMI
    res_AMMI_Q1 = organise_classical_AMMI(classical_AMMI_Q1, data, Q=1)
    res_AMMI_Q2 = organise_classical_AMMI(classical_AMMI_Q2, data, Q=2)
    res_AMMI_Q3 = organise_classical_AMMI(classical_AMMI_Q3, data, Q=3)

    # Get parameter estimates from the Bayesian AMMI
    res_bAMMI_with_post_Q1 = organise_bayesian_AMMI_WITH_postprocessing(bayesian_AMMI_Q1, data, Q=1)
    res_bAMMI_with_post_Q2 = organise_bayesian_AMMI_WITH_postprocessing(bayesian_AMMI_Q2, data, Q=2)
    res_bAMMI_with_post_Q3 = organise_bayesian_AMMI_WITH_postprocessing(bayesian_AMMI_Q3, data, Q=3)

    res_bAMMI_without_post_Q1 = organise_bayesian_AMMI_WITHOUT_postprocessing(bayesian_AMMI_Q1, data, Q=1)
    res_bAMMI_without_post_Q2 = organise_bayesian_AMMI_WITHOUT_postprocessing(bayesian_AMMI_Q2, data, Q=2)
    res_bAMMI_without_post_Q3 = organise_bayesian_AMMI_WITHOUT_postprocessing(bayesian_AMMI_Q3, data, Q=3)

    # Get parameter estimates from AMBARTI
    res_ambarti = organise_AMBARTI(ambarti, data)

    # Compute the RMSE and RRMSE
    metrics_AMMI_Q1 = get_metrics(res_AMMI_Q1, data, j)
    metrics_AMMI_Q2 = get_metrics(res_AMMI_Q2, data, j)
    metrics_AMMI_Q3 = get_metrics(res_AMMI_Q3, data, j)

    metrics_bAMMI_post_Q1 = get_metrics(res_bAMMI_with_post_Q1, data, j)
    metrics_bAMMI_post_Q2 = get_metrics(res_bAMMI_with_post_Q2, data, j)
    metrics_bAMMI_post_Q3 = get_metrics(res_bAMMI_with_post_Q3, data, j)

    metrics_bAMMI_Q1 = get_metrics(res_bAMMI_without_post_Q1, data, j)
    metrics_bAMMI_Q2 = get_metrics(res_bAMMI_without_post_Q2, data, j)
    metrics_bAMMI_Q3 = get_metrics(res_bAMMI_without_post_Q3, data, j)

    metrics_ambarti = get_metrics(res_ambarti, data, j)

    save_results = rbind(save_results,
                         metrics_AMMI_Q1,
                         metrics_AMMI_Q2,
                         metrics_AMMI_Q3,
                         metrics_bAMMI_post_Q1,
                         metrics_bAMMI_post_Q2,
                         metrics_bAMMI_post_Q3,
                         metrics_bAMMI_Q1,
                         metrics_bAMMI_Q2,
                         metrics_bAMMI_Q3,
                         metrics_ambarti)

    print(paste('comb = ', i, ' out of ', n_comb , ' (rep = ', j, ')', sep=''))
  }
}

save(save_results, file = paste(save_file, '00_AMB_results_consolidated.RData',    sep=''))
