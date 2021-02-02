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
nseed  = 0 # start seed

for (i in 1:n_comb){

  comb    = all_comb[i,] # Get the row of the combination i
  I       = comb$I # Number of genotypes
  J       = comb$J # Number of environments
  s_alpha = comb$s_alpha # standard deviation of alpha
  s_beta  = comb$s_beta # standard deviation of alpha
  s_y     = comb$s_y # standard deviation of y

  for (j in 1:n_rep){
    # Set a seed to make it reproducible
    set.seed(nseed)

    data = generate_data_AMBARTI(I, J, s_alpha, s_beta, s_y, ntrees = 200)

    # run classical AMMI
    classical_AMMI_Q1 = run_classical_AMMI(data, Q=1)
    classical_AMMI_Q2 = run_classical_AMMI(data, Q=2)
    classical_AMMI_Q3 = run_classical_AMMI(data, Q=3)

    # run Bayesian AMMI
    bayesian_AMMI_Q1 = run_bayesian_AMMI(data,Q=1)
    bayesian_AMMI_Q2 = run_bayesian_AMMI(data,Q=2)
    bayesian_AMMI_Q3 = run_bayesian_AMMI(data,Q=3)

    # run AMBARTI
    ambarti = run_AMBARTI(data, ntrees = 50, nburn = 2000, npost = 1000)

    # Increment the seed number by 1
    nseed = nseed + 1

    print(paste('comb = ', i, ' out of ', n_comb , ' (rep = ', j, ')', sep=''))

    # save files
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

    save(data,              file = paste(save_file, data_filename,    sep=''))
    save(classical_AMMI_Q1, file = paste(save_file, ammi_filenameQ1,  sep=''))
    save(classical_AMMI_Q2, file = paste(save_file, ammi_filenameQ2,  sep=''))
    save(classical_AMMI_Q3, file = paste(save_file, ammi_filenameQ3,  sep=''))
    save(bayesian_AMMI_Q1,  file = paste(save_file, bammi_filenameQ1, sep=''))
    save(bayesian_AMMI_Q2,  file = paste(save_file, bammi_filenameQ2, sep=''))
    save(bayesian_AMMI_Q3,  file = paste(save_file, bammi_filenameQ3, sep=''))
    save(ambarti,           file = paste(save_file, ambarti_filename, sep=''))

  }
}
