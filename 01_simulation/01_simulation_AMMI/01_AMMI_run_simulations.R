library(devtools)
install_github("ebprado/AMBARTI/R package",
               ref = 'main')

library(AMBARTI)

save_file = "~/R/AMBARTI/01_simulation_AMMI/results/"

I = c(10) # c(10, 25, 50) # Number of genotypes
J = c(10) # c(10, 25, 50) # Number of environments
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
nseed = 0 # start seed

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
      # nseed      = 0 # start seed
      
      # Set a seed to make it reproducible
      set.seed(nseed)
  
      data = generate_data_AMMI(I, J, s_alpha, s_beta, s_y, lambda)
  
      # run classical AMMI
      classical_AMMI = run_classical_AMMI(data)
  
      # run Bayesian AMMI
      bayesian_AMMI = run_bayesian_AMMI(data)
  
      # run AMBARTI
      ambarti = run_AMBARTI(data, ntrees = 200, nburn = 1000, npost = 1000)
  
      # Increment the seed number by 1
      nseed = nseed + 1
  
      print(paste('comb = ', i, ' out of ', n_comb , ' (rep = ', j, ')', sep=''))
  
      # save files
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
  
      save(data,           file = paste(save_file, data_filename,    sep=''))
      save(classical_AMMI, file = paste(save_file, ammi_filename,    sep=''))
      save(bayesian_AMMI,  file = paste(save_file, bammi_filename,   sep=''))
      save(ambarti,        file = paste(save_file, ambarti_filename, sep=''))

  }
}
