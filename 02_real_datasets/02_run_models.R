library(dplyr)
library(devtools)

install_github("ebprado/AMBARTI/R package")

library(AMBARTI)

# Set up where to save the results ----
save_path = "~/R/AMBARTI/04_real_data_sets/results/"

run_real_data_sets <- function(year){

  set.seed(001)

  filename = paste('Ireland_VCU_',
                   year,
                   sep='')

  load(paste(filename, '_data.RData',              sep=''))

  # run classical AMMI
  classical_AMMI_Q1 = run_classical_AMMI(data, Q=1)
  classical_AMMI_Q2 = run_classical_AMMI(data, Q=2)
  classical_AMMI_Q3 = run_classical_AMMI(data, Q=3)

  # run AMBARTI
  ambarti = run_AMBARTI(data, ntrees = 200, nburn = 2000, npost = 1000)

  # save files
  filename = paste('Ireland_VCU_',
                   year,
                   sep='')

  data_filename    = paste(filename, '_data.RData',              sep='')
  ammi_filenameQ1  = paste(filename, '_classical_AMMI_Q1.RData', sep='')
  ammi_filenameQ2  = paste(filename, '_classical_AMMI_Q2.RData', sep='')
  ammi_filenameQ3  = paste(filename, '_classical_AMMI_Q3.RData', sep='')
  ambarti_filename = paste(filename, '_AMBARTI.RData',           sep='')

  # save(data,              file = paste(save_path, data_filename,    sep=''))
  save(classical_AMMI_Q1, file = paste(save_path, ammi_filenameQ1,  sep=''))
  save(classical_AMMI_Q2, file = paste(save_path, ammi_filenameQ2,  sep=''))
  save(classical_AMMI_Q3, file = paste(save_path, ammi_filenameQ3,  sep=''))
  save(ambarti,           file = paste(save_path, ambarti_filename, sep=''))
}

all_years = unique(ireland$Year)

for (k in all_years){
  print(paste('Year: ', k, sep=''))
  run_real_data_sets(k)
}
