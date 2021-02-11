library(dplyr)
library(devtools)

install_github("ebprado/AMBARTI/R package",
               ref = 'main',
               auth_token = '363d4ad84eaa25d3eb26752732cc208f7e698086')

library(AMBARTI)

file_path = "~/Documents/GitHub/AMBARTI/03_real_datasets/01_data_sets/"
file_name = 'Historical_data_Ireland_VCU.csv'

ireland = read.csv(paste(file_path, file_name, sep=''), sep = ';', dec = ',')

ireland_year = ireland %>% 
  filter(Year==2010) %>% 
  mutate(g = as.character(Genotype),
         e = as.character(Location),
         y = yld_ton_ha) %>% 
  select(g,e,y)

ireland_year = ireland_year %>% 
  group_by(g, e) %>% 
  summarise(y = mean(y)) %>% 
  data.frame()

data = list()
data[['x']] = ireland_year[,c('g','e')]
data[['y']] = ireland_year[,'y']
data[['I']] = length(unique(ireland_year[,'g']))
data[['J']] = length(unique(ireland_year[,'e']))
    
    # run classical AMMI
    classical_AMMI_Q1 = run_classical_AMMI(data, Q=1)
    classical_AMMI_Q2 = run_classical_AMMI(data, Q=2)
    classical_AMMI_Q3 = run_classical_AMMI(data, Q=3)
    
    # run Bayesian AMMI
    # bayesian_AMMI_Q1 = run_bayesian_AMMI(data,Q=1)
    # bayesian_AMMI_Q2 = run_bayesian_AMMI(data,Q=2)
    # bayesian_AMMI_Q3 = run_bayesian_AMMI(data,Q=3)
    
    # run AMBARTI
    ambarti = run_AMBARTI(data, ntrees = 200, nburn = 1000, npost = 1000)
    
    # Increment the seed number by 1
    nseed = nseed + 1
    
    print(paste('comb = ', i, ' out of ', n_comb , ' (rep = ', j, ')', sep=''))
    
    # save files
    filename = paste('AMBARTI',
                     year,
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
    # save(bayesian_AMMI_Q1,  file = paste(save_file, bammi_filenameQ1, sep=''))
    # save(bayesian_AMMI_Q2,  file = paste(save_file, bammi_filenameQ2, sep=''))
    # save(bayesian_AMMI_Q3,  file = paste(save_file, bammi_filenameQ3, sep=''))
    save(ambarti,           file = paste(save_file, ambarti_filename, sep=''))
