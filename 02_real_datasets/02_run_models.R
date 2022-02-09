library(dplyr)
library(devtools)
library(gss)
library(BASS)
library(earth)
install_github("ebprado/AMBARTI/R package",
               ref = 'main')

library(AMBARTI)
# Get real data sets ------ 
file_path = "~/R/AMBARTI/04_real_data_sets/01_data_sets/"
file_name = 'Historical_data_Ireland_VCU.csv'
ireland = read.csv(paste(file_path, file_name, sep=''), sep = ';', dec = ',')

# Set up where to save the results ---- 
save_path = "~/R/AMBARTI/04_real_data_sets/results/"

run_mars <- function(data){
  mars <- earth(data$y~., data = data$x, degree = 2)
  return(mars = mars)
}

run_bmars <- function(data, nmcmc = 10000, nburn = 9000, thin = 1){
  mod <- bass(data$x, data$y, nmcmc = nmcmc, nburn = nburn, thin = thin)
  return(mod)
}

run_gss <- function(data){
  data <- data.frame(y = data$y, g = data$x$g, e = data$x$e)
  ss <- ssanova(y~ g + e + g*e, data = data)
  return(ss)
}

run_real_data_sets <- function(year){
  
  set.seed(001)
    
  ireland_year = ireland %>% 
    filter(Year==year) %>% 
    mutate(g = factor(Genotype, levels = unique(Genotype), labels=paste('g', 1:length(unique(Genotype)), sep='')),
           e = factor(Location, levels = unique(Location), labels=paste('e', 1:length(unique(Location)), sep='')),
           y = yld_ton_ha) %>% 
    dplyr::select(g, Genotype, e, Location, y)
  
  ireland_year = ireland_year %>% 
    group_by(g,e) %>% 
    mutate(id = 1:n())
    
  ireland_year = ireland_year %>% 
    filter(id %in% c(1,3)) %>% 
    group_by(g, e) %>% 
    summarise(y = mean(y)) %>% 
    dplyr::select(g,e,y) %>% 
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
  bayesian_AMMI_Q1 = run_bayesian_AMMI(data,Q=1)
  bayesian_AMMI_Q2 = run_bayesian_AMMI(data,Q=2)
  bayesian_AMMI_Q3 = run_bayesian_AMMI(data,Q=3)
  
  # run GSS
  gss <- run_gss(data)
  
  # run BASS
  bmars <- run_bmars(data)
  
  # run MARS
  mars <- run_mars(data)
  
  # run AMBARTI
  ambarti = run_AMBARTI(data, ntrees = 200, nburn = 2000, npost = 1000)

  # save files
  filename = paste('Ireland_VCU_',
                   year,
                   sep='')
  
  data_filename     = paste(filename, '_data.RData',              sep='')
  ammi_filenameQ1   = paste(filename, '_classical_AMMI_Q1.RData', sep='')
  ammi_filenameQ2   = paste(filename, '_classical_AMMI_Q2.RData', sep='')
  ammi_filenameQ3   = paste(filename, '_classical_AMMI_Q3.RData', sep='')
  bammi_filenameQ1  = paste(filename, '_bayesian_AMMI_Q1.RData',  sep='')
  bammi_filenameQ2  = paste(filename, '_bayesian_AMMI_Q2.RData',  sep='')
  bammi_filenameQ3  = paste(filename, '_bayesian_AMMI_Q3.RData',  sep='')
  gss_filename      = paste(filename, '_gss.RData',               sep='')
  bmars_filename    = paste(filename, '_bmars.RData',             sep='')
  mars_filename     = paste(filename, '_mars.RData',              sep='')
  ambarti_filename  = paste(filename, '_AMBARTI.RData',           sep='')
  
  save(data,              file = paste(save_path, data_filename,    sep=''))
  save(classical_AMMI_Q1, file = paste(save_path, ammi_filenameQ1,  sep=''))
  save(classical_AMMI_Q2, file = paste(save_path, ammi_filenameQ2,  sep=''))
  save(classical_AMMI_Q3, file = paste(save_path, ammi_filenameQ3,  sep=''))
  save(bayesian_AMMI_Q1, file = paste(save_path, bammi_filenameQ1,  sep=''))
  save(bayesian_AMMI_Q2, file = paste(save_path, bammi_filenameQ2,  sep=''))
  save(bayesian_AMMI_Q3, file = paste(save_path, bammi_filenameQ3,  sep=''))
  save(gss,              file = paste(save_path, gss_filename,      sep=''))
  save(bmars,            file = paste(save_path, bmars_filename,    sep=''))
  save(mars,             file = paste(save_path, mars_filename,     sep=''))
  save(ambarti,          file = paste(save_path, ambarti_filename,  sep=''))
}

all_years = unique(ireland$Year)

for (k in all_years){
  print(paste('Year: ', k, sep=''))
  run_real_data_sets(k)
}
