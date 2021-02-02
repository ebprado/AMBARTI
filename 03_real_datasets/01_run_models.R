library(ggplot2)
library(tidyverse)

file_path = "~/Documents/GitHub/AMBARTI/03_real_datasets/01_data_sets/"
file_name = 'Historical_data_Ireland_VCU.csv'

ireland = read.csv(paste(file_path, file_name, sep=''), sep = ';', dec = ',')
ireland$Year = as.factor(ireland$Year)

# g = Genotype
# e = Location
# y = yld_ton_ha

# Some descriptive stats

aux = ireland %>% select(Plot_ID,
                         Genotype,
                         Bloc,
                         Location,
                         Year,
                         yld_ton_ha)

tab1 = aux %>%
  group_by(Year,
           Location,
           Genotype) %>%
  summarise(mean_y=mean(yld_ton_ha))

write.csv(tab1, file='test.csv')

tab2 = aux %>%
  group_by(Year) %>%
  summarise(n=n())

tab3 = aux %>%
  group_by(Location) %>%
  summarise(n=n())

tab4 = aux %>%
  group_by(Year,
           Location) %>%
  summarise(n=n())

tab5 = aux %>%
  group_by(Year,
           Genotype) %>%
  summarise(n=n())

## Filter some sites and genotypes to ease the visualisation
data = ireland %>% filter(Year %in% c(2010, 2011, 2012),
                          Location %in% c('BN', 'CK', 'CKC'),
                          Genotype %in% c('Alchemy', 'Alves', 'Avatar', 'Cordiale'))

## plot 1
data %>%
  ggplot(aes(x = Location , y = yld_ton_ha, colour=Genotype)) +
  geom_boxplot() +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        legend.position = 'bottom',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_grid(~ Year)

# plot 2
data %>%
  ggplot(aes(x = Year , y = yld_ton_ha, colour=Genotype)) +
  geom_boxplot() +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        legend.position = 'bottom',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# plot 3
data %>%
  ggplot(aes(x = Year , y = yld_ton_ha, colour=Location)) +
  geom_boxplot() +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        legend.position = 'bottom',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
