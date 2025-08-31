rm(list = ls()); gc()
#necessary packages

packs<-c("dplyr","tidyverse","readxl","purrr","foreach","utils","stringr")
for (package in packs) { #Installs packages if not yet installed
  if (!requireNamespace(package, quietly = TRUE))
    install.packages(package)
}

library(dplyr) #for general data manipulation
library(tidyverse) #for general data manipulation
library(readxl) #to read in xlsx files
library(foreach) #to do parallel processing in loop
library(purrr)  #for general data manipulation
library(utils)
library(arrow)
library(stringr)
#setwd("~/shared/OAR/OAP/CCD/CSIB/Methane-Ozone/MOMM-RFT")

#path    <- file.path("//iec.local/jen/CIRA_Health_Analyses/Global PM Analysis/11_RFT/3_Final code")
Inputs  <- file.path("input")
Inputs.RFT <- file.path("output","rft")
Outputs <- file.path("output","npd")

allRFF = '1' #(set to 1 if reading in all RFF scenario, otherwise read in mean)
PULSEYEAR = 2030

countries <- read.csv(file.path(Inputs,"Final Country Grid Col_Row Index.csv"))[,c(1,3:4,6,8)]
countries <- countries %>% filter(Region != "")


read_rff_rft = 
  function(x){
    ttemp <- 
      read_parquet(x) %>%
      filter(year >= PULSEYEAR) %>%
      select(!contains("Deaths")) %>% #subset for monetary results only
      rename_with(~str_replace(.x,pattern="Annual_impacts",replacement = "damages")) %>%
      group_by_at(c('year','LocID','Pollutant', 'Model','trial')) %>% # sum across all countries
      pivot_wider(id_cols = c(year,LocID,Pollutant,Model,trial,pop,gdp,Country,Region),
                  names_from  = damageType, 
                  values_from = c("damages",'Temperature'),
                  values_fill = NA) %>%
      summarise(damages.baseline = sum(damages_Baseline),
                damages.perturbed = sum(damages_Perturbed),
                pop = sum(pop),
                gdp = sum(gdp), 
                temp_perturbed = mean(Temperature_Perturbed),
                temp_baseline = mean(Temperature_Baseline),
                .groups='keep') %>%
      drop_na() %>%
      ungroup() #%>%
  }

names_list = c('1_2','501_1000','1001_2000','2001_3000','3001_4000','4001_5000','5001_6000','6001_7000','7001_8000','8001_9000','9001_10000')
file_list = list(1:2,501:1000,1001:2000,2001:3000,3001:4000,4001:5000,5001:6000,6001:7000,7001:8000,8001:9000,9001:10000)

for (i in 1:length(names_list)) {
if (allRFF==1){
  # To read all files:
  print(paste0('Reading in files: ',list.files(Inputs.RFT, pattern = "damages_\\d+\\_rft.*", full.names = T)[unlist(file_list[i])]))
  damages = 
    list.files(Inputs.RFT, pattern = "damages_\\d+\\_rft.*", full.names = T)[unlist(file_list[i])] %>% 
    map_df(~read_rff_rft(.))
}

countries <- damages %>% ungroup() %>% distinct(LocID)

#damages are calculated as a function of GDP (growth is corrected for damages)
# calculate percent of GDP that are damages
vars.damages <- names(damages)[grep("damages",names(damages))]

if (allRFF==1){
  damages      <- damages %>%
    group_by(trial,LocID,Pollutant,Model) %>%
    mutate(across(all_of(vars.damages),  ~ 1 - (1/(1+(./gdp))), .names="{.col}_pct"), # fraction of GDP that are impacts
           across(all_of(vars.damages),  ~ (get(paste0(cur_column(),"_pct")) * gdp), .names="{.col}_cor"),   #calculate the corrected level of damages
           damages.marginal              = (damages.perturbed_cor-damages.baseline_cor) * (12/44) * (1e-9), #convert FaIR's GtC to tCO2
           ypc                           = ((1-damages.baseline_pct) * gdp)/pop,
           base.ypc                      = case_when(year==PULSEYEAR ~ypc, T ~0),
           base.ypc                      = max(base.ypc)
    ) %>%
    select(-c(contains('pct'), contains('cor')))
}

  damages %>% 
    write_parquet(file.path(Outputs, paste0('damages_',names_list[i],'.parquet')))
  print(paste0('Done ',i,' of ',length(file_list)))
}