#######################################################################################
### Title:        4_npv_calculation.R
### Purpose:      This file will read all output files from the Global AQ Reduced Form
##                Tool and calculates the net present value per ton of GHG emissions
##                across all trials. (for each model & country & pollutant)
##                The calculations are done for each model and country and combined at
##                the end to help reduce the memory usage when processing all 10,000 RFF scenarios
### Written by:   US EPA, Climate Change Division (OAP) 
### Date Created: 3/15/2023      
### Last Edited:  8/29/2025 - Melanie Jackson, IEc, E. McDuffie (OAP)
##      Inputs:                                                                       ##
##        -damages_mean_X_rft.parquet                                                 ##
##        -Final Country Grid Col_Row Index.csv                                       ##
##     Outputs:
##        -npd_full_streams_rff_X_COUNTRY.parquet                                     ##
##        - npd_country_rff_means_COUNTRY.parquet                                     ##
##        - npd_global_rff_means.csv                                                  ##
## Units: All $ values for $2020                                                      ##
########################################################################################

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

read_rft = 
  function(x){
    ttemp <- 
    read_parquet(x) %>%
      select(!contains("Deaths")) %>% #subset for monetary results only
      rename_with(~str_replace(.x,pattern="Annual_impacts",replacement = "damages")) %>%
      group_by_at(c('year','LocID','Pollutant', 'Model','trial')) #%>% # sum across all countries
      pivot_wider(id_cols = c(year,LocID,Pollutant,Model,trial,pop,gdp,Country,Region),
                  names_from  = damageType, 
                  values_from = c("damages", "damages_2_5", "damages_97_5", "damages_wlag","damages_wlag_2_5", "damages_wlag_97_5",'Temperature'),
                  values_fill = NA) %>%
      #group_by(year,LocID, Pollutant, Model, trial) %>%
      summarise(damages.baseline = sum(damages_Baseline),
                damages_wlag.baseline = sum(damages_wlag_Baseline),
                damages_2_5.baseline = sum(damages_2_5_Baseline),
                damages_wlag_2_5.baseline = sum(damages_wlag_2_5_Baseline),
                damages_97_5.baseline = sum(damages_97_5_Baseline),
                damages_wlag_97_5.baseline = sum(damages_wlag_97_5_Baseline),
                damages.perturbed = sum(damages_Perturbed),
                damages_wlag.perturbed = sum(damages_wlag_Perturbed),
                damages_2_5.perturbed = sum(damages_2_5_Perturbed),
                damages_wlag_2_5.perturbed = sum(damages_wlag_2_5_Perturbed),
                damages_97_5.perturbed = sum(damages_97_5_Perturbed),
                damages_wlag_97_5.perturbed = sum(damages_wlag_97_5_Perturbed),
                pop = sum(pop),
                gdp = sum(gdp), 
                temp_perturbed = mean(Temperature_Perturbed),
                temp_baseline = mean(Temperature_Baseline),
                .groups='keep') %>%
           ungroup() #%>%
  }

read_rff_rft = 
  function(x){
    ttemp <- 
      read_parquet(x) %>%
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


if (allRFF==1){
    # To read all files:
    NAME = '1_500'
    damages = 
      list.files(Inputs.RFT, pattern = "damages_\\d+\\_rft.*", full.names = T)[1:500] %>% 
      map_df(~read_rff_rft(.))
} else {
    # to read specific file:
    damages =
      read_rft(file.path(Inputs.RFT,paste0('damages_mean_rft.parquet')))
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
} else {
  damages      <- damages %>%
    group_by(trial,LocID,Pollutant,Model) %>%
    mutate(across(all_of(vars.damages),  ~ 1 - (1/(1+(./gdp))), .names="{.col}_pct"), # fraction of GDP that are impacts
           across(all_of(vars.damages),  ~ (get(paste0(cur_column(),"_pct")) * gdp), .names="{.col}_cor"),   #calculate the corrected level of damages
           damages.marginal              = (damages.perturbed_cor-damages.baseline_cor) * (12/44) * (1e-9), #convert FaIR's GtC to tCO2
           damages_wlag.marginal         = (damages_wlag.perturbed_cor-damages_wlag.baseline_cor) * (12/44) * (1e-9),
           damages_2_5.marginal          = (damages_2_5.perturbed_cor-damages_2_5.baseline_cor) * (12/44) * (1e-9),
           damages_wlag_2_5.marginal     = (damages_wlag_2_5.perturbed_cor-damages_wlag_2_5.baseline_cor) * (12/44) * (1e-9),
           damages_97_5.marginal         = (damages_97_5.perturbed_cor-damages_97_5.baseline_cor) * (12/44) * (1e-9),
           damages_wlag_97_5.marginal    = (damages_wlag_97_5.perturbed_cor-damages_wlag_97_5.baseline_cor) * (12/44) * (1e-9),
           ypc                           = ((1-damages.baseline_pct) * gdp)/pop,
           base.ypc                      = case_when(year==PULSEYEAR ~ypc, T ~0),
           base.ypc                      = max(base.ypc),
           ypc_wlag                      = ((1-damages_wlag.baseline_pct) * gdp)/pop,
           base.ypc_wlag                 = case_when(year==PULSEYEAR ~ypc_wlag, T ~0),
           base.ypc_wlag                 = max(base.ypc_wlag),
           #these are the 95% CRF confidence intervals (with and without lag)
           ypc_2.5                       = ((1-damages_2_5.baseline_pct) * gdp)/pop,
           base_2.5.ypc                  = case_when(year==PULSEYEAR ~ypc_2.5, T ~0),
           base_2.5.ypc                  = max(base_2.5.ypc),
           ypc_2.5_wlag                  = ((1-damages_wlag_2_5.baseline_pct) * gdp)/pop,
           base_2.5.ypc_wlag             = case_when(year==PULSEYEAR ~ypc_2.5_wlag, T ~0),
           base_2.5.ypc_wlag             = max(base_2.5.ypc_wlag),
           ypc_97.5                       = ((1-damages_97_5.baseline_pct) * gdp)/pop,
           base_97.5.ypc                  = case_when(year==PULSEYEAR ~ypc_97.5, T ~0),
           base_97.5.ypc                  = max(base_97.5.ypc),
           ypc_97.5_wlag                  = ((1-damages_wlag_97_5.baseline_pct) * gdp)/pop,
           base_97.5.ypc_wlag             = case_when(year==PULSEYEAR ~ypc_97.5_wlag, T ~0),
           base_97.5.ypc_wlag             = max(base_97.5.ypc_wlag)
             ) %>%
    select(-c(contains('pct'), contains('cor')))
}
  
###SET
NAME = '501_1000'
damages <- read_parquet(file.path(Outputs,'damages_501_1000.parquet'))

  ## discount rates
  ## Using stochastic Ramsey discount rate calculation, following Rennert et al., 2022
rates = tibble(rate = c('1.5% Ramsey', '2.0% Ramsey', '2.5% Ramsey', '3.0% Ramsey', '7.0% Ramsey', '2.0% CDR', '3.0% CDR','7.0% CDR'),
               rho  = c(exp(0.000091496)-1, exp(0.001972641)-1, exp(0.004618785)-1, exp(0.007702710)-1, exp(0.040238690)-1, 0.02, 0.03,0.07), ## under discrete time, need to transform the rho that was calibrated using continuous time 
               #rho = c(0.000091496, 0.001972641, 0.004618785, 0.007702711, 0.02, 0.03),
               eta  = c(1.016010261, 1.244459020, 1.421158057, 1.567899403, 2.167368813, 0, 0, 0))

  
for (COUNTRY in 1:nrow(countries)) {
    
  ## object to store data
  data = tibble()
  # object to store data
  means = tibble()
  net_means = tibble()
    
  print(COUNTRY)
    
  for (RATE in 1:length(rates$rate)){
  
    ## get damage parameters
    rate = rates$rate[[RATE]]
    rho  = rates$rho[[RATE]]
    eta  = rates$eta[[RATE]]
  
    ## code calculates a stochastic Ramsey discount factor due to the discrete time nature of the results
    ## More info in Rennert et al., 2022
    ## get streams of discounted damages and net present damages
    if (allRFF==1){
      data = 
        bind_rows(
          data,
          damages %>%
            filter(year >= PULSEYEAR) %>%
            filter(LocID == countries$LocID[COUNTRY]) %>%
            group_by(trial,LocID,Pollutant,Model) %>% 
            mutate(discount.rate                    = rate,
                   discount.factor                  = case_when(grepl("Ramsey",rate) ~ (base.ypc/ypc)^eta/(1+rho)^(year-PULSEYEAR),
                                                                T ~ 1/(1+rho)^(year-PULSEYEAR)),
                   damages.marginal.discounted          = damages.marginal * discount.factor,
                   npd                                  = sum(damages.marginal.discounted, na.rm=F),
                   cert.eq.adj                          = case_when(grepl("Ramsey", rate) ~ (base.ypc^-eta)/mean(base.ypc^-eta, na.rm = T),
                                                                    T ~ 1)
            ) %>%
            ungroup()
        )
    } else {
      data = 
        bind_rows(
        data,
        damages %>%
          filter(year >= PULSEYEAR) %>%
          filter(LocID == countries$LocID[COUNTRY]) %>%
          group_by(trial,LocID,Pollutant,Model) %>% 
          mutate(discount.rate                    = rate,
                 discount.factor                  = case_when(grepl("Ramsey",rate) ~ (base.ypc/ypc)^eta/(1+rho)^(year-PULSEYEAR),
                                                             T ~ 1/(1+rho)^(year-PULSEYEAR)),
                 discount.factor_wlag             = case_when(grepl("Ramsey",rate) ~ (base.ypc_wlag/ypc_wlag)^eta/(1+rho)^(year-PULSEYEAR),
                                                             T ~ 1/(1+rho)^(year-PULSEYEAR)),
                 discount.factor_2.5                  = case_when(grepl("Ramsey",rate) ~ (base_2.5.ypc/ypc_2.5)^eta/(1+rho)^(year-PULSEYEAR),
                                                              T ~ 1/(1+rho)^(year-PULSEYEAR)),
                 discount.factor_2.5_wlag             = case_when(grepl("Ramsey",rate) ~ (base_2.5.ypc_wlag/ypc_2.5_wlag)^eta/(1+rho)^(year-PULSEYEAR),
                                                              T ~ 1/(1+rho)^(year-PULSEYEAR)),
                 discount.factor_97.5                  = case_when(grepl("Ramsey",rate) ~ (base_97.5.ypc/ypc_97.5)^eta/(1+rho)^(year-PULSEYEAR),
                                                              T ~ 1/(1+rho)^(year-PULSEYEAR)),
                 discount.factor_97.5_wlag             = case_when(grepl("Ramsey",rate) ~ (base_97.5.ypc_wlag/ypc_97.5_wlag)^eta/(1+rho)^(year-PULSEYEAR),
                                                              T ~ 1/(1+rho)^(year-PULSEYEAR)),
                 damages.marginal.discounted          = damages.marginal * discount.factor,
                 npd                                  = sum(damages.marginal.discounted, na.rm=F),
                 damages.marginal_wlag.discounted     = damages_wlag.marginal * discount.factor_wlag,
                 npd_wlag                             = sum(damages.marginal_wlag.discounted, na.rm = F),
                 damages_2.5.marginal.discounted      = damages_2_5.marginal * discount.factor_2.5,
                 npd_crf_2.5                          = sum(damages_2.5.marginal.discounted, na.rm=F),
                 damages_2.5.marginal_wlag.discounted = damages_wlag_2_5.marginal * discount.factor_2.5_wlag,
                 npd_crf_2.5_wlag                     = sum(damages_2.5.marginal_wlag.discounted, na.rm = F),
                 damages_97.5.marginal.discounted     = damages_97_5.marginal * discount.factor_97.5,
                 npd_crf_97.5                         = sum(damages_97.5.marginal.discounted, na.rm=F),
                 damages_97.5.marginal_wlag.discounted = damages_wlag_97_5.marginal * discount.factor_97.5_wlag,
                 npd_crf_97.5_wlag                    = sum(damages_97.5.marginal_wlag.discounted, na.rm = F),
                 cert.eq.adj                          = case_when(grepl("Ramsey", rate) ~ (base.ypc^-eta)/mean(base.ypc^-eta, na.rm = T),
                                                              T ~ 1),
                 cert.eq.adj.wlag                     = case_when(grepl("Ramsey", rate) ~ (base.ypc_wlag^-eta)/mean(base.ypc_wlag^-eta, na.rm = T),
                                                              T ~ 1),
                 cert.eq.adj_2.5                      = case_when(grepl("Ramsey", rate) ~ (base_2.5.ypc^-eta)/mean(base_2.5.ypc^-eta, na.rm = T),
                                                                  T ~ 1),
                 cert.eq.adj_2.5.wlag                 = case_when(grepl("Ramsey", rate) ~ (base_2.5.ypc_wlag^-eta)/mean(base_2.5.ypc_wlag^-eta, na.rm = T),
                                                                  T ~ 1),
                 cert.eq.adj_97.5                     = case_when(grepl("Ramsey", rate) ~ (base_97.5.ypc^-eta)/mean(base_97.5.ypc^-eta, na.rm = T),
                                                                      T ~ 1),
                 cert.eq.adj_97.5.wlag                = case_when(grepl("Ramsey", rate) ~ (base_97.5.ypc_wlag^-eta)/mean(base_97.5.ypc_wlag^-eta, na.rm = T),
                                                                      T ~ 1),
                 ) %>%
          ungroup()
    )
    }
  }
  
  ## export full streams
  if (allRFF ==1){
    data %>%
      write_parquet(file.path(Outputs, paste0('npd_full_streams_rff_',countries$LocID[COUNTRY],'.parquet')))
  } else {
    data %>%
     write_parquet(file.path(Outputs, paste0('npd_full_streams_',countries$LocID[COUNTRY],'.parquet')))
  }


  # recover summary statistics across all trials
  # Note that the stats and certainty equivalent corrections are only relevant when model is run for all trials (not the mean RFF-SP)
  if (allRFF ==1){
    means =
      bind_rows(
        means,
        data %>%
          filter(year == PULSEYEAR) %>%
          group_by(LocID,discount.rate,Pollutant,Model) %>% #calculate statistics for each country, discount rate, pollutant and model
          summarise(mean_npd         = mean(npd),
                    mean.npd.cert.eq = mean(npd * cert.eq.adj, na.rm = T),
                    npd_2.5          = quantile(npd, .025, na.rm = T), #these uncertainties are socioeconomic uncertainties (not BenMAP) - but only when results are run for all RFF-SPs (not just the mean)
                    npd_97.5         = quantile(npd, .975, na.rm = T),
                    median           = median(npd, na.rm = T),
                    .groups = 'drop'))
    net_means = 
      bind_rows(
        net_means,
        data %>%
          filter(year == PULSEYEAR) %>% #npds are the same for each year
          group_by(LocID,discount.rate,trial,Model) %>% #sum across pollutants
          summarise(net_npd         = sum(npd),
                    cert.eq.adj     = mean(cert.eq.adj),
                    .groups = 'drop') %>%
          ungroup() %>%
          group_by(LocID,discount.rate,Model) %>% #calculate statistics for each country, discount rate, and model
          summarise(mean_npd         = mean(net_npd),
                    mean.npd.cert.eq = mean(net_npd * cert.eq.adj, na.rm = T),
                    npd_2.5          = quantile(net_npd, .025, na.rm = T), #these uncertainties are socioeconomic uncertainties (not BenMAP) - but only when results are run for all RFF-SPs (not just the mean)
                    npd_97.5         = quantile(net_npd, .975, na.rm = T),
                    median           = median(net_npd, na.rm = T),
                    .groups = 'drop')
        )
  } else {
  means =
    bind_rows(
      means,
      data %>%
        group_by(LocID,discount.rate,Pollutant,Model) %>% #calculate statistics for each country, discount rate, pollutant and model
        summarise(mean_npd         = mean(npd),
                  mean.npd.cert.eq = mean(npd * cert.eq.adj, na.rm = T),
                  npd_2.5          = quantile(npd, .025, na.rm = T), #these uncertainties are socioeconomic uncertainties (not BenMAP) - but only when results are run for all RFF-SPs (not just the mean)
                  npd_97.5         = quantile(npd, .975, na.rm = T),
                  median           = median(npd, na.rm = T),
                  mean_npd_wlag    = mean(npd_wlag),
                  mean.npd.wlag.cert.eq = mean(npd_wlag * cert.eq.adj.wlag, na.rm = T),
                  npd_2.5_wlag     = quantile(npd_wlag, .025, na.rm = T),
                  npd_97.5_wlag    = quantile(npd_wlag, .975, na.rm = T),
                  median_wlag      = median(npd_wlag, na.rm = T),
                  #these are the CRF 95% confidence intervals (with and without lag)
                  mean_npd_crf_2.5         = mean(npd_crf_2.5),
                  mean.npd.cert.eq_crf_2.5 = mean(npd_crf_2.5 * cert.eq.adj_2.5, na.rm = T),
                  npd_2.5_crf_2.5          = quantile(npd_crf_2.5, .025, na.rm = T), #these uncertainties are socioeconomic uncertainties (not BenMAP) - but only when results are run for all RFF-SPs (not just the mean)
                  npd_97.5_crf_2.5         = quantile(npd_crf_2.5, .975, na.rm = T),
                  median_crf_2.5           = median(npd_crf_2.5, na.rm = T),
                  mean_npd_wlag_crf_2.5    = mean(npd_crf_2.5_wlag),
                  mean.npd.wlag.cert.eq_crf_2.5 = mean(npd_crf_2.5_wlag * cert.eq.adj_2.5.wlag, na.rm = T),
                  npd_2.5_wlag_crf_2.5     = quantile(npd_crf_2.5_wlag, .025, na.rm = T),
                  npd_97.5_wlag_crf_2.5    = quantile(npd_crf_2.5_wlag, .975, na.rm = T),
                  median_wlag_crf_2.5      = median(npd_crf_2.5_wlag, na.rm = T),
                  mean_npd_crf_97.5         = mean(npd_crf_97.5),
                  mean.npd.cert.eq_crf_97.5 = mean(npd_crf_97.5 * cert.eq.adj_97.5, na.rm = T),
                  npd_2.5_crf_97.5          = quantile(npd_crf_97.5, .025, na.rm = T), #these uncertainties are socioeconomic uncertainties (not BenMAP) - but only when results are run for all RFF-SPs (not just the mean)
                  npd_97.5_crf_97.5         = quantile(npd_crf_97.5, .975, na.rm = T),
                  median_crf_97.5           = median(npd_crf_97.5, na.rm = T),
                  mean_npd_wlag_crf_97.5    = mean(npd_crf_97.5_wlag),
                  mean.npd.wlag.cert.eq_crf_97.5 = mean(npd_crf_97.5_wlag * cert.eq.adj_97.5.wlag, na.rm = T),
                  npd_2.5_wlag_crf_97.5     = quantile(npd_crf_97.5_wlag, .025, na.rm = T),
                  npd_97.5_wlag_crf_97.5    = quantile(npd_crf_97.5_wlag, .975, na.rm = T),
                  median_wlag_crf_97.5      = median(npd_crf_97.5_wlag, na.rm = T),
                  .groups = 'drop'))

  }
  ## export summary stats
  if (allRFF ==1){
    means %>%
      write_csv(file.path(Outputs,paste0('npd_country_rff_means_',countries$LocID[COUNTRY],'.csv')))
    net_means %>%
      write_csv(file.path(Outputs,paste0('npd_country_rff_net_means_',countries$LocID[COUNTRY],'.csv')))
  } else {
    means %>%
      write_csv(file.path(Outputs,paste0('npd_country_means_',countries$LocID[COUNTRY],'.csv')))
  }

}


#Do global calculation

read_results = function(x){
    ttemp <- 
      read_parquet(x) %>%
      filter(year == PULSEYEAR) #%>% #all npd years are the same
}

read_net_results = function(x){
  ttemp <- 
    read_parquet(x) %>%
    filter(year == PULSEYEAR) %>% #all npd years are the same
    group_by(LocID,discount.rate,trial,Model) %>% #sum across pollutants
    summarise(net_npd         = sum(npd),
              cert.eq.adj     = mean(cert.eq.adj),
              .groups = 'drop') 
}


if (allRFF ==1){
  Results_comb = 
    list.files(Outputs, pattern = paste0("npd_full_streams_rff_"),full.names = T) %>% 
    map_df(~read_results(.))
  Results_comb_net = 
    list.files(Outputs, pattern = paste0("npd_full_streams_rff_"),full.names = T) %>% 
    map_df(~read_net_results(.))
} else {
  Results_comb = 
    list.files(Outputs, pattern = paste0("npd_full_streams_"),full.names = T) %>% 
    map_df(~read_results(.))
}

#  #write sum of country data
if (allRFF ==1){
  Results_comb %>%
    write_csv(file.path(Outputs,paste0('npd_all_country_rff_means_',NAME,'.csv')))
  Results_comb_net %>%
    write_csv(file.path(Outputs,paste0('npd_all_country_net_rff_means_',NAME,'.csv')))
} else {
  Results_comb %>%
    write_csv(file.path(Outputs,paste0('npd_all_country_means.csv')))
}


#export global stats
glob_means = tibble()
glob_net_means = tibble()

# Note that the stats and certainty equivalent corrections are only relevant when model is run for all trials (not the mean RFF-SP)
if (allRFF ==1){  
  glob_means =
    bind_rows(
      glob_means,
      Results_comb %>%
        group_by(discount.rate,trial,Pollutant,Model) %>% #group by rate, trial, pollutant, and model, and then sum across all countries
        summarise(
          npd.cert.eq = sum(npd * cert.eq.adj, na.rm = T),
          npd      = sum(npd),
          .groups = 'keep') %>%
        ungroup() %>%
        group_by(discount.rate,Pollutant,Model) %>% #next, group by discount rate, pollutant, and model to calculate stats across all trials
        summarise(mean_npd      = mean(npd),
                  mean.npd.cert.eq = mean(npd.cert.eq),
                  npd_2.5       = quantile(npd, .025, na.rm = T), #these are the socioeconomic stats (not BenMAP)
                  npd_97.5      = quantile(npd, .975, na.rm = T),
                  median        = median(npd, na.rm = T),
                  .groups = 'keep')
    )
  glob_net_means =
    bind_rows(
      glob_net_means,
      Results_comb_net %>%
        group_by(discount.rate,trial,Model) %>% #group by rate, trial, and model, and then sum across all countries
        summarise(
          npd.cert.eq = sum(net_npd * cert.eq.adj, na.rm = T),
          npd      = sum(net_npd),
          .groups = 'keep') %>%
        ungroup() %>%
        group_by(discount.rate,Model) %>% #next, group by discount rate and model to calculate stats across all trials
        summarise(mean_npd      = mean(npd),
                  mean.npd.cert.eq = mean(npd.cert.eq),
                  npd_2.5       = quantile(npd, .025, na.rm = T), #these are the socioeconomic stats (not BenMAP)
                  npd_97.5      = quantile(npd, .975, na.rm = T),
                  median        = median(npd, na.rm = T),
                  .groups = 'keep')
    )
  
  } else {
  glob_means =
      bind_rows(
        glob_means,
        Results_comb %>%
        group_by(discount.rate,trial,Pollutant,Model) %>% #group by rate, trial, pollutant, and model, and then sum across all countries
        summarise(
            npd.cert.eq = sum(npd * cert.eq.adj, na.rm = T),
            npd.wlag.cert.eq = sum(npd_wlag * cert.eq.adj.wlag, na.rm = T),
            npd_wlag = sum(npd_wlag),
            npd      = sum(npd),
            #these are crf uncertainties
            npd_crf_2.5.cert.eq = sum(npd_crf_2.5 * cert.eq.adj_2.5, na.rm = T),
            npd_crf_2.5.wlag.cert.eq = sum(npd_crf_2.5_wlag * cert.eq.adj_2.5.wlag, na.rm = T),
            npd_crf_2.5_wlag = sum(npd_crf_2.5_wlag),
            npd_crf_2.5      = sum(npd_crf_2.5),
            npd_crf_97.5.cert.eq = sum(npd_crf_97.5 * cert.eq.adj_97.5, na.rm = T),
            npd_crf_97.5.wlag.cert.eq = sum(npd_crf_97.5_wlag * cert.eq.adj_97.5.wlag, na.rm = T),
            npd_crf_97.5_wlag = sum(npd_crf_97.5_wlag),
            npd_crf_97.5      = sum(npd_crf_97.5),
            .groups = 'keep') %>%
        ungroup() %>%
        group_by(discount.rate,Pollutant,Model) %>% #next, group by discount rate, pollutant, and model to calculate stats across all trials
        summarise(mean_npd      = mean(npd),
                    mean.npd.cert.eq = mean(npd.cert.eq),
                    npd_2.5       = quantile(npd, .025, na.rm = T), #these are the socioeconomic stats (not BenMAP)
                    npd_97.5      = quantile(npd, .975, na.rm = T),
                    median        = median(npd, na.rm = T),
                    mean_npd_wlag = mean(npd_wlag),
                    mean.npd.wlag.cert.eq = mean(npd.wlag.cert.eq),
                    npd_2.5_wlag  = quantile(npd_wlag, .025, na.rm = T),
                    npd_97.5_wlag = quantile(npd_wlag, .975, na.rm = T),
                    median_wlag   = median(npd_wlag, na.rm = T),
                  #these are the crf uncertainty results (with and without lag)
                  mean_npd_crf_2.5      = mean(npd_crf_2.5),
                  mean.npd.cert.eq_crf_2.5 = mean(npd_crf_2.5.cert.eq),
                  npd_2.5_crf_2.5       = quantile(npd_crf_2.5, .025, na.rm = T), #these are the socioeconomic stats (not BenMAP) + BenMAP uncertainty
                  npd_97.5_crf_2.5      = quantile(npd_crf_2.5, .975, na.rm = T),
                  median_crf_2.5        = median(npd_crf_2.5, na.rm = T),
                  mean_npd_wlag_crf_2.5 = mean(npd_crf_2.5_wlag),
                  mean.npd.wlag.cert.eq_2.5 = mean(npd_crf_2.5.wlag.cert.eq),
                  npd_2.5_wlag_2.5  = quantile(npd_crf_2.5_wlag, .025, na.rm = T),
                  npd_97.5_wlag_2.5 = quantile(npd_crf_2.5_wlag, .975, na.rm = T),
                  median_wlag_2.5   = median(npd_crf_2.5_wlag, na.rm = T),
                  mean_npd_crf_97.5      = mean(npd_crf_97.5),
                  mean.npd.cert.eq_crf_97.5 = mean(npd_crf_97.5.cert.eq),
                  npd_2.5_crf_97.5       = quantile(npd_crf_97.5, .025, na.rm = T), #these are the socioeconomic stats (not BenMAP)
                  npd_97.5_crf_97.5      = quantile(npd_crf_97.5, .975, na.rm = T),
                  median_crf_97.5        = median(npd_crf_97.5, na.rm = T),
                  mean_npd_wlag_crf_97.5 = mean(npd_crf_97.5_wlag),
                  mean.npd.wlag.cert.eq_97.5 = mean(npd_crf_97.5.wlag.cert.eq),
                  npd_2.5_wlag_97.5  = quantile(npd_crf_97.5_wlag, .025, na.rm = T),
                  npd_97.5_wlag_97.5 = quantile(npd_crf_97.5_wlag, .975, na.rm = T),
                  median_wlag_97.5   = median(npd_crf_97.5_wlag, na.rm = T),
                    .groups = 'keep')
  )
  }
# #
#  #write data
if (allRFF ==1){
  glob_means %>%
      write_csv(file.path(Outputs,paste0('npd_global_rff_means_',NAME,'.csv')))
  glob_net_means %>%
    write_csv(file.path(Outputs,paste0('npd_global_rff_net_means_',NAME,'.csv')))
} else {
    glob_means %>%
      write_csv(file.path(Outputs,paste0('npd_global_means.csv')))
}

glob_avg_spc <- glob_means %>% 
  group_by(discount.rate,Pollutant) %>% 
  summarise(mean_npd = mean(mean_npd),
            mean.npd.cert.eq = mean(mean.npd.cert.eq),
            npd_2.5 = quantile(npd_2.5, .025, na.rm = T),
            npd_97.5 = quantile(npd_97.5, .975, na.rm = T)) 
glob_avg_net <- glob_means %>% 
  group_by(discount.rate,Model) %>% 
  summarise(mean_npd = sum(mean_npd),
            mean.npd.cert.eq = sum(mean.npd.cert.eq),
            npd_2.5 = quantile(npd_2.5, .025, na.rm = T),
            npd_97.5 = quantile(npd_97.5, .975, na.rm = T))%>% 
  ungroup %>% 
  group_by(discount.rate) %>% 
  summarise(mean_npd = mean(mean_npd),
            mean.npd.cert.eq = mean(mean.npd.cert.eq),
            npd_2.5 = mean(npd_2.5),npd_97.5 = mean(npd_97.5))
  #CODE END  
