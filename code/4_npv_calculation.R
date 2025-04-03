#######################################################################################
### Title:        4_npv_calculation.R
### Purpose:      This file will read all output files from the Global AQ Reduced Form
##                Tool and calculates the net present value per ton of GHG emissions
##                across all trials. (for each model & country)
##                The calculations are done for each model and country and combined at
##                the end to help reduce the memory usage when processing all 10,000 RFF scenarios
### Written by:   US EPA, Climate Change Division (OAP) 
### Date Created: 3/15/2023      
### Last Edited:  4/4/2025 - Melanie Jackson, IEc, E. McDuffie (OAP)
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
Inputs  <- file.path(path,"inputs")
Inputs.RFT <- file.path(path,"outputs","rft")
Outputs <- file.path(path,"outputs","npd")

allRFF = 'mean' #(set to 1 if reading in all RFF scenario, otherwise read in mean)

countries <- read.csv(file.path(Inputs,"Final Country Grid Col_Row Index.csv"))[,c(1,3:4,6,8)]
countries <- countries %>% filter(Region != "")

read_rft = 
  function(x){
    ttemp <- 
    read_parquet(x) %>%
      select(!contains("Deaths")) %>% #subset for monetary results only
      rename_with(~str_replace(.x,pattern="Annual_impacts",replacement = "Damages")) %>%
      ungroup()
      #mutate(damages = annual_impacts,
      #     damages_wlag = annual_impacts_wlag) %>%
      #group_by(Year,LocID,trial) %>% # sum across all countries
      #summarise(damages = sum(damages),
      #          damages_wlag = sum(damages_wlag),
      #          pop = sum(pop),
      #          gdp = sum(gdp),
      #          .groups='keep') %>%
      #ungroup()
  }


  if (allRFF==1){
    # To read all files:
    damages = 
      list.files(Inputs.RFT, pattern = "damages_\\d+\\_*", full.names = T) %>% 
      map_df(~read_rft(.))
  } else {
    # to read specific file:
    damages =
      read_rft(file.path(Inputs.RFT,paste0('damages_mean_rft.parquet')))
  }
  
  countries <- damages %>% distinct(LocID)

  #damages are calculated as a function of GDP (growth is corrected for damages)
  # calculate percent of GDP that are damages
  vars.damages <- names(damages)[grep("Damages",names(damages))]
  damages      <- damages %>%
    group_by(trial,LocID,Pollutant,Model) %>%
    mutate(across(all_of(vars.damages), ~ 1 - 1/(1+(./gdp)), .names="{.col}_pct"), # calculate gdp per capita, with gdp reduced by % of damages
           across(all_of(vars.damages), ~ ((1 - get(paste0(cur_column(),"_pct"))) * gdp)/pop, .names="{.col}_ypc"),   
           across(all_of(vars.damages), ~ case_when(Year == 2020 ~ get(paste0(cur_column(),"_ypc")) , T~0), .names="{.col}_ypcBase"), # calculate base year gdp/cap
           across(all_of(vars.damages), ~ get(paste0(cur_column(),"_pct")) * gdp, .names="{.col}_marg"))   #calculate the damages as a % of GDP
  vars.damages.ypcBase <- names(damages)[grep("_ypcBase",names(damages))]  
  damages      <- damages %>%
    mutate(across(all_of(vars.damages.ypcBase), ~ max(.))) %>% #clunky way to ensure that the base ypc for discrete time discount factor can be used in the discounting function, fix later
    select(-c(contains('pct'))) %>%  # remove pct columns, is intermediary calc
    ungroup()

  ## discount rates
  ## Using stochastic Ramsey discount rate calculation, following Rennert et al., 2022
  rates = tibble(rate = c('1.5% Ramsey', '2.0% Ramsey', '2.5% Ramsey', '3.0% Ramsey', '2.0% CDR', '3.0% CDR'),
               rho  = c(exp(0.000091496)-1, exp(0.001972641)-1, exp(0.004618785)-1, exp(0.007702711)-1, 0.02, 0.03), ## under discrete time, need to transform the rho that was calibrated using continuous time 
               #rho = c(0.000091496, 0.001972641, 0.004618785, 0.007702711, 0.02, 0.03),
               eta  = c(1.016010261, 1.244459020, 1.421158057, 1.567899395, 0, 0))

  
  for (COUNTRY in 1:nrow(countries)) {
    
    ## object to store data
    data = tibble()
    ## object to store data
    means = tibble()
    
    print(COUNTRY)
    
    for (RATE in 1:length(rates$rate)){
  
    ## get damage parameters
    rate = rates$rate[[RATE]]
    rho  = rates$rho[[RATE]]
    eta  = rates$eta[[RATE]]
  
    ## code calculates a stochastic Ramsey discount factor due to the discrete time nature of the results
    ## More info in Rennert et al., 2022
    ## get streams of discounted damages and net present damages
    data = 
      bind_rows(
        data,
        damages %>%
          filter(LocID == countries$LocID[COUNTRY]) %>%
          group_by(trial,LocID,Pollutant,Model) %>% #LocID
          mutate(discount.rate                   = rate) %>%
          mutate(
            across(all_of(vars.damages), ~ case_when(grepl("Ramsey", rate) ~ 
                                                     (get(paste0(cur_column(),"_ypcBase"))/get(paste0(cur_column(),"_ypc")))^eta/(1+rho)^(Year-2020),
                                                     T ~ 1/(1+rho)^(Year-2020)),.names = "{.col}_DiscFact"),
            across(all_of(vars.damages), ~ get(paste0(cur_column(),"_marg")) * get(paste0(cur_column(),"_DiscFact")),.names="{.col}_discounted"),
            across(all_of(vars.damages), ~ sum(get(paste0(cur_column(),"_discounted")), na.rm=F),.names="{.col}_npd"),
            across(all_of(vars.damages), ~  case_when(grepl("Ramsey", rate) ~
                                                      (get(paste0(cur_column(),"_ypcBase"))^-eta)/mean(get(paste0(cur_column(),"_ypcBase"))^-eta,na.rm=T),
                                                       T ~ 1),.names = "{.col}_cert.eq.adj")) %>%
        select(-contains(c("ypc","DiscFact"))) %>%
        ungroup()
    )
  

  ## export full streams
 if (allRFF ==1){
  data %>%
    write_parquet(file.path(Outputs, paste0('npd_full_streams_rff_vsl10_',countries$LocID[COUNTRY],'.parquet')))
 } else {
   data %>%
     write_parquet(file.path(Outputs, paste0('npd_full_streams_vsl10_',countries$LocID[COUNTRY],'.parquet')))
 }
}
  # recover summary statistics across all trials
  # Note that the stats and certainty equivalent corrections are only relevant when model is run for all trials (not the mean RFF-SP)
  means =
    bind_rows(
      means,
      data %>%
        group_by(LocID,discount.rate,Pollutant,Model) %>% #calculate statistics for each country, discount rate, pollutant and model
        summarise(
          across(all_of(vars.damages),~mean(get(paste0(cur_column(),"_npd"))), .names="{.col}_mean.NPD"),
          across(all_of(vars.damages),~mean(get(paste0(cur_column(),"_npd")) * get(paste0(cur_column(),"_cert.eq.adj"))), .names="{.col}_mean.NPD.cert.eq"),
          across(all_of(vars.damages),~quantile(get(paste0(cur_column(),"_npd")), 0.025, na.rm=T), .names="{.col}_NPD.2.5"),
          across(all_of(vars.damages),~quantile(get(paste0(cur_column(),"_npd")), 0.975, na.rm=T), .names="{.col}_NPD.97.5"),
          across(all_of(vars.damages),~median(get(paste0(cur_column(),"_npd"))), .names="{.col}_median.NPD"),
      .groups = 'drop')
  )


  ## export summary stats
  if (allRFF ==1){
    means %>%
        write_csv(file.path(Outputs,paste0('npd_country_rff_means_vsl10_',countries$LocID[COUNTRY],'.csv')))
  } else {
    means %>%
      write_csv(file.path(Outputs,paste0('npd_country_means_vsl10_',countries$LocID[COUNTRY],'.csv')))
  }

}


#Do global calculation

read_results = 
  function(x){
    ttemp <- 
      read_parquet(x) %>%
      filter(Year == 2020) #%>% #all npd years are the same
      #select(LocID,trial,discount.rate,npd,npd_wlag)
  }


if (allRFF ==1){
  Results_comb = 
    list.files(Outputs, pattern = paste0("full_streams_rff_vsl10_"),full.names = T) %>% 
    map_df(~read_results(.))
} else {
  Results_comb = 
    list.files(Outputs, pattern = paste0("full_streams_vsl10_"),full.names = T) %>% 
    map_df(~read_results(.))
}


#export global stats
glob_means = tibble()
# Note that the stats and certainty equivalent corrections are only relevant when model is run for all trials (not the mean RFF-SP)
  glob_means =
      bind_rows(
        glob_means,
        Results_comb %>%
        group_by(discount.rate, trial,Pollutant,Model) %>% #group by rate, trial, pollutant, and model, and then sum across all countries
        summarise(
          across(all_of(vars.damages),~sum(get(paste0(cur_column(),"_npd"))), .names="{.col}_NPD"),
          across(all_of(vars.damages),~sum(get(paste0(cur_column(),"_npd")) * get(paste0(cur_column(),"_cert.eq.adj")),na.rm=T), .names="{.col}_NPD.cert.eq"),
          .groups = 'keep') %>%
          ungroup() %>%
        group_by(discount.rate,Pollutant,Model) %>% #next, group by discount rate, pollutant, and model to calculate stats across all trials
        summarise(
          across(all_of(paste0(vars.damages,"_NPD")),~mean(get(cur_column())), .names="{.col}_mean"),
          across(all_of(paste0(vars.damages,"_NPD.cert.eq")),~mean(get(cur_column())), .names="{.col}_mean"),
          across(all_of(paste0(vars.damages,"_NPD")),~quantile(get(cur_column()), 0.025, na.rm=T), .names="{.col}_2.5"),
          across(all_of(paste0(vars.damages,"_NPD")),~quantile(get(cur_column()), 0.975, na.rm=T), .names="{.col}_97.5"),
          across(all_of(paste0(vars.damages,"_NPD")),~median(get(cur_column())), .names="{.col}_median"),
          .groups = 'keep')
  )
# #
#  #write data
  if (allRFF ==1){
    glob_means %>%
      write_csv(file.path(Outputs,paste0('npd_global_rff_means_vsl10.csv')))
  } else {
    glob_means %>%
      write_csv(file.path(Outputs,paste0('npd_global_means_vsl10.csv')))
  }

  #CODE END  