#####
## Title:   2_build_inputs_for_rtf.R
## Purpose: This file is meant to read in raw RFF scenarios (from Zenodo) and 
##          builds complete GDP (in $2020) and national population 
##          timeseries (2010-2300) for input into the ozone and PM mortality model 
##          - RFF sp population is in thousands and is converted to # of people.
##          - RFF sp is in 2011$ and is converted to 2020$
##          This file also formats the MimiGIVE temperature files for use in the rft (from sc-ghg work)
##          This file also creates a look-up table of the RFF scenarios' baseline mortality
##           baseline mortality = mortality rate * population
## Inputs:  inputs/RFF/pop_income/rffsp_pop_income_run_TRIAL.feather
##          inputs/IFs/mx_trajectories.csv
##          inputs/IFs/pop_trajectories.csv
##          inputs/IFs/adjustment_factors.csv
##          inputs/IFs/country_crosswalk_2024.csv
##          inputs/IFs/Crosswalk RFF Mort Year.xlsx
##          inputs/IFs/ifs_country_codes_2024.csv
##          inputs/RFF/temperature/path_to_temperature_file
##          inputs/RFF/temperature/rffsp_fair_sequence.csv
## Outputs: inputs/RFF/rft_inputs/global_mean_surface_temperature_baseline.parquet
##          inputs/RFF/rft_inputs/global_mean_surface_temperature_perturbed_GAS_YEAR.parquet
##          inputs/RFF/rft_inputs/rffsp_pop_gdp_all_trials.parquet
##          inputs/RFF/rft_inputs/All Trajectories Baseline Mortality_Cause Specific.parquet
## Written by: US EPA, Climate Change Division; March 2023 & Melanie Jackson, IEc
## Last updated: 4/3/2025 by E. McDuffie, EPA
#####


##########################
#################  library
##########################

## Clear worksace
rm(list = ls())
gc()

## This function will check if a package is installed and, if not, install it
list.of.packages <- c('magrittr','tidyverse',
                      'readxl', 'arrow',
                      'foreach','doParallel','pbapply')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.rstudio.com/")
lapply(list.of.packages, library, character.only = TRUE)

##########################
##############  data paths
##########################
inputpath_rff = file.path('input','RFF','pop_income')
input_path_temp = file.path('input','RFF','temperature')
input_path_mort <- file.path('input','IFs')
inputpath_ifs = file.path('input','IFs')
rft_input_path = file.path('input','RFF','rft_inputs')

##########################
####################  data
##########################

### Read-in RFF-SP Population and GSP data, format to correct units
#rff-sps are in units of 1000's of people and 2011 dollars. Covert to people and 2020$
gdp_2011_to_2020 = 105.361/91.481 # GDP Implicit Price Deflators (https://apps.bea.gov/iTable/?reqid=19&step=3&isuri=1&select_all_years=0&nipa_table_list=13&series=a&first_year=2006&last_year=2020&scale=-99&categories=survey&thetable= )
#last access: April 1, 2025

#### read rffsps
#collect file names
##### Collect Files ####
c_iteration  <- inputpath_rff %>% list.files(pattern = ".feather") %>%
  (function(x){sub(".feather", "", x)}) %>%
  (function(x){sub("rffsp_pop_income_run_", "", x)}) %>%
  as.numeric %>% sort; 
c_iteration %>% length


rffsp_pop_gdp <- 
  pblapply(1:length(c_iteration), function(i){
    ### File name
    infile_i  <- inputpath_rff %>%
      file.path("rffsp_pop_income_run") %>%
      paste(c_iteration[i], sep="_") %>%
      paste0(".", "feather")
    ### Read in data and return
    data_i    <- infile_i %>% read_feather
    ### Filter data for model type, national total, desired sectors, baseline scenario
    data_i    <- data_i   %>% 
      filter(Year <= 2300) %>% #filter for pre-2300
      mutate(gdp = GDP * gdp_2011_to_2020 * 1e6) %>% #RFF data is in millions of 2011$, convert to 2020$
      mutate(pop = Pop * 1e3) %>%                    #RFF population data is in 1000s, convert to individual ppl
      mutate(dollar.year = 2020) %>%
      select(-c('GDP','Pop')) %>%
      mutate(trial = i) #%>%
    ## export scenario-specific population files to read into the RFT
    #write_feather(rft_input_path %>% file.path("RFF",paste0('rffsp_pop_gdp_', i, '.feather')))
  })  %>%
  
  #also bind all together
  (function(x){
    do.call(rbind, x)
  }); rffsp_pop_gdp %>% glimpse

### Save single file with data from all trials
rffsp_pop_gdp %>%
  write_parquet(rft_input_path %>% file.path("rffsp_pop_gdp_all_trials.parquet"))

############################
#### format temperature files  ####
############################
### these are baseline and perturbed temperature output from MimiGIVE (as prepared for EPA sc-ghg 2025 guidance)
### MimiGIVE was run with rff-sp GHG emission inputs in FaIR v1.6.2, with carbon feedbacks
# temperatures have been re-baselined so that they represent the temperature change relative to the 
# 2005-2014 average period. This baseline can be changed here in future analyses. 

#read in rff-sp trial sample numbers
## final sample randomly selected in MimiGIVE
rffsp_sample = 
  input_path_temp %>%
  file.path('rffsp_fair_sequence.csv') %>%
  read_csv(show_col_types = FALSE) %>% 
  select(-fair.id)

# baseline temperature data
temp = 
  input_path_temp %>%
  file.path('CO2-2030-climate_variables_with_carbon_feedbacks','model_1','TempNorm_1850to1900_global_temperature_norm.csv') %>%
  read_csv(show_col_types = FALSE) %>% 
  rename(year=time, temp_C_global=2, trial=trialnum) %>% 
  filter(year > 2005)

## recover 2005-2014 baseline scale
temp.relative.baseline = 
  temp %>% 
  filter(year %in% seq(2005, 2014, 5)) %>% 
  group_by(trial) %>% 
  summarise(base.2005.2014 = mean(temp_C_global)) %>% 
  ungroup()

## rescale temp and trim for input into the rft
temp %<>% 
  left_join(temp.relative.baseline) %>% 
  mutate(temp_C_global = temp_C_global - base.2005.2014) %>% 
  filter(year %in% seq(2020, 2300, 5)) %>%
  select(-base.2005.2014) %>% 
  left_join(rffsp_sample) %>% 
  relocate(trial, rffsp.id, year)

## export as inputs to fredi
temp %>% 
  write_parquet(rft_input_path %>% file.path('global_mean_surface_temperature_baseline.parquet'))

## export scenario-specific files to read into fredi
#foreach(j = jList, .packages=c('tidyverse')) %dopar% {
#  temp %>% 
#    filter(trial==j) %>% 
#    select(year, temp_C_global) %>% 
#    write_csv(fredi_input_path %>% file.path('temp_baseline', paste0('temp_baseline_', j, '.csv')))
#}

#### process data for 2030 co2 emissions pulse (perturbed year) ########
# can upload additional temperature/gas scenarios from DMAP
GAS = 'co2'
YEAR = 2030

temp =
  input_path_temp %>%
  file.path(file.path('CO2-2030-climate_variables_with_carbon_feedbacks','model_2','TempNorm_1850to1900_global_temperature_norm.csv')) %>%
  read_csv(show_col_types = FALSE) %>% 
  rename(year          = time, 
         temp_C_global = 2, 
         trial         = trialnum) %>% 
  filter(year > 2005)

## rescale temp and trim for fredi
temp %<>% 
  left_join(temp.relative.baseline,
            by = 'trial') %>% 
  mutate(temp_C_global = temp_C_global - base.2005.2014) %>% 
  filter(year %in% seq(2020, 2300, 5)) %>%
  select(-base.2005.2014) %>% 
  left_join(rffsp_sample,
            by = 'trial') %>% 
  relocate(trial, rffsp.id, year)

## export as inputs to fredi
temp %>% 
  write_parquet(rft_input_path %>% 
                  file.path(paste0('global_mean_surface_temperature_perturbed_', GAS, '_', YEAR, '.parquet')))

###############################
#### Process baseline mortality projections ####
##################################
### These are also from the rff-sps, they are the age-specific background mortality
### rate data used to develop the rff-sp population trajectories

### Read in inputs
#input data (takes a few minutes to load, these files are large)
mx  <- read.csv(file.path(input_path_mort,"mx_trajectories.csv"))
pop <- read.csv(file.path(input_path_mort,"pop_trajectories.csv"))
# must crosswalk mortality year ranges to years first
Cross.year <- read_xlsx(file.path(input_path_mort,"Crosswalk RFF Mort Year.xlsx"))
mx  <- left_join(mx,Cross.year,by=c("Year"="RFF year"))
# multiply population by 1000, RFF stores population in 1,000s
pop$Pop <- pop$Pop*1000

#International Futures (IFs) respiratory-related mortality to all-cause mortality
# adjustments - requires crosswalking IFs country ID to RFF country ID
IFsFactors   <- read.csv(file.path(inputpath_ifs,"adjustment_factors.csv"))
IFsCountry   <- read.csv(file.path(inputpath_ifs,"ifs_country_codes_2024.csv"))
IFsCtryCross <- read.csv(file.path(inputpath_ifs,"country_crosswalk_2024.csv"))
IFsFactors   <-left_join(IFsFactors,IFsCountry,by=c("Region"="RegionNum"))
IFsFactors   <- left_join(IFsFactors,IFsCtryCross,by=c("RegionName"="IFsName"))
# remove NAs, two IFs countries not in RFF list
IFsFactors   <- IFsFactors[!is.na(IFsFactors$RFFName),]
IFsFactors   < -IFsFactors[,1:8]
names(IFsFactors)[7] <- "LocID"
#extend to 2300 (hold 2100 rates constant)
IFsFactors<- IFsFactors %>%
  group_by(LocID,RFF_Age_range,Region,RegionName,RFFName)%>%
  complete(Years=2019:2300)%>%
  fill(Resp_AF,NCD_LRI_AF)%>%
  ungroup()
rm(IFsCountry,IFsCtryCross)

### Calculate baseline mortality
# join mortality by age-with population
mx <- left_join(mx,pop,by=c("LocID","Min.Year"="Year","Age","Trajectory"))
# remove pop object to create space
rm(pop)
# join with IFs factors
mx <- left_join(mx,IFsFactors,by=c("LocID","Age"="RFF_Age_range","Year"="Years"))
# calculate cause-specific rates
mx$Resp_Rate     <-mx$mx * mx$Resp_AF
mx$NCD_LRI_Rate  <- mx$mx * mx$NCD_LRI_AF
# calculate cause-specific baseline mortality counts
mx$Resp_BaseMort <- mx$Resp_Rate * mx$Pop
mx$NCD_LRI_BaseMort <- mx$NCD_LRI_Rate * mx$Pop

# summarize across age groups for each endpoint (Resp is 0-99, NCD+LRI is 25-99)
NCDAges <- c("25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90+")
# create respiratory summary (all ages, 0-99)      
mx.resp <- mx %>%
  group_by(LocID,Period,Min.Year,Trajectory) %>%
  summarize(Resp_BaseMort=sum(Resp_BaseMort))
names(mx.resp)[3] <- "Year"
mx.NCD.LRI <- mx[mx$Age %in% NCDAges,] %>%
  group_by(LocID,Period,Min.Year,Trajectory) %>%
  summarize(NCD_LRI_BaseMort=sum(NCD_LRI_BaseMort))      
names(mx.NCD.LRI)[3] <- "Year"

# Join cause-specific summaries into single dataframe
mx.mort <- left_join(mx.resp,mx.NCD.LRI,by=c("LocID","Period","Year","Trajectory"))
rm(mx.resp,mx.NCD.LRI)

# interpolate for years between 5-year bins
mx.mort <- mx.mort %>%
  group_by(LocID,Trajectory) %>%
  complete(Year=2020:2295)%>%
  mutate(Resp_BaseMort = approx(x=Year,y=Resp_BaseMort,xout=2020:2295)$y,
         NCD_LRI_BaseMort = approx(x=Year,y=NCD_LRI_BaseMort,xout=2020:2295)$y)
# extend 2095 to 2300
mx.mort <- mx.mort %>%
  group_by(LocID,Trajectory) %>%
  complete(Year=2095:2300) %>%
  fill(Resp_BaseMort,NCD_LRI_BaseMort)

# export
# remove  period column to reduce file size
mx.mort <- mx.mort[,-4]
write_parquet(mx.mort,file.path(rft_input_path,"All Trajectories Baseline Mortality_Cause Specific.parquet")) #RFT pulls full list of all trials
