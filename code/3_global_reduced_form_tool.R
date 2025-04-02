### Title: Global Air Quality Reduced Form Tool.R                                     ##
### Purpose: To adjust the global health effects impact generated in BenCloud  ####
### Created by: Melanie Jackson, IEc 
### Adapted by: E. McDuffie (EPA) & Melanie Jackson, IEc
### Date Created: 3/15/2023      
### Last Edited: 4/1/2025
##      Inputs:                                                                       ##
##        -2095 Country Results 2095 by Model, SSP & Temperature                      ##
##        -All Trajectories Baseline Respiratory Mortality_Ages0to99.rds (6 parts)    ##
##        -All Trajectories Baseline NCD+LRI Mortality_Ages25to99.rds (6 parts)       ##
##        -sampled_pop_trajectory_numbers.csv                                         ##
##        -Final Country Grid Col_Row Index.csv                                       ##
##     Outputs:
##        - output/rft/
## Units: Code adjusts all $ values for $2020                                         ##
########################################################################################

#### Read in necessary packages, set working paths, set constant variables  ####
#clean up space before loop
rm(list = ls()); gc()

#necessary packages
packs<-c("dplyr","tidyverse","readxl","purrr","foreach","utils","pbmcapply","arrow")
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

#working directories
path    <- file.path("//iec.local/jen/CIRA_Health_Analyses/Global PM Analysis/11_RFT/3_Final code")
Inputs  <- file.path(path,"Inputs")
Outputs <- file.path(path,"Outputs","RFT")


########
##### Define constants 
#######

#constant variables (These are the conditions run through the BenMAP WebTool)
#gdp_2011_to_2020    = 113.784/98.164 # GDP Implicit Price Deflators (https://apps.bea.gov/iTable/?reqid=19&step=3&isuri=1&select_all_years=0&nipa_table_list=13&series=a&first_year=2006&last_year=2020&scale=-99&categories=survey&thetable= )
#last access: March 30, 2023
gdp_2011_to_2020    = 105.361/91.481 # GDP Implicit Price Deflators (https://apps.bea.gov/iTable/?reqid=19&step=3&isuri=1&select_all_years=0&nipa_table_list=13&series=a&first_year=2006&last_year=2020&scale=-99&categories=survey&thetable= )
#last access: February 10, 2025
base_vsl            = 10.05e6 # $2020 VSL from Rennert et al., 2022 (and EPA SC-GHG report) 9.33e6  # = USD VSL in 2006$, inflated to 2020$ (from EPA 2010)
Elasticity          = 1  # income elasticity
endYear             = 2300
Years               <- c(2020:endYear) # set simulation years
lags                = c(6, 2.5, 2.5, 2.5, 2.5, 4/15, 4/15, 4/15, 4/15,4/15, 4/15, 4/15, 4/15, 
                        4/15, 4/15, 4/15, 4/15, 4/15, 4/15, 4/15)/20
# These are EPA standard 20 year cessation lags (%) for particulate matter and
# represent the time delay between the year of exposure and distribution of attributable
# deaths over the next 20 years

########
##### Read in Original BenMAP Simulation Data ####
########


## 1) Original BenMAP results 
CloudFxns <- readRDS(file.path(Inputs,"BenMAP","Country Damage Functions.rds"))
# Description: 
#  change in deaths per baseline mortality (referred to as function or fxn) 
#   following change in temperature (SSP scenario) of each pollutant (PM or Ozone) and model (GISS or CESM2)

## Country Index Key ##
countries <- read.csv(file.path(Inputs,"Final Country Grid Col_Row Index.csv"))[,c(1,3:4,6,8,9)]
# Description: The country name and LocID crosswalk (for BenMAP 0.5x0.5 grid)

## global percent difference between 2.5th and 97.5th BenMAP results and mean BenMAP results
Global.percentiles <- read.csv(file.path(Inputs,"BenMAP","Global CI Percent Diff from PE Mean.csv"))


########
##### Set User-Defined Scenario Inputs
#######

# Specify Population, Mortality, GDP data
# The tool is currently set up to use data from the RFF-SPs
# The user is asked to specify a specific RFF-SP trajectory Number, OR
# specify 'ALL' to calculate the mortality estimates for all 10,000 trajectories
RFF_TrajNumber = 'mean'  #Options: numerical value 1-10000, or 'All', or 'mean'

# Implement cessation lags (1 = yes, 0 = no (default=1))
cessation_flag = 1 #if set to 1, with calculate results with and without lag


##########
###### Pre-Processing Steps
#########

# 1) Load in cross walks: a) between public RFF and underlying mortality trajectories
#                         b) between public RFF and underlying temperatures of projections, subset for years selected
SampleIDs <- read.csv(file.path(Inputs,"sampled_pop_trajectory_numbers.csv"))
# Description: Cross-walk between mortality and RFF trajectory & trial #s
RFFTemps <- read.csv(file.path(Inputs,"TempNorm_1850to1900_global_temperature_norm.csv")) %>% filter(time %in% Years)
names(RFFTemps)[1]<-"Year"
# Description: annual temperature associated with the RFF trial number

# 2) Set the trajectory numbers for the user-defined scenario
if (RFF_TrajNumber == 'All') {
  Trajectory = SampleIDs$x
  message("Running all 10,000 projections")
} else if (RFF_TrajNumber == 'mean') {
  Trajectory =1 #will trigger mean data in main loop
  message("Running mean of 10,000 projections")
} else {
  Trajectory = SampleIDs$x[as.integer(RFF_TrajNumber)]
  message(paste0("Running trajectory #",Trajectory))
}

# 3) Read in pre-processed baseline mortality values and GDP data
BaselineMort <- read_parquet(file.path(Inputs,"RFF","All Trajectories Baseline Mortality_Cause Specific.parquet"))
# Description: Interpolated timeseries of baseline respiratory mortality (ages 0-99) and NCD+LRI mortality (ages 25-99) for
# all countries in 1000 RFF scenarios
# These data are already corrected with Int'l Futures respiratory / NCD+LRI to all-cause mortality ratio
pop_gdp_file <- read_parquet(file.path(Inputs,"RFF",'rffsp_pop_gdp_all_trials.parquet'))
# Description: GDP (in $2020) and national population of 10000 RFF trials


########
##### Begin loop through Specified SocioEconomic Scenarios ####
########    

###### Run RFT ######
# Start the clock!
ptm   <- proc.time(); time1 <- Sys.time()

## start parallel (use this to specify the start file [if simulation is disrupted])
startFile = 1

####
##Begin Analysis ##
####
#loop through the trajectories (or scenarios)
Results <- foreach(itrial = startFile:length(Trajectory),
                   .packages=c('purrr','dplyr','utils','arrow','tidyverse')) %do% { #dopar takes longer to run than do
                     
                     # 1) Subset mortality, population, gdp, and temperature data for the given trial
                     if (RFF_TrajNumber == 'mean') {
                       pop_gdp_data <- pop_gdp_file %>%
                         mutate(gdp_per_cap = gdp/pop) %>%
                         group_by_at(.vars = c('Year','Country')) %>% #sum across trials for each country
                         summarise_at(.vars = c('gdp_per_cap','pop','gdp'), mean) %>%
                         ungroup() %>%
                         #interpolate between 5 year intervals
                         group_by(Country) %>%
                         complete(Year=min(Year):endYear)%>%
                         mutate(gdp_per_cap = approx(x=Year,y=gdp_per_cap,xout=2020:endYear, rule =2)$y,
                                pop = approx(x=Year,y=pop,xout=2020:endYear, rule=2)$y,
                                gdp = approx(x=Year,y=gdp,xout=2020:endYear, rule=2)$y)
                       pop_gdp_data <- right_join(countries,pop_gdp_data, by= c("RFF_iso_code"="Country"),multiple='all') %>%
                         select(COL,Year,pop,gdp,gdp_per_cap)
                       
                       BaselineMort_data <- BaselineMort %>%
                         group_by(LocID,Year) %>% #average across each country and year
                         summarize(Resp_BaseMort=mean(Resp_BaseMort),
                                   NCD_LRI_BaseMort=mean(NCD_LRI_BaseMort)) %>%
                         ungroup()
                       
                       Temps <- RFFTemps %>%
                         group_by(Year) %>%  #average across country and year
                         summarize(Temperature=mean(global_temperature_norm)) %>%
                         ungroup()
                       
                     } else {
                       # Important reminder, gdp and temp data have 10,000 unique projections (trials) 
                       # while mortality has 10,000 projects from 1,000 unique samples.
                       # Therefore, gdp and temp are subset from the 10,000 trials while mortality is subset from the 1,000 samples
                       if(RFF_TrajNumber == 'All'){
                         pop_gdp_data <- pop_gdp_file[pop_gdp_file$trial==itrial,]
                         Temps <- RFFTemps[RFFTemps$trialnum==itrial,] %>%
                           rename("Temperature" = "global_temperature_norm") %>%
                           select(!trialnum)
                       }else{
                         pop_gdp_data <- pop_gdp_file[pop_gdp_file$trial==RFF_TrajNumber,] #RFF_TrajNumber can be 1-10,000 projections/trials
                         Temps <- RFFTemps[RFFTemps$trialnum==RFF_TrajNumber,] %>%
                           rename("Temperature" = "global_temperature_norm") %>%
                           select(!trialnum)
                       }
                       pop_gdp_data <- right_join(countries,pop_gdp_data, by= c("RFF_iso_code"="Country"),multiple='all')
                       pop_gdp_data <- pop_gdp_data %>% 
                         select(COL,Year,pop,gdp) %>%
                         mutate(gdp_per_cap = gdp/pop) %>%
                         group_by_at(.vars = c('COL')) %>% 
                         #interpolate between 5 year intervals
                         complete(Year=min(Year):endYear)%>%
                         mutate(gdp_per_cap = approx(x=Year,y=gdp_per_cap,xout=2020:endYear, rule=2)$y,
                                pop = approx(x=Year,y=pop,xout=2020:endYear, rule=2)$y,
                                gdp = approx(x=Year,y=gdp,xout=2020:endYear, rule=2)$y)
                       
                       BaselineMort_data <- BaselineMort[BaselineMort$Trajectory==Trajectory[itrial],] #Trajectory is the sampleID that aligns with whichever trial/projection (RFF_TrajNumber) is selected
                     }
                     
                     # 2) Identify the BenMAP damage function(s) (change in deaths per baseline mortality following change in temperature)
                     # corresponding to the modeled change in temperature of the RFF trial, for each pollutant/model.
                     # RFF temperature change will (likely) fall between two temperature changes run through BenMAP,
                     # Therefore, identify the BenMAP runs with the nearest temperature modeled above and below the RFF temperature.
                     # We will then create a linear model to calculate the damage function of the RFF temperature change.
                     # If the temperature falls above or below the max/min temperatures modeled in BenMAP, use the max/min damage function.
                     
                     # Drop contextual columns from cloud fxns table
                     Fxns <- CloudFxns %>% 
                       rename("LocID" = "Column") %>%
                       select(-LM)
                     
                     # Join the damage function dataset of the closest matching temperature (above or below) to the RFF temperatre 
                     by_low_1  <- join_by(closest(Temperature>=min_temp)) # first join doesn't need to be matched by countries
                     by_low_2  <- join_by(LocID,Country,closest(Temperature>=min_temp)) # join matching countries
                     # set negative temperature changes to 0, not modeled in this tool
                     Analysis <- Temps %>% mutate(Temperature = replace(Temperature, Temperature < 0, 0))
                     # Must add suffix after each 2 joins
                     Analysis <- left_join(Analysis,Fxns[Fxns$Pollutant=="Ozone"&Fxns$Model=="GISS",],by_low_1,multiple="all") %>% 
                       select(!c(Pollutant,Model)) %>%
                       left_join(Fxns[Fxns$Pollutant=="Ozone"&Fxns$Model=="CESM2",],by_low_2,suffix=c("_O3_GISS","_O3_CESM2")) %>% 
                       select(!c(Pollutant,Model)) %>%
                       left_join(Fxns[Fxns$Pollutant=="PM"&Fxns$Model=="GISS",],by_low_2) %>% 
                       select(!c(Pollutant,Model)) %>%
                       left_join(Fxns[Fxns$Pollutant=="PM"&Fxns$Model=="CESM2",],by_low_2,suffix=c("_PM_GISS","_PM_CESM2")) %>% 
                       select(!c(Pollutant,Model))
                     
                     # pivot results longer so there's a column for model/pollutant
                     Analysis <- Analysis %>%
                       pivot_longer(.,-c(Year,Temperature,LocID,Country),  # convert annual results from wide to long
                                    names_pattern = '^(.*)_(.*)_(.*)$',
                                    names_to = c(".value","Pollutant","Model"))
                     
                     # calculate damage function from linear model slopes and intercepts (this is faster than using the predict function on lm objects)
                     Analysis <- Analysis %>%
                       mutate(Fxn_Int = (LM_Slope * Temperature) + LM_Int) %>%
                       select(c(Year,LocID,Temperature,Pollutant,Model,Fxn_Int))
                     
                     # Join damage functions with RFF trial baseline mortality data (pop/gdp info not needed yet).
                     Analysis <- left_join(Analysis,BaselineMort_data,by=c("LocID","Year"))
                     
                     # Calculate the change in deaths for each pollutant/model and country/year
                     # Baseline mortality is cause-specific, Respiratory (Resp) is a result of Ozone (ages 0-99)
                     # Non-Communicable Diseases + Lower Respiratory Infection (NCD+LRI) is a result of PM (ages 25-99)
                     Analysis <- Analysis %>%
                       mutate(Deaths = case_when(Pollutant == "O3" ~ Fxn_Int * Resp_BaseMort,
                                                 Pollutant == "PM" ~ Fxn_Int * NCD_LRI_BaseMort))
                     
                     #This spreads the annual mortality data across 20 years into the future using
                     # EPA standard cessation lags
                     if (cessation_flag ==1) {
                       
                       # calculate results with 20-yr cessation lag
                       # Arrange data by year and group by country, pollutant, and model, then apply cessation lags to the previous [n]years results
                       # note: more concise code to do this utilizes complete and fill, however, this is too slow when looping through 'all' scenarios
                       Analysis_cess <- Analysis %>%
                         select(c('Year','LocID','Pollutant','Model','Deaths')) %>%
                         arrange(Year) %>%
                         group_by(LocID,Pollutant,Model) %>%
                         mutate(Deaths_wlag = Deaths * lags[1] +
                                  ifelse(Year-1 < 2020, 0, lag(Deaths, 1, order_by = Year) * lags[2]) +
                                  ifelse(Year-2 < 2020, 0, lag(Deaths, 2, order_by = Year) * lags[3]) +
                                  ifelse(Year-3 < 2020, 0, lag(Deaths, 3, order_by = Year) * lags[4]) +
                                  ifelse(Year-4 < 2020, 0, lag(Deaths, 4, order_by = Year) * lags[5]) +
                                  ifelse(Year-5 < 2020, 0, lag(Deaths, 5, order_by = Year) * lags[6]) +
                                  ifelse(Year-6 < 2020, 0, lag(Deaths, 6, order_by = Year) * lags[7]) +
                                  ifelse(Year-7 < 2020, 0, lag(Deaths, 7, order_by = Year) * lags[8]) +
                                  ifelse(Year-8 < 2020, 0, lag(Deaths, 8, order_by = Year) * lags[9]) +
                                  ifelse(Year-9 < 2020, 0, lag(Deaths, 9, order_by = Year) * lags[10]) +
                                  ifelse(Year-10 < 2020, 0, lag(Deaths,10, order_by = Year) * lags[11]) +
                                  ifelse(Year-11 < 2020, 0, lag(Deaths, 11, order_by = Year) * lags[12]) +
                                  ifelse(Year-12 < 2020, 0, lag(Deaths, 12, order_by = Year) * lags[13]) +
                                  ifelse(Year-13 < 2020, 0, lag(Deaths, 13, order_by = Year) * lags[14]) +
                                  ifelse(Year-14 < 2020, 0, lag(Deaths, 14, order_by = Year) * lags[15]) +
                                  ifelse(Year-15 < 2020, 0, lag(Deaths, 15, order_by = Year) * lags[16]) +
                                  ifelse(Year-16 < 2020, 0, lag(Deaths, 16, order_by = Year) * lags[17]) +
                                  ifelse(Year-17 < 2020, 0, lag(Deaths, 17, order_by = Year) * lags[18]) +
                                  ifelse(Year-18 < 2020, 0, lag(Deaths, 18, order_by = Year) * lags[19]) +
                                  ifelse(Year-19 < 2020, 0, lag(Deaths, 19, order_by = Year) * lags[20]))
                       
                       Analysis <- Analysis_cess %>%
                         left_join(pop_gdp_data,by=c("LocID"="COL","Year"="Year")) %>%        
                         left_join(countries[,c("COUNTRY","COL","Region")],by=c("LocID"="COL")) %>%
                         rename("Country"="COUNTRY") %>%
                         ungroup()
                       
                       rm(Analysis_cess)
                     }else{
                       # join pop/gdp data and region names
                       Analysis <- Analysis %>%
                         select(c(LocID,Country,Year,contains("Deaths"))) %>%
                         left_join(pop_gdp_data,by=c("LocID"="COL","Year"="Year")) %>%
                         left_join(countries[,c("COUNTRY","COL","Region")],by=c("LocID"="COL")) %>%
                         rename("Country"="COUNTRY") %>%
                         ungroup()
                     }
                     
                     ######
                     ####VALUATION STEP ###
                     #####
                     #Valuation of annual (non-discounted impacts)
                     #1. Calculate reference vsl
                     usa_base_income <- Analysis %>%
                       filter(LocID == 840,           #USA 
                              Year == 2020,
                              Pollutant == "O3", Model == "GISS") %>% # this filter doesn't change results, just produces 1 output instead of 4
                       select(gdp_per_cap) %>%
                       rename(base_income = gdp_per_cap)
                     
                     # 2. Calculate country-specific VSL
                     Analysis <- Analysis %>%
                       mutate(vsl = base_vsl*((gdp_per_cap/usa_base_income$base_income)^Elasticity))
                     
                     # 3. calculate regional VSL & assign to countries with missing vsl data
                     VSLResults.Reg <- Analysis %>%
                       filter(Pollutant == "O3", Model == "GISS") %>% # this filter doesn't change results, just produces 1 output instead of 4 (2 polls * 2 models of results)
                       group_by_at(.vars = c('Region','Year')) %>%
                       summarise_at(.vars= c('gdp','pop'), sum, na.rm=TRUE) %>%
                       mutate(region_gdp_per_cap = gdp/pop) %>%
                       #select(-c('gdp','pop')) %>%
                       mutate(vsl_reg = base_vsl*((region_gdp_per_cap/usa_base_income$base_income)^Elasticity))
                     
                     # Assign regional vsl when there is no country vsl available
                     # join regional data, then use coalesce (fills in NAs)
                     Analysis <- left_join(Analysis,VSLResults.Reg,by=c("Region","Year"),suffix=c("","_reg"))
                     Analysis <- Analysis %>%
                       mutate(vsl = coalesce(vsl,vsl_reg),
                              pop = coalesce(pop,pop_reg),
                              gdp = coalesce(gdp,gdp_reg))
                     
                     # 4. Calculate annual monetary (undiscounted damages)
                     PhysicalImpacts <- names(Analysis)[grep("Deaths",names(Analysis))] # this method includes wlag if it exists
                     Analysis <- Analysis %>%
                       mutate(across(all_of(PhysicalImpacts), ~ . * vsl,.names = "Ann_{.col}"),
                              trial = as.numeric(itrial))
                     colnames(Analysis) <- gsub("Ann_Deaths","Annual_$",colnames(Analysis))  
                     
                     Analysis$Pollutant[Analysis$Pollutant == "O3"] <- "Ozone"
                     
                     # For mean scenario, calculate 2.5th and 97.5th percentiles using global average percent difference from mean 
                     # Due to time constraints, unable to add full monte carlo
                     ResultNames <- names(Analysis)[grep("Deaths|Annual",names(Analysis))]
                     if (RFF_TrajNumber == 'mean' ) {
                       Analysis <- left_join(Analysis,Global.percentiles,by=c("Pollutant"="pollutant","Model"="model")) %>%
                         mutate(across(all_of(ResultNames), ~ . * (1-PE_2_5_PctDiff) ,.names="{.col}_2_5"),
                                across(all_of(ResultNames), ~ . * (1+PE_97_5_PctDiff),.names="{.col}_97_5")) %>%
                         select(-c(PE_2_5_PctDiff, PE_97_5_PctDiff))
                     } 
                     
                     #cleanup/format results and column names
                     Analysis <- Analysis %>%
                       select(-c(gdp_per_cap,vsl,gdp_reg,pop_reg,vsl_reg,region_gdp_per_cap,vsl_reg))
                     ResultNames <- names(Analysis)[grep("Deaths|Annual",names(Analysis))]
                     Analysis <- Analysis[,c("Year","LocID","Country","Region","pop","gdp","Pollutant","Model",ResultNames,"trial")]
                     
                     
                     if (RFF_TrajNumber == 'mean' ) {
                       Analysis %>%
                         write_parquet(file.path(Outputs,paste0('damages_mean_vsl10_rft.parquet'))) %>%
                         write.csv(file.path(Outputs,paste0('damages_mean_vsl10_rft.csv')),row.names=F)
                     } else {
                       Analysis %>% 
                         write_parquet(file.path(Outputs,paste0('damages_',itrial,'_vsl10_rft.parquet')))
                     }
                     
                     rm()
                   }#End Scenario Loop

### stop the clock\
time2 <- Sys.time(); time2 - time1
proc.time() - ptm
###### Finish #####
### stop cluster
#  parallel::stopCluster(cl = my.cluster)

# Code Complete. Discounting of annual damages done in Code file #3    
