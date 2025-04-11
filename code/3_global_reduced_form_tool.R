### Title: Global Air Quality Reduced Form Tool.R                                     ##
### Purpose: To use the global air quality temperature impact functions to calculate  ####
###           the change in deaths for a user defined scenario
### Created by: Melanie Jackson, IEc 
### Adapted by: E. McDuffie (EPA) & Melanie Jackson, IEc
### Date Created: 3/15/2023      
### Last Edited: 4/3/2025
##      Inputs:    
##        input/RFF/rft_inputs/rff_pop_gdp_all_trials.parquet                         ##
##        input/RFF/rft_inputs/Trajectories Baseline Mortality_Cause Specific.parquet ##
##        input/Final Country Grid Col_Row Index.csv                                  ##
##        input/Country Damage Functions.rds                                          ##
##        input/BenMAP/Global CI Percent Diff from PE Mean.csv                        ##
##        input/RFF/rft_inputs/sampled_pop_trajectory_numbers.csv                     ##
##        input/RFF/rft_inputs/global_mean_surface_temperature_baseline.csv           ##
##     Outputs:
##        output/rft/damages_mean_rft.csv (and parquet)
##        output/rft/damages_itrial_rft.csv
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
#path    <- file.path("//iec.local/jen/CIRA_Health_Analyses/Global PM Analysis/11_RFT/3_Final code")
Inputs  <- file.path("input")
Outputs <- file.path("output","rft")


########
##### Define constants 
#######
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
##### Set User-Defined Scenario Inputs
#######

# Specify Population, Mortality, GDP data
# The tool is currently set up to use data from the RFF-SPs
# The user is asked to specify a specific RFF-SP trajectory Number, OR
# specify 'ALL' to calculate the mortality estimates for all 10,000 trajectories
# OR specify 'mean' to run the model for the mean across the temperature and population
RFF_TrajNumber = 'mean'  #Options: numerical value 1-10000, or 'All', or 'mean'

# Implement cessation lags (1 = yes, 0 = no (default=1))
cessation_flag = 1 #if set to 1, with calculate results with and without lag


########
##### Read in Impact Damages Functions and Supporting Data ####
########

## 1) By degree damage functions (created in code step 1)
ImpactFxns <- readRDS(file.path(Inputs,"Country Damage Functions.rds")) %>%
  rename("LocID" = "Column") %>%
  select(-LM)
# Description: 
#  change in deaths per baseline mortality (referred to as function or fxn) 
#   following change of each pollutant (PM or Ozone) for each temperature scenario (SSP scenario) and model (GISS or CESM2)
#   Remove the columns that aren't needed

## 2) Country Index Key ##
countries <- read.csv(file.path(Inputs,"Final Country Grid Col_Row Index.csv"))[,c(1,3:4,6,8,9)]
# Description: The country name and LocID crosswalk (for analysis 0.5x0.5 grid)

## 3) global percent difference between 2.5th and 97.5th impact results and mean impact results
Global.percentiles <- read.csv(file.path(Inputs,"BenMAP","Global CI Percent Diff from PE Mean.csv"))





##########
###### Pre-Processing Steps
#########
# 1) Load in cross walks: a) between public RFF and underlying mortality trajectories
#                         b) between public RFF and underlying temperature projections, subset for years selected
SampleIDs <- read.csv(file.path(Inputs,"RFF","rft_inputs","sampled_pop_trajectory_numbers.csv"))
# Description: Cross-walk between mortality and RFF trajectory & trial #s
RFFTemps_baseline <- read_parquet(file.path(Inputs,"RFF","rft_inputs","global_mean_surface_temperature_baseline.parquet")) %>% filter(year %in% Years)
RFFTemps_perturbed <- read_parquet(file.path(Inputs,"RFF","rft_inputs","global_mean_surface_temperature_perturbed_co2_2030.parquet")) %>% filter(year %in% Years)

#names(RFFTemps)[1]<-"Year"
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
BaselineMort <- read_parquet(file.path(Inputs,"RFF","rft_inputs", "All Trajectories Baseline Mortality_Cause Specific.parquet"))
# Description: Interpolated timeseries of baseline respiratory mortality (ages 0-99) and NCD+LRI mortality (ages 25-99) for
# all countries in 1000 RFF scenarios
# These respiratory mortality data are already corrected with Int'l Futures respiratory / NCD+LRI to all-cause mortality ratio
pop_gdp_file <- read_parquet(file.path(Inputs,"RFF","rft_inputs",'rffsp_pop_gdp_all_trials.parquet'))
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
#loop through the desired rff-sp trials (mean, all, or specific trial number)
Results <- foreach(itrial = startFile:length(Trajectory),
                   .packages=c('purrr','dplyr','utils','arrow','tidyverse')) %do% { #dopar takes longer to run than do
                     
                     # 1) Subset mortality, population, gdp, and temperature data for the given trial
                     if (RFF_TrajNumber == 'mean') {
                       pop_gdp_data <- pop_gdp_file %>%
                         mutate(gdp_per_cap = gdp/pop) %>%
                         group_by_at(.vars = c('Year','Country')) %>% #mean across trials for each country and year
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
                         group_by(LocID,Year) %>% #average across all trials for each country and year
                         summarize(Resp_BaseMort=mean(Resp_BaseMort),
                                   NCD_LRI_BaseMort=mean(NCD_LRI_BaseMort)) %>%
                         ungroup()
                       
                       Temps.b <- RFFTemps_baseline %>%
                         group_by(year) %>%  #average across all trials for each country and year
                         summarize(Temperature=mean(temp_C_global)) %>%
                         ungroup() 
                       Temps.p <- RFFTemps_perturbed %>%
                         group_by(year) %>%  #average across all trials for each country and year
                         summarize(Temperature=mean(temp_C_global)) %>%
                         ungroup()
                       
                     } else {
                       # Important reminder, gdp and temp data have 10,000 unique projections (trials) 
                       # while mortality has 10,000 projects from 1,000 unique samples.
                       # Therefore, gdp and temp are subset from the 10,000 trials while mortality is subset from the 1,000 samples
                       if(RFF_TrajNumber == 'All'){
                         pop_gdp_data <- pop_gdp_file[pop_gdp_file$trial==itrial,]
                         Temps.b <- RFFTemps_baseline[RFFTemps_baseline$trial==itrial,] %>%
                           rename("Temperature" = "temp_C_global") %>%
                           select(!trial)
                         Temps.p <- RFFTemps_pertrubed[RFFTemps_pertrubed$trial==itrial,] %>%
                           rename("Temperature" = "temp_C_global") %>%
                           select(!trial)
                       } else {
                         pop_gdp_data <- pop_gdp_file[pop_gdp_file$trial==RFF_TrajNumber,] #RFF_TrajNumber can be 1-10,000 projections/trials
                         Temps.b <- RFFTemps_baseline[RFFTemps_baseline$trial==RFF_TrajNumber,] %>%
                           rename("Temperature" = "temp_C_global") %>%
                           select(!trial)
                         Temps.p <- RFFTemps_pertrubed[RFFTemps_pertrubed$trial==RFF_TrajNumber,] %>%
                           rename("Temperature" = "temp_C_global") %>%
                           select(!trial)
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
                     
                     # 2) Use the Impact functions (change in deaths per baseline mortality per temperature change) to calculate
                     #    the number of deaths in each country in each year in the selected trials
                     #  The temperature change of the RFF scenario will (likely) fall between two temperature changes (or SSP scenario results) 
                     ##   run through BenMAP. Therefore, this chunk of code identifies the two nearest temperature points of the
                     ##   impact function (for each model) that are above and below the RFF temperature in each year. 
                     ##  The LM slope and intercepts associated with the particular temperature range for each year will then 
                     ##   be applied to the RFF temperature in that year to calculate the resulting impact (each year, for each country)
                     # If the temperature falls above or below the max/min temperatures of the impact functions, use the max/min impact function values.
                     
                     
                     # Join the damage function data set of the closest matching temperature (above or below) to the RFF temperature 
                     by_low_1  <- join_by(closest(Temperature>=min_temp)) # first join doesn't need to be matched by countries
                     by_low_2  <- join_by(LocID,Country,closest(Temperature>=min_temp)) # join matching countries
                     # set negative temperature changes to 0, not modeled in this tool
                     Impact_Results.b <- Temps.b %>% mutate(Temperature = replace(Temperature, Temperature < 0, 0))
                     Impact_Results.p <- Temps.p %>% mutate(Temperature = replace(Temperature, Temperature < 0, 0))
                     # Must add suffix after each 2 joins
                     Impact_Results.b <- left_join(Impact_Results.b,ImpactFxns[ImpactFxns$Pollutant=="Ozone"&ImpactFxns$Model=="GISS",],by_low_1,multiple="all") %>% 
                       select(!c(Pollutant,Model)) %>%
                       left_join(ImpactFxns[ImpactFxns$Pollutant=="Ozone"&ImpactFxns$Model=="CESM2",],by_low_2,suffix=c("_O3_GISS","_O3_CESM2")) %>% 
                       select(!c(Pollutant,Model)) %>%
                       left_join(ImpactFxns[ImpactFxns$Pollutant=="PM"&ImpactFxns$Model=="GISS",],by_low_2) %>% 
                       select(!c(Pollutant,Model)) %>%
                       left_join(ImpactFxns[ImpactFxns$Pollutant=="PM"&ImpactFxns$Model=="CESM2",],by_low_2,suffix=c("_PM_GISS","_PM_CESM2")) %>% 
                       select(!c(Pollutant,Model))
                     Impact_Results.p <- left_join(Impact_Results.p,ImpactFxns[ImpactFxns$Pollutant=="Ozone"&ImpactFxns$Model=="GISS",],by_low_1,multiple="all") %>% 
                       select(!c(Pollutant,Model)) %>%
                       left_join(ImpactFxns[ImpactFxns$Pollutant=="Ozone"&ImpactFxns$Model=="CESM2",],by_low_2,suffix=c("_O3_GISS","_O3_CESM2")) %>% 
                       select(!c(Pollutant,Model)) %>%
                       left_join(ImpactFxns[ImpactFxns$Pollutant=="PM"&ImpactFxns$Model=="GISS",],by_low_2) %>% 
                       select(!c(Pollutant,Model)) %>%
                       left_join(ImpactFxns[ImpactFxns$Pollutant=="PM"&ImpactFxns$Model=="CESM2",],by_low_2,suffix=c("_PM_GISS","_PM_CESM2")) %>% 
                       select(!c(Pollutant,Model))
                     
                     # pivot results longer so there's a column for model/pollutant
                     Impact_Results.b <- Impact_Results.b %>%
                       pivot_longer(.,-c(year,Temperature,LocID,Country),  # convert annual results from wide to long
                                    names_pattern = '^(.*)_(.*)_(.*)$',
                                    names_to = c(".value","Pollutant","Model"))
                     Impact_Results.p <- Impact_Results.p %>%
                       pivot_longer(.,-c(year,Temperature,LocID,Country),  # convert annual results from wide to long
                                    names_pattern = '^(.*)_(.*)_(.*)$',
                                    names_to = c(".value","Pollutant","Model"))
                     
                     # calculate damage function from linear model slopes and intercepts (this is faster than using the predict function on lm objects)
                     Impact_Results.b <- Impact_Results.b %>%
                       mutate(y_DeathsPerBaseMort = (LM_Slope * Temperature) + LM_Int) %>%
                       select(c(year,LocID,Temperature,Pollutant,Model,y_DeathsPerBaseMort))
                     Impact_Results.p <- Impact_Results.p %>%
                       mutate(y_DeathsPerBaseMort = (LM_Slope * Temperature) + LM_Int) %>%
                       select(c(year,LocID,Temperature,Pollutant,Model,y_DeathsPerBaseMort))
                     
                     # Join damage functions with RFF trial baseline mortality data (pop/gdp info not needed yet).
                     Impact_Results.b <- left_join(Impact_Results.b,BaselineMort_data,by=c("LocID","year"="Year"))
                     Impact_Results.p <- left_join(Impact_Results.p,BaselineMort_data,by=c("LocID","year"="Year"))
                     
                     # Calculate the change in deaths for each pollutant/model and country/year
                     # Baseline mortality is cause-specific, Respiratory (Resp) is a result of Ozone (ages 0-99)
                     # Non-Communicable Diseases + Lower Respiratory Infection (NCD+LRI) is a result of PM (ages 25-99)
                     Impact_Results.b <- Impact_Results.b %>%
                       mutate(Deaths = case_when(Pollutant == "O3" ~ y_DeathsPerBaseMort * Resp_BaseMort,
                                                 Pollutant == "PM" ~ y_DeathsPerBaseMort * NCD_LRI_BaseMort))
                     Impact_Results.p <- Impact_Results.p %>%
                       mutate(Deaths = case_when(Pollutant == "O3" ~ y_DeathsPerBaseMort * Resp_BaseMort,
                                                 Pollutant == "PM" ~ y_DeathsPerBaseMort * NCD_LRI_BaseMort))
                     
                     #This spreads the annual mortality data across 20 years into the future using
                     # EPA standard cessation lags (for O3 only?)
                     if (cessation_flag ==1) {
                       
                       # calculate results with 20-yr cessation lag
                       # Arrange data by year and group by country, pollutant, and model, then apply cessation lags to the previous [n]years results
                       # note: more concise code to do this utilizes complete and fill, however, this is too slow when looping through 'all' scenarios
                       Impact_Results_wlags.b <- Impact_Results.b %>%
                         select(c('year','LocID','Temperature','Pollutant','Model','Deaths')) %>%
                         arrange(year) %>%
                         group_by(LocID,Pollutant,Model) %>%
                         mutate(Deaths_wlag = Deaths * lags[1] +
                                  ifelse(year-1 < 2020, 0, lag(Deaths, 1, order_by = year) * lags[2]) +
                                  ifelse(year-2 < 2020, 0, lag(Deaths, 2, order_by = year) * lags[3]) +
                                  ifelse(year-3 < 2020, 0, lag(Deaths, 3, order_by = year) * lags[4]) +
                                  ifelse(year-4 < 2020, 0, lag(Deaths, 4, order_by = year) * lags[5]) +
                                  ifelse(year-5 < 2020, 0, lag(Deaths, 5, order_by = year) * lags[6]) +
                                  ifelse(year-6 < 2020, 0, lag(Deaths, 6, order_by = year) * lags[7]) +
                                  ifelse(year-7 < 2020, 0, lag(Deaths, 7, order_by = year) * lags[8]) +
                                  ifelse(year-8 < 2020, 0, lag(Deaths, 8, order_by = year) * lags[9]) +
                                  ifelse(year-9 < 2020, 0, lag(Deaths, 9, order_by = year) * lags[10]) +
                                  ifelse(year-10 < 2020, 0, lag(Deaths,10, order_by = year) * lags[11]) +
                                  ifelse(year-11 < 2020, 0, lag(Deaths, 11, order_by = year) * lags[12]) +
                                  ifelse(year-12 < 2020, 0, lag(Deaths, 12, order_by = year) * lags[13]) +
                                  ifelse(year-13 < 2020, 0, lag(Deaths, 13, order_by = year) * lags[14]) +
                                  ifelse(year-14 < 2020, 0, lag(Deaths, 14, order_by = year) * lags[15]) +
                                  ifelse(year-15 < 2020, 0, lag(Deaths, 15, order_by = year) * lags[16]) +
                                  ifelse(year-16 < 2020, 0, lag(Deaths, 16, order_by = year) * lags[17]) +
                                  ifelse(year-17 < 2020, 0, lag(Deaths, 17, order_by = year) * lags[18]) +
                                  ifelse(year-18 < 2020, 0, lag(Deaths, 18, order_by = year) * lags[19]) +
                                  ifelse(year-19 < 2020, 0, lag(Deaths, 19, order_by = year) * lags[20]))
                       
                       Impact_Results.b <- Impact_Results_wlags.b %>%
                         left_join(pop_gdp_data,by=c("LocID"="COL","year"="Year")) %>%        
                         left_join(countries[,c("COUNTRY","COL","Region")],by=c("LocID"="COL")) %>%
                         rename("Country"="COUNTRY") %>%
                         ungroup()
                       rm(Impact_Results_wlags.b)
                       Impact_Results_wlags.p <- Impact_Results.p %>%
                         select(c('year','LocID','Temperature','Pollutant','Model','Deaths')) %>%
                         arrange(year) %>%
                         group_by(LocID,Pollutant,Model) %>%
                         mutate(Deaths_wlag = Deaths * lags[1] +
                                  ifelse(year-1 < 2020, 0, lag(Deaths, 1, order_by = year) * lags[2]) +
                                  ifelse(year-2 < 2020, 0, lag(Deaths, 2, order_by = year) * lags[3]) +
                                  ifelse(year-3 < 2020, 0, lag(Deaths, 3, order_by = year) * lags[4]) +
                                  ifelse(year-4 < 2020, 0, lag(Deaths, 4, order_by = year) * lags[5]) +
                                  ifelse(year-5 < 2020, 0, lag(Deaths, 5, order_by = year) * lags[6]) +
                                  ifelse(year-6 < 2020, 0, lag(Deaths, 6, order_by = year) * lags[7]) +
                                  ifelse(year-7 < 2020, 0, lag(Deaths, 7, order_by = year) * lags[8]) +
                                  ifelse(year-8 < 2020, 0, lag(Deaths, 8, order_by = year) * lags[9]) +
                                  ifelse(year-9 < 2020, 0, lag(Deaths, 9, order_by = year) * lags[10]) +
                                  ifelse(year-10 < 2020, 0, lag(Deaths,10, order_by = year) * lags[11]) +
                                  ifelse(year-11 < 2020, 0, lag(Deaths, 11, order_by = year) * lags[12]) +
                                  ifelse(year-12 < 2020, 0, lag(Deaths, 12, order_by = year) * lags[13]) +
                                  ifelse(year-13 < 2020, 0, lag(Deaths, 13, order_by = year) * lags[14]) +
                                  ifelse(year-14 < 2020, 0, lag(Deaths, 14, order_by = year) * lags[15]) +
                                  ifelse(year-15 < 2020, 0, lag(Deaths, 15, order_by = year) * lags[16]) +
                                  ifelse(year-16 < 2020, 0, lag(Deaths, 16, order_by = year) * lags[17]) +
                                  ifelse(year-17 < 2020, 0, lag(Deaths, 17, order_by = year) * lags[18]) +
                                  ifelse(year-18 < 2020, 0, lag(Deaths, 18, order_by = year) * lags[19]) +
                                  ifelse(year-19 < 2020, 0, lag(Deaths, 19, order_by = year) * lags[20]))
                       Impact_Results.p <- Impact_Results_wlags.p %>%
                         left_join(pop_gdp_data,by=c("LocID"="COL","year"="Year")) %>%        
                         left_join(countries[,c("COUNTRY","COL","Region")],by=c("LocID"="COL")) %>%
                         rename("Country"="COUNTRY") %>%
                         ungroup()
                       rm(Impact_Results_wlags.p)
                     } else {
                       # join pop/gdp data and region names
                       Impact_Results.b <- Impact_Results.b %>%
                         select(c(LocID,Country,Year,Temperature,contains("Deaths"))) %>%
                         left_join(pop_gdp_data,by=c("LocID"="COL","year"="Year")) %>%
                         left_join(countries[,c("COUNTRY","COL","Region")],by=c("LocID"="COL")) %>%
                         rename("Country"="COUNTRY") %>%
                         ungroup()
                       Impact_Results.p <- Impact_Results.p %>%
                         select(c(LocID,Country,Year,Temperature,contains("Deaths"))) %>%
                         left_join(pop_gdp_data,by=c("LocID"="COL","year"="Year")) %>%
                         left_join(countries[,c("COUNTRY","COL","Region")],by=c("LocID"="COL")) %>%
                         rename("Country"="COUNTRY") %>%
                         ungroup()
                     }
                     
                     ######
                     ####VALUATION STEP ###
                     #####
                     #Valuation of annual (non-discounted impacts)
                     #1. Calculate reference vsl
                     usa_base_income <- Impact_Results.b %>%
                       filter(LocID == 840,           #USA 
                              year == 2020,
                              Pollutant == "O3", Model == "GISS") %>% # this filter doesn't change results, just produces 1 output instead of 4
                       select(gdp_per_cap) %>%
                       rename(base_income = gdp_per_cap)
                     
                     # 2. Calculate country-specific VSL in each year
                     Impact_Results.b <- Impact_Results.b %>%
                       mutate(vsl = base_vsl*((gdp_per_cap/usa_base_income$base_income)^Elasticity))
                     Impact_Results.p <- Impact_Results.p %>%
                       mutate(vsl = base_vsl*((gdp_per_cap/usa_base_income$base_income)^Elasticity))
                     
                     # 3. calculate regional VSL & assign to countries with missing vsl data
                     VSLResults.Reg.b <- Impact_Results.b %>%
                       filter(Pollutant == "O3", Model == "GISS") %>% # this filter doesn't change results, just produces 1 output instead of 4 (2 polls * 2 models of results)
                       group_by_at(.vars = c('Region','year')) %>%
                       summarise_at(.vars= c('gdp','pop'), sum, na.rm=TRUE) %>%
                       mutate(region_gdp_per_cap = gdp/pop) %>%
                       mutate(vsl_reg = base_vsl*((region_gdp_per_cap/usa_base_income$base_income)^Elasticity))
                     VSLResults.Reg.p <- Impact_Results.p %>%
                       filter(Pollutant == "O3", Model == "GISS") %>% # this filter doesn't change results, just produces 1 output instead of 4 (2 polls * 2 models of results)
                       group_by_at(.vars = c('Region','year')) %>%
                       summarise_at(.vars= c('gdp','pop'), sum, na.rm=TRUE) %>%
                       mutate(region_gdp_per_cap = gdp/pop) %>%
                       mutate(vsl_reg = base_vsl*((region_gdp_per_cap/usa_base_income$base_income)^Elasticity))
                     
                     # Assign regional vsl when there is no country vsl available
                     # join regional data, then use coalesce (fills in NAs)
                     Impact_Results.b <- left_join(Impact_Results.b,VSLResults.Reg.b,by=c("Region","year"="year"),suffix=c("","_reg"))
                     Impact_Results.b <- Impact_Results.b %>%
                       mutate(vsl = coalesce(vsl,vsl_reg),
                              pop = coalesce(pop,pop_reg),
                              gdp = coalesce(gdp,gdp_reg))
                     Impact_Results.p <- left_join(Impact_Results.p,VSLResults.Reg.p,by=c("Region","year"="year"),suffix=c("","_reg"))
                     Impact_Results.p <- Impact_Results.p %>%
                       mutate(vsl = coalesce(vsl,vsl_reg),
                              pop = coalesce(pop,pop_reg),
                              gdp = coalesce(gdp,gdp_reg))
                     
                     
                     # 4. Calculate annual monetary (un-discounted) damages
                     PhysicalImpacts.b <- names(Impact_Results.b)[grep("Deaths",names(Impact_Results.b))] # this method includes wlag if it exists
                     Impact_Results.b <- Impact_Results.b %>%
                       mutate(across(all_of(PhysicalImpacts.b), ~ . * vsl,.names = "Ann_{.col}"),
                              trial = as.numeric(itrial))
                     colnames(Impact_Results.b) <- gsub("Ann_Deaths","Annual_impacts",colnames(Impact_Results.b))  
                     Impact_Results.b$Pollutant[Impact_Results.b$Pollutant == "O3"] <- "Ozone"
                     PhysicalImpacts.p <- names(Impact_Results.p)[grep("Deaths",names(Impact_Results.p))] # this method includes wlag if it exists
                     Impact_Results.p <- Impact_Results.p %>%
                       mutate(across(all_of(PhysicalImpacts.p), ~ . * vsl,.names = "Ann_{.col}"),
                              trial = as.numeric(itrial))
                     colnames(Impact_Results.p) <- gsub("Ann_Deaths","Annual_impacts",colnames(Impact_Results.p))  
                     Impact_Results.p$Pollutant[Impact_Results.p$Pollutant == "O3"] <- "Ozone"
                     
                     
                     # For mean scenario, calculate 2.5th and 97.5th percentiles using global average percent difference from mean 
                     # Due to time constraints, unable to add full monte carlo
                     ResultNames.b <- names(Impact_Results.b)[grep("Deaths|Annual",names(Impact_Results.b))]
                     ResultNames.p <- names(Impact_Results.p)[grep("Deaths|Annual",names(Impact_Results.p))]
                     if (RFF_TrajNumber == 'mean' ) {
                       Impact_Results.b <- left_join(Impact_Results.b,Global.percentiles,by=c("Pollutant"="pollutant","Model"="model")) %>%
                         mutate(across(all_of(ResultNames.b), ~ . * (1-PE_2_5_PctDiff) ,.names="{.col}_2_5"),
                                across(all_of(ResultNames.b), ~ . * (1+PE_97_5_PctDiff),.names="{.col}_97_5")) %>%
                         select(-c(PE_2_5_PctDiff, PE_97_5_PctDiff))
                       Impact_Results.p <- left_join(Impact_Results.p,Global.percentiles,by=c("Pollutant"="pollutant","Model"="model")) %>%
                         mutate(across(all_of(ResultNames.p), ~ . * (1-PE_2_5_PctDiff) ,.names="{.col}_2_5"),
                                across(all_of(ResultNames.p), ~ . * (1+PE_97_5_PctDiff),.names="{.col}_97_5")) %>%
                         select(-c(PE_2_5_PctDiff, PE_97_5_PctDiff))
                     } 
                     
                     #cleanup/format results and column names
                     Impact_Results.b <- Impact_Results.b %>%
                       select(-c(gdp_per_cap,vsl,gdp_reg,pop_reg,vsl_reg,region_gdp_per_cap,vsl_reg))
                     ResultNames.b <- names(Impact_Results.b)[grep("Deaths|Annual",names(Impact_Results.b))]
                     Impact_Results.b <- Impact_Results.b[,c("year","LocID","Country","Region","Temperature","pop","gdp","Pollutant","Model",ResultNames.b,"trial")]
                     Impact_Results.p <- Impact_Results.p %>%
                       select(-c(gdp_per_cap,vsl,gdp_reg,pop_reg,vsl_reg,region_gdp_per_cap,vsl_reg))
                     ResultNames.p <- names(Impact_Results.p)[grep("Deaths|Annual",names(Impact_Results.p))]
                     Impact_Results.p <- Impact_Results.p[,c("year","LocID","Country","Region","Temperature","pop","gdp","Pollutant","Model",ResultNames.p,"trial")]
                     
                     ##COMBINE 
                     results_i <- Impact_Results.b %>% mutate(damageType = "Baseline") %>% 
                       rbind( Impact_Results.p %>% mutate(damageType = "Perturbed") )
                     ###
                     
                     if (RFF_TrajNumber == 'mean' ) {
                       results_i %>%
                         write_parquet(file.path(Outputs,paste0('damages_mean_rft.parquet'))) %>%
                         write.csv(file.path(Outputs,paste0('damages_mean_rft.csv')),row.names=F)
                     } else {
                       results_i %>% 
                         write_parquet(file.path(Outputs,paste0('damages_',itrial,'_rft.parquet')))
                     }
                     
                     rm()
                   }#End Scenario Loop

### stop the clock\
time2 <- Sys.time(); time2 - time1
proc.time() - ptm
###### Finish #####
### stop cluster
#  parallel::stopCluster(cl = my.cluster)


## print mean baseline 2100 results
filter_results <- Impact_Results.b %>% 
  filter(year ==2100) %>%
  group_by_at(.vars = c('Country','Pollutant')) %>% #mean across models
  summarise_at(.vars = c('Deaths','Deaths_2_5','Deaths_97_5','Annual_impacts','Annual_impacts_2_5','Annual_impacts_97_5'), mean) %>%
  ungroup()
print('Global PM Results: ')
print(filter_results %>% filter(Pollutant == 'PM') %>% summarise_at(.vars = c('Deaths','Deaths_2_5','Deaths_97_5','Annual_impacts','Annual_impacts_2_5','Annual_impacts_97_5'), sum))
print('Global Ozone Results: ')
print(filter_results %>% filter(Pollutant == 'Ozone') %>% summarise_at(.vars = c('Deaths','Deaths_2_5','Deaths_97_5','Annual_impacts','Annual_impacts_2_5','Annual_impacts_97_5'), sum))
print('Global Net Results: ')
print(filter_results %>% summarise_at(.vars = c('Deaths','Deaths_2_5','Deaths_97_5','Annual_impacts','Annual_impacts_2_5','Annual_impacts_97_5'), sum))


# Code Complete. Discounting of annual damages done in Code file #3    
