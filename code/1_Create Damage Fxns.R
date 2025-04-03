### Title: 1_Create Damage Functions.R
### Purpose: To convert the webtool outputs (deaths per baseline mortality per scenario 
###           into temperature-specific impact functions, by country
### Created by: Melanie Jackson, IEc, edited by E. McDuffie, OAP/EPA
### Date Created: 2/4/2025
### Last Edited: 4/4/2025
##      Inputs:
##        - Country level results for each model and pollutant of 2095
##        - SSP temperature crosswalk
######################################################################################

#clear environment
rm(list = ls()); gc()

#### Step 0. Set data paths ####
#setwd("~/CCD-CSIB/GitHub/Code/globalAQ_rft")
Outpath<-"input"
benmap_data <- "analysis/data"

library(dplyr)
library(tidyverse)

#####################
### Step 1. Read in summary results data from BenMAP ###
#### Clean up country summaries, subset for relevant data ####

# read in summary BenMAP results
CountryResults   <- read.csv(file.path(benmap_data,"result_summary_by_country.csv"))
CountryResults <- CountryResults[CountryResults$year==2095,]

# rename columns, re-order, re-format, add origin point
names(CountryResults) <- c("Column","ChangeInDeaths","Population","BaseMort","Pollutant","Model",
                             "Scenario","Year","ChangePerCapita","DeathsPerBaseMort","Temp","Country")
CountryResults<-CountryResults[,c("Column","Country","Pollutant","Model","Scenario","ChangeInDeaths","DeathsPerBaseMort","Temp")]
CountryResults$Model[CountryResults$Model=="NCAR"]<-"CESM2"
origins <- CountryResults[CountryResults$Scenario==126,]
origins[,c("Scenario","ChangeInDeaths",
               "DeathsPerBaseMort","Temp")] <- 0
CountryResults <- rbind(CountryResults,origins)
write.csv(CountryResults,file.path(Outpath,"Intermediate_Country Results.csv"),row.names=FALSE)
 
######################
### Step 2. Create impact functions 

## Create linear regressions, join next highest temperature data as column to table
# create a table with the current value and the next highest value (to create regression between)
join <- join_by(Column,Country,Pollutant,Model,closest(Temp<Temp))
# subset dataframe for only the needed columns (i.e., do not need SPP name and change in deaths to create linear regression)
CountryResults   <- CountryResults[,-c(5)]
CountryFunctions <- left_join(CountryResults,CountryResults,by=join,suffix=c("_L","_H"))
# where _H column is NA, the low temperature is the highest modeled temperature, set _H results = _L (use coalesce)
# find the first non-missing value and use that to replace missing high values (i.e., keep the functions constant after highest temperature)
CountryFunctions <- CountryFunctions %>%
   mutate(DeathsPerBaseMort_H = coalesce(DeathsPerBaseMort_H,DeathsPerBaseMort_L),
          Temp_H = coalesce(Temp_H,Temp_L))

# Transform high/low data points from wide to long (necessary to create linear functions by country/year)   
Models <- CountryFunctions %>%
   mutate(min_temp = Temp_L , max_temp = Temp_H) %>%
   pivot_longer(cols = -c("Column","Country","Pollutant","Model","min_temp","max_temp"),
                names_to =c('.value','Direction'),
                names_pattern = '(.*)_(.)$') # this pattern creates a single column for each variable that has a high/low variation (e.g., Fxn_O3_GISS_L and Fxn_O3_GISS_H will go into Fxn_O3_GISS)
# Create linear regressions between deaths per base mort and global temperature for each country, pollutant, and model
Models <- Models %>%
   group_by(Column,Country,Pollutant,Model,min_temp,max_temp) %>%
   do(LM = lm(DeathsPerBaseMort~Temp,data=.))
  
# add coefficients as fields to dataframe
## the RFT works faster using the coefs in formula rather than using predict function
# round results to 10, otherwise get decimals out to 19th place and intersects that should be 0 have an on-zero value
Models$LM_Int        <- as.numeric(lapply(Models$LM, function(x) coef(x)[[1]])) %>% round(.,digits=10)
Models$LM_Slope      <- as.numeric(lapply(Models$LM, function(x) coef(x)[[2]])) %>% round(.,digits=10) %>% 
   replace(is.na(.),0) #NA's occur at highest temperature modeled, therefore slope is 0
 
# export as rds to retain linear model objects
saveRDS(Models,file.path(Outpath,"Country Damage Functions.rds"))
  
  