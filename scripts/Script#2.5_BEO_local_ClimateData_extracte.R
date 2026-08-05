### ------------------------------------------------------------------------------------------------------------

### 10/01/25 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to use the query_timeseries() function of the 'rTubeDB' package of Stephan Wöllauer (Instrumentation Core Team)
### to download the hourly temperature and precipitation and soil moisture (and wind maybe) data from the Explortaories' weather stations

### R script to: 
###  - Test the different options of the 'rTubeDB' package; Identify the best parameters to download the gap-filled
###    hourly data
### - Download all variables of interest at daily scale (mean,min,max) and store on home dir of 'climcal3' server

### Resources: 
# https://environmentalinformatics-marburg.github.io/tubedb/
# http://137.248.186.133:61036/content/visualisation_meta/visualisation_meta.html for list of sensors and what they measure 

### Last update: 24/09/25 (Download daily PAR and radiation veraibles for as many EPs as possible -> help interpret air temperature offsets daily variations)

### ------------------------------------------------------------------------------------------------------------

# Libraries install.packages("tidyverse")
library("devtools")
library("dplyr")
library("reshape2")
library("data.table")
library("lubridate")
library("viridis")
library("scales") 
library("RColorBrewer")

### ------------------------------------------------------------------------------------------------------------

### 1°) Install the rTubeDB package and automatically install updated versions.
# In some cases a restart of R is needed to work with a updated version of rTubeDB package (in RStudio - Session - Terminate R).
if(!require('remotes')) install.packages('remotes')
if(!require('httr')) install.packages('httr')
remotes::install_github('environmentalinformatics-marburg/tubedb/rTubeDB')
library("rTubeDB")
# ?TubeDB


### 2°) Open connexion to server and explore options/sensors/variables
# Open TubeDB server connection
tubedb <- rTubeDB::TubeDB(url = 'http://137.248.186.133:61036', user = 'fabio.benedetti', password = 'Dkxn7QiE')

# Get regions/projects as data.frame
regionDF <- rTubeDB::query_regions(tubedb)

# Get plots of Exploratories as data.frame
plotDF <- rTubeDB::query_region_plots(tubedb, 'BE')
# a data.frame containing: plot info and metadata
# str(plotDF)
# dim(plotDF) # 635 x 7
# head(plotDF)

# Get all sensors of all plots of Exploratories as data.frame. Not all that sensors are from the main plots.
sensorDF <- rTubeDB::query_region_sensors(tubedb, 'BE') # Also a ddf
# head(sensorDF)
# rownames(sensorDF) # vector of variables measured in the field?

# Get sensors (and derived sensors) from one plot
sensorList <- rTubeDB::query_region_plot_sensors(tubedb,'BE','AEG01')
sensorList # yes, list of variables avaiable for plot HEG01

# Get sensor list for all 300 EPs (based on 'plotDF') and identify those sensors common to ALL EP
# Get vector of EP names 
setwd("/home/fbenedetti/Exploratories/EP metadata/1000_10_Dataset")
plots <- read.csv("1000_10_data.csv", h = T)
# dim(plots); str(plots)
# Keep those that are not 'na'
plots <- plots[plots$EP_Plot_ID != "na",]
unique(plots$EP_Plot_ID) # OK, still need to add some 0 though
# FUN from Marc Beringer
BEplotZeros <- function(dat, column, plotnam = "PlotSTD"){
      dat <- as.data.frame(dat)
      funz <- function(x) ifelse((nchar(as.character(x))==4), gsub("(.)$", "0\\1", x), as.character(x)) # eo funz
      dat[,plotnam] <- sapply(dat[,column],funz)
      return(dat)
} # eo FUN
# Apply
plots <- BEplotZeros(plots, "EP_Plot_ID", plotnam = "EP")
plots <- unique(plots$EP); plots

# With lapply, extract sensorList for all EPs
p <- "HEG40" # for testing
require("parallel")
res <- lapply(X = plots, FUN = function(p) {
          message(p)
          s <- rTubeDB::query_region_plot_sensors(tubedb,'BE',p)
          d <- data.frame(EP = p, Sensors = s)
          return(d)
      }
) # eo lapply
# Rbind
table.sensors <- dplyr::bind_rows(res)
# dim(table.sensors) ; head(table.sensors)
# unique(table.sensors$EP); unique(table.sensors$Sensors)

# Identify those sensors that are available for all 300 EPs
# Find sensors that are present in all 300 locations
common_sensors <- data.frame(
      table.sensors %>%
      group_by(factor(Sensors)) %>%
      summarise(n_locations = n_distinct(factor(EP)))
)
colnames(common_sensors) <- c("Sensors","N_plots")
common_sensors <- common_sensors[order(common_sensors$N_plots, decreasing = TRUE),]
common_sensors[common_sensors$N_plots >= 299,"Sensors"]
# [1] CycleCounter                         precipitation_dwd_year              
# [3] precipitation_radolan                precipitation_radolan_acc           
# [5] precipitation_radolan_rain_days      precipitation_radolan_rain_days_RR1 
# [7] precipitation_radolan_rain_days_RR10 precipitation_radolan_rain_days_RR20
# [9] precipitation_radolan_rain_days_RR5  rH_200                              
#[11] rH_200_DMR                           rH_200_max                          
#[13] rH_200_min                           SM_10                               
#[15] SM_10_max                            SM_10_min                           
#[17] SM_10_U                              Ta_10                               
#[19] Ta_10_max                            Ta_10_min                           
#[21] Ta_200                               Ta_200_cold_days                    
#[23] Ta_200_cold_sum                      Ta_200_cool_days                    
#[25] Ta_200_dew_point                     Ta_200_dew_point_approx             
#[27] Ta_200_DTR                           Ta_200_extremely hot days           
#[29] Ta_200_extremely_cold_days           Ta_200_growing_degree_days_10       
#[31] Ta_200_growth_sum                    Ta_200_gruenlandtemperatur          
#[33] Ta_200_gruenlandtemperatursumme      Ta_200_heat_index                   
#[35] Ta_200_heating_degree_days           Ta_200_heating_degree_sum           
#[37] Ta_200_humidex                       Ta_200_ice_days                     
#[39] Ta_200_max                           Ta_200_min                          
#[41] Ta_200_ref                           Ta_200_ref_diff                     
#[43] Ta_200_SSI                           Ta_200_summer_days                  
#[45] Ta_200_tropical_days                 Ta_200_tropical_nights              
#[47] Ta_200_vapour_pressure               Ta_200_warm_days                    
#[49] Ta_200_warm10_sum                    Ta_200_warm20_sum                   
#[51] Ta_200_wet_bulb_temperature          Ts_05                               
#[53] Ts_05_max                            Ts_05_min                           
#[55] Ts_10                                Ts_10_max                           
#[57] Ts_10_min                            Ts_20                               
#[59] Ts_20_max                            Ts_20_min                           
#[61] Ts_50                                Ts_50_max                           
#[63] Ts_50_min                            UB    

### NOTE:   http://137.248.186.133:61036/content/visualisation_meta/visualisation_meta.html
###         for list of sensors and what they measure

### NOTE: rH = relative air humidity? -> saturation vapor pressure at T (air temperature) and at T_d (dew point)
### Magnus formula: T = air temperature, Td = dew point temperature
### Requires temperatures in °C (Celsius), not Kelvin
calc_rh <- function(T, Td) {  
    # Magnus formula constants
    a <- 17.62
    b <- 243.12
    # Saturation vapor pressures (hPa)
    es_T  <- 6.112 * exp((a * T)  / (b + T))
    es_Td <- 6.112 * exp((a * Td) / (b + Td))
    # Relative humidity (%)
    RH <- 100 * (es_Td / es_T)
    return(RH)
} # eo FUN - calc_rh
# Example
calc_rh(T = 20, Td = 10)

# In ERA5-Land: https://codes.ecmwf.int/grib/param-db/168 

### 3°) Query TS for some plots and variables
# Get climate time series hours of two sensors over 2 plots
# For testing, start with few plots as it takes several minutes to process and interpolate many plots.
# ?query_timeseries

test <- rTubeDB::query_timeseries(
  tubedb,
  plot=c('HEG01','HEG02','HEG01'),
  sensor=c('Ta_200','rH_200'),
  datetimeFormat='POSIXlt',
  start=2010,
  end=2017,
  aggregation='hour',
  quality='empirical',
  interpolated=TRUE,
  colPlot=TRUE,
  colYear=TRUE,
  colDayOfYear=TRUE,
  colMonth=TRUE,
  colDay=TRUE,
  colHour=TRUE,
  casted=FALSE, # column style order or row style order
)
# Check output of test
# class(test) # a ddf, good
# str(test) # very good structure
dim(test); head(test)

# Get climate time series days of over some plots as data.frame
test2 <- rTubeDB::query_timeseries(
      tubedb,
      plot = c('HEG01','HEG02','HEG03'),
      sensor = c('Ta_200','Ta_200_min','Ta_200_max'),
      datetimeFormat = 'POSIXlt',
      start = 2008,
      end = 2025,
      aggregation = 'day',
      quality = 'empirical',
      interpolated = TRUE,
      quality_counter = TRUE, # count interpolated hour values (quality flag)
      colPlot = TRUE,
      colYear = TRUE,
      colDayOfYear = TRUE,
      colMonth = TRUE,
      colDay = TRUE,
      casted = FALSE,
) # Quite fast! :)
# Check output of test2
str(test2)
dim(test2)
head(test2)

### NOTES about 'quality' argument --> quality checks, one of: no, physical, step, empirical.
### "The physical check implemented in TubeDB is an ab- solute method, step checking is 
### relative to values within one time-series, and empirical checking is relative to values
### between time-series" (Wöllauer et al., 2021)

### - Physical checking filters measurements based on expectable value range
###   (such as temperatures between − 40 and 60 ◦C). This basic check ensures
###   that all values are in the specified range, and coarse outliers are excluded

### - Step checking filters the data by the strength of a fluctuation in a temporal value.
###   Most climate measurements do change slowly over time. Abrupt value changes above a
###   threshold can be used to identify outliers beyond the expected value range

### - Empirical checking compares measurements to reference time-series. This check can reveal
###   subtle errors that may be missed by visual time- series inspection, such as slow sensor
###   degradation or unintentionally changed measurement conditions at a climate station

### According to the discussion I had with Stephan W. and Paul M. in Wernigerode, I should 
### use interpolation == TRUE and quality == 'empirical'

# To write to file
# data.table::fwrite(subset(tsDF_hours, select=-c(datetime)), 'climate_data_hours.csv')
# data.table::fwrite(subset(tsDF_days, select=-c(datetime)), 'climate_data_days.csv')

### ------------------------------------------------------------------------------------------------------------

### 11/03/25: Write FUN to download all interpolated daily variables needed

# In the query_timeseries() FUN, download all EPs --> Need vector with all EP names
setwd("/home/fbenedetti/Exploratories/EP metadata/1000_10_Dataset")
plots <- read.csv("1000_10_data.csv", h = T)
plots <- plots[plots$EP_Plot_ID != "na",]
# Need to add the 0 in the EP id though
BEplotZeros <- function(dat, column, plotnam = "PlotSTD"){
      dat <- as.data.frame(dat)
      funz <- function(x) ifelse((nchar(as.character(x))==4), gsub("(.)$", "0\\1", x), as.character(x)) # eo funz
      dat[,plotnam] <- sapply(dat[,column],funz)
      return(dat)
} # eo FUN
# Apply
plots <- BEplotZeros(plots, "EP_Plot_ID", plotnam = "EP")
plots <- unique(plots$EP)
# plots # tut gut

# Go back to dir
setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/")

### 19/03/25: Write FUN to query the local climate measurements of the BE with 'rTubeDB'
# sensor <- "Ts_10"
# system <- "grasslands"
# region <- "SCH"
# time <- "month"
# start <- "2008"
# end <- "2025"

query_local_climate <- function(sensor, system, region, time, start, end) {

      #' This function takes four arguments and returns a formatted data.frame:
      #'
      #' @param sensor the sensor to query from 'TubeDB' (character)
      #' @param system System name (character): 'grasslands' or 'forests'
      #' @param region Region name (character): 'SCH' or 'HND' or 'SWA'
      #' @param time the scale at which the sensor data should be aggregated (character): 'hour' 'day' 'week' 'month' or 'year'
      #' @param start the starting year of the time series (character)
      #' @param the last year of the time series (character)
      #' 
      #' @return A formatted data.frame combining the inputs.

      # Message
      message(paste("Querying ",sensor," measurements for the ",system," of the ",region,
                  " at the ",time," scale from ",start," to ",end, sep = ""))

      # Identify plots of interest
      if(system == "grasslands") {
            # Keep plots with 'G' inside
            plots_sub <- plots[grepl("G",plots)]
      } else {
            # Keep plots with 'W' inside
            plots_sub <- plots[grepl("W",plots)]
      } # eo if else loop - systems

      # Identify plots of interest based on region of interest
      if(region == "HND") {
            # Keep plots with 'H' inside
            plots_sub2 <- plots_sub[grepl("H",plots_sub)]
      } else if (region == "SCH") {
            # Keep plots with 'S' inside
            plots_sub2 <- plots_sub[grepl("S",plots_sub)]
      } else {
            # Keep plots with 'A' inside
            plots_sub2 <- plots_sub[grepl("A",plots_sub)]
      } # eo if else loop - regions

      # Put vector of plots in order
      plots_sub2 <- plots_sub2[order(as.numeric(gsub("HEG", "", plots_sub2)))]

      # Got to dir of interest
      if(time == "day") {
            setwd(paste("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily/",system,"/",region,"/", sep = ""))
      } # eo if loop

      if(time == "month") {
            setwd(paste("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/monthly/", sep = ""))
      } # eo if loop

      # Extracting the data
      clim_data <- rTubeDB::query_timeseries(
                  tubedb,
                  plot = plots_sub2,
                  sensor = sensor,
                  datetimeFormat = 'character',
                  start = start,
                  end = end,
                  aggregation = time,
                  quality = 'empirical',
                  interpolated = TRUE,
                  quality_counter = TRUE, # count interpolated hour values (quality flag)
                  colPlot = TRUE,
                  colYear = TRUE,
                  colDayOfYear = TRUE,
                  colMonth = TRUE,
                  colDay = TRUE,
                  casted = FALSE,
      ) 
      # Check 'clim_data'
      # head(clim_data)
      # summary(clim_data)
      # str(clim_data)

      # Adjust some colnames to your liking
      colnames(clim_data)[c(1,2)] <- c("EP","Date")
      # Convert to actual 'Date' vector if time is day
      if(time == "day") {
            clim_data$Date <- as.Date(clim_data$Date)
      } # eo if loop

      # Looks good, can be saved
      # Message
      message(paste("Saving ",sensor," measurements for the ",system," of the ",region,
                  " at the ",time," scale from ",start," to ",end,"\n", sep = ""))

      if(time == "day") {
             data.table::fwrite(x = clim_data, file = paste("table_all_EPs_interp_daily_",sensor,"_",system,"_",region,"_",start,"-",end,"_",Sys.Date(),".csv", sep = "") )
      } else if(time == "month") {
            data.table::fwrite(x = clim_data, file = paste("table_all_EPs_interp_monthly_",sensor,"_",system,"_",region,"_",start,"-",end,"_",Sys.Date(),".csv", sep = "") )
      } # eo if else loop
      
      # Go back to dir and clean
      setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/")
      rm(clim_data,plots_sub2); gc()

} # eo FUN - query_local_climate

# Apply FUN above
# query_local_climate(sensor = "precipitation_radolan", system = "grasslands", region = "HND", time = "day", start = "2008", end = "2025") # works well
sensors <- c("precipitation_radolan",
      "Ta_200","Ta_200_min","Ta_200_max",
      "Ta_10","Ta_10_min","Ta_10_max",
      "SM_10","SM_10_min","SM_10_max",
      "Ts_05","Ts_05_min","Ts_05_max",
      "Ts_10","Ts_10_min","Ts_10_max",
      "Ts_20","Ts_20_min","Ts_20_max",
      "Ts_50","Ts_50_min","Ts_50_max",
      "rH_200","rH_200_min","rH_200_max"
) # eo sensors

# Loop over systems and regions
systems <- c("grasslands","forests")
regions <- c("HND","SCH","SWA")


### 20/03/25: Or use the parallel version on climcal3 to go a lot faster:
library("parallel")
for(s in systems) {
      for(r in regions) {
            mclapply(sensors, function(sens) { 
                        query_local_climate(sensor = sens, system = s, region = r, time = "day", start = "2008", end = "2025") 
                  }, mc.cores = length(sensors)
            ) # eo mclapply
      } # eo 2nd for loop - regions 
} # eo first for loop - systems 


### 05/06/2025: Downloading monthly Ts data for Thu Zar Nwe
library("parallel")
sensors <- c("Ts_10","Ts_20")
regions <- c("HND","SCH","SWA")
for(r in regions) {
      mclapply(sensors, function(sens) { 
                  query_local_climate(sensor = sens, system = "grasslands", region = r, time = "month", start = "2008", end = "2025") 
            }, mc.cores = length(sensors)
      ) # eo mclapply
} # eo 2nd for loop - regions 


### ------------------------------------------------------------------------------------------------------------

### 24/09/2025: Identify grassland and forest EPs that have sensors measuring:
# - radiation (PAR, Net radiation - Rn)
# - shortwave upward radiation = how much incoming shortwave radiation is being reflected away through albedo)
#     --> Compare SW (or albedo) of forests vs grasslands in early spring: 
# Grasslands: darker, already green → absorb more sunlight (low SW). But then they lose a large share of this extra energy as latent heat (LE),
# which cools the surface
# Forests: reflect more (higher SW), so they start with less net energy. But since they don’t spend much on evaporation yet,
# most of the remaining energy goes into sensible heat (H), which directly warms the air inside the canopy.
# - evaporation 

### -> These should help you interpret patterns in air temperature offsets between the grasslands and the forests

# PAR [µmol/(m^2*s)]
table.sensors[table.sensors$Sensors == "PAR","EP"] # 299 EPs - gut

# evaporation [mm/h]
table.sensors[table.sensors$Sensors == "evaporation","EP"] # 21 EPs - and ONLY in the grasslands...OK can still work

# Net radiation = Rn [kWh/m^2]
table.sensors[table.sensors$Sensors == "Rn","EP"] # same as "evaporation" 

# shortwave upward radiation = SWUR [kWh/m^2]
table.sensors[table.sensors$Sensors == "SWUR","EP"] # same as "evaporation" 
table.sensors[table.sensors$Sensors == "SWUR_300","EP"] # same but at 3m height - # same EPs as "evaporation" -> only grasslands EPs then


### Extracting & saving the data
s <- "2009"
e <- "2025"
t <- "day"

## A°) PAR
setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily/"); dir()

clim_data <- rTubeDB::query_timeseries(
      tubedb,
      plot = unique(table.sensors[table.sensors$Sensors == "PAR","EP"]),
      sensor = "PAR",
      datetimeFormat = 'character',
      start = s,
      end = e,
      aggregation = t,
      quality = 'empirical',
      interpolated = TRUE,
      quality_counter = TRUE,
      colPlot = TRUE,
      colYear = TRUE,
      colDayOfYear = TRUE,
      colMonth = TRUE,
      colDay = TRUE,
      casted = FALSE,
) 
# Check 'clim_data'
head(clim_data)
summary(clim_data)
str(clim_data)

# Adjust some colnames to your liking
colnames(clim_data)[c(1,2)] <- c("EP","Date")
# Convert to actual 'Date' vector if time is day
clim_data$Date <- as.Date(clim_data$Date)

# Save 
data.table::fwrite(x = clim_data, file = paste("table_all_EPs_interp_daily_PAR_",s,"-",e,"_",Sys.Date(),".csv", sep = "") )


## B°) Evaporation
clim_data <- rTubeDB::query_timeseries(
      tubedb,
      plot = unique(table.sensors[table.sensors$Sensors == "evaporation","EP"]),
      sensor = "evaporation",
      datetimeFormat = 'character',
      start = s,
      end = e,
      aggregation = t,
      quality = 'empirical',
      interpolated = TRUE,
      quality_counter = TRUE,
      colPlot = TRUE,
      colYear = TRUE,
      colDayOfYear = TRUE,
      colMonth = TRUE,
      colDay = TRUE,
      casted = FALSE,
) 
# Check 'clim_data'
head(clim_data)
summary(clim_data)
str(clim_data)

# Adjust some colnames to your liking
colnames(clim_data)[c(1,2)] <- c("EP","Date")
# Convert to actual 'Date' vector if time is day
clim_data$Date <- as.Date(clim_data$Date)

# Save 
data.table::fwrite(x = clim_data, file = paste("table_subset_grasslands_EPs_interp_daily_evaporation_",s,"-",e,"_",Sys.Date(),".csv", sep = "") )


## C°) Net Radiation - Rn
clim_data <- rTubeDB::query_timeseries(
      tubedb,
      plot = unique(table.sensors[table.sensors$Sensors == "Rn","EP"]),
      sensor = "Rn",
      datetimeFormat = 'character',
      start = s,
      end = e,
      aggregation = t,
      quality = 'empirical',
      interpolated = TRUE,
      quality_counter = TRUE,
      colPlot = TRUE,
      colYear = TRUE,
      colDayOfYear = TRUE,
      colMonth = TRUE,
      colDay = TRUE,
      casted = FALSE,
) 
# Check 'clim_data'
head(clim_data)
summary(clim_data)
str(clim_data)

# Adjust some colnames to your liking
colnames(clim_data)[c(1,2)] <- c("EP","Date")
clim_data$Date <- as.Date(clim_data$Date)

# Save 
data.table::fwrite(x = clim_data, file = paste("table_subset_grasslands_EPs_interp_daily_Rn_",s,"-",e,"_",Sys.Date(),".csv", sep = "") )


## D°) SWUR
clim_data <- rTubeDB::query_timeseries(
      tubedb,
      plot = unique(table.sensors[table.sensors$Sensors == "SWUR","EP"]),
      sensor = "SWUR",
      datetimeFormat = 'character',
      start = s,
      end = e,
      aggregation = t,
      quality = 'empirical',
      interpolated = TRUE,
      quality_counter = TRUE,
      colPlot = TRUE,
      colYear = TRUE,
      colDayOfYear = TRUE,
      colMonth = TRUE,
      colDay = TRUE,
      casted = FALSE,
) 
# Check 'clim_data'
head(clim_data)
summary(clim_data)
str(clim_data)

# Adjust some colnames to your liking
colnames(clim_data)[c(1,2)] <- c("EP","Date")
clim_data$Date <- as.Date(clim_data$Date)

# Save 
data.table::fwrite(x = clim_data, file = paste("table_subset_grasslands_EPs_interp_daily_SWUR_",s,"-",e,"_",Sys.Date(),".csv", sep = "") )


### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------