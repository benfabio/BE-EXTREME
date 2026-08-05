### ------------------------------------------------------------------------------------------------------------

### 31/03/25 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to compare the daily temperature and precipitation data from E-OBS and ERA5 Land measured 
### at the EP's location against the locally-measured climate data.  
### Will help determine how close E-OBS/ERA5 Land data are to our local conditions and how much they need to be corrected
### to be used in a more precise definition of local ECEs.

### R script to create a function that:
### - source() previous R Script#6.3.2 and the functions within
### - apply these functions (compute eval metrics and make plots)
### - analyze outputs

### Last update: 14/04/26 (Re-running evaluate_daily_stat() for writing Methods of the ESSD paper)

### ------------------------------------------------------------------------------------------------------------

# Libraries 
library("dplyr")
library("ggplot2")
library("data.table")
library("purrr")
library("reshape2")
require("scales")
library("ggpubr")
library("lubridate")
library("parallel")

# Source Script#6.3.2 from GitHub
setwd("/home/fbenedetti/R/scripts/")
source("Script#6.3.2_BEO_E-OBS_daily_data_E-OBS+ERA5Land_compare.r") # Have to do it this way for now because my repo is still set to 'private'

### ------------------------------------------------------------------------------------------------------------

### 1°) Use evaluate_daily_stat(). Arguments: var, stat, region, anoms

### 1.A. With anoms == FALSE - based on regular values
regions <- c("HND","SCH","SWA")

# Run evaluate_daily_stat() for all 3 regions and rbind
res <- lapply(regions, function(r) {
        message(paste("Running evaluate_daily_stat for the ",r, sep = ""))
        d <- evaluate_daily_stat(var = "SM_10", stat = "min", region = r, anoms = TRUE, std.anoms = FALSE)
        return(d)
    } # eo FUN
) # eo lapply - regions
# Rbind
table <- dplyr::bind_rows(res)
# dim(table) ; head(table)
rm(res); gc()

# Examine distibution of metrics
# summary(table)
summary(table[table$region == "HND",])
summary(table[table$region == "SCH",])
summary(table[table$region == "SWA",])


### 13-14/04/26: Re-running evaluate_daily_stat() for writing Methods of the ESSD paper
# E-OBS
mean(table$mae_EOBS) ; sd(table$mae_EOBS)
mean(table$mbe_EOBS) ; sd(table$mbe_EOBS)
mean(table$corr_EOBS) ; sd(table$corr_EOBS)
mean(table$rmse_EOBS) ; sd(table$rmse_EOBS)

# ERA5-Land 
mean(table$mae_ERA5) ; sd(table$mae_ERA5)
mean(table$mbe_ERA5) ; sd(table$mbe_ERA5)
mean(table$corr_ERA5) ; sd(table$corr_ERA5)
mean(table$rmse_ERA5) ; sd(table$rmse_ERA5)


## max Ta_200: Both ERA5-Land & E-OBS have very high correlation (0.99 on median for both).
## But ERA5-Land has higher error estimates than E-OBS in all 3 regions: 
## median ERA5-Land MAE: 1.17; 0.63 for E-OBS
## median ERA5-Land RMSE: 1.50; 0.85 for E-OBS

## min Ta_200: Less clear differences with min Ta_200, but E-OBS still has better evaluation metrics.
## E-OBS and ERA5-Land closer in the HND than in the SCH/SWA

## total precipitation: corr coeff much lower than for temperature (not surprising)
## but E-OBS still much better (corr ~ 0.8 instead of 0.6-0.65).
## Errors MUCH higher for ERA5-Land than E-OBS

## max Ta_10: corr coeff between 0.97 & 0.98; MAE between 1.73 and 2.10; RMSE between 2.27 & 2.81
## errors higher in the SWA

## min Ta_10: corr coeff between 0.91 & 0.95; MAE between 1.6 and 2.2; RMSE between 2.5 & 3.3
## errors higher in the SWA too

## max Ts_05: corr coeff between 0.96 & 0.98; MAE between 1.6 and 2.1; RMSE between 2.0 & 2.8
## errors higher in the...SWA! 

## min Ts_05: corr coeff between 0.97 & 0.98; MAE between 1.7 and 2.55; RMSE between 2.0 & 2.9
## errors higher in the SWA again

## max Ts_10: Should be the same as max Ts_05

## max Ts_20: corr coeff between 0.97 & 0.98; MAE between 1.27 & 1.61; RMSE between 1.59 & 2.13
## errors lower in the HND relatve to the other 2 regions

## min Ts_20: corr coeff between 0.97 & 0.983; MAE between 1.35 & 1.72; RMSE between 1.70 & 2.07
## errors lower in the SCH relatve to the other 2 regions this time.
# Very good for deeper soil temperatures I'd say.

### 07/04/25
## max SM_10: corr coeff between 0.53 and 0.71; MAE between 5.8% and 10.8%; RMSE between 7.55 and 13.1
## Errors highest in the SCH and lowest in the HND. corr coeff highest in the HND (0.70), lowest in the SWA (0.51)

## min SM_10: same patterns as max SM_10


### 23/04/25: SM_10 against level_2 ERA5-Land data
## max SM_10: corr coeff between 0.63 and 0.79 (median = 0.72); MAE between 5.04% and 9.42%; RMSE between 6.48 and 11.17
## Again: Errors highest in the SCH (median RMSE = 11.4) and lowest in the HND (6.43).
## corr coeff highest in the HND (0.784), lowest in the SWA (0.65)

### --> level_2 shows lower errors and higher corr coeff against local BEO data!!!


## min SM_10: corr coeff between 0.64 and 0.80; MAE between 4.9% and 9.8%; RMSE between 6.35 and 11.41
## corr coeff highest in the HND (0.79), lowest in the SWA (0.66)

### -------------------------------------------------

### 1.B. With anoms == TRUE - based on anomalies to the monthly mean

res <- lapply(regions, function(r) {
        message(paste("Running evaluate_daily_stat for the ",r, sep = ""))
        d <- evaluate_daily_stat(var = "SM_10", stat = "max", region = r, anoms = TRUE, std.anoms = FALSE)
        return(d)
    } # eo FUN
) # eo lapply - regions
# Rbind
table <- dplyr::bind_rows(res)
# dim(table); summary(table)
head(table)
rm(res); gc()

# Examine distibution of metrics
summary(table)
summary(table[table$region == "HND",])
summary(table[table$region == "SCH",])
summary(table[table$region == "SWA",])


## max Ta_200: Same as before, E-OBS much better than ERA5-land
## min Ta_200: Same as with normal data above

## total precip: Same as with regular data

### --> E-OBS better than ERA5-Land in all cases applicable; anomalies or not

### -------------------------------------------------

### 03/04/25: Compute evaluation metrics based on std.anomalies: Only changes corr coefficients!
res <- lapply(regions, function(r) {
        message(paste("Running evaluate_daily_stat for the ",r, sep = ""))
        d <- evaluate_daily_stat(var = "SM_10", stat = "min", region = r, anoms = FALSE, std.anoms = TRUE)
        return(d)
    } # eo FUN
) # eo lapply - regions
# Rbind
table <- dplyr::bind_rows(res)
# dim(table); summary(table)
head(table)
rm(res); gc()

# Examine distribution of metrics
summary(table)
summary(table[table$region == "HND",])
summary(table[table$region == "SCH",])
summary(table[table$region == "SWA",])

### NOTE: Remember that corr coeff based on std.anomalies allow you to evaluate how 
###       ERA5-Land and E-OBS reproduce temporal variability if observations, not absolute values 
##        (see section 1.A for this)

## max Ta_200: E-OBS > ERA5-Land (0.9786 > 0.931) # very good
## min Ta_200: E-OBS > ERA5-Land (0.9168 > 0.850) # very good

## total precipitation: E-OBS > ERA5-Land (0.712 > 0.611) # good

## max Ta_10: 0.694 to 0.721 in general. Much lower in the SWA (0.547-0.721)
## compared to the other 2 regions (~0.81) --> influence of mountains?
## min Ta_10: 0.45 to 0.72. Again much lower in the SWA (0.37-0.45) than
## in the HND (0.71) and SCH (0.67)
### Overall good fit of ERA5-Land, not amazing but still good.
### Minimum seems harder to model than maximum; Also true for Ta_200! 

## max Ts_05: 0.55 to 0.68 in general. Lower in the SCH (0.50) this time compared to
## the SWA (0.63) and HND (0.65)
## min Ts_05: 0.49 to 0.67 in general. Higher in the HND (0.64) compared to
### the SCH (0.51) and the SWA (0.55)

## max Ts_10: 0.485 to 0.60 in general. Lower in the SCH again (0.44) relative to
## the HND (0.58) and the SWA (0.56)
## min Ts_10: 0.435 to 0.60 in general. Lower in the SCH again (0.44)

### --> corr coeff based on std.anomalies quite lower than for raw values 
### --> temporal variability harder to model than absolute values

### 07/04/25
## max Ts_10: corr coeff between and in general 0.05 to 0.15 - VERY low compared to raw values.
##            temporal variability a LOT harder to model apparently. Lowest in the SCH.

## min SM_10: Same as max SM_10

### ------------------------------------------------------------------------------------------------------------

### 2°) Make plots wth plot_daily_stat_comparison() (2nd FUN)
# arguments are: region, var, stat, mon.anoms, biases, boxp, histo, scatt, time_series

### 03/04/25: Making distribution plots based on std.anomalies with plot_daily_stat_comparison()
# c("precipitation","Ta_200","Ta_10","Ts_05","Ts_10")
for(v in c("SM_10")) {
    for(s in c("max","min")) {
        for(r in c("HND","SCH","SWA")) {

            plot_daily_stat_comparison(region = r, var = v, stat = s,
                    std.anoms = FALSE, mon.anoms = FALSE, biases = FALSE,
                    boxp = TRUE, histo = TRUE, scatt = TRUE, time_series = TRUE
            )

        } # eo for loop - r
    } # eo for loop - s
} # eo for loop - v

### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------