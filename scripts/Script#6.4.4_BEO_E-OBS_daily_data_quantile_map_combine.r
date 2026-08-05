### ------------------------------------------------------------------------------------------------------------

### 05/05/25 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to combine the EP-level and variable-level TS of each of the four QM strategies with 
### the local measurements.
### Will use these 'final' TS for ECE detection and quantification. Will then compare the ECE metrics  
### stemming from various QM strategies to investigate how these strategies impact ECE detection.

### NOTE: If most QM strategies do not seem to provide drastically different results, why not 
###       use an ensemble approach?

### Recyling script 6.4.3 to write a master FUN that will:
###  - Load the corrected daily TS data from R script 6.4.2 (quantile_mapper outputs)
###  - Load the corresponding local measurements at EP level
###  - Combine them into one final TS (all < 2009-01-01 -> QM outputs; all > 2009-01-01 -> local obs)

### Last update: 21/01/26 (Re-running combine_quantile_maps() for forest precipitation reconstruction)

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

# Dir to save the combined TS in:
file.dir <- "/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/grasslands"

### ------------------------------------------------------------------------------------------------------------

### Master FUN - combine_quantile_maps()
# To test combine_quantile_maps while you're writing it:
# var <- "precipitation"
# stat <- "max"
# region <- "HND"
# method <- "mw"

combine_quantile_maps <- function(var, stat, region, method) {

        #' This function takes four arguments and returns a formatted data.frame:
        #' @param var the climate variable to process (character) - one of the following: 
        #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20", "Ts_50", "precipitation" or "SM_10"
        #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character): 'mean', 'max' or 'min'
        #' @param region Region name (character): 'SCH' or 'HND' or 'SWA'
        #' @param method Which quantile mapping approach to run (character): 'global' or 'monthly' or 'mw' or 'anomalies'
        #' @return A formatted data.frame combining the daily statistics

        message(paste("Loading the modelled and observational TS of ",paste(stat,var, sep = " ")," for the ",region," based on the ",method," QM strategy", sep = ""))

        ## Go to local obs directory and extract TS of interest
        setwd(paste("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily/grasslands/",region, sep = ""))

        ### Overwrite 'stat' if var == "precipitation" 
        if( var == "precipitation" ) {
          stat <- "total"
          files <- dir()[grepl(paste(var, sep = "_"),dir())]
        } else {
          # Check that that the daily stats of your variable are available 
          files <- dir()[grepl(paste(var,stat, sep = "_"),dir())]
          if( length(files) < 1 ) {
              stop(
                paste("!!! ERROR: Missing observed ",paste(var,stat, sep = "_")," file for the ",region, sep = "")
              )
          } # eo if loop 
        } # eo if else loop - var == "precipitation"

        # Read the file containing your observed daily stats
        obs_daily_stat <- read.csv(files, h = T, sep = ",", dec = ".")
        # Adjust colnames
        colnames(obs_daily_stat)[3] <- paste(var,stat, sep = "_")

        # Sanity check
        if( exists("obs_daily_stat") == FALSE ) {
            stop(
              paste("!!! ERROR: Could not load observed ",paste(var,stat, sep = "_")," file for the ",region, sep = "")
            )
        } # eo if loop - sanity check

        ## Go to quantile mapping outputs dir and load corresponding data 
        setwd(paste("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/",method, sep = ""))
        file <- dir()[grepl(paste(stat,var,region, sep = "_"),dir())] # file
        corr_data <- get(load(file))
        # dim(corr_data) ; str(corr_data) ; colnames(corr_data)
        
        # If loop: if method == anomalies, use the PENULTIMATE column, not the last
        if( method == "anomalies" ) {
            colnames(corr_data)[length(corr_data) - 1] <- paste(var,stat, sep = "_")
        } else {
            colnames(corr_data)[length(corr_data)] <- paste(var,stat, sep = "_")
        } # eo if else loop - anomalies

        # To make sure Date format is homogeneous across all 3 tables
        if( class(obs_daily_stat$Date) != "Date" ) {
            obs_daily_stat$Date <- as.Date(obs_daily_stat$Date)
        } # eo if loop

        if( class(corr_data$Date) != "Date" ) {
            corr_data$Date <- as.Date(corr_data$Date)
        } # eo if loop

        # Define names to merge by
        names <- c("EP","Date",paste(var,stat, sep = "_"))

        # Merge using full_join iteratively
        merged_df <- reduce(
                      list(obs_daily_stat[,names],
                           corr_data[,names]
                      ), full_join, by = c("EP","Date")
        ) # eo reduce
        # dim(merged_df) ; str(merged_df)
        # summary(merged_df) # OK
        rm(obs_daily_stat,corr_data)
        gc()

        # Adjust colnames
        colnames(merged_df)[c(3,4)] <- c("obs","corr_model")
        
        # Arrange by EP and Date
        merged_df <- merged_df %>% arrange(EP,Date)

        # Restrict dates only before "2025-01-01"
        merged_df <- merged_df[merged_df$Date < as.Date("2025-01-01"),] # We do not care about 2025 data (for now)
        # summary(merged_df)
        # merged_df[which(merged_df$Date >= "2008-11-01"),][1:500,]

        # Combine 
        # Use mutate() to keep the values of 'corr_model' going beyond 2009-01-01
        # and thus also detect recent ECEs based on corr_model too and validate against 'obs' again
        ### '2009-05-01' because we miss Ta data in the SWA until May 2009
        combined_df <- merged_df %>% 
                        mutate(
                            EP,
                            Date, 
                            final_value = ifelse(Date >= as.Date("2009-05-01"), obs, corr_model) 
                        ) # eo mutate 

        # Check 'combined_df'
        # summary(merged_df) 
        # summary(combined_df)

        # Save in appropriate dir 
        if( exists("combined_df") == FALSE ) {
            stop(
              paste("!!! ERROR: Could not find final TS of ",paste(stat,var, sep = " ")," file for the ",region," based on the ",method," QM strategy", sep = "")
            )
        } else {
            setwd(file.dir)
            paste("Saving the final TS of ",paste(stat,var, sep = " ")," file for the ",region," based on the ",method," QM strategy \n", sep = "")
            save(x = combined_df, file = paste("table_combined_obs+corr_",method,"_",stat,"_",var,"_grasslands_",region,".Rdata", sep = "") )
        } # eo if else loop 

} # eo master FUN - combine_quantile_maps


### Apply combine_quantile_maps() to all cases possible
### takes ~ 40min
for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","precipitation","SM_10")) {
    for(s in c("max","min")) {
        for(r in c("HND","SCH","SWA")) {
            for(m in c('global','monthly','mw','anomalies')) {

                combine_quantile_maps(var = v, stat = s, region = r, method = m)

            } # eo for loop - m
        } # eo for loop - r
    } # eo for loop - s
} # eo for loop - v


### ------------------------------------------------------------------------------------------------------------

### 21/01/26: Re-running combine_quantile_maps() for forest precipitation reconstruction

file.dir <- "/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/forests"

# To test combine_forest_quantile_maps() while you're writing it:
var <- "precipitation"
stat <- "total"
region <- "HND"
method <- "mw"

combine_forest_quantile_maps <- function(var = "precipitation", stat = "total", region, method) {

        #' This function takes four arguments and returns a formatted data.frame:
        #' @param var the climate variable to process (character) - Default = "precipitation"
        #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character) - Default = "total"
        #' @param region Region name (character): 'SCH' or 'HND' or 'SWA'
        #' @param method Which quantile mapping approach to run (character): 'global' or 'monthly' or 'mw' or 'anomalies'
        #' @return A formatted data.frame combining the daily statistics

        message(paste("Loading the modelled and observational TS of ",paste(stat,var, sep = " ")," for the ",region," based on the ",method," QM strategy", sep = ""))

        ## Go to local obs directory and extract TS of interest
        setwd(paste("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily/forests/",region, sep = ""))

        ### Overwrite 'stat' if var == "precipitation" 
        if( var == "precipitation" ) {
          stat <- "total"
          files <- dir()[grepl(paste(var, sep = "_"),dir())]
        } else {
          # Check that that the daily stats of your variable are available 
          files <- dir()[grepl(paste(var,stat, sep = "_"),dir())]
          if( length(files) < 1 ) {
              stop(
                paste("!!! ERROR: Missing observed ",paste(var,stat, sep = "_")," file for the ",region, sep = "")
              )
          } # eo if loop 
        } # eo if else loop - var == "precipitation"

        # Read the file containing your observed daily stats
        obs_daily_stat <- read.csv(files, h = T, sep = ",", dec = ".")
        # Adjust colnames
        colnames(obs_daily_stat)[3] <- paste(var,stat, sep = "_")

        # Sanity check
        if( exists("obs_daily_stat") == FALSE ) {
            stop(
              paste("!!! ERROR: Could not load observed ",paste(var,stat, sep = "_")," file for the ",region, sep = "")
            )
        } # eo if loop - sanity check

        ## Go to quantile mapping outputs dir and load corresponding data 
        setwd(paste("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/",method, sep = ""))
        file <- dir()[grepl("forest",dir())] # file
        file <- file[grepl(region,file)] # file
        corr_data <- get(load(file))
        # dim(corr_data) ; str(corr_data) ; colnames(corr_data)
        
        # If loop: if method == anomalies, use the PENULTIMATE column, not the last
        if( method == "anomalies" ) {
            colnames(corr_data)[length(corr_data) - 1] <- paste(var,stat, sep = "_")
        } else {
            colnames(corr_data)[length(corr_data)] <- paste(var,stat, sep = "_")
        } # eo if else loop - anomalies

        # To make sure Date format is homogeneous across all 3 tables
        if( class(obs_daily_stat$Date) != "Date" ) {
            obs_daily_stat$Date <- as.Date(obs_daily_stat$Date)
        } # eo if loop

        if( class(corr_data$Date) != "Date" ) {
            corr_data$Date <- as.Date(corr_data$Date)
        } # eo if loop

        # Define names to merge by
        names <- c("EP","Date",paste(var,stat, sep = "_"))

        # Merge using full_join iteratively
        merged_df <- reduce(
                      list(obs_daily_stat[,names],
                           corr_data[,names]
                      ), full_join, by = c("EP","Date")
        ) # eo reduce
        # dim(merged_df) ; str(merged_df)
        # summary(merged_df) # OK
        rm(obs_daily_stat,corr_data)
        gc()

        # Adjust colnames
        colnames(merged_df)[c(3,4)] <- c("obs","corr_model")
        
        # Arrange by EP and Date
        merged_df <- merged_df %>% arrange(EP,Date)

        # Restrict dates only before "2025-01-01"
        merged_df <- merged_df[merged_df$Date < as.Date("2025-01-01"),] # We do not care about 2025 data (for now)
        # summary(merged_df)
        # merged_df[which(merged_df$Date >= "2008-11-01"),][1:500,]

        # Combine 
        # Use mutate() to keep the values of 'corr_model' going beyond the date of the first precipitation value
        # (obs != NA) and thus also detect recent ECEs based on corr_model too and validate against 'obs' again
        combined_df <- merged_df %>% mutate(EP, Date, final_value = ifelse(Date >= "2009-01-01", obs, corr_model) ) # eo mutate 
        # Check 'combined_df'
        # summary(merged_df) 
        # summary(combined_df)
        # Check
        # combined_df[is.na(combined_df$final_value),]

        # Save in appropriate dir 
        if( exists("combined_df") == FALSE ) {
            stop(
              paste("!!! ERROR: Could not find final TS of ",paste(stat,var, sep = " ")," file for the ",region," based on the ",method," QM strategy", sep = "")
            )
        } else {
            setwd(file.dir)
            paste("Saving the final TS of ",paste(stat,var, sep = " ")," file for the forests of the ",region," based on the ",method," QM strategy \n", sep = "")
            save(x = combined_df, file = paste("table_combined_obs+corr_",method,"_",stat,"_",var,"_forests_",region,".Rdata", sep = "") )
        } # eo if else loop 

} # eo master FUN - combine_forest_quantile_maps

### Apply combine_forest_quantile_maps()
for(r in c("HND","SCH","SWA")) {
    for(m in c("global","monthly","mw","anomalies")) {
        combine_forest_quantile_maps(var = "precipitation", stat = "total", region = r, method = m)
    } # eo for loop - m
} # eo for loop - r

### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------