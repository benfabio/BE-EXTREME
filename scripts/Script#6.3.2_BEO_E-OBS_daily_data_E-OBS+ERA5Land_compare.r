### ------------------------------------------------------------------------------------------------------------

### 25/03/25 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to compare the daily temperature and precipitation data from E-OBS and ERA5 Land measured 
### at the EP's location against the locally-measured climate data.  
### Will help determine how close E-OBS/ERA5 Land data are to our local conditions and how much they need to be corrected
### to be use din a more precise definition of local ECEs.

### R script to create a function that:
###  - Loads daily data from E-OBS (Script#6.0) & ERA Land (Script#5.2) and the daily local Exploratories
###    climate data (Script#2.5)
###  - Joins the three datasets by date, EP and variables
###  - Compares our observed daily data to the E-OBS and ERA Land daily data
###    (i.e., calculates error metrics, corr coeff, etc. see Script#6.3.1)``
###  - Returns a table containing the evaluation metrics etc. for each EP and variable and each 
###    daily statistic (one metric for obs vs. E-OBS and one for obs vs. ERA5 Land)
###  - Develop another FUN for returning plots (histograms, time series etc. see Script#6.3.1)

### !!! NOTE 1: Advice from Adrian Huerta (from the 05/03/25): 
### Assess corr coef and errors with ANOMALIES too (z-scores, or anomalies to the monthly means)
### In the function, maybe think about adding a switch that computes the metrics based on anomalies 
### as a function argument (i.e., anoms == TRUE/FALSE)

### Assess corr coef and errors with STANDARDIZED ANOMALIES too!
## Keeps Physical Meaning of Anomalies
# Standardized anomalies are based on the observed temperature’s mean (μ) and standard deviation (σ),
# ensuring that the anomalies reflect deviations relative to the observed climate.
# Z-scores, on the other hand, would be computed separately for each dataset (observed and modeled),
# meaning the models would be standardized independently, making direct comparisons harder.
# --> Only applies to corr coeff computation though

### !!! NOTE 2: Focus on grassland data until you find a proper way to acounnt for the forests'
### microclimatic processes (imapcting local observations but not E-OBS or ERA5 Land data)

### !!! NOTE 3: E-OBS can only be used for daily min/max air temperature (Ta_200) and precipitation

### Last update: 23/04/25 (Evaluating ERA5-Land volumetric_soil_water_layer_2 data aginst SM_10)

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

# 28/03/25: Dir to save plots: 
plot.dir <- "/home/fbenedetti/plots/climate_reconstructions_evaluation"

### ------------------------------------------------------------------------------------------------------------

### Write master FUN to return table containing the evaluation metrics. What argument should the FUN comprise?
### - variable (use labels from local obs: Ta_200, Ta_10, precipitation, Ts_10 etc.)
### - stat (min/max/mean)
### - region (SCH/HND/SWA)
### - anoms (TRUE/FALSE) - should evaluation metrics be computed on anomalies to the monthly mean

### Master FUN - evaluate_daily_stat
# To test evaluate_daily_stat while you're writing it: 
#var <- "SM_10"
#stat <- "max"
#region <- "HND"
#anoms <- FALSE
#std.anoms <- FALSE

evaluate_daily_stat <- function(var, stat, region, anoms, std.anoms) {

        #' This function takes four arguments and returns a formatted data.frame:
        #' @param var the climate variable to process (character) - one of the following: 
        #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20", "Ts_50", "precipitation" or "SM_10"
        #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character): 'mean', 'max' or 'min'
        #' @param region Region name (character): 'SCH' or 'HND' or 'SWA'
        #' @param anoms Switch - Whether evaluation should be computed based on anomalies to the monly mean (BOOLEAN)
        #' @param std.anoms Switch - Whether correlation coeffs should be computed based on standardized anomalies (BOOLEAN)
        #' @return A formatted data.frame combining the daily statistics

        ## Go to local obs directory and extract TS of interes
        setwd(paste("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily/grasslands/",region, sep = ""))

        ### WARNING: precipitation is cumulative (sum per day), so it cannot have mean/max/min
        ### overwrite 'stat' if var == "precipitation" 
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

        ### If std.anoms == TRUE, then anoms shoudl always be == FALSE
        if( std.anoms ) {
          anoms <- FALSE
        } # eo if else loop - zscores vs. anoms

        # Read the file containing your observed daily stats
        obs_daily_stat <- read.csv(files, h = T, sep = ",", dec = ".")
        # dim(obs_daily_stat); str(obs_daily_stat); summary(obs_daily_stat)

        # Sanity check
        if( exists("obs_daily_stat") == FALSE ) {
            stop(
              paste("!!! ERROR: Could not load observed ",paste(var,stat, sep = "_")," file for the ",region, sep = "")
            )
        } # eo if loop - sanity check

        # Message - depends on variable
        if( var %in% c("Ta_200","precipitation") ) {

            ## Go to E-OBS and ERA5-Land dirs and load their daily data too
            message(paste("Comparing observed daily ",stat," ",var," for the ",region," against E-OBS and ERA5-Land data", sep = ""))
            
            ## Go to E-OBS dir and load corresponding dataset
            setwd("/home/fbenedetti/E-OBS/Explos")
            if( var == "Ta_200" & stat == "max" ) {
                EOBS_daily_stat <- get(load("table_daily_max_temp_all_EPs_E-OBS_1950-2024_09.01.25.Rdata"))
            } else if( var == "Ta_200" & stat == "min" ) {
                EOBS_daily_stat <- get(load("table_daily_min_temp_all_EPs_E-OBS_1950-2024_09.01.25.Rdata"))
            } else if( var == "precipitation" ) {
                EOBS_daily_stat <- get(load("table_daily_precip_all_EPs_E-OBS_1950-2024_09.01.25.Rdata"))
            } # eo if else loop - var

            # Sanity check
            if( exists("EOBS_daily_stat") == FALSE ) {
              stop(
                paste("!!! ERROR: Could not load E-OBS ",paste(var,stat, sep = "_")," data", sep = "")
              )
            } # eo if loop - sanity check

            # Subset region fo interest in 'EOBS_daily_stat' (contains all EPs be default - my bad)
            if( region == "SCH" ) {
                # identify EPs to subset & subt them from 'EOBS_daily_stat'
                plots2subset <- unique(EOBS_daily_stat$EP)
                plots2subset <- plots2subset[grepl("S",plots2subset)] # these are G+W
                plots2subset <- plots2subset[grepl("G",plots2subset)] # these are G only
                # subset
                sub_EOBS_daily_stat <- EOBS_daily_stat[EOBS_daily_stat$EP %in% plots2subset,] # unique(sub_EOBS_daily_stat$EP)
                # delete 'EOBS_daily_stat' - no use anymore
                rm(EOBS_daily_stat,plots2subset); gc()
            } else if( region == "HND" ) {
                # identify EPs to subset & subt them from 'EOBS_daily_stat'
                plots2subset <- unique(EOBS_daily_stat$EP)
                plots2subset <- plots2subset[grepl("H",plots2subset)] # these are G+W
                plots2subset <- plots2subset[grepl("G",plots2subset)] # these are G only
                # subset
                sub_EOBS_daily_stat <- EOBS_daily_stat[EOBS_daily_stat$EP %in% plots2subset,] # unique(sub_EOBS_daily_stat$EP)
                # delete 'EOBS_daily_stat' - no use anymore
                rm(EOBS_daily_stat,plots2subset); gc()
            } else if( region == "SWA" ) {
                # identify EPs to subset & subt them from 'EOBS_daily_stat'
                plots2subset <- unique(EOBS_daily_stat$EP)
                plots2subset <- plots2subset[grepl("A",plots2subset)] # these are G+W
                plots2subset <- plots2subset[grepl("G",plots2subset)] # these are G only
                # subset
                sub_EOBS_daily_stat <- EOBS_daily_stat[EOBS_daily_stat$EP %in% plots2subset,] # unique(sub_EOBS_daily_stat$EP)
                # delete 'EOBS_daily_stat' - no use anymore
                rm(EOBS_daily_stat,plots2subset); gc()
            } # eo if else loop - var
            # dim(sub_EOBS_daily_stat); summary(sub_EOBS_daily_stat)
            
            # Re-name to avoid more of else loops later
            colnames(sub_EOBS_daily_stat) <- c("EP","Date","value")

            ## Go to ERA5-Land dir and load corresponding dataset
            setwd(paste("/home/fbenedetti/ERA5-Land-DEU-processed/daily/grassland/",region, sep = ""))
            # NOTE: Beware ERA5-Land variable do not follow the same convention as your observed data
            # --> Need to create a key to associate both nomenclatures (e.g., Ta_200 = 2m_temperature, etc.)
            # ECMWF's soil layers: 
            # Layer 1: 0 - 7 cm
            # Layer 2: 7 - 28 cm
            # Layer 3: 28 - 100 cm
            # Layer 4: 100 - 289 cm
            if(var == "Ta_200") {
                era5_var <- "2m_temperature"
            } else if(var == "Ta_10") {
                era5_var <- "skin_temperature" # in ECMWF, skin temperature is temp. at the surface - Ta_10 is the closest to it
            } else if(var == "precipitation") {
                era5_var <- "precipitation"
            } else if(var == "Ts_05") {
                era5_var <- "soil_temperature_level_1"
            } else if(var == "Ts_10") {
                era5_var <- "soil_temperature_level_1"
            } else if(var == "Ts_20") {
                era5_var <- "soil_temperature_level_2"
            } else if(var == "Ts_50") {
                era5_var <- "soil_temperature_level_3"
            } else if(var == "SM_10") {
                era5_var <- "volumetric_soil_water_layer_2" # test with volumetric_soil_water_layer_2 too
            } # eo if else loop

            # Sanity check
            if( exists("era5_var") == FALSE ) {
              stop(
                paste("!!! ERROR: Could not find matching ERA5-Land ",paste(stat,var, sep = " ")," data", sep = "")
              )
            } # eo if loop - sanity check
            
            # Identify files of interest and load them
            era5_files <- dir()[grepl(era5_var,dir())] # era5_files
            era5_data <- lapply(era5_files, function(f) {
                      d <- get(load(f))
                      return(d)
                } # eo fun in lapply
            ) # eo lapply - era5_files
            # Bind in a data.frame/tibble
            era5_ddf <- dplyr::bind_rows(era5_data)
            rm(era5_data,era5_files); gc()

            # Keep 'stat' of interest (discard the others)
            era5_ddf <- era5_ddf[,c("region","system","plot","date","day","month","year",stat)]

            # For each EP: match the daily data that overlap across all 3 sources and combine in a single ddf
            # Rename columns so names match across 3 data.frames
            colnames(obs_daily_stat)[3] <- paste(var,stat, sep = "_")
            colnames(sub_EOBS_daily_stat)[3] <- paste(var,stat, sep = "_")
            colnames(era5_ddf)[c(3,4,length(era5_ddf))] <- c("EP","Date",paste(var,stat, sep = "_"))
            names <-  c("EP","Date",paste(var,stat, sep = "_")) # vector of colnames to join the 3 ddf by

            # To make sure Date format is homogeneous across all 3 tables
            if( class(obs_daily_stat$Date) != "Date" ) {
                obs_daily_stat$Date <- as.Date(obs_daily_stat$Date)
            } # eo if loop

            if( class(era5_ddf$Date) != "Date" ) {
                era5_ddf$Date <- as.Date(era5_ddf$Date)
            } # eo if loop

            if( class(sub_EOBS_daily_stat$Date) != "Date" ) {
                sub_EOBS_daily_stat$Date <- as.Date(sub_EOBS_daily_stat$Date)
            } # eo if loop

            # Merge using full_join iteratively
            merged_df <- reduce(
                            list(obs_daily_stat[,names],
                                 sub_EOBS_daily_stat[,names],
                                 era5_ddf[,names]
                            ), full_join, by = c("EP","Date")
            ) # eo reduce
            # dim(merged_df) ; str(merged_df)
            # summary(merged_df) # OK
            rm(obs_daily_stat,sub_EOBS_daily_stat,era5_ddf)
            gc()

            # Adjust colnames
            colnames(merged_df)[c(3:5)] <- c("obs","E_OBS","ERA5_Land") 
            # Remove NAs for evaluation
            merged_df <- na.omit(merged_df)

            # Vector of EP IDs - use to subset 'merged_df'
            plots <- unique(merged_df$EP)

            ### SWITCH: anoms == TRUE/FALSE
            ### If anoms == TRUE --> compute metrics based on anomalies to the monthly mean
            ### If anoms == FALSE --> compute metrics on normal data

            ### Calculate evaluation metrics: mean bias, corr coef, RMSE, mean absolute bias (MAE)
            ### for each EP separately!
            # 1 = obs vs. ERA5-Land
            # 2 = obs vs. E-OBS

            if( anoms ) {

                # Compute monthly means for each data source
                merged_df <- merged_df %>%
                  mutate(month = format(Date, "%Y-%m")) %>%
                  group_by(EP,month) %>%
                  mutate(
                    mon_mean_obs = mean(obs, na.rm = TRUE),
                    mon_mean_EOBS = mean(E_OBS, na.rm = TRUE),
                    mon_mean_ERA5 = mean(ERA5_Land, na.rm = TRUE)
                  ) %>% 
                  ungroup()
                # summary(merged_df)

                # Compute monthly anomalies
                merged_df <- merged_df %>%
                  mutate(
                    mon_anom_obs = obs - mon_mean_obs,
                    mon_anom_EOBS = E_OBS - mon_mean_EOBS,
                    mon_anom_ERA5 = ERA5_Land - mon_mean_ERA5
                  )
                # summary(merged_df)
            
                # Subset merged_df per EP and return evaluation metrics
                eval_metrics <- lapply(plots, function(p) {
                      # p <- plots[13]
                      message(paste("Computing evaluation metrics for EP: ",p," based on monthly anomalies", sep = ""))
                      sub_merged_df <- merged_df[merged_df$EP == p,]
                      ## Mean bias error (MBE)
                      mbe_1 <- mean(sub_merged_df$mon_anom_ERA5 - sub_merged_df$mon_anom_obs) # mbe_1
                      mbe_2 <- mean(sub_merged_df$mon_anom_EOBS - sub_merged_df$mon_anom_obs) # mbe_2
                      ## Mean absolute error (MAE)
                      mae_1 <- mean(abs(sub_merged_df$mon_anom_ERA5 - sub_merged_df$mon_anom_obs)) # mae_1
                      mae_2 <- mean(abs(sub_merged_df$mon_anom_EOBS - sub_merged_df$mon_anom_obs)) # mae_2
                      ## Root mean square error (RMSE)
                      rmse_1 <- sqrt(mean((sub_merged_df$mon_anom_ERA5 - sub_merged_df$mon_anom_obs)^2, na.rm = TRUE)) # rmse_1
                      rmse_2 <- sqrt(mean((sub_merged_df$mon_anom_EOBS - sub_merged_df$mon_anom_obs)^2, na.rm = TRUE)) # rmse_2
                      ## Corr coeff
                      corr_coeff_1 <- cor(sub_merged_df$mon_anom_obs, sub_merged_df$mon_anom_ERA5, use = "complete.obs") # corr_coeff_1
                      corr_coeff_2 <- cor(sub_merged_df$mon_anom_obs, sub_merged_df$mon_anom_EOBS, use = "complete.obs") # corr_coeff_2
                      
                      # Return evaluation metrics in a data.frame
                      eval_metrics <- data.frame(EP = p,
                            mbe_ERA5 = mbe_1, mbe_EOBS = mbe_2,
                            mae_ERA5 = mae_1, mae_EOBS = mae_2,
                            rmse_ERA5 = rmse_1, rmse_EOBS = rmse_2,
                            corr_ERA5 = corr_coeff_1, corr_EOBS = corr_coeff_2
                      ) # eo ddf

                      # Discard subset
                      rm(sub_merged_df)
                      # Return
                      return(eval_metrics)

                  } # eo FUN - plots

                ) # eo lapply - eval_metrics

                # Rbind
                table_eval_metrics <- dplyr::bind_rows(eval_metrics)
                rm(eval_metrics); gc()
                # summary(table_eval_metrics) # quick check

                # Return and clean
                table_eval_metrics$region <- region
                # table_eval_metrics$system <-  ### only grasslands for now
                table_eval_metrics$variable <- var
                table_eval_metrics$stat <- stat
            
                message(paste("Returning evaluation metrics for observed daily ",stat," ",var,
                  " from E-OBS and ERA5-Land data for the ",region," - based on monthly anomalies", sep = ""))
                return(table_eval_metrics)
                rm(merged_df); gc()

            } else if( std.anoms ) {

                # Compute standardized anomalies
                merged_df <- merged_df %>%
                  mutate(month = format(Date, "%Y-%m")) %>%
                  group_by(EP,month) %>%
                  mutate(
                    SA_obs = (obs - mean(obs)) / sd(obs),
                    SA_EOBS = (E_OBS - mean(obs)) / sd(obs),
                    SA_ERA5_Land = (ERA5_Land - mean(obs)) / sd(obs)
                  ) %>% 
                  ungroup()
                # summary(merged_df)

                # Subset merged_df per EP and return evaluation metrics based on z-scores
                eval_metrics <- lapply(plots, function(p) {
      
                      message(paste("Computing evaluation metrics for EP: ",p," based on Z-SCORES", sep = ""))
                      sub_merged_df <- merged_df[merged_df$EP == p,]
                      ## Mean bias error (MBE)
                      mbe_1 <- mean(sub_merged_df$ERA5_Land - sub_merged_df$obs) # mbe_1
                      mbe_2 <- mean(sub_merged_df$E_OBS - sub_merged_df$obs) # mbe_2
                      ## Mean absolute error (MAE)
                      mae_1 <- mean(abs(sub_merged_df$ERA5_Land - sub_merged_df$obs)) # mae_1
                      mae_2 <- mean(abs(sub_merged_df$E_OBS - sub_merged_df$obs)) # mae_2
                      ## Root mean square error (RMSE)
                      rmse_1 <- sqrt(mean((sub_merged_df$ERA5_Land - sub_merged_df$obs)^2, na.rm = TRUE)) # rmse_1
                      rmse_2 <- sqrt(mean((sub_merged_df$E_OBS - sub_merged_df$obs)^2, na.rm = TRUE)) # rmse_2
                      ## Corr coeff based on standardized anomalies 
                      corr_coeff_1 <- cor(sub_merged_df$SA_obs, sub_merged_df$SA_ERA5_Land, use = "complete.obs") # corr_coeff_1
                      corr_coeff_2 <- cor(sub_merged_df$SA_obs, sub_merged_df$SA_EOBS, use = "complete.obs") # corr_coeff_2
                      
                      # Return evaluation metrics in a data.frame
                      eval_metrics <- data.frame(EP = p,
                          mbe_ERA5 = mbe_1, mbe_EOBS = mbe_2,
                          mae_ERA5 = mae_1, mae_EOBS = mae_2,
                          rmse_ERA5 = rmse_1, rmse_EOBS = rmse_2,
                          corr_ERA5 = corr_coeff_1, corr_EOBS = corr_coeff_2
                      ) # eo ddf

                      # Discard subset
                      rm(sub_merged_df)
                      # Return
                      return(eval_metrics)

                  } # eo FUN - plots

                ) # eo lapply - eval_metrics

                # Rbind
                table_eval_metrics <- dplyr::bind_rows(eval_metrics)
                rm(eval_metrics); gc()
                # summary(table_eval_metrics)

                # Return and clean
                table_eval_metrics$region <- region
                # table_eval_metrics$system <-  ### only grasslands for now
                table_eval_metrics$variable <- var
                table_eval_metrics$stat <- stat
            
                message(paste("Returning evaluation metrics for observed daily ",stat," ",var,
                  " from E-OBS and ERA5-Land data for the ",region," - based on standardized anomalies", sep = ""))
                return(table_eval_metrics)
                rm(merged_df); gc()

            } else {

                # Subset merged_df per EP and return evaluation metrics
                eval_metrics <- lapply(plots, function(p) {

                      message(paste("Computing evaluation metrics for EP: ",p, sep = ""))
                      # subset
                      sub_merged_df <- merged_df[merged_df$EP == p,]
                      ## Mean bias error (MBE)
                      mbe_1 <- mean(sub_merged_df$ERA5_Land - sub_merged_df$obs) # mbe_1
                      mbe_2 <- mean(sub_merged_df$E_OBS - sub_merged_df$obs) # mbe_2
                      ## Mean absolute error (MAE)
                      mae_1 <- mean(abs(sub_merged_df$ERA5_Land - sub_merged_df$obs)) # mae_1
                      mae_2 <- mean(abs(sub_merged_df$E_OBS - sub_merged_df$obs)) # mae_2
                      ## Root mean square error (RMSE)
                      rmse_1 <- sqrt(mean((sub_merged_df$ERA5_Land - sub_merged_df$obs)^2, na.rm = TRUE)) # rmse_1
                      rmse_2 <- sqrt(mean((sub_merged_df$E_OBS - sub_merged_df$obs)^2, na.rm = TRUE)) # rmse_2
                      ## Corr coeff
                      corr_coeff_1 <- cor(sub_merged_df$obs, sub_merged_df$ERA5_Land, use = "complete.obs") # corr_coeff_1
                      corr_coeff_2 <- cor(sub_merged_df$obs, sub_merged_df$E_OBS, use = "complete.obs") # corr_coeff_2
                      
                      # Return evaluation metrics in a data.frame
                      eval_metrics <- data.frame(EP = p,
                            mbe_ERA5 = mbe_1, mbe_EOBS = mbe_2,
                            mae_ERA5 = mae_1, mae_EOBS = mae_2,
                            rmse_ERA5 = rmse_1, rmse_EOBS = rmse_2,
                            corr_ERA5 = corr_coeff_1, corr_EOBS = corr_coeff_2
                      ) # eo ddf

                      # Discard subset
                      rm(sub_merged_df)
                      # Return
                      return(eval_metrics)

                  } # eo FUN - plots

                ) # eo lapply - eval_metrics

                # Rbind
                table_eval_metrics <- dplyr::bind_rows(eval_metrics)
                rm(eval_metrics); gc()
                # summary(table_eval_metrics) # quick check

                # Return and clean
                table_eval_metrics$region <- region
                # table_eval_metrics$system <-  ### only grasslands for now
                table_eval_metrics$variable <- var
                table_eval_metrics$stat <- stat
            
                message(paste("Returning evaluation metrics for observed daily ",stat," ",var,
                  " from E-OBS and ERA5-Land data for the ",region," - based on the raw data", sep = ""))
                return(table_eval_metrics)
                rm(merged_df); gc()

            } # eo if else loop - anoms == T

        } else {

            message(paste("Comparing observed daily ",stat," ",var," for the ",region," against ERA5-Land data", sep = ""))

            ## Go to ERA5-Land dir and load corresponding dataset
            setwd(paste("/home/fbenedetti/ERA5-Land-DEU-processed/daily/grassland/",region, sep = ""))
            # see line above for this part
            if(var == "Ta_200") {
                era5_var <- "2m_temperature"
            } else if(var == "Ta_10") {
                era5_var <- "skin_temperature" 
            } else if(var == "precipitation") {
                era5_var <- "precipitation"
            } else if(var == "Ts_05") {
               era5_var <- "soil_temperature_level_1"
            } else if(var == "Ts_10") {
              era5_var <- "soil_temperature_level_1"
            } else if(var == "Ts_20") {
              era5_var <- "soil_temperature_level_2"
            } else if(var == "Ts_50") {
              era5_var <- "soil_temperature_level_3"
            } else if(var == "SM_10") {
              era5_var <- "volumetric_soil_water_layer_2"
            } # eo if else loop

            # Sanity check
            if( exists("era5_var") == FALSE ) {
              stop(
                paste("!!! ERROR: Could not find matching ERA5-Land ",paste(stat,var, sep = " ")," data", sep = "")
              )
            } # eo if loop - sanity check
            
            # Identify files of interest and load them
            era5_files <- dir()[grepl(era5_var,dir())] # era5_files
            era5_data <- lapply(era5_files, function(f) {
                      d <- get(load(f))
                      return(d)
                } # eo fun in lapply
            ) # eo lapply - era5_files
            era5_ddf <- dplyr::bind_rows(era5_data)
            rm(era5_data,era5_files); gc()

            era5_ddf <- era5_ddf[,c("region","system","plot","date","day","month","year",stat)]

            # For each EP: match the daily data that overlap across all 3 sources and combine in a single ddf
            # Rename columns so names match across 3 data.frames
            colnames(obs_daily_stat)[3] <- paste(var,stat, sep = "_")
            colnames(era5_ddf)[c(3,4,length(era5_ddf))] <- c("EP","Date",paste(var,stat, sep = "_"))
            names <-  c("EP","Date",paste(var,stat, sep = "_")) # vector of colnames to join the 3 ddf by

            # To make sure Date format is homogeneous across all 3 tables
            if( class(obs_daily_stat$Date) != "Date" ) {
                obs_daily_stat$Date <- as.Date(obs_daily_stat$Date)
            } # eo if loop

            if( class(era5_ddf$Date) != "Date" ) {
                era5_ddf$Date <- as.Date(era5_ddf$Date)
            } # eo if loop

            # Merge using full_join iteratively
            merged_df <- reduce(
                            list(
                              obs_daily_stat[,names],
                              era5_ddf[,names]
                            ), full_join, by = c("EP","Date")
            ) # eo reduce

            # Adjust colnames
            colnames(merged_df)[c(3,4)] <- c("obs","ERA5_Land") 
            # Remove NAs for evaluation
            merged_df <- na.omit(merged_df)

            ### Calculate evaluation metrics: mean bias, corr coef, RMSE, mean absolute bias (MAE)
            ### for each EP separately
            # Vector of EP IDs - use to subset 'merged_df'
            plots <- unique(merged_df$EP)

            ### SWITCH: anoms == TRUE/FALSE
            ### If anoms == TRUE --> compute metrics based on anomalies to the monthly mean
            ### If anoms == FALSE --> compute metrics on normal data

            ### Calculate evaluation metrics: mean bias, corr coef, RMSE, mean absolute bias (MAE)
            ### for each EP separately!
            # 1 = obs vs. ERA5-Land
            # 2 = obs vs. E-OBS

            if( anoms ) { 

                # Compute monthly means for each data source
                merged_df <- merged_df %>%
                  mutate(month = format(Date, "%Y-%m")) %>%
                  group_by(EP,month) %>%
                  mutate(
                    mon_mean_obs = mean(obs, na.rm = TRUE),
                    mon_mean_ERA5 = mean(ERA5_Land, na.rm = TRUE)
                  ) %>% 
                  ungroup()
                # summary(merged_df)

                # Compute monthly anomalies
                merged_df <- merged_df %>%
                  mutate(
                    mon_anom_obs = obs - mon_mean_obs,
                    mon_anom_ERA5 = ERA5_Land - mon_mean_ERA5
                  )
                # summary(merged_df)

                # Subset merged_df per EP and return evaluation metrics based on monthly anomalies
                eval_metrics <- lapply(plots, function(p) {
                      message(paste("Computing evaluation metrics for EP: ",p, sep = ""))
                      # subset
                      sub_merged_df <- merged_df[merged_df$EP == p,]
                      ## Mean bias error (MBE)
                      mbe_1 <- mean(sub_merged_df$mon_anom_ERA5 - sub_merged_df$mon_anom_obs) # mbe_1
                      ## Mean absolute error (MAE)
                      mae_1 <- mean(abs(sub_merged_df$mon_anom_ERA5 - sub_merged_df$mon_anom_obs)) # mae_1
                      ## Corr coeff
                      corr_coeff_1 <- cor(sub_merged_df$mon_anom_obs, sub_merged_df$mon_anom_ERA5, use = "complete.obs") # corr_coeff_1
                      ## Root mean square error (RMSE)
                      rmse_1 <- sqrt(mean((sub_merged_df$mon_anom_obs - sub_merged_df$mon_anom_ERA5)^2, na.rm = TRUE)) # rmse_1
                      # Return evaluation metrics in a data.frame
                      eval_metrics <- data.frame(EP = p,
                            mbe_ERA5 = mbe_1, mae_ERA5 = mae_1,
                            rmse_ERA5 = rmse_1, corr_ERA5 = corr_coeff_1
                      ) # eo ddf
                      # Discard subset
                      rm(sub_merged_df)
                      # Return
                      return(eval_metrics)
                  } # eo FUN - plots
                ) # eo lapply - eval_metrics
                # Rbind
                table_eval_metrics <- dplyr::bind_rows(eval_metrics)
                rm(eval_metrics); gc()
                # summary(table_eval_metrics) # quick check

                # Return and clean
                table_eval_metrics$region <- region
                # table_eval_metrics$system <-  ### only grasslands for now
                table_eval_metrics$variable <- var
                table_eval_metrics$stat <- stat

                message(paste("Returning evaluation metrics for observed daily ",stat," ",var,
                  " from ERA5-Land data for the ",region," based on monthly anomalies", sep = ""))
                return(table_eval_metrics)
                rm(merged_df); gc()
            
            } else if( std.anoms ) {

                # Compute standardized anomalies
                merged_df <- merged_df %>%
                  mutate(month = format(Date, "%Y-%m")) %>%
                  group_by(EP,month) %>%
                  mutate(
                    SA_obs = (obs - mean(obs)) / sd(obs),
                    SA_ERA5_Land = (ERA5_Land - mean(obs)) / sd(obs)
                  ) %>% 
                  ungroup()
                # summary(merged_df)

                # Subset merged_df per EP and return evaluation metrics based on z-scores
                eval_metrics <- lapply(plots, function(p) {

                      message(paste("Computing evaluation metrics for EP: ",p," based on Z-SCORES", sep = ""))
                      sub_merged_df <- merged_df[merged_df$EP == p,]
                      ## Mean bias error (MBE)
                      mbe_1 <- mean(sub_merged_df$ERA5_Land - sub_merged_df$obs)
                      ## Mean absolute error (MAE)
                      mae_1 <- mean(abs(sub_merged_df$ERA5_Land - sub_merged_df$obs))
                      ## Root mean square error (RMSE)
                      rmse_1 <- sqrt(mean((sub_merged_df$ERA5_Land - sub_merged_df$obs)^2, na.rm = TRUE)) # rmse_1
                      ## Corr coeff based on standardized anomalies 
                      corr_coeff_1 <- cor(sub_merged_df$SA_obs, sub_merged_df$SA_ERA5_Land, use = "complete.obs") # corr_coeff_1
                      
                      # Return evaluation metrics in a data.frame
                      eval_metrics <- data.frame(EP = p,
                            mbe_ERA5 = mbe_1, mae_ERA5 = mae_1,
                            rmse_ERA5 = rmse_1, corr_ERA5 = corr_coeff_1
                      ) # eo ddf

                      # Discard subset
                      rm(sub_merged_df)
                      # Return
                      return(eval_metrics)

                  } # eo FUN - plots

                ) # eo lapply - eval_metrics

                # Rbind
                table_eval_metrics <- dplyr::bind_rows(eval_metrics)
                rm(eval_metrics); gc()
                # summary(table_eval_metrics)

                # Return and clean
                table_eval_metrics$region <- region
                # table_eval_metrics$system <-  ### only grasslands for now
                table_eval_metrics$variable <- var
                table_eval_metrics$stat <- stat
            
                message(paste("Returning evaluation metrics for observed daily ",stat," ",var,
                  " from E-OBS and ERA5-Land data for the ",region," - based on standardized anomalies", sep = ""))
                return(table_eval_metrics)
                rm(merged_df); gc()
            
            } else {

                # Subset merged_df per EP and return evaluation metrics
                eval_metrics <- lapply(plots, function(p) {
                      message(paste("Computing evaluation metrics for EP: ",p, sep = ""))
                      # subset
                      sub_merged_df <- merged_df[merged_df$EP == p,]
                      ## Mean bias error (MBE)
                      mbe_1 <- mean(sub_merged_df$ERA5_Land - sub_merged_df$obs)
                      ## Mean absolute error (MAE)
                      mae_1 <- mean(abs(sub_merged_df$ERA5_Land - sub_merged_df$obs))
                      ## Root mean square error (RMSE)
                      rmse_1 <- sqrt(mean((sub_merged_df$ERA5_Land - sub_merged_df$obs)^2, na.rm = TRUE)) 
                      ## Corr coeff
                      corr_coeff_1 <- cor(sub_merged_df$obs, sub_merged_df$ERA5_Land, use = "complete.obs")
                      # Return evaluation metrics in a data.frame
                      eval_metrics <- data.frame(EP = p,
                            mbe_ERA5 = mbe_1, mae_ERA5 = mae_1,
                            rmse_ERA5 = rmse_1, corr_ERA5 = corr_coeff_1
                      ) # eo ddf
                      # Discard subset
                      rm(sub_merged_df)
                      # Return
                      return(eval_metrics)
                  } # eo FUN - plots
                ) # eo lapply - eval_metrics
                # Rbind
                table_eval_metrics <- dplyr::bind_rows(eval_metrics)
                rm(eval_metrics); gc()
                # summary(table_eval_metrics) # quick check

                # Return and clean
                table_eval_metrics$region <- region
                # table_eval_metrics$system <-  ### only grasslands for now
                table_eval_metrics$variable <- var
                table_eval_metrics$stat <- stat

                message(paste("Returning evaluation metrics for observed daily ",stat," ",var,
                  " from ERA5-Land data for the ",region," based on raw data", sep = ""))
                return(table_eval_metrics)
                rm(merged_df); gc()

            } # eo if else loop - anoms

        } # eo if else loop

} # eo master FUN - evaluate_daily_stat


### 27/03/25: Write another FUN to make plots such as histograms, time series or scatter plots
### The FUN will understand its arguments as 'switches' to display certain type of plots 
### It will use the same plots arguments to make them based on biases (model - obs) and
### monthly anomalies too thanks to the 'mon.anoms' and 'biases' arguments

### Master plotting FUN
# To test while you're writing it:
#region <- "SCH"
#var <- "SM_10"
#stat <- "min"
#mon.anoms <- TRUE
#boxp <- TRUE
#histo <- TRUE
#scatt <- TRUE
#time_series <- TRUE
#biases <- TRUE 
#std.anoms <- TRUE

plot_daily_stat_comparison <- function(region, var, stat, mon.anoms, std.anoms, biases, boxp, histo, scatt, time_series) {

    #' This function takes 9 arguments and returns a formatted data.frame:
    #' @param var the climate variable to process (character) - one of the following: 
    #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20", "Ts_50", "precipitation" or "SM_10"
    #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character): 'mean', 'max' or 'min' 
    #' @param region Region name (character): 'SCH' or 'HND' or 'SWA'
    #' @param mon.anoms Switch (BOOLEAN) - Whether plots should be made based on anomalies to the monly mean 
    #' @param std.anoms Switch (BOOLEAN) - Whether plots should be made based on standardized anomalies
    #' Better when aiming to evaluate variability
    #' @param biases Switch (BOOLEAN) - Whether plots should be made for the biases (model - obs) too
    #' @param boxp Switch (BOOLEAN) - Whether boxplots should be made
    #' @param histo Switch (BOOLEAN) - Whether gghistograms should be made
    #' @param scatt Switch (BOOLEAN) - Whether scaterr plots should be made
    #' @param time_series Switch (BOOLEAN) - Whether time series plots should be made
    #' @return the plots selected through the arguments above

      ## Go to local obs directory and extract TS of interest
      setwd(paste("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily/grasslands/",region, sep = ""))

      ### WARNING: precipitation is cumulative (sum per day), so it cannot have mean/max/min
      ### --> overwrite 'stat' if var == "precipitation" 
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

      ### WARNING: Do not use monthly anomalies when 'biases == TRUE'
      ### (only plot biases between obs and models based on normal values)
      ### --> overwrite 'mon.anoms' if 'biases == TRUE'
      if( biases ) {
        mon.anoms <- FALSE
      } # eo if loop - biases

      # Same with std.anoms
      if( std.anoms ) {
        biases <- FALSE
        mon.anoms <- FALSE
      } # eo if loop - biases

      # Read the file containing your observed daily stats
      obs_daily_stat <- read.csv(files, h = T, sep = ",", dec = ".")

      # Sanity check
      if( exists("obs_daily_stat") == FALSE ) {
          stop(
            paste("!!! ERROR: Could not load observed ",paste(var,stat, sep = "_")," file for the ",region, sep = "")
          )
      } # eo if loop - sanity check

      # Message - depends on variable
      if( var %in% c("Ta_200","precipitation") ) {

          ## Go to E-OBS and ERA5-Land dirs and load their daily data too
          message(paste("Plotting distribution of observed daily ",stat," ",var," for the ",region," against E-OBS and ERA5-Land data", sep = ""))
            
          ## Go to E-OBS dir and load corresponding dataset
          setwd("/home/fbenedetti/E-OBS/Explos")
          if( var == "Ta_200" & stat == "max" ) {
              EOBS_daily_stat <- get(load("table_daily_max_temp_all_EPs_E-OBS_1950-2024_09.01.25.Rdata"))
          } else if( var == "Ta_200" & stat == "min" ) {
              EOBS_daily_stat <- get(load("table_daily_min_temp_all_EPs_E-OBS_1950-2024_09.01.25.Rdata"))
          } else if( var == "precipitation" ) {
              EOBS_daily_stat <- get(load("table_daily_precip_all_EPs_E-OBS_1950-2024_09.01.25.Rdata"))
          } # eo if else loop - var

          # Sanity check
          if( exists("EOBS_daily_stat") == FALSE ) {
            stop(
              paste("!!! ERROR: Could not load E-OBS ",paste(var,stat, sep = "_")," data", sep = "")
            )
          } # eo if loop - sanity check

          # Subset region fo interest in 'EOBS_daily_stat' (contains all EPs be default - my bad)
          if( region == "SCH" ) {
                # identify EPs to subset & subt them from 'EOBS_daily_stat'
                plots2subset <- unique(EOBS_daily_stat$EP)
                plots2subset <- plots2subset[grepl("S",plots2subset)] # these are G+W
                plots2subset <- plots2subset[grepl("G",plots2subset)] # these are G only
                # subset
                sub_EOBS_daily_stat <- EOBS_daily_stat[EOBS_daily_stat$EP %in% plots2subset,] # unique(sub_EOBS_daily_stat$EP)
                # delete 'EOBS_daily_stat' - no use anymore
                rm(EOBS_daily_stat,plots2subset); gc()
          } else if( region == "HND" ) {
                # identify EPs to subset & subt them from 'EOBS_daily_stat'
                plots2subset <- unique(EOBS_daily_stat$EP)
                plots2subset <- plots2subset[grepl("H",plots2subset)] # these are G+W
                plots2subset <- plots2subset[grepl("G",plots2subset)] # these are G only
                # subset
                sub_EOBS_daily_stat <- EOBS_daily_stat[EOBS_daily_stat$EP %in% plots2subset,] # unique(sub_EOBS_daily_stat$EP)
                # delete 'EOBS_daily_stat' - no use anymore
                rm(EOBS_daily_stat,plots2subset); gc()
          } else if( region == "SWA" ) {
                # identify EPs to subset & subt them from 'EOBS_daily_stat'
                plots2subset <- unique(EOBS_daily_stat$EP)
                plots2subset <- plots2subset[grepl("A",plots2subset)] # these are G+W
                plots2subset <- plots2subset[grepl("G",plots2subset)] # these are G only
                # subset
                sub_EOBS_daily_stat <- EOBS_daily_stat[EOBS_daily_stat$EP %in% plots2subset,] # unique(sub_EOBS_daily_stat$EP)
                # delete 'EOBS_daily_stat' - no use anymore
                rm(EOBS_daily_stat,plots2subset); gc()
          } # eo if else loop - var
          # dim(sub_EOBS_daily_stat); summary(sub_EOBS_daily_stat)
            
          # Re-name to avoid more of else loops later
          colnames(sub_EOBS_daily_stat) <- c("EP","Date","value")

          ## Go to ERA5-Land dir and load corresponding dataset
          setwd(paste("/home/fbenedetti/ERA5-Land-DEU-processed/daily/grassland/",region, sep = ""))
          # NOTE: Beware ERA5-Land variable do not follow the same convention as your observed data
          # --> Need to create a key to associate both nomenclatures (e.g., Ta_200 = 2m_temperature, etc.)
          if(var == "Ta_200") {
                era5_var <- "2m_temperature"
          } else if(var == "Ta_10") {
                era5_var <- "skin_temperature" # in ECMWF, skin temperature is temp. at the surface - Ta_10 is the closest to it
          } else if(var == "precipitation") {
                era5_var <- "precipitation"
          } else if(var == "Ts_05") {
               era5_var <- "soil_temperature_level_1"
          } else if(var == "Ts_10") {
              era5_var <- "soil_temperature_level_1"
          } else if(var == "Ts_20") {
              era5_var <- "soil_temperature_level_2"
          } else if(var == "Ts_50") {
              era5_var <- "soil_temperature_level_3"
          } else if(var == "SM_10") {
              era5_var <- "volumetric_soil_water_layer_2"
          } # eo if else loop

          # Sanity check
          if( exists("era5_var") == FALSE ) {
            stop(
              paste("!!! ERROR: Could not find matching ERA5-Land ",paste(stat,var, sep = " ")," data", sep = "")
            )
          } # eo if loop - sanity check
            
          # Identify files of interest and load them
          era5_files <- dir()[grepl(era5_var,dir())] # era5_files
          era5_data <- lapply(era5_files, function(f) {
                  d <- get(load(f))
                  return(d)
            } # eo fun in lapply
          ) # eo lapply - era5_files
          era5_ddf <- dplyr::bind_rows(era5_data)
          rm(era5_data,era5_files); gc()

          # Keep 'stat' of interest (discard the others)
          era5_ddf <- era5_ddf[,c("region","system","plot","date","day","month","year",stat)]

          # For each EP: match the daily data that overlap across all 3 sources and combine in a single ddf
          # Rename columns so names match across 3 data.frames
          colnames(obs_daily_stat)[3] <- paste(var,stat, sep = "_")
          colnames(sub_EOBS_daily_stat)[3] <- paste(var,stat, sep = "_")
          colnames(era5_ddf)[c(3,4,length(era5_ddf))] <- c("EP","Date",paste(var,stat, sep = "_"))
          names <-  c("EP","Date",paste(var,stat, sep = "_")) # vector of colnames to join the 3 ddf by

          # To make sure Date format is homogeneous across all 3 tables
          if( class(obs_daily_stat$Date) != "Date" ) {
              obs_daily_stat$Date <- as.Date(obs_daily_stat$Date)
          } # eo if loop

          if( class(era5_ddf$Date) != "Date" ) {
              era5_ddf$Date <- as.Date(era5_ddf$Date)
          } # eo if loop

          if( class(sub_EOBS_daily_stat$Date) != "Date" ) {
              sub_EOBS_daily_stat$Date <- as.Date(sub_EOBS_daily_stat$Date)
          } # eo if loop

          # Merge using full_join iteratively
          merged_df <- reduce(
                        list(obs_daily_stat[,names],
                              sub_EOBS_daily_stat[,names],
                              era5_ddf[,names]
                        ), full_join, by = c("EP","Date")
          ) # eo reduce
          # summary(merged_df) # OK
          rm(obs_daily_stat,sub_EOBS_daily_stat,era5_ddf)
          gc()

          # Adjust colnames
          colnames(merged_df)[c(3:5)] <- c("obs","E_OBS","ERA5_Land") 
          # Remove NAs for plotting
          merged_df <- na.omit(merged_df)

          ### SWITCHES:
          ### If anoms == TRUE --> make plots based on anomalies to the monthly mean
          ### If std.anoms == TRUE --> make plots on standardized anomalies
          ### else --> make plots on normal data

          if( mon.anoms ) {

            # Compute monthly means for each data source
            merged_df <- merged_df %>%
                mutate(month = format(Date, "%Y-%m")) %>%
                group_by(EP,month) %>%
                mutate(
                    mon_mean_obs = mean(obs, na.rm = TRUE),
                    mon_mean_EOBS = mean(E_OBS, na.rm = TRUE),
                    mon_mean_ERA5 = mean(ERA5_Land, na.rm = TRUE)
                ) %>% 
                ungroup()

            # Compute monthly anomalies
            merged_df <- merged_df %>%
                mutate(
                    mon_anom_obs = obs - mon_mean_obs,
                    mon_anom_EOBS = E_OBS - mon_mean_EOBS,
                    mon_anom_ERA5 = ERA5_Land - mon_mean_ERA5
                )
            
            ### Save plots in plot.dir
            setwd(plot.dir)
            colors_sources <- c("obs" = "#00AFBB", "E_OBS" = "#e31a1c", "ERA5_Land" = "#E7B800")

            # 1st: boxplots 
            if( boxp ) {
                # Boxplot of stat distribution - needs melting
                m_merged_df <- reshape2::melt(merged_df[,c("EP","Date","mon_anom_obs","mon_anom_EOBS","mon_anom_ERA5")],
                                        id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Re-name source levels so they match the 'colors_sources' palette above; unique(m_merged_df$source)
                m_merged_df$source <- as.character(m_merged_df$source)
                m_merged_df[m_merged_df$source == "mon_anom_obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "mon_anom_EOBS","source"] <- "E_OBS"
                m_merged_df[m_merged_df$source == "mon_anom_ERA5","source"] <- "ERA5_Land"
                # Make the plot
                plot <- ggplot(data = m_merged_df, aes(x = factor(source), y = value, fill = factor(source))) +
                  geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = .2) + 
                  geom_hline(yintercept = 0, linetype = "dashed") + 
                  scale_fill_manual(values = colors_sources) + labs(y = paste(stat,var,sep = " "), x = "") + 
                  theme_bw() + theme(legend.position = "none") +
                  ggtitle(paste("Distribution of ",stat," ",var," monthly anomalies in the ",region, sep = ""))
                # Save plot in dir  
                ggsave(plot = plot, filename = paste("boxplot_evaluation_anomalies_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 5, height = 5)  
            } # eo if loop - boxp

            # 2nd: histograms - also need melting
            if( histo ) {
                # Needs melting too
                m_merged_df <- reshape2::melt(merged_df[,c("EP","Date","mon_anom_obs","mon_anom_EOBS","mon_anom_ERA5")],
                                        id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Re-name source levels so they match the 'colors_sources' palette above; unique(m_merged_df$source)
                m_merged_df$source <- as.character(m_merged_df$source)
                m_merged_df[m_merged_df$source == "mon_anom_obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "mon_anom_EOBS","source"] <- "E_OBS"
                m_merged_df[m_merged_df$source == "mon_anom_ERA5","source"] <- "ERA5_Land"
                # Make the plot
                plot <- gghistogram(m_merged_df, x = "value", add = "mean", rug = FALSE,
                          fill = "source", palette = colors_sources,
                          ylab = "Count", xlab = paste(stat,var,sep = " "),
                          title = paste("Distribution of ",stat," ",var," monthly anomalies in the ",region, sep = "")
                )
                # Save plot in dir  
                ggsave(plot = plot, filename = paste("gghistogram_evaluation_anomalies_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 7, height = 5)  
            } # eo if loop - histo

            # 3rd: scatt - make 2: obs vs E-OBS & obs vs ERA5-Land
            if( scatt ) {
                # Make plot 1
                p1 <- ggplot(merged_df, aes(x = mon_anom_obs, y = mon_anom_EOBS)) + geom_point(alpha = 0.2, colour = "grey50") +
                        geom_abline(intercept = 0, slope = 1, color = "#d53e4f", linetype = "dashed") + 
                        geom_smooth(method = "lm", se = TRUE, colour = "#00AFBB") + 
                        geom_hline(yintercept = 0, linetype = "dashed") + 
                        geom_vline(xintercept = 0, linetype = "dashed") + 
                        labs(title = paste("Local vs E-OBS daily ",stat," ",var,"\nmonthly anomalies", sep = ""),
                            x = "Local measurements anomalies", y = "E-OBS anomalies") +
                        theme_bw()
                # Make plot 2
                p2 <- ggplot(merged_df, aes(x = mon_anom_obs, y = mon_anom_ERA5)) + geom_point(alpha = 0.2, colour = "grey50") +
                        geom_abline(intercept = 0, slope = 1, color = "#d53e4f", linetype = "dashed") + 
                        geom_smooth(method = "lm", se = TRUE, colour = "#00AFBB") + 
                        geom_hline(yintercept = 0, linetype = "dashed") + 
                        geom_vline(xintercept = 0, linetype = "dashed") + 
                        labs(title = paste("Local vs ERA5-Land daily ",stat," ",var,"\nmonthly anomalies", sep = ""),
                            x = "Local measurements anomalies", y = "ERA5-Land anomalies") +
                        theme_bw()
                # Save as panel in dir  
                panel <- ggarrange(p1, p2, align = "hv", nrow = 1, ncol = 2)
                ggsave(plot = panel, filename = paste("scatter_evaluation_anomalies_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 10, height = 5)  
            } # eo if loop - scatt

            # 4th: time_series
            if( time_series ) {
                p <- ggplot(merged_df, aes(x = Date)) + 
                  geom_line(aes(y = mon_anom_obs), linewidth = 1, colour = "#00AFBB", alpha = .25) +
                  geom_line(aes(y = mon_anom_EOBS), linewidth = 1, colour = "#e31a1c", alpha = .25) +
                  geom_line(aes(y = mon_anom_ERA5), linewidth = 1, colour = "#E7B800", alpha = .25) +
                  geom_hline(yintercept = 0, linetype = "dashed") + 
                  labs(title = paste("Time series of ",stat," ",var," monthly anomalies in the ",region," (all EPs)", sep = ""),
                      x = "Date", y = paste(stat,var,sep = " ")) +
                  theme_minimal()
                # Save
                 ggsave(plot = p, filename = paste("time_series_evaluation_anomalies_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 10, height = 4)  
            } # eo if loop - time_series

          } else if( std.anoms ) {
          
            # Compute monthly anomalies
            merged_df <- merged_df %>%
                mutate(
                    SA_obs = (obs - mean(obs)) / sd(obs),
                    SA_EOBS = (E_OBS - mean(obs)) / sd(obs),
                    SA_ERA5_Land = (ERA5_Land - mean(obs)) / sd(obs)
                )
            
            ### Save plots in plot.dir
            setwd(plot.dir)
            colors_sources <- c("obs" = "#00AFBB", "E_OBS" = "#e31a1c", "ERA5_Land" = "#E7B800")

            # 1st: boxplots 
            if( boxp ) {
                # Boxplot of stat distribution - needs melting
                m_merged_df <- reshape2::melt(merged_df[,c("EP","Date","SA_obs","SA_EOBS","SA_ERA5_Land")], id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Re-name to match color palette's labels
                m_merged_df$source <- as.character(m_merged_df$source)
                m_merged_df[m_merged_df$source == "mon_anom_obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "mon_anom_EOBS","source"] <- "E_OBS"
                m_merged_df[m_merged_df$source == "mon_anom_ERA5","source"] <- "ERA5_Land"
                # Make the plot
                plot <- ggplot(data = m_merged_df, aes(x = factor(source), y = value, fill = factor(source))) +
                  geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = .2) + 
                  scale_fill_manual(values = colors_sources) + labs(y = paste(stat,var,sep = " "), x = "") + 
                  theme_bw() + theme(legend.position = "none") +
                  ggtitle(paste("Distribution of ",stat," ",var," in the ",region, sep = ""))
                # Save plot in dir  
                ggsave(plot = plot, filename = paste("boxplot_evaluation_std.anoms_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 5, height = 5)  
            } # eo if loop - boxp

            # 2nd: histograms - also need melting
            if( histo ) {
                # Boxplot of stat distribution - needs melting
                m_merged_df <- reshape2::melt(merged_df[,c("EP","Date","SA_obs","SA_EOBS","SA_ERA5_Land")], id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Re-name to match color palette's labels
                m_merged_df$source <- as.character(m_merged_df$source)
                m_merged_df[m_merged_df$source == "mon_anom_obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "mon_anom_EOBS","source"] <- "E_OBS"
                m_merged_df[m_merged_df$source == "mon_anom_ERA5","source"] <- "ERA5_Land"
                # Make the plot
                plot <- gghistogram(m_merged_df, x = "value", add = "mean", rug = FALSE,
                          fill = "source", palette = colors_sources,
                          ylab = "Count", xlab = paste(stat,var,sep = " "),
                          title = paste("Distribution of ",stat," ",var,"\nstandardized anomalies in the ",region, sep = "")
                )
                # Save plot in dir  
                ggsave(plot = plot, filename = paste("gghistogram_evaluation_std.anoms_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 7, height = 5)  
            } # eo if loop - histo

            # 3rd: scatt - make 2: obs vs E-OBS & obs vs ERA5-Land
            if( scatt ) {
                # Make plot 1
                p1 <- ggplot(merged_df, aes(x = SA_obs, y = SA_EOBS)) + geom_point(alpha = 0.2, colour = "grey50") +
                        geom_abline(intercept = 0, slope = 1, color = "#d53e4f", linetype = "dashed") + 
                        geom_smooth(method = "lm", se = TRUE, colour = "#00AFBB") + 
                        labs(title = paste("Local vs E-OBS daily ",stat," ",var,"\nbased on standardized anomalies", sep = ""),
                            x = "Local measurements", y = "E-OBS") +
                        theme_bw()
                # Make plot 2
                p2 <- ggplot(merged_df, aes(x = SA_obs, y = SA_ERA5_Land)) + geom_point(alpha = 0.2, colour = "grey50") +
                        geom_abline(intercept = 0, slope = 1, color = "#d53e4f", linetype = "dashed") + 
                        geom_smooth(method = "lm", se = TRUE, colour = "#00AFBB") + 
                        labs(title = paste("Local vs ERA5-Land daily ",stat," ",var,"\nbased on standardized anomalies", sep = ""),
                            x = "Local measurements", y = "ERA5-Land") +
                        theme_bw()
                # Save as panel in dir  
                panel <- ggarrange(p1, p2, align = "hv", nrow = 1, ncol = 2)
                ggsave(plot = panel, filename = paste("scatter_evaluation_std.anoms_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 10, height = 5)  
            } # eo if loop - histo

            # 4th: time_series
            if( time_series ) {
                p <- ggplot(merged_df, aes(x = Date)) + 
                  geom_line(aes(y = SA_obs), linewidth = 1, colour = "#00AFBB", alpha = .25) +
                  geom_line(aes(y = SA_EOBS), linewidth = 1, colour = "#e31a1c", alpha = .25) +
                  geom_line(aes(y = SA_ERA5_Land), linewidth = 1, colour = "#E7B800", alpha = .25) +
                  labs(title = paste("Time series of ",stat," ",var,"\nstandardized anomalies in the ",region," (all EPs)", sep = ""),
                      x = "Date", y = paste(stat,var,sep = " ")) +
                  theme_minimal()
                # Save
                 ggsave(plot = p, filename = paste("time_series_evaluation_std.anoms_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 10, height = 4)  
            } # eo if loop - time_series
          
          } else {

            ### Save plots depending on the switches activated as function's arguments: boxp, histo, scatt, time_series
            setwd(plot.dir)
            colors_sources <- c("obs" = "#00AFBB", "E_OBS" = "#e31a1c", "ERA5_Land" = "#E7B800")

            # 1st: boxplots 
            if( boxp ) {
                # Boxplot of stat distribution - needs melting
                m_merged_df <- reshape2::melt(merged_df, id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Re-name source levels so they match the 'colors_sources' palette above; unique(m_merged_df$source)
                m_merged_df$source <- as.character(m_merged_df$source)
                m_merged_df[m_merged_df$source == "SA_obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "SA_ERA5_Land","source"] <- "ERA5_Land"
                m_merged_df[m_merged_df$source == "SA_EOBS","source"] <- "E_OBS"
                # Make the plot
                plot <- ggplot(data = m_merged_df, aes(x = factor(source), y = value, fill = factor(source))) +
                  geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = .2) + 
                  scale_fill_manual(values = colors_sources) + labs(y = paste(stat,var,sep = " "), x = "") + 
                  theme_bw() + theme(legend.position = "none") +
                  ggtitle(paste("Distribution of ",stat," ",var," in the ",region, sep = ""))
                # Save plot in dir  
                ggsave(plot = plot, filename = paste("boxplot_evaluation_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 5, height = 5)  
            } # eo if loop - boxp

            # 2nd: histograms - also need melting
            if( histo ) {
                # Boxplot of stat distribution - needs melting
                m_merged_df <- reshape2::melt(merged_df, id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Re-name source levels so they match the 'colors_sources' palette above; unique(m_merged_df$source)
                m_merged_df$source <- as.character(m_merged_df$source)
                m_merged_df[m_merged_df$source == "SA_obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "SA_ERA5_Land","source"] <- "ERA5_Land"
                m_merged_df[m_merged_df$source == "SA_EOBS","source"] <- "E_OBS"
                # Make the plot
                plot <- gghistogram(m_merged_df, x = "value", add = "mean", rug = FALSE,
                          fill = "source", palette = colors_sources,
                          ylab = "Count", xlab = paste(stat,var,sep = " "),
                          title = paste("Distribution of ",stat," ",var," in the ",region, sep = "")
                )
                # Save plot in dir  
                ggsave(plot = plot, filename = paste("gghistogram_evaluation_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 7, height = 5)  
            } # eo if loop - histo

            # 3rd: scatt - make 2: obs vs E-OBS & obs vs ERA5-Land
            if( scatt ) {
                # Make plot 1
                p1 <- ggplot(merged_df, aes(x = obs, y = E_OBS)) + geom_point(alpha = 0.2, colour = "grey50") +
                        geom_abline(intercept = 0, slope = 1, color = "#d53e4f", linetype = "dashed") + 
                        geom_smooth(method = "lm", se = TRUE, colour = "#00AFBB") + 
                        labs(title = paste("Local vs E-OBS daily ",stat," ",var, sep = ""),
                            x = "Local measurements", y = "E-OBS") +
                        theme_bw()
                # Make plot 2
                p2 <- ggplot(merged_df, aes(x = obs, y = ERA5_Land)) + geom_point(alpha = 0.2, colour = "grey50") +
                        geom_abline(intercept = 0, slope = 1, color = "#d53e4f", linetype = "dashed") + 
                        geom_smooth(method = "lm", se = TRUE, colour = "#00AFBB") + 
                        labs(title = paste("Local vs ERA5-Land daily ",stat," ",var, sep = ""),
                            x = "Local measurements", y = "ERA5-Land") +
                        theme_bw()
                # Save as panel in dir  
                panel <- ggarrange(p1, p2, align = "hv", nrow = 1, ncol = 2)
                ggsave(plot = panel, filename = paste("scatter_evaluation_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 10, height = 5)  
            } # eo if loop - histo

            # 4th: time_series
            if( time_series ) {
                p <- ggplot(merged_df, aes(x = Date)) + 
                  geom_line(aes(y = obs), linewidth = 1, colour = "#00AFBB", alpha = .25) +
                  geom_line(aes(y = E_OBS), linewidth = 1, colour = "#e31a1c", alpha = .25) +
                  geom_line(aes(y = ERA5_Land), linewidth = 1, colour = "#E7B800", alpha = .25) +
                  labs(title = paste("Time series of ",stat," ",var," in the ",region," (all EPs)", sep = ""),
                      x = "Date", y = paste(stat,var,sep = " ")) +
                  theme_minimal()
                # Save
                 ggsave(plot = p, filename = paste("time_series_evaluation_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 10, height = 4)  
            } # eo if loop - time_series

          } # eo if else loop - anoms == T

          if( biases ) {

                message(paste("Plotting distribution of biases in ",stat," ",var," for the ",region," against E-OBS and ERA5-Land data", sep = ""))
                
                ## Compute biases in 'merged_df'
                merged_df$biases_EOBS <- merged_df$E_OBS - merged_df$obs
                merged_df$biases_ERA5 <- merged_df$ERA5_Land - merged_df$obs

                ## Use plots' switches again to tell the function which type of plot to save to map biases distribution
                # 1st: boxplots 
                if( boxp ) {
                    # Boxplot of stat distribution - needs melting
                    m_merged_df <- reshape2::melt(merged_df[,c("EP","Date","biases_EOBS","biases_ERA5")], id.vars = c("EP","Date"))
                    colnames(m_merged_df)[c(3)] <- c("source")
                    # Re-name before ggplotting
                    m_merged_df$source <- as.character(m_merged_df$source)
                    m_merged_df[m_merged_df$source == "biases_EOBS","source"] <- "E_OBS"
                    m_merged_df[m_merged_df$source == "biases_ERA5","source"] <- "ERA5_Land"
                    # Plot
                    plot <- ggplot(data = m_merged_df, aes(x = factor(source), y = value, fill = factor(source))) +
                        geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = .2) + 
                        scale_fill_manual(values = colors_sources) + labs(y = paste(stat,var,sep = " "), x = "") + 
                        geom_hline(yintercept = 0, linetype = "dashed") + 
                        theme_bw() + theme(legend.position = "none") +
                        ggtitle(paste("Distribution of ",stat," ",var," biases\nin the ",region, sep = ""))
                    # Save 
                    ggsave(plot = plot, filename = paste("boxplot_evaluation_biases_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 5, height = 5)  
                } # eo if loop - boxp

                # 2nd: histograms - also need melting
                if( histo ) {
                    # Boxplot of stat distribution - needs melting
                    m_merged_df <- reshape2::melt(merged_df[,c("EP","Date","biases_EOBS","biases_ERA5")], id.vars = c("EP","Date"))
                    colnames(m_merged_df)[c(3)] <- c("source")
                    # Re-name before ggplotting
                    m_merged_df$source <- as.character(m_merged_df$source)
                    m_merged_df[m_merged_df$source == "biases_EOBS","source"] <- "E_OBS"
                    m_merged_df[m_merged_df$source == "biases_ERA5","source"] <- "ERA5_Land"
                    # Make the plot
                    plot <- gghistogram(m_merged_df, x = "value", add = "mean", rug = FALSE,
                          fill = "source", palette = colors_sources,
                          ylab = "Count", xlab = paste(stat,var,sep = " "),
                          title = paste("Distribution of ",stat," ",var," biases in the ",region, sep = "")
                        ) # eo gghistogram
                    # Save 
                    ggsave(plot = plot, filename = paste("gghistogram_evaluation_biases_",stat,"_",var,"_",region,".jpg",sep = ""),
                        dpi = 300, width = 7, height = 5)  
                } # eo if loop - histo

                # 3rd: scatt 
                if( scatt ) {
                    p1 <- ggplot(merged_df, aes(x = biases_EOBS, y = biases_ERA5)) + geom_point(alpha = 0.2, colour = "grey50") +
                        geom_abline(intercept = 0, slope = 1, color = "#d53e4f", linetype = "dashed") + 
                        geom_smooth(method = "lm", se = TRUE, colour = "#00AFBB") + 
                        geom_hline(yintercept = 0, linetype = "dashed") + 
                        geom_vline(xintercept = 0, linetype = "dashed") + 
                        labs(title = paste("E-OBS vs ERA5-Land biases in\ndaily ",stat," ",var, sep = ""),
                            x = "E-OBS biases", y = "ERA5-Land biases") + theme_bw()
                    # Save
                    ggsave(plot = p1, filename = paste("scatter_evaluation_biases_",stat,"_",var,"_",region,".jpg",sep = ""),
                        dpi = 300, width = 5, height = 5)  
                } # eo if loop - histo

                # 4th: time_series
                if( time_series ) {
                    p <- ggplot(merged_df, aes(x = Date)) + 
                      geom_line(aes(y = biases_EOBS), linewidth = 1, colour = "#e31a1c", alpha = .25) +
                      geom_line(aes(y = biases_ERA5), linewidth = 1, colour = "#E7B800", alpha = .25) +
                      geom_hline(yintercept = 0, linetype = "dashed") + 
                      labs(title = paste("Time series of ",stat," ",var," biases in the ",region," (all EPs)", sep = ""),
                        x = "Date", y = paste(stat,var,sep = " ")) +
                      theme_minimal()
                    # Save
                    ggsave(plot = p, filename = paste("time_series_evaluation_biases_",stat,"_",var,"_",region,".jpg",sep = ""),
                        dpi = 300, width = 10, height = 4)  
                } # eo if loop - time_series

          } # eo if loop - biases

        } else {

          message(paste("Plotting observed daily ",stat," ",var," for the ",region," against ERA5-Land data", sep = ""))

          ## Go to ERA5-Land dir and load corresponding dataset
          setwd(paste("/home/fbenedetti/ERA5-Land-DEU-processed/daily/grassland/",region, sep = ""))

          if(var == "Ta_200") {
              era5_var <- "2m_temperature"
          } else if(var == "Ta_10") {
              era5_var <- "skin_temperature" 
          } else if(var == "precipitation") {
              era5_var <- "precipitation"
          } else if(var == "Ts_05") {
              era5_var <- "soil_temperature_level_1"
          } else if(var == "Ts_10") {
              era5_var <- "soil_temperature_level_1"
          } else if(var == "Ts_20") {
              era5_var <- "soil_temperature_level_2"
          } else if(var == "Ts_50") {
              era5_var <- "soil_temperature_level_3"
          } else if(var == "SM_10") {
              era5_var <- "volumetric_soil_water_layer_2"
          } # eo if else loop

          # Sanity check
          if( exists("era5_var") == FALSE ) {
            stop(
              paste("!!! ERROR: Could not find matching ERA5-Land ",paste(stat,var, sep = " ")," data", sep = "")
            )
          } # eo if loop - sanity check
            
          # Identify files of interest and load them
          era5_files <- dir()[grepl(era5_var,dir())] # era5_files
          era5_data <- lapply(era5_files, function(f) {
                    d <- get(load(f))
                    return(d)
              } # eo fun in lapply
          ) # eo lapply - era5_files
          era5_ddf <- dplyr::bind_rows(era5_data)
          rm(era5_data,era5_files); gc()

          era5_ddf <- era5_ddf[,c("region","system","plot","date","day","month","year",stat)]

          # For each EP: match the daily data that overlap across all 3 sources and combine in a single ddf
          # Rename columns so names match across 3 data.frames
          colnames(obs_daily_stat)[3] <- paste(var,stat, sep = "_")
          colnames(era5_ddf)[c(3,4,length(era5_ddf))] <- c("EP","Date",paste(var,stat, sep = "_"))
          names <-  c("EP","Date",paste(var,stat, sep = "_")) # vector of colnames to join the 2 ddfs by

          # To make sure Date format is homogeneous across all 3 tables
          if( class(obs_daily_stat$Date) != "Date" ) {
              obs_daily_stat$Date <- as.Date(obs_daily_stat$Date)
          } # eo if loop

          if( class(era5_ddf$Date) != "Date" ) {
              era5_ddf$Date <- as.Date(era5_ddf$Date)
          } # eo if loop

          # Merge using full_join iteratively
          merged_df <- reduce(
                list(obs_daily_stat[,names],
                    era5_ddf[,names]
                ), full_join, by = c("EP","Date")
          ) # eo reduce

          # Adjust colnames
          colnames(merged_df)[c(3,4)] <- c("obs","ERA5_Land") 
          
          # Remove NAs for evaluation
          merged_df <- na.omit(merged_df)

          ### SWITCHES:
          ### If anoms == TRUE --> make plots based on anomalies to the monthly mean
          ### If std.anoms == TRUE --> make plots based on standardized anomalies
          ### else --> make plots on normal data

          if( mon.anoms ) { 

            # Compute monthly means for each data source
            merged_df <- merged_df %>%
                mutate(month = format(Date, "%Y-%m")) %>%
                group_by(EP,month) %>%
                mutate(
                  mon_mean_obs = mean(obs, na.rm = TRUE),
                  mon_mean_ERA5 = mean(ERA5_Land, na.rm = TRUE)
                ) %>% 
                ungroup()
            # summary(merged_df)

            # Compute monthly anomalies
            merged_df <- merged_df %>%
                mutate(
                  mon_anom_obs = obs - mon_mean_obs,
                  mon_anom_ERA5 = ERA5_Land - mon_mean_ERA5
                )

            ### Save plots depending on the switches activated as function's arguments: boxp, histo, scatt, time_series
            setwd(plot.dir)
            colors_sources <- c("obs" = "#00AFBB", "ERA5_Land" = "#E7B800")

            # 1st: boxplots 
            if( boxp ) {
                # Boxplot of stat distribution - needs melting
                m_merged_df <- reshape2::melt(merged_df[,c("EP","Date","mon_anom_obs","mon_anom_ERA5")],
                                        id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Re-name source levels so they match the 'colors_sources' palette above; unique(m_merged_df$source)
                m_merged_df$source <- as.character(m_merged_df$source)
                m_merged_df[m_merged_df$source == "mon_anom_obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "mon_anom_ERA5","source"] <- "ERA5_Land"
                # Make the plot
                plot <- ggplot(data = m_merged_df, aes(x = factor(source), y = value, fill = factor(source))) +
                  geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = .2) + 
                  geom_hline(yintercept = 0, linetype = "dashed") + 
                  scale_fill_manual(values = colors_sources) + labs(y = paste(stat,var,sep = " "), x = "") + 
                  theme_bw() + theme(legend.position = "none") +
                  ggtitle(paste("Distribution of ",stat," ",var,"\nmonthly anomalies in the ",region, sep = ""))
                # Save plot in dir  
                ggsave(plot = plot, filename = paste("boxplot_evaluation_anomalies_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 5, height = 5)  
            } # eo if loop - boxp

            # 2nd: histograms - also need melting
            if( histo ) {

                m_merged_df <- reshape2::melt(merged_df[,c("EP","Date","mon_anom_obs","mon_anom_ERA5")],
                                        id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Re-name source levels so they match the 'colors_sources' palette above; unique(m_merged_df$source)
                m_merged_df$source <- as.character(m_merged_df$source)
                m_merged_df[m_merged_df$source == "mon_anom_obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "mon_anom_ERA5","source"] <- "ERA5_Land"
                # Make the plot
                plot <- gghistogram(m_merged_df, x = "value", add = "mean", rug = FALSE,
                          fill = "source", palette = colors_sources,
                          ylab = "Count", xlab = paste(stat,var,sep = " "),
                          title = paste("Distribution of ",stat," ",var," monthly anomalies in the ",region, sep = "")
                )
                # Save plot in dir  
                ggsave(plot = plot, filename = paste("gghistogram_evaluation_anomalies_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 7, height = 5)  
            } # eo if loop - histo

            # 3rd: scatt - make 2: obs vs E-OBS & obs vs ERA5-Land
            if( scatt ) {
                p1 <- ggplot(merged_df, aes(x = mon_anom_obs, y = mon_anom_ERA5)) + geom_point(alpha = 0.2, colour = "grey50") +
                        geom_abline(intercept = 0, slope = 1, color = "#d53e4f", linetype = "dashed") + 
                        geom_smooth(method = "lm", se = TRUE, colour = "#00AFBB") + 
                        geom_hline(yintercept = 0, linetype = "dashed") + 
                        geom_vline(xintercept = 0, linetype = "dashed") + 
                        labs(title = paste("Local vs ERA5-Land daily ",stat," ",var,"\nmonthly anomalies", sep = ""),
                            x = "Local measurements anomalies", y = "ERA5-Land anomalies") +
                        theme_bw()
                # Save
                ggsave(plot = p1, filename = paste("scatter_evaluation_anomalies_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 5, height = 5)  
            } # eo if loop - scatt

            # 4th: time_series
            if( time_series ) {
                p <- ggplot(merged_df, aes(x = Date)) + 
                  geom_line(aes(y = mon_anom_obs), linewidth = 1, colour = "#00AFBB", alpha = .25) +
                  geom_line(aes(y = mon_anom_ERA5), linewidth = 1, colour = "#E7B800", alpha = .25) +
                  geom_hline(yintercept = 0, linetype = "dashed") + 
                  labs(title = paste("Time series of ",stat," ",var," monthly anomalies in the ",region," (all EPs)", sep = ""),
                      x = "Date", y = paste(stat,var,sep = " ")) +
                  theme_minimal()
                # Save
                 ggsave(plot = p, filename = paste("time_series_evaluation_anomalies_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 10, height = 4)  
            } # eo if loop - time_series

          
           } else if( std.anoms ) {
          
            # Compute monthly anomalies
            merged_df <- merged_df %>%
                mutate(
                    SA_obs = (obs - mean(obs)) / sd(obs),
                    SA_ERA5_Land = (ERA5_Land - mean(obs)) / sd(obs)
                )
            
            ### Save plots depending on the switches activated as function's arguments: boxp, histo, scatt, time_series
            setwd(plot.dir)
            colors_sources <- c("obs" = "#00AFBB", "ERA5_Land" = "#E7B800")

            # 1st: boxplots 
            if( boxp ) {
 
                m_merged_df <- reshape2::melt(merged_df[,c("EP","Date","SA_obs","SA_ERA5_Land")], id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Re-name to match color palette labels
                m_merged_df$source <- as.character(m_merged_df$source)
                m_merged_df[m_merged_df$source == "SA_obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "SA_ERA5_Land","source"] <- "ERA5_Land"
                # Make the plot
                plot <- ggplot(data = m_merged_df, aes(x = factor(source), y = value, fill = factor(source))) +
                  geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = .2) + 
                  scale_fill_manual(values = colors_sources) + labs(y = paste(stat,var,sep = " "), x = "") + 
                  theme_bw() + theme(legend.position = "none") +
                  ggtitle(paste("Distribution of ",stat," ",var,"\nstandardized anomalies in the ",region, sep = ""))
                # Save plot in dir  
                ggsave(plot = plot, filename = paste("boxplot_evaluation_std.anoms_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 5, height = 5)  
            } # eo if loop - boxp

            # 2nd: histograms - also need melting
            if( histo ) {
 
                m_merged_df <- reshape2::melt(merged_df[,c("EP","Date","SA_obs","SA_ERA5_Land")], id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Re-name to match color palette labels
                m_merged_df$source <- as.character(m_merged_df$source)
                m_merged_df[m_merged_df$source == "SA_obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "SA_ERA5_Land","source"] <- "ERA5_Land"
                # Make the plot
                plot <- gghistogram(m_merged_df, x = "value", add = "mean", rug = FALSE,
                          fill = "source", palette = colors_sources,
                          ylab = "Count", xlab = paste(stat,var,sep = " "),
                          title = paste("Distribution of ",stat," ",var,"\nstandardized anomalies in the ",region, sep = "")
                )
                # Save plot in dir  
                ggsave(plot = plot, filename = paste("gghistogram_evaluation_std.anoms_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 7, height = 5)  
            } # eo if loop - histo

            # 3rd: scatt - make 2: obs vs E-OBS & obs vs ERA5-Land
            if( scatt ) {
                p <- ggplot(merged_df, aes(x = SA_obs, y = SA_ERA5_Land)) + geom_point(alpha = 0.2, colour = "grey50") +
                        geom_abline(intercept = 0, slope = 1, color = "#d53e4f", linetype = "dashed") + 
                        geom_smooth(method = "lm", se = TRUE, colour = "#00AFBB") + 
                        labs(title = paste("Local vs ERA5-Land daily ",stat," ",var,"\nbased on standardized anomalies", sep = ""),
                            x = "Local measurements", y = "ERA5-Land") +
                        theme_bw()
                # Save plot  
                ggsave(plot = p, filename = paste("scatter_evaluation_std.anoms_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 5, height = 5)  
            } # eo if loop - histo

            # 4th: time_series
            if( time_series ) {
                p <- ggplot(merged_df, aes(x = Date)) + 
                  geom_line(aes(y = SA_obs), linewidth = 1, colour = "#00AFBB", alpha = .25) +
                  geom_line(aes(y = SA_ERA5_Land), linewidth = 1, colour = "#E7B800", alpha = .25) +
                  labs(title = paste("Time series of ",stat," ",var,"\nstandardized anomalies in the ",region," (all EPs)", sep = ""),
                      x = "Date", y = paste(stat,var,sep = " ")) +
                  theme_minimal()
                # Save
                 ggsave(plot = p, filename = paste("time_series_evaluation_std.anoms_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 10, height = 4)  
            } # eo if loop - time_series
          
          } else {

            ### Save plots depending on the switches activated as function's arguments: boxp, histo, scatt, time_series
            setwd(plot.dir)
            colors_sources <- c("obs" = "#00AFBB", "ERA5_Land" = "#E7B800")

            # 1st: boxplots 
            if( boxp ) {
                # Boxplot of stat distribution - needs melting
                m_merged_df <- reshape2::melt(merged_df, id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Make the plot
                plot <- ggplot(data = m_merged_df, aes(x = factor(source), y = value, fill = factor(source))) +
                  geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = .2) + 
                  scale_fill_manual(values = colors_sources) + labs(y = paste(stat,var,sep = " "), x = "") + 
                  theme_bw() + theme(legend.position = "none") +
                  ggtitle(paste("Distribution of ",stat," ",var," in the ",region, sep = ""))
                # Save 
                ggsave(plot = plot, filename = paste("boxplot_evaluation_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 5, height = 5)  
            } # eo if loop - boxp

            # 2nd: histograms - also need melting
            if( histo ) {
                # Boxplot of stat distribution - needs melting
                m_merged_df <- reshape2::melt(merged_df, id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Make the plot
                plot <- gghistogram(m_merged_df, x = "value", add = "mean", rug = FALSE,
                          fill = "source", palette = colors_sources,
                          ylab = "Count", xlab = paste(stat,var,sep = " "),
                          title = paste("Distribution of ",stat," ",var," in the ",region, sep = "")
                )
                # Save 
                ggsave(plot = plot, filename = paste("gghistogram_evaluation_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 7, height = 5)  
            } # eo if loop - histo

            # 3rd: scatt 
            if( scatt ) {
                p1 <- ggplot(merged_df, aes(x = obs, y = ERA5_Land)) + geom_point(alpha = 0.2, colour = "grey50") +
                        geom_abline(intercept = 0, slope = 1, color = "#d53e4f", linetype = "dashed") + 
                        geom_smooth(method = "lm", se = TRUE, colour = "#00AFBB") + 
                        labs(title = paste("Local vs ERA5-Land daily ",stat," ",var, sep = ""),
                            x = "Local measurements", y = "ERA5-Land") +
                        theme_bw()
                # Save
                ggsave(plot = p1, filename = paste("scatter_evaluation_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 5, height = 5)  
            } # eo if loop - histo

            # 4th: time_series
            if( time_series ) {
                p <- ggplot(merged_df, aes(x = Date)) + 
                  geom_line(aes(y = obs), linewidth = 1, colour = "#00AFBB", alpha = .25) +
                  geom_line(aes(y = ERA5_Land), linewidth = 1, colour = "#E7B800", alpha = .25) +
                  labs(title = paste("Time series of ",stat," ",var," in the ",region," (all EPs)", sep = ""),
                      x = "Date", y = paste(stat,var,sep = " ")) +
                  theme_minimal()
                # Save
                 ggsave(plot = p, filename = paste("time_series_evaluation_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 10, height = 4)  
            } # eo if loop - time_series

            ### Plot biases distribution if "biases = TRUE"
            if( biases ) {

                message(paste("Plotting distribution of ",stat," ",var," biases for the ",region," against ERA5-Land data", sep = ""))
                # Re-use the previous plts switches to choose how you want to plot biases distributions

                setwd(plot.dir)
                colors_sources <- c("obs" = "#00AFBB", "ERA5_Land" = "#E7B800")

                ## Compute biases in 'merged_df'
                merged_df$biases <- merged_df$ERA5_Land - merged_df$obs

                ### SKIP BOXPLOTS SINCE ONLY ONE SOURCE OF BIASES: ERA5-Land

                # 1st: histograms - also need melting
                if( histo ) {
                    # Make the plot
                    plot <- gghistogram(merged_df, x = "biases", add = "mean", rug = FALSE,
                          fill = "#E7B800", ylab = "Count", xlab = paste(stat,var,sep = " "),
                          title = paste("Distribution of ",stat," ",var," biases in the ",region, sep = "")
                        ) # eo gghistogram
                    # Save 
                    ggsave(plot = plot, filename = paste("gghistogram_evaluation_biases_",stat,"_",var,"_",region,".jpg",sep = ""),
                        dpi = 300, width = 7, height = 5)  
                } # eo if loop - histo

                ### SKIP SCATTER PLOTS TOO

                # 2nd: time_series
                if( time_series ) {
                    p <- ggplot(merged_df, aes(x = Date)) + 
                      geom_line(aes(y = biases), linewidth = 1, colour = "#E7B800", alpha = .25) +
                      geom_hline(yintercept = 0, linetype = "dashed") + 
                      labs(title = paste("Time series of ",stat," ",var," biases in the ",region," (all EPs)", sep = ""),
                        x = "Date", y = paste(stat,var,sep = " ")) +
                      theme_minimal()
                    # Save
                    ggsave(plot = p, filename = paste("time_series_evaluation_biases_",stat,"_",var,"_",region,".jpg",sep = ""),
                        dpi = 300, width = 10, height = 4)  
                } # eo if loop - time_series

            } # eo if loop - biases

          } # eo if else loop - anoms

        } # eo if else loop

} # eo plotting FUN - plot_daily_stat_comparison

### --> Script#6.3.3 to compute biases (models - obs) like in the plotting FUN above and save as ddf
### --> Script#6.3.4 to run these two functions and analyze outputs

### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------