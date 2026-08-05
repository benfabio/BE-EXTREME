### ------------------------------------------------------------------------------------------------------------

### 08/04/25 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)
###           with input from Adrian Huerta (Klimatologie group, GIUB, Uni Bern)

### R script to perform quantile mapping (QM) with the 'qmap' R package to correct the time series
### of daily climate data (i.e., daily minimum temperature) from E-OBS (air temperature, precipitation) and 
### ERA5-Land (all other variables).
### Previous scripts (Script#6.3) showed that E-OBS matched the Exploratories local climate data better
### than ERA5-Land and that systemic temporal biases exist (biases tend to be higher in summer).
### Furthermore, this bias varies between seasons and decreases through time --> need to perform QM within each EP!  

### However, various QM could be implemented in diverse ways and it is unclear which strategy is the best: 
### - local, based on the full time series (TS) of each EP
### - local & monthly, decompose full TS into 12 monthly TS and apply 12 QM models per EP, 
###   and then reconstruct full TS
### - local with moving window: for each day of the year (1-365), take the 15 days before 
###   and after and run a QM model on that TS, then reconstruct full TS

### R script to write a master FUN that will:
### - Load both the local daily data and the reconstructed data (E-OBS or ERA5-Land)
### - Join the two datasets by overlapping 
### - Perform QM (fitQmapQUANT()) according to the 3 strategies mentioned above
###   (i.e., add a switch argument in the master FUN to choose QM option: 'global', 'monthly', 'mw')
### - Return fully corrected TS and save it in appropriate directory

### --> Script#6.4.3 to evaluate corrected TS and decide which strategy is the best

### NOTE: Depending on the results, consider performing QM on the ANOMALIES to the monthly mean?

### Last update: 21/01/26 (Testing quantile_forest_mapper())

### ------------------------------------------------------------------------------------------------------------

# Libraries 
library("dplyr")
library("tidyr")
library("data.table")
library("purrr")
library("tibble")
library("reshape2")
require("scales")
library("lubridate")
require("RColorBrewer")
library("viridis")
library("ggpubr")
library("parallel")
library("qmap") # to perform quantile mapping

### ------------------------------------------------------------------------------------------------------------

### NOTE
## With 'qmap', you can do both 'fitQmapQUANT' (basic) or 'fitQmapRQUANT'
## fitQmapRQUANT() adds regularization/smoothing to the quantile mapping function.
## It internally uses monotonic smoothing splines instead of linear interpolation between quantiles.
## It avoids abrupt changes and smooths the quantile–quantile relationship.
## Overall, it's better at handling extreme values (tails) and smooths mapping function and 
## avoids artificial discontinuities.

### Write master FUN to return the corrected E-OBS and ERA5-Land TS. What argument should the FUN comprise?
### - variable (use labels from local obs: Ta_200, Ta_10, precipitation, Ts_10 etc.)
### - stat (min/max/mean)
### - region (SCH/HND/SWA) - but not EP, use mclapply() to run QM on all EP of the region
### - method ('global','monthly','mw' or 'anoms')

### Master FUN - evaluate_daily_stat
# To test quantile_mapper() while you're writing it: 
# var <- "precipitation"
# stat <- "total"
# region <- "SWA"
# method <- "monthly"
# w_size = 15


quantile_mapper <- function(var, stat, region, method, w_size) {

        #' This function takes four arguments and returns a formatted data.frame:
        #' @param var the climate variable to process (character) - one of the following: 
        #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20", "Ts_50", "precipitation" or "SM_10"
        #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character): 'mean', 'max' or 'min'
        #' @param region Region name (character): 'SCH' or 'HND' or 'SWA'
        #' @param method Which quantile mapping approach to run (character): 'global' or 'monthly' or 'mw' or 'anoms'
        #' @param w_size The size (in days) of the moving window (integer)

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

        # Read the file containing your observed daily stats
        obs_daily_stat <- read.csv(files, h = T, sep = ",", dec = ".")

        # Sanity check
        if( exists("obs_daily_stat") == FALSE ) {
            stop(
              paste("!!! ERROR: Could not load observed ",paste(var,stat, sep = "_")," file for the ",region, sep = "")
            )
        } # eo if loop - sanity check

        # For Ta-200 and precipitation --> E-OBS data
        if( var %in% c("Ta_200","precipitation") ) { 
            
            ## Go to E-OBS dir and load their daily data too
            message(paste("\nCorrecting E-OBS with observed daily ",stat," ",var," for the ",region,
                    " based on the ",method," quantile mapping strategy\n", sep = ""))
            
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
            
            # Re-name 
            colnames(sub_EOBS_daily_stat) <- c("EP","Date","value")

            # For each EP: match the daily data that overlap across all 3 sources and combine in a single ddf
            # Rename columns so names match across 3 data.frames
            colnames(obs_daily_stat)[3] <- paste(var,stat, sep = "_")
            colnames(sub_EOBS_daily_stat)[3] <- paste(var,stat, sep = "_")
            names <-  c("EP","Date",paste(var,stat, sep = "_")) # vector of colnames to join the 3 ddf by

            # To make sure Date format is homogeneous across all 3 tables
            if( class(obs_daily_stat$Date) != "Date" ) {
                obs_daily_stat$Date <- as.Date(obs_daily_stat$Date)
            } # eo if loop

            if( class(sub_EOBS_daily_stat$Date) != "Date" ) {
                sub_EOBS_daily_stat$Date <- as.Date(sub_EOBS_daily_stat$Date)
            } # eo if loop

            # Merge using full_join iteratively
            merged_df <- reduce(
                            list(obs_daily_stat[,names],
                                 sub_EOBS_daily_stat[,names]
                            ), full_join, by = c("EP","Date")
            ) # eo reduce

            # Adjust colnames
            colnames(merged_df)[c(3:4)] <- c("obs","E_OBS")
            
            # Remove NAs for evaluation
            merged_df <- na.omit(merged_df)
            # Base data.frame with EP/Date/obs & model data

            # Vector of EP IDs - use to subset 'merged_df'
            plots <- unique(merged_df$EP) # plots

            ### Start if else loop based on QM method: 
            ## if method == "global" --> Perform basic QM for each EP separately based on overlapping data (May 2009-2024)
            ## if method == "monthly" --> Split and join the data into 12 monthly TS, perform QM on those
            ## if method == "mw" --> Using mowing window to smooth the threshold effects you may encounter
            ## if method == "anoms" --> Correcting the anomalies to the mean instead of the raw values like in the other 3 methods

            if( method == "global" ) {

                pred_QM_all_plots <- mclapply(plots, function(p) {
                        
                        # p <- plots[27]
                        message(paste("Performing quantile mapping for ",p," based on the ",method," approach", sep = ""))
                        sub_merged_df <- merged_df[merged_df$EP == p,]
                        
                        # Subset 'sub_EOBS_daily_stat' to have the full TS to project the QM on below
                        sub_EOBS <- sub_EOBS_daily_stat[sub_EOBS_daily_stat$EP == p,]
                        
                        # Remove 2008 because unreliable year for observations
                        # "2009-05-01" instead of "2009-01-01" because of tempertaure data in the SWA
                        sub_merged_df2 <- sub_merged_df[which(sub_merged_df$Date >= "2009-05-01"),]
                        
                        ### In the case of precipitation, only train the QM model by excluding the days 
                        ### that are observed and modelled as dry (precip. = 0).
                        ### Then, correct the model data that are also non dry
                        if( var == "precipitation" ) {

                            # Fit and apply quantile mapping by excluding the days that are dry in both model and obs (wet.day = 0)
                            # Setting wet.days to 0 disables the wet day correction process.
                            # This means the function will not try to match the fraction of wet days between the modeled and observed data.
                            # Instead, any values below the specified threshold (which can be a numeric value or derived from the observed data)
                            # will be treated as zero.
                            # The result is that the modeled data is not adjusted to match the observed data's wet day frequency.
                            # This can be useful in situations where the wet day correction is not appropriate or desired.
                            # For example, it might be used to analyze data where the frequency of wet days is not of primary concern. 
                            qm_model <- fitQmapRQUANT(obs = sub_merged_df2$obs, mod = sub_merged_df2$E_OBS, qstep = 0.01, wet.day = 0)
                            
                            # NOTE: Here, I am using the mod and obs data that overlap in time of course (2009-2023 period).
                            # However, the daily model data I want to predict the model on (with doQmapRQUANT()) span a much longer period 
                            # (1950-2023).
                            # Yet, we do not want to do the doQmapRQUANT() on those model data correctly equal to 0 in the
                            # 2009-2023 period. Since the obs data do not go back before 2008, there is no way I can tell if the dry days
                            # (precipitation = 0) before 2009 are correct or not. We need to make sure to use doQmapRQUANT() on
                            # all mod data < 2009 plus only those days > 2009 that are wrongly modelled as dry days. 

                            # Indices for dates before and after 2009
                            idx_EOBS_pre2009 <- which(sub_EOBS[,"Date"] < as.Date("2009-05-01"))

                            # Get the indices of the pre 2009 ('idx_EOBS_pre2009') + the days 
                            # of the post-2009 model data that were positively modelled as 0 (dry days)
                            # To do so, identify the dates where obs and E-OBS are both equal to 0 and remove them with '!'.
                            # Then, use those to get the IDs from 'idx_EOBS_post2009'
                            dates2keep <- sub_merged_df2[which(!(sub_merged_df2$obs == 0 & sub_merged_df2$E_OBS == 0)),"Date"]
                            
                            # Subset sub_EOBS to retain those days where 0 where NOT correctly modelled (= dates2keep) by their indices
                            sub_EOBS_dates2keep <- which(sub_EOBS$Date %in% dates2keep)

                            # Combine all pre-2009 + the selected post-2009 days to be corrected (idx_EOBS_pre2009 + sub_EOBS_dates2keep)
                            correct_indices <- c(idx_EOBS_pre2009, sub_EOBS_dates2keep)
                            
                            # Predict full TS based on 'qm_model' object using 'correct_indices'
                            pred <- doQmapRQUANT(x = sub_EOBS[correct_indices,3], fobj = qm_model, type = "linear")
                            sub_EOBS2 <- sub_EOBS # retain uncorrected E-OBS data
                            sub_EOBS2[correct_indices,3] <- pred # add corrected data at the right indices

                            # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                            ddf <- data.frame(
                                    EP = p,
                                    Date = sub_EOBS[,"Date"],
                                    raw_value = sub_EOBS[,3],
                                    corrected_value = sub_EOBS2[,3]
                            ) # eo ddf
                            ## Check corrected data against initial obs data quickly
                            # summary(sub_merged_df2$obs)
                            # summary(ddf[ddf$Date >= as.Date("2009-01-01"),"corrected_value"]) # EXTREMELY good agreement

                        } else {

                            # Fit and apply quantile mapping
                            qm_model <- fitQmapRQUANT(
                                    obs = sub_merged_df2$obs,
                                    mod = sub_merged_df2$E_OBS,
                                    qstep = 0.01,
                                    wet.day = FALSE
                            ) # eo fitQmapRQUANT
                        
                            # Predict full TS based on 'qm_model' object
                            pred <- doQmapRQUANT(
                                    x = sub_EOBS[,3],
                                    fobj = qm_model,
                                    type = "linear"
                            ) # eo doQmapRQUANT
                        
                            # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                            ddf <- data.frame(
                                    EP = p,
                                    Date = sub_EOBS[,"Date"],
                                    raw_value = sub_EOBS[,3],
                                    corrected_value = pred
                            ) # eo ddf

                        } # eo if else loop - var == precipitation
                        
                        # Clean & return ddf
                        rm(pred,qm_model,sub_merged_df2,sub_EOBS,sub_merged_df)
                        return(ddf)

                    }, mc.cores = 25

                ) # eo mclapply - plots

                # Rbind and save in proper dir ('/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/global')
                table <- dplyr::bind_rows(pred_QM_all_plots)
                # dim(table) ; summary(table) ; table[134050:136050,]
                rm(pred_QM_all_plots); gc()

                # Add stat, vars, region before saving 
                table <- add_column(table, region = region, .after = "EP")
                table <- add_column(table, var = var, .after = "region")
                table <- add_column(table, stat = stat, .after = "var")

                # Save
                setwd("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/global") # dir()
                save(x = table, file = paste("table_pred_QM_",method,"_",stat,"_",var,"_",region,".Rdata", sep = ""))
                # Remove some stuff
                rm(table); gc()

            } else if ( method == "monthly" ) {

                ### 21/05: FIX ISSUE WITH PRECIPITATION DATA HERE

                pred_QM_all_plots <- mclapply(plots, function(p) {
                        
                        # p <- plots[25]; p
                        message(paste("Performing quantile mapping for ",p," based on the ",method," approach", sep = ""))
                        sub_merged_df <- merged_df[merged_df$EP == p,]
                        # dim(sub_merged_df)
                        
                        # Subset 'sub_EOBS_daily_stat' to have the full TS to priject the QM on below
                        sub_EOBS <- sub_EOBS_daily_stat[sub_EOBS_daily_stat$EP == p,]
                        
                        # Split this one into 12 monthly series to use as objects for QM prediction
                        monthly_EOBS_for_pred <- split(sub_EOBS, format(sub_EOBS$Date, "%m")) 
                        # monthly_EOBS_for_pred[[4]][1:30,] # very good
                        
                        # Remove 2008 because unreliable year for observations
                        sub_merged_df2 <- sub_merged_df[which(sub_merged_df$Date >= "2009-05-01"),]

                        # Split the data series frame by month
                        monthly_series <- split(sub_merged_df2, format(sub_merged_df2$Date, "%m"))
                        # str(monthly_series); monthly_series[[10]][1:50,] # Worked! 
                        
                        # For each element m of 'monthly_series' (= month), perform QM
                        # m <- 3
                        res_monthly_preds <- lapply(c(1:12), function(m) {
       
                                if( var == "precipitation" ) {

                                    # Fit and apply quantile mapping
                                    qm_model <- fitQmapRQUANT(
                                                obs = monthly_series[[m]][,3],
                                                mod = monthly_series[[m]][,4],
                                                qstep = 0.01,
                                                wet.day = 0
                                    ) # eo fitQmapRQUANT
                            
                                    # Apply only on wet days in the model: get their index
                                    # idx_EOBS_pre2009 <- which(sub_EOBS[,"Date"] < as.Date("2009-01-01"))
                            
                                    # Get the indices of the pre 2009 ('idx_EOBS_pre2009') + the days 
                                    # of the post-2009 model data that were positively modelled as 0 (dry days)
                                    # To do so, identify the dates where obs and E-OBS are both equal to 0 and remove them with '!'.
                                    # Then, use those to get the IDs from 'idx_EOBS_post2009'
                                    # dates2keep <- monthly_series[[m]][which(!(monthly_series[[m]]$obs == 0 & monthly_series[[m]]$E_OBS == 0)),"Date"]
                            
                                    # Use those to identify the IDs 'sub_EOBS'
                                    # sub_EOBS_dates2keep <- which(sub_EOBS$Date %in% dates2keep)
                                    # duplicated_rows <- duplicated(sub_EOBS_dates2keep) ; unique(duplicated_rows)

                                    # Combine: all pre-2009 + selective post-2009
                                    # correct_indices <- c(idx_EOBS_pre2009, sub_EOBS_dates2keep)

                                    ### 21/05/25: Correcting the following issue: 
                                    ### I applied each month-specific quantile model to ALL relevant dates (not filtered by month)
                                    ### -> Ended up with 12 monthly corrections but for each single 'Date'
                                    ### -> Ensure each month’s model is only applied to dates in that month (format(Date,"%m") == m)

                                    # Limit everything to month `m`
                                    # Format the current month index (m) as a two-digit string (e.g., "01", "02", ..., "12")
                                    month_str <- sprintf("%02d", m)
                                    
                                    # Identify the row indices in sub_EOBS corresponding to the current month
                                    sub_EOBS_month_idx <- which(format(sub_EOBS$Date, "%m") == month_str)
                                    
                                    # Identify indices in sub_EOBS that are:
                                    #   (1) before 2009-01-01 (pre-2009), &
                                    #   (2) within the current month
                                    idx_EOBS_pre2009_m <- intersect(which(sub_EOBS$Date < as.Date("2009-01-01")), sub_EOBS_month_idx)
                                    
                                    # For post-2009 data, identify the dates where neither obs nor E-OBS is zero
                                    # These are the "wet days" that should be corrected
                                    dates2keep <- monthly_series[[m]][which(!(monthly_series[[m]]$obs == 0 & monthly_series[[m]]$E_OBS == 0)),"Date"]
                                    
                                    # From sub_EOBS, identify indices of the dates that: are in dates2keep & belong to the current month 
                                    sub_EOBS_dates2keep_m <- intersect(which(sub_EOBS$Date %in% dates2keep), sub_EOBS_month_idx)

                                    # Combine the pre-2009 dates (from the same month) and the selected post-2009 wet days
                                    # These are the rows on which quantile mapping will be applied for month `m`
                                    correct_indices <- c(idx_EOBS_pre2009_m, sub_EOBS_dates2keep_m)
                            
                                    # Predict full TS based on 'qm_model' object using 'correct_indices'
                                    pred <- doQmapRQUANT(
                                                x = sub_EOBS[correct_indices,3],
                                                fobj = qm_model,
                                                type = "linear"
                                    ) # eo - doQmapRQUANT
                                    # sub_EOBS2 <- sub_EOBS
                                    # sub_EOBS2[correct_indices,3] <- pred

                                    # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                                    # ddf <- data.frame(
                                    #            EP = p,
                                    #            Date = sub_EOBS[,"Date"],
                                    #            raw_value = sub_EOBS[,3],
                                    #            corrected_value = sub_EOBS2[,3]
                                    # ) # eo ddf

                                    # data.frame above replaced by this:
                                    ddf <- data.frame(
                                                EP = p,
                                                Date = sub_EOBS[correct_indices, "Date"],
                                                raw_value = sub_EOBS[correct_indices,3],
                                                corrected_value = pred
                                    ) # eo ddf
                                    
                                    ## Check corrected data against initial obs data quickly
                                    # summary( monthly_series[[m]][,3] )
                                    # summary(ddf[ddf$Date >= as.Date("2009-01-01"),"corrected_value"])
                                    # sum(duplicated(ddf[,c("EP","Date")]))

                                } else {

                                    # Fit QM model on TS of corresponding month m
                                    qm_model <- fitQmapRQUANT(
                                                obs = monthly_series[[m]][,3],
                                                mod = monthly_series[[m]][,4],
                                                qstep = 0.01,
                                                wet.day = FALSE
                                    ) # eo fitQmapRQUANT

                                    # Predict on monthly_EOBS_for_pred[[m]]
                                    pred <- doQmapRQUANT(
                                                x = monthly_EOBS_for_pred[[m]][,3],
                                                fobj = qm_model,
                                                type = "linear"
                                    ) # eo doQmapRQUANT

                                    # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                                    ddf <- data.frame(
                                                EP = p,
                                                Date = monthly_EOBS_for_pred[[m]][,"Date"],
                                                raw_value = monthly_EOBS_for_pred[[m]][,3],
                                                corrected_value = pred
                                    ) # eo ddf

                                } # eo if else loop - var = precipitation 
                        
                                # Clean & return ddf
                                rm(pred,qm_model)
                                return(ddf)

                            } # eo FUN 

                        ) # eo lapply - m

                        # Combine 'res_monthly_preds' back into one data frame & restore chronological order
                        tible <- dplyr::bind_rows(res_monthly_preds)
                        tible <- tible[order(tible$Date),]

                        # Check for duplicates?
                        # duplicated_rows <- duplicated(tible) ; unique(duplicated_rows)
                        # sum(duplicated(tible[,c("EP","Date")])) # should return 0

                        rm(monthly_series,sub_merged_df2,monthly_EOBS_for_pred,sub_EOBS,sub_merged_df)
                        gc()

                        # Return
                        return(tible)

                    }, mc.cores = 25

                ) # eo mclapply - plots

                # Rbind and save in proper dir
                table <- dplyr::bind_rows(pred_QM_all_plots)
                # dim(table) ; summary(table) ; table[1:300,]
                # unique(table$EP) ; unique(table$Date)
                rm(pred_QM_all_plots); gc()

                # Add stat, vars, region before saving 
                table <- add_column(table, region = region, .after = "EP")
                table <- add_column(table, var = var, .after = "region")
                table <- add_column(table, stat = stat, .after = "var")

                # duplicated_rows <- duplicated(table) ; unique(duplicated_rows)
                # summary(factor(duplicated_rows)) # Still some duplicates somehow
                # sum(duplicated(table[,c("EP","Date")])) # should return 0
                # table[1220322:1231322,]

                ### 22/05/25: Check for rows with only NAs
                # only_na_rows <- apply(table, 1, function(x) all(is.na(x)))
                # any(only_na_rows) # OK, should be good
                
                # Save
                setwd("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/monthly") # dir()
                save(x = table, file = paste("table_pred_QM_",method,"_",stat,"_",var,"_",region,".Rdata", sep = "") )
                # Remove some stuff
                rm(table); gc()

            } else if ( method == "mw" ) {
                
                # Define window size for the mowing window
                window_size <- w_size

                pred_QM_all_plots <- mclapply(plots, function(p) {

                        # p <- plots[13]
                        message(paste("Performing quantile mapping for ",p," based on the ",method," approach", sep = ""))
                        sub_merged_df <- merged_df[merged_df$EP == p,]
                        
                        # Add vector of unique days of year (1–366) to subset
                        sub_merged_df$doy <- as.integer(format(sub_merged_df$Date, "%j"))
                        doys <- sort(unique(sub_merged_df$doy))
                        
                        # Subset 'sub_EOBS_daily_stat' to have the full TS to priject the QM on below
                        sub_EOBS <- sub_EOBS_daily_stat[sub_EOBS_daily_stat$EP == p,]

                        # Add vector of unique days of year (1–366) to EOBS data (full TS needed for prediction)
                        sub_EOBS$doy <- as.integer(format(sub_EOBS$Date, "%j"))

                        # Use lapply to process each day-of-year (doy)
                        ### WARNING: moving window option takes a lot longer !
                        res_mw_preds <- lapply(doys, function(d) {
                                
                                    # Moving window: handle wrap-around at start/end of year
                                    window_days <- ((d - window_size):(d + window_size)) %% 366
                                    
                                    # if there is a 0 in the moving window, replace it 0 by 366 (loop around the year)
                                    if( 0 %in% window_days ) {
                                        window_days[window_days == 0] <- 366
                                    } # eo if loop - 0
  
                                    # Select those rows of 'sub_merged_df' that fall within the moving window
                                    in_window <- sub_merged_df$doy %in% window_days
  
                                    # Fit quantile mapping model
                                    # if loop: if var == "precipitation" use wet.day = 0, otherwise wet.day = FALSE
                                    if( var == "precipitation" ) {
                                        
                                        qm_model <- fitQmapRQUANT(
                                            obs = sub_merged_df$obs[in_window],
                                            mod = sub_merged_df$E_OBS[in_window],
                                            qstep = 0.01,
                                            wet.day = 0
                                        ) # eo - fitQmapRQUANT

                                    } else {

                                        qm_model <- fitQmapRQUANT(
                                            obs = sub_merged_df$obs[in_window],
                                            mod = sub_merged_df$E_OBS[in_window],
                                            qstep = 0.01,
                                            wet.day = FALSE
                                        ) # eo - fitQmapRQUANT

                                    } # eo if else loop 
                                    
                                    # Apply model to the full DOY values from 'sub_EOBS'
                                    # First, identify rows of full EP's TS that correspond to DOY
                                    apply_idx <- which(sub_EOBS$doy == d)
                                    
                                    # Second, predict based on the QM model for these
                                    pred <- doQmapRQUANT(
                                            x = sub_EOBS[apply_idx,3],
                                            fobj = qm_model,
                                            type = "linear"
                                    ) # eo - doQmapRQUANT

                                    # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                                    ddf <- data.frame(
                                            EP = p,
                                            Date = sub_EOBS[apply_idx,2],
                                            DOY = d,
                                            raw_value = sub_EOBS[apply_idx,3],
                                            corrected_value = pred
                                    ) # eo ddf
                        
                                    # Clean & return ddf
                                    rm(pred,qm_model,in_window,window_days,apply_idx)
                                    return(ddf)
    
                                } # eo FUN

                        ) # eo lapply - day of the year (doy)

                        # Combine 'res_monthly_preds' back into one data frame & restore chronological order
                        table_plot <- dplyr::bind_rows(res_mw_preds)
                        table_plot <- table_plot[order(table_plot$Date),]
                        # table_plot[2721:2921,]
                        
                        # Clean and return
                        rm(res_mw_preds,sub_merged_df)
                        gc()

                        # Return
                        return(table_plot)

                    }, mc.cores = 25
                
                ) # eo mclapply - plots
                
                # Rbind and save in proper dir
                table <- dplyr::bind_rows(pred_QM_all_plots)
                rm(pred_QM_all_plots); gc()
                # table[13205:13405,]

                # Add stat, vars, region before saving 
                table <- add_column(table, region = region, .after = "EP")
                table <- add_column(table, var = var, .after = "region")
                table <- add_column(table, stat = stat, .after = "var")

                # Save
                setwd("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/mw") # dir()
                save(x = table, file = paste("table_pred_QM_",method,"_",stat,"_",var,"_",region,".Rdata", sep = "") )
                # Remove some stuff
                rm(table); gc()

            } else if ( method == "anoms" ) {
                
                # Set training period
                train_start <- as.Date("2009-01-01")
                train_end <- max(merged_df$Date)

                # Run anomaly-based QM in parallel with mclapply()
                pred_QM_all_plots <- mclapply(plots, function(p) {
  
                        # p <- plots[sample(c(1:50),1,1)] # sample random EP for testing script below
                        message(paste("Performing quantile mapping for ",p," based on the ",method," approach", sep = ""))
                        # EP-level ddf with Date/obs & model data 
                        sub_merged_df <- merged_df[merged_df$EP == p,]
                        # summary(sub_merged_df)

                        # Subset 'sub_EOBS_daily_stat' to have the full TS to priject the QM on below
                        sub_EOBS <- sub_EOBS_daily_stat[sub_EOBS_daily_stat$EP == p,]
                        
                        # Compute monthly means from 'sub_EOBS' to get the full TS of daily anomalies to the monthly mean
                        # Add month column for grouping
                        sub_EOBS$Month <- format(sub_EOBS$Date, "%m")
                        monthly_means_mod <- tapply(sub_EOBS[,3], sub_EOBS$Month, mean, na.rm = TRUE)
                        # monthly_means_mod

                        # Compute anomalies to the month for the FULL model TS (for QM prediction later)
                        sub_EOBS$anom <- sub_EOBS[,3] - monthly_means_mod[sub_EOBS$Month]
                        # summary(sub_EOBS)
                        ### NOTE: sub_EOBS is not to be used to train the QM model, just for the prediction

                        # Remove 2008 because unreliable year for observations - not to be used to train QM model
                        sub_merged_df2 <- sub_merged_df[which(sub_merged_df$Date >= "2009-05-01"),]
  
                        # Add month column for grouping
                        sub_merged_df2$Month <- format(sub_merged_df2$Date, "%m")
  
                        # Get training data (overlap)
                        sub_merged_df_train <- sub_merged_df2[sub_merged_df2$Date >= train_start & sub_merged_df2$Date <= train_end,]
  
                        # Compute monthly means from OBS in training period
                        monthly_means_obs <- tapply(sub_merged_df_train$obs, sub_merged_df_train$Month, mean, na.rm = TRUE)
                        # monthly_means_obs ; monthly_means_mod
  
                        # Compute anomalies for the post 2009 period
                        sub_merged_df_train$mod_anom <- sub_merged_df_train$E_OBS - monthly_means_mod[sub_merged_df_train$Month]
                        sub_merged_df_train$obs_anom <- sub_merged_df_train$obs - monthly_means_obs[sub_merged_df_train$Month]
  
                        # Loop over months (can later parallelize if needed)
                        # m <- "10"
                        pred_QM_mon_anoms <- lapply(sprintf("%02d",1:12), function(m) {

                                # Subset obs anomalies for month 'm'
                                obs_anom_m <- sub_merged_df_train$obs_anom[sub_merged_df_train$Month == m]
                                
                                # Subset obs anomalies for month 'm'
                                mod_anom_m <- sub_merged_df_train$mod_anom[sub_merged_df_train$Month == m]
                                # summary(obs_anom_m) ; summary(mod_anom_m)
                                
                                # Dates to be corrected
                                dates <- sub_merged_df_train$Date[sub_merged_df_train$Month == m]
                                
                                # IDs of the point corresponding to month 'm'
                                idx_month <- which(sub_merged_df2$Month == m)
                                
                                # Fit QM model        
                                qmodel <- fitQmapRQUANT(
                                                obs = obs_anom_m,
                                                mod = mod_anom_m,
                                                qstep = 0.01,
                                                wet.day = FALSE
                                ) # eo fitQmapRQUANT
                                
                                # Correct model-based anomalies with QM model
                                pred <- doQmapRQUANT(
                                                x = sub_EOBS[sub_EOBS$Month == m,"anom"],
                                                fobj = qmodel,
                                                type = "linear"
                                ) # eo doQmapRQUANT
                                # summary(mod_anom_m); summary(pred)
                                
                                # Add anomaly back to monthly mean from obs or mod
                                # Here, we want to match observations as closely as possible so we add 'monthly_means_obs[m]'
                                corrected_value_means_obs <- pred + monthly_means_obs[m]
                                # Or monthly_means_mod[m] if wan to preserve model's climatology while correcting variability
                                corrected_value_means_mod <- pred + monthly_means_mod[m]	

                                # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                                ddf <- data.frame(
                                                EP = p,
                                                Date = sub_EOBS[sub_EOBS$Month == m,"Date"],
                                                raw_value = sub_EOBS[sub_EOBS$Month == m,3],
                                                corrected_value_means_obs = corrected_value_means_obs,
                                                corrected_value_means_mod = corrected_value_means_mod
                                ) # eo ddf
                                # summary(ddf)

                                ### 25/04/25: In the case of total precipitation with the 'anoms' approach, 
                                ### you may ignore the 'wet.day' argument (leave it to FALSE) since you are 
                                ### working with anomalies.
                                ### HOWEVER, after adding the 'monthly_means_obs', you may never find values == 0 
                                ### (i.e., dry days) or even negative precip. values! 
                                ### -> Need to correct that here by concerting the smallest values > 0 (i.e., 0.02766) to 0
                                
                                ## Find said min value & replace by 0
                                if( var == "precipitation" ) {
                                    
                                    min_val <- min(ddf$corrected_value_means_obs[ddf$corrected_value_means_obs > 0]) # min_val
                                    # Convert them to 0 in 'ddf'
                                    ddf[ddf$corrected_value_means_obs == min_val,"corrected_value_means_obs"] <- 0

                                    ### 29/04/25: Same for v2 (model-based monthly means used to get final corrected values)
                                    min_val <- min(ddf$corrected_value_means_mod[ddf$corrected_value_means_mod > 0]) # min_val
                                    ddf[ddf$corrected_value_means_mod == min_val,"corrected_value_means_mod"] <- 0
                                
                                    ### If there are any negative values...convert to 0 too
                                    ddf[ddf$corrected_value_means_obs < 0,"corrected_value_means_obs"] <- 0
                                    ddf[ddf$corrected_value_means_mod < 0,"corrected_value_means_mod"] <- 0

                                } # eo var == "precipitation"
                                
                                # Return & clean
                                return(ddf)
                                rm(min_val,pred,qmodel,idx_month,obs_anom_m,mod_anom_m,
                                    corrected_value_means_obs,corrected_value_means_mod)

                            } # eo FUN
                            
                        ) # eo lapply - m

                        # Combine 'res_monthly_preds' back into one data frame & restore chronological order
                        table <- dplyr::bind_rows(pred_QM_mon_anoms)
                        table <- table[order(table$Date),]
                        # dim(table); str(table); summary(table)

                        rm(sub_merged_df2,sub_EOBS,sub_merged_df)
                        gc()

                        # Return
                        return(table)
  
                    }, mc.cores = 25
                    
                ) # eo mclapply - p in plots 

                # Rbind and save in proper dir
                table <- dplyr::bind_rows(pred_QM_all_plots)
                # dim(table) ; summary(table) ; table[136050:136150,]
                rm(pred_QM_all_plots); gc()

                # Add stat, vars, region before saving 
                table <- add_column(table, region = region, .after = "EP")
                table <- add_column(table, var = var, .after = "region")
                table <- add_column(table, stat = stat, .after = "var")

                # Save
                setwd("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/anomalies") # dir()
                save(x = table, file = paste("table_pred_QM_",method,"_",stat,"_",var,"_",region,".Rdata", sep = "") )
                # Remove some stuff
                rm(table); gc()

            } # eo if else loop


        } else {


            ### ERA5-compatible variables only from here (Ta_10, Ts & SM)
             message(paste("Correcting ERA5-Land data with observed daily ",stat," ",var," for the ",region,
                    " based on the ",method," quantile mapping strategy\n", sep = ""))
            
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
              era5_var <- "volumetric_soil_water_layer_2" # volumetric_soil_water_layer_2 was slightly better than volumetric_soil_water_layer_1
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
                ),  full_join, by = c("EP","Date")
            ) # eo reduce

            # Adjust colnames
            colnames(merged_df)[c(3,4)] <- c("obs","ERA5_Land")
          
            # Remove NAs for evaluation
            merged_df <- na.omit(merged_df)

            # Vector of EP IDs - use to subset 'merged_df'
            plots <- unique(merged_df$EP)

            ### Start if else loop based on QM method: 
            ## if method == "global" --> Perform basic QM for each EP based on overlapping data
            ## if method == "monthly" --> Split and join the data into 12 monthly TS, perform QM on those
            ## if method == "mw" --> Using mowing window to smooth the threshold effects you may encounter 
            ## with the other 2 methods

            if( method == "global" ) {

                pred_QM_all_plots <- mclapply(plots, function(p) {
                        
                        # p <- plots[33]
                        message(paste("Performing quantile mapping for ",p," based on the ",method," approach", sep = ""))
                        sub_merged_df <- merged_df[merged_df$EP == p,]
                        
                        # Subset 'era5_ddf' to have the full TS to project the QM on below
                        sub_ERA5 <- era5_ddf[era5_ddf$EP == p,]
                        
                        # Remove 2008 because unreliable year for observations
                        sub_merged_df2 <- sub_merged_df[which(sub_merged_df$Date >= "2009-05-01"),]
                        
                        # Fit and apply quantile mapping
                        qm_model <- fitQmapRQUANT(
                                    obs = sub_merged_df2$obs,
                                    mod = sub_merged_df2$ERA5_Land,
                                    qstep = 0.01,
                                    wet.day = FALSE
                        ) # eo qm_model
                        
                        # Predict full TS based on 'qm_model' object
                        pred <- doQmapRQUANT(x = sub_ERA5[,8], fobj = qm_model, type = "linear")
                        
                        # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                        ddf <- data.frame(
                                    EP = p,
                                    Date = sub_ERA5[,4],
                                    raw_value = sub_ERA5[,8],
                                    corrected_value = pred
                        ) # eo ddf 
                        colnames(ddf) <- c("EP","Date","raw_value","corrected_value")
                        
                        # Clean & return ddf
                        rm(pred,qm_model,sub_merged_df2,sub_ERA5,sub_merged_df)
                        return(ddf)

                    }, mc.cores = 25

                ) # eo mclapply - plots

                # Rbind and save in proper dir ('/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/global')
                table <- dplyr::bind_rows(pred_QM_all_plots)
                # summary(table)
                rm(pred_QM_all_plots); gc()

                # Add stat, vars, region before saving 
                table <- add_column(table, region = region, .after = "EP")
                table <- add_column(table, var = var, .after = "region")
                table <- add_column(table, stat = stat, .after = "var")

                # Save
                setwd("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/global") # dir()
                save(x = table, file = paste("table_pred_QM_",method,"_",stat,"_",var,"_",region,".Rdata", sep = "") )
                # Remove some stuff
                rm(table); gc()

            } else if ( method == "monthly" ) {

                pred_QM_all_plots <- mclapply(plots, function(p) {
                        
                        # p <- plots[13]
                        message(paste("Performing quantile mapping for ",p," based on the ",method," approach", sep = ""))
                        sub_merged_df <- merged_df[merged_df$EP == p,]
                        
                        # Subset 'sub_EOBS_daily_stat' to have the full TS to project the QM on below
                        sub_ERA5 <- era5_ddf[era5_ddf$EP == p,]
                        
                        # Split this one into 12 monthly series to use as objects for QM prediction
                        monthly_ERA5_for_pred <- split(sub_ERA5, format(sub_ERA5$Date, "%m")) 
                        # monthly_ERA5_for_pred[[4]][1:30,] # very good
                        
                        # Remove 2008 because unreliable year for observations
                        sub_merged_df2 <- sub_merged_df[which(sub_merged_df$Date >= "2009-05-01"),]

                        # Split the data series frame by month
                        monthly_series <- split(sub_merged_df2, format(sub_merged_df2$Date, "%m"))
                        # monthly_series[[3]][1:20,] # very good
                        
                        # For each element m of 'monthly_series' (= month), perform QM
                        res_monthly_preds <- lapply(c(1:12), function(m) {
       
                                # Fit QM model on TS of corresponding month m
                                qm_model <- fitQmapRQUANT(
                                            obs = monthly_series[[m]][,3],
                                            mod = monthly_series[[m]][,4],
                                            qstep = 0.01,
                                            wet.day = FALSE
                                ) # eo - fitQmapRQUANT

                                # Predict on monthly_EOBS_for_pred[[m]]
                                pred <- doQmapRQUANT(x = monthly_ERA5_for_pred[[m]][,8], fobj = qm_model, type = "linear")

                                # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                                ddf <- data.frame(
                                            EP = p,
                                            Date = monthly_ERA5_for_pred[[m]][,"Date"],
                                            raw_value = monthly_ERA5_for_pred[[m]][,8],
                                            corrected_value = pred
                                ) # eo ddf
                                colnames(ddf) <- c("EP","Date","raw_value","corrected_value")
                        
                                # Clean & return ddf
                                rm(pred,qm_model)
                                return(ddf)

                            } # eo FUN 

                        ) # eo lapply - m

                        # Combine 'res_monthly_preds' back into one data frame & restore chronological order
                        table.mon <- dplyr::bind_rows(res_monthly_preds)
                        table.mon <- table.mon[order(table.mon$Date),]

                        rm(monthly_series,sub_merged_df2,monthly_ERA5_for_pred,sub_ERA5,sub_merged_df)
                        gc()

                        # Return
                        return(table.mon)

                    }, mc.cores = 25

                ) # eo mclapply - plots

                # Rbind and save in proper dir
                table <- dplyr::bind_rows(pred_QM_all_plots)
                rm(pred_QM_all_plots); gc()

                # Add stat, vars, region before saving 
                table <- add_column(table, region = region, .after = "EP")
                table <- add_column(table, var = var, .after = "region")
                table <- add_column(table, stat = stat, .after = "var")

                # Save
                setwd("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/monthly") # dir()
                save(x = table, file = paste("table_pred_QM_",method,"_",stat,"_",var,"_",region,".Rdata", sep = "") )
                # Remove some stuff
                rm(table); gc()

            } else if ( method == "mw" ) {
                
                # Define window size for the mowing window
                window_size <- w_size

                pred_QM_all_plots <- mclapply(plots, function(p) {

                        # p <- plots[1] # for testing
                        message(paste("Performing quantile mapping for ",p," based on the ",method," approach", sep = ""))
                        sub_merged_df <- merged_df[merged_df$EP == p,]

                        # Add vector of unique days of year (1–366) to subset
                        sub_merged_df$doy <- as.integer(format(sub_merged_df$Date, "%j"))
                        doys <- sort(unique(sub_merged_df$doy))

                        # Subset 'sub_EOBS_daily_stat' to have the full TS to project the QM on below
                        sub_ERA5 <- era5_ddf[era5_ddf$EP == p,]

                        # Add vector of unique days of year (1–366) to ERA5-Land data (full TS needed for prediction)
                        sub_ERA5$doy <- as.integer(format(sub_ERA5$Date, "%j")) # summary(sub_EOBS$doy)

                        # Use lapply to process each day-of-year (doy)
                        ### WARNING: moving window option takes a lot longer !
                        res_mw_preds <- lapply(doys, function(d) {
                                
                                    # Moving window: handle wrap-around at start/end of year
                                    window_days <- ((d - window_size):(d + window_size)) %% 366
                                    # if there is a 0 in the moving window, replace it 0 by 366 (loop around the year)
                                    if( 0 %in% window_days ) {
                                        window_days[window_days == 0] <- 366
                                    } # eo if loop - 0
  
                                    # Select those rows of 'sub_merged_df' that fall within the moving window
                                    in_window <- sub_merged_df$doy %in% window_days
  
                                    # Fit quantile mapping model
                                    qm_model <- fitQmapRQUANT(
                                            obs = sub_merged_df$obs[in_window],
                                            mod = sub_merged_df$ERA5_Land[in_window],
                                            qstep = 0.01,
                                            wet.day = FALSE
                                    ) # eo - fitQmapRQUANT

                                    # Apply model to the full DOY values from 'sub_EOBS'
                                    # First, identify rows of full EP's TS that correspond to DOY
                                    apply_idx <- which(sub_ERA5$doy == d)
                                    # Second, predict based on the QM model for these
                                    pred <- doQmapRQUANT(
                                            x = sub_ERA5[apply_idx,8],
                                            fobj = qm_model,
                                            type = "linear"
                                    ) # eo - doQmapRQUANT

                                    # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                                    ddf <- data.frame(
                                            EP = p,
                                            Date = sub_ERA5[apply_idx,4],
                                            raw_value = sub_ERA5[apply_idx,8],
                                            corrected_value = pred
                                    ) # eo ddf
                                    colnames(ddf) <- c("EP","Date","raw_value","corrected_value")
                        
                                    # Clean & return ddf
                                    rm(pred,qm_model,in_window,window_days,apply_idx)
                                    return(ddf)
    
                                } # eo FUN

                        ) # eo lapply - day of the year (doy)

                        # Combine 'res_monthly_preds' back into one data frame & restore chronological order
                        table_plot <- dplyr::bind_rows(res_mw_preds)
                        table_plot <- table_plot[order(table_plot$Date),]

                        rm(res_mw_preds,sub_merged_df)
                        gc()

                        # Return
                        return(table_plot)

                    }, mc.cores = 25
                
                ) # eo mclapply - plots
                
                # Rbind and save in proper dir
                table <- dplyr::bind_rows(pred_QM_all_plots)
                rm(pred_QM_all_plots); gc()

                # Add stat, vars, region before saving 
                table <- add_column(table, region = region, .after = "EP")
                table <- add_column(table, var = var, .after = "region")
                table <- add_column(table, stat = stat, .after = "var")

                # Save
                setwd("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/mw") # dir()
                save(x = table, file = paste("table_pred_QM_",method,"_",stat,"_",var,"_",region,".Rdata", sep = "") )
                # Remove some stuff
                rm(table); gc()

            } else if ( method == "anoms" ) {
            
                # Set training period
                train_start <- as.Date("2009-05-01")
                train_end <- max(merged_df$Date)

                # Run anomaly-based QM in parallel with mclapply()
                pred_QM_all_plots <- mclapply(plots, function(p) {
  
                        # p <- plots[sample(c(1:50),1,1)] # sample random EP for testing script below
                        message(paste("Performing quantile mapping for ",p," based on the ",method," approach", sep = ""))
                        
                        # EP-level ddf with Date/obs & model data 
                        sub_merged_df <- merged_df[merged_df$EP == p,]
                        # summary(sub_merged_df)

                        # Subset 'sub_EOBS_daily_stat' to have the full TS to priject the QM on below
                        sub_ERA5 <- era5_ddf[era5_ddf$EP == p,] # summary(sub_ERA5)
                        sub_ERA5 <- as.data.frame(sub_ERA5) # otherwise gets annoying a bit later on
                        
                        # Compute monthly means from 'sub_EOBS' to get the full TS of daily anomalies to the monthly mean
                        sub_ERA5$Month <- format(sub_ERA5$Date, "%m") # str(sub_ERA5)
                        monthly_means_mod <- tapply(X = sub_ERA5[,8], INDEX = sub_ERA5[,"Month"], FUN = mean, na.rm = TRUE)
                        # monthly_means_mod

                        # Compute anomalies to the month for the FULL model TS (for QM prediction later)
                        sub_ERA5$anom <- sub_ERA5[,8] - monthly_means_mod[sub_ERA5$Month]
                        # summary(sub_ERA5)
                        ### NOTE: sub_ERA5 is not to be used to train the QM model, just for the prediction

                        # Remove 2008 because unreliable year for observations - not to be used to train QM model
                        sub_merged_df2 <- sub_merged_df[which(sub_merged_df$Date >= "2009-05-01"),]
  
                        # Add month column for grouping
                        sub_merged_df2$Month <- format(sub_merged_df2$Date, "%m")
  
                        # Get training data (overlap)
                        sub_merged_df_train <- sub_merged_df2[sub_merged_df2$Date >= train_start & sub_merged_df2$Date <= train_end,]
  
                        # Compute monthly means from OBS in training period
                        monthly_means_obs <- tapply(sub_merged_df_train$obs, sub_merged_df_train$Month, mean, na.rm = TRUE)
                        # monthly_means_obs ; monthly_means_mod
  
                        # Compute anomalies for the post 2009 period
                        sub_merged_df_train$mod_anom <- sub_merged_df_train$ERA5_Land - monthly_means_mod[sub_merged_df_train$Month]
                        sub_merged_df_train$obs_anom <- sub_merged_df_train$obs - monthly_means_obs[sub_merged_df_train$Month]
  
                        # lapply over months
                        # m <- "06"
                        pred_QM_mon_anoms <- lapply(sprintf("%02d",1:12), function(m) {

                                # Subset obs anomalies for month 'm'
                                obs_anom_m <- sub_merged_df_train$obs_anom[sub_merged_df_train$Month == m]
                                # Subset obs anomalies for month 'm'
                                mod_anom_m <- sub_merged_df_train$mod_anom[sub_merged_df_train$Month == m]
                                # summary(obs_anom_m) ; summary(mod_anom_m)
                                # Dates to be corrected
                                dates <- sub_merged_df_train$Date[sub_merged_df_train$Month == m]
                                # IDs of the point corresponding to month 'm'
                                idx_month <- which(sub_merged_df2$Month == m)
                                
                                # Fit QM model        
                                qmodel <- fitQmapRQUANT(
                                                obs = obs_anom_m,
                                                mod = mod_anom_m,
                                                qstep = 0.01,
                                                wet.day = FALSE
                                ) # eo fitQmapRQUANT
                                
                                # Correct model-based anomalies with QM model
                                pred <- doQmapRQUANT(
                                                x = sub_ERA5[sub_ERA5$Month == m,"anom"],
                                                fobj = qmodel,
                                                type = "linear"
                                ) # eo doQmapRQUANT
                                # summary(mod_anom_m); summary(pred)
                                
                                # Add anomaly back to monthly mean
                                corrected_value_means_obs <- pred + monthly_means_obs[m]
                                # Or based on monthly_means_mod[m] - save both to test
                                corrected_value_means_mod <- pred + monthly_means_mod[m]

                                # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                                ddf <- data.frame(
                                                EP = p,
                                                Date = sub_ERA5[sub_ERA5$Month == m,"Date"],
                                                raw_value = sub_ERA5[sub_ERA5$Month == m,8],
                                                corrected_value_means_obs = corrected_value_means_obs,
                                                corrected_value_means_mod = corrected_value_means_mod
                                ) # eo ddf
                                # summary(ddf)
                                
                                # Return & clean
                                return(ddf)
                                rm(pred,qmodel,idx_month,obs_anom_m,mod_anom_m,corrected_value)

                            } # eo FUN
                            
                        ) # eo lapply - m

                        # Combine 'res_monthly_preds' back into one data frame & restore chronological order
                        table <- dplyr::bind_rows(pred_QM_mon_anoms)
                        table <- table[order(table$Date),]
                        # dim(table); str(table); summary(table)

                        rm(sub_merged_df2,sub_ERA5,sub_merged_df)
                        gc()

                        # Return
                        return(table)
  
                    }, mc.cores = 25
                    
                ) # eo mclapply - p in plots 

                # Rbind and save in proper dir
                table <- dplyr::bind_rows(pred_QM_all_plots)
                # dim(table) ; summary(table) ; table[13514:13714,]
                rm(pred_QM_all_plots); gc()

                # Add stat, vars, region before saving 
                table <- add_column(table, region = region, .after = "EP")
                table <- add_column(table, var = var, .after = "region")
                table <- add_column(table, stat = stat, .after = "var")

                # Save
                setwd("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/anomalies") # dir()
                save(x = table, file = paste("table_pred_QM_",method,"_",stat,"_",var,"_",region,".Rdata", sep = "") )
                # Remove some stuff
                rm(table); gc()

            } # eo if else loop - method

        } # eo if else loop - var %in% c("Ta_200","precipitation")

} # eo master FUN - quantile_mapper


### 29/04/25: Re-apply quantile_mapper() and save corrected climate reconstructions in their respective directories
# c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","precipitation","SM_10")
for(v in c("precipitation","SM_10")) {
    for(s in c("max","min")) {
        for(r in c("HND","SCH","SWA")) {
            for(m in c('global','monthly','mw','anoms')) {

                quantile_mapper(var = v, stat = s, region = r, method = m, w_size = 15)

            } # eo for loop - m
        } # eo for loop - r
    } # eo for loop - s
} # eo for loop - v


### ------------------------------------------------------------------------------------------------------------

### 21/01/26: Re-write quantile_mapper() to adapt it to forets precipitation 

# To test quantile_mapper() while you're writing it: 
#var <- "precipitation"
#stat <- "total"
#region <- "SCH"
#method <- "global"
#w_size = 15


quantile_forest_mapper <- function(var = "precipitation", stat = "total", region, method, w_size) {

        #' This function takes four arguments and returns a formatted data.frame:
        #' @param var the climate variable to process (character) - Default is "precipitation"
        #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character) - Default is "total"
        #' @param region Region name (character): 'SCH' or 'HND' or 'SWA'
        #' @param method Which quantile mapping approach to run (character): 'global' or 'monthly' or 'mw' or 'anoms'
        #' @param w_size The size (in days) of the moving window (integer)

        ## Go to local obs directory and extract TS of interes
        setwd(paste("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily/forests/",region, sep = ""))

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

        # Read the file containing your observed daily stats
        obs_daily_stat <- read.csv(files, h = T, sep = ",", dec = ".")

        # Sanity check
        if( exists("obs_daily_stat") == FALSE ) {
            stop(
              paste("!!! ERROR: Could not load observed ",paste(var,stat, sep = "_")," file for the ",region, sep = "")
            )
        } # eo if loop - sanity check
            
        ## Go to E-OBS dir and load their daily data too
        message(paste("\nCorrecting E-OBS with observed daily ",stat," ",var," for the ",region,
                    " based on the ",method," quantile mapping strategy\n", sep = ""))
            
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
            
        ### Retain forest EPs only 
        ep2keep <- unique( EOBS_daily_stat$EP[grepl("W",EOBS_daily_stat$EP)] ) # ep2keep
        EOBS_daily_stat <- EOBS_daily_stat[EOBS_daily_stat$EP %in% ep2keep,]
        # summary(EOBS_daily_stat)

        # Subset region fo interest in 'EOBS_daily_stat' (contains all EPs be default - my bad)
        if( region == "SCH" ) {
                # identify EPs to subset & subt them from 'EOBS_daily_stat'
                plots2subset <- unique(EOBS_daily_stat$EP)
                plots2subset <- plots2subset[grepl("S",plots2subset)]
                # subset
                sub_EOBS_daily_stat <- EOBS_daily_stat[EOBS_daily_stat$EP %in% plots2subset,] # unique(sub_EOBS_daily_stat$EP)
                # delete 'EOBS_daily_stat' - no use anymore
                rm(EOBS_daily_stat,plots2subset); gc()
        } else if( region == "HND" ) {
                # identify EPs to subset & subt them from 'EOBS_daily_stat'
                plots2subset <- unique(EOBS_daily_stat$EP)
                plots2subset <- plots2subset[grepl("H",plots2subset)]
                # subset
                sub_EOBS_daily_stat <- EOBS_daily_stat[EOBS_daily_stat$EP %in% plots2subset,] # unique(sub_EOBS_daily_stat$EP)
                # delete 'EOBS_daily_stat' - no use anymore
                rm(EOBS_daily_stat,plots2subset); gc()
        } else if( region == "SWA" ) {
                # identify EPs to subset & subt them from 'EOBS_daily_stat'
                plots2subset <- unique(EOBS_daily_stat$EP)
                plots2subset <- plots2subset[grepl("A",plots2subset)]
                # subset
                sub_EOBS_daily_stat <- EOBS_daily_stat[EOBS_daily_stat$EP %in% plots2subset,] # unique(sub_EOBS_daily_stat$EP)
                # delete 'EOBS_daily_stat' - no use anymore
                rm(EOBS_daily_stat,plots2subset); gc()
        } # eo if else loop - var
            
        # Re-name 
        colnames(sub_EOBS_daily_stat) <- c("EP","Date","value")

        # For each EP: match the daily data that overlap across all 3 sources and combine in a single ddf
        # Rename columns so names match across 3 data.frames
        colnames(obs_daily_stat)[3] <- paste(var,stat, sep = "_")
        colnames(sub_EOBS_daily_stat)[3] <- paste(var,stat, sep = "_")
        names <- c("EP","Date",paste(var,stat, sep = "_")) # vector of colnames to join the 3 ddf by

        # To make sure Date format is homogeneous across all 3 tables
        if( class(obs_daily_stat$Date) != "Date" ) {
            obs_daily_stat$Date <- as.Date(obs_daily_stat$Date)
        } # eo if loop

        if( class(sub_EOBS_daily_stat$Date) != "Date" ) {
            sub_EOBS_daily_stat$Date <- as.Date(sub_EOBS_daily_stat$Date)
        } # eo if loop
        # summary(obs_daily_stat); summary(sub_EOBS_daily_stat)
        # head(obs_daily_stat); head(sub_EOBS_daily_stat)

        # Merge using full_join iteratively
        merged_df <- reduce(
                list(obs_daily_stat[,names],
                    sub_EOBS_daily_stat[,names]
                ), full_join, by = c("EP","Date")
        ) # eo reduce

        # Adjust colnames
        colnames(merged_df)[c(3:4)] <- c("obs","E_OBS")
        # summary(merged_df)
            
        # Remove NAs for evaluation
        merged_df <- na.omit(merged_df)
        # Base data.frame with EP/Date/obs & model data

        # Vector of EP IDs - use to subset 'merged_df'
        plots <- unique(merged_df$EP) # plots

        ### Start if else loop based on QM method: 
        # if method == "global" --> Perform basic QM for each EP based on overlapping data
        # if method == "monthly" --> Split and join the data into 12 monthly TS, perform QM on those
        # if method == "mw" --> Using mowing window to smooth the threshold effects you may encounter
        # if method == "anoms" --> Correcting the anomalies to the mean instead of the raw values

        if( method == "global" ) {

            pred_QM_all_plots <- mclapply(plots, function(p) {
                        
                        # p <- plots[27]
                        message(paste("Performing quantile mapping for ",p," based on the ",method," approach", sep = ""))
                        sub_merged_df <- merged_df[merged_df$EP == p,]
                        
                        # Subset 'sub_EOBS_daily_stat' to have the full TS to project the QM on below
                        sub_EOBS <- sub_EOBS_daily_stat[sub_EOBS_daily_stat$EP == p,]
                        
                        # Remove 2008 because unreliable year for observations
                        # "2009-05-01" instead of "2009-01-01" because of temperature data in the SWA
                        # sub_merged_df2 <- sub_merged_df[which(sub_merged_df$Date >= "2009-05-01"),]
                        
                        ### In the case of precipitation, only train the QM model by excluding the days 
                        ### that are observed and modelled as dry (precip. = 0).
                        ### Then, correct the model data that are also non dry
                        if( var == "precipitation" ) {

                            # Fit and apply quantile mapping by excluding the days that are dry in both model and obs (wet.day = 0)
                            # Setting wet.days to 0 disables the wet day correction process.
                            # This means the function will not try to match the fraction of wet days between the modeled and observed data.
                            # Instead, any values below the specified threshold (which can be a numeric value or derived from the observed data)
                            # will be treated as zero.
                            # The result is that the modeled data is not adjusted to match the observed data's wet day frequency.
                            # This can be useful in situations where the wet day correction is not appropriate or desired.
                            # For example, it might be used to analyze data where the frequency of wet days is not of primary concern. 
                            qm_model <- fitQmapRQUANT(obs = sub_merged_df$obs, mod = sub_merged_df$E_OBS, qstep = 0.01, wet.day = 0)
                            
                            # NOTE: Here, I am using the mod and obs data that overlap in time of course (2009-2023 period).
                            # However, the daily model data I want to predict the model on (with doQmapRQUANT()) span a much longer period 
                            # (1950-2023).
                            # Yet, we do not want to do the doQmapRQUANT() on those model data correctly equal to 0 in the
                            # 2009-2023 period. Since the obs data do not go back before 2008, there is no way I can tell if the dry days
                            # (precipitation = 0) before 2009 are correct or not. We need to make sure to use doQmapRQUANT() on
                            # all mod data < 2009 plus only those days > 2009 that are wrongly modelled as dry days. 

                            # Indices for dates before and after 2009
                            idx_EOBS_pre2009 <- which(sub_EOBS[,"Date"] < as.Date("2009-05-01"))

                            # Get the indices of the pre 2009 ('idx_EOBS_pre2009') + the days 
                            # of the post-2009 model data that were positively modelled as 0 (dry days)
                            # To do so, identify the dates where obs and E-OBS are both equal to 0 and remove them with '!'.
                            # Then, use those to get the IDs from 'idx_EOBS_post2009'
                            dates2keep <- sub_merged_df[which(!(sub_merged_df$obs == 0 & sub_merged_df$E_OBS == 0)),"Date"]
                            
                            # Subset sub_EOBS to retain those days where 0 where NOT correctly modelled (= dates2keep) by their indices
                            sub_EOBS_dates2keep <- which(sub_EOBS$Date %in% dates2keep)

                            # Combine all pre-2009 + the selected post-2009 days to be corrected (idx_EOBS_pre2009 + sub_EOBS_dates2keep)
                            correct_indices <- c(idx_EOBS_pre2009, sub_EOBS_dates2keep)
                            
                            # Predict full TS based on 'qm_model' object using 'correct_indices'
                            pred <- doQmapRQUANT(x = sub_EOBS[correct_indices,3], fobj = qm_model, type = "linear")
                            sub_EOBS2 <- sub_EOBS # retain uncorrected E-OBS data
                            sub_EOBS2[correct_indices,3] <- pred # add corrected data at the right indices

                            # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                            ddf <- data.frame(
                                    EP = p,
                                    Date = sub_EOBS[,"Date"],
                                    raw_value = sub_EOBS[,3],
                                    corrected_value = sub_EOBS2[,3]
                            ) # eo ddf
                            ## Check corrected data against initial obs data quickly
                            # summary(sub_merged_df2$obs)
                            # summary(ddf[ddf$Date >= as.Date("2009-01-01"),"corrected_value"]) # EXTREMELY good agreement

                        } else {

                            # Fit and apply quantile mapping
                            qm_model <- fitQmapRQUANT(
                                    obs = sub_merged_df$obs,
                                    mod = sub_merged_df$E_OBS,
                                    qstep = 0.01,
                                    wet.day = FALSE
                            ) # eo fitQmapRQUANT
                        
                            # Predict full TS based on 'qm_model' object
                            pred <- doQmapRQUANT(
                                    x = sub_EOBS[,3],
                                    fobj = qm_model,
                                    type = "linear"
                            ) # eo doQmapRQUANT
                        
                            # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                            ddf <- data.frame(
                                    EP = p,
                                    Date = sub_EOBS[,"Date"],
                                    raw_value = sub_EOBS[,3],
                                    corrected_value = pred
                            ) # eo ddf

                        } # eo if else loop - var == precipitation
                        
                        # Clean & return ddf
                        rm(pred,qm_model,sub_EOBS,sub_merged_df)
                        return(ddf)

                }, mc.cores = 25

            ) # eo mclapply - plots

            # Rbind and save in proper dir ('/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/global')
            table <- dplyr::bind_rows(pred_QM_all_plots)
            # dim(table) ; summary(table) ; table[134050:134100,]
            rm(pred_QM_all_plots); gc()

            # Add stat, vars, region before saving 
            table <- add_column(table, region = region, .after = "EP")
            table <- add_column(table, var = var, .after = "region")
            table <- add_column(table, stat = stat, .after = "var")

            # Save
            setwd("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/global") # dir()
            save(x = table, file = paste("table_pred_QM_forest_",method,"_",stat,"_",var,"_",region,".Rdata", sep = ""))
            # Remove some stuff
            rm(table); gc()

        } else if ( method == "monthly" ) {

            pred_QM_all_plots <- mclapply(plots, function(p) {
                        
                        # p <- plots[3]; p
                        message(paste("Performing quantile mapping for ",p," based on the ",method," approach", sep = ""))
                        sub_merged_df <- merged_df[merged_df$EP == p,]
                        # dim(sub_merged_df); summary(sub_merged_df)
                        
                        # Subset 'sub_EOBS_daily_stat' to have the full TS to priject the QM on below
                        sub_EOBS <- sub_EOBS_daily_stat[sub_EOBS_daily_stat$EP == p,]
                        
                        # Split this one into 12 monthly series to use as objects for QM prediction
                        monthly_EOBS_for_pred <- split(sub_EOBS, format(sub_EOBS$Date, "%m")) 
                        # monthly_EOBS_for_pred[[4]][1:30,] # very good
                        
                        # Remove 2008 because unreliable year for observations - NOT NECESSARY FOR PRECIPITATION
                        # sub_merged_df2 <- sub_merged_df[which(sub_merged_df$Date >= "2009-05-01"),]

                        # Split the data series frame by month
                        monthly_series <- split(sub_merged_df, format(sub_merged_df$Date, "%m"))
                        # str(monthly_series); monthly_series[[10]][1:50,]
                        
                        # For each element m of 'monthly_series' (= month), perform QM
                        # m <- 4
                        res_monthly_preds <- lapply(c(1:12), function(m) {
       
                                if( var == "precipitation" ) {

                                    # Fit and apply quantile mapping
                                    qm_model <- fitQmapRQUANT(
                                                obs = monthly_series[[m]][,3],
                                                mod = monthly_series[[m]][,4],
                                                qstep = 0.01,
                                                wet.day = 0
                                    ) # eo fitQmapRQUANT
                            
                                    # Apply only on wet days in the model: get their index
                                    # idx_EOBS_pre2009 <- which(sub_EOBS[,"Date"] < as.Date("2009-01-01"))
                            
                                    # Get the indices of the pre 2009 ('idx_EOBS_pre2009') + the days 
                                    # of the post-2009 model data that were positively modelled as 0 (dry days)
                                    # To do so, identify the dates where obs and E-OBS are both equal to 0 and remove them with '!'.
                                    # Then, use those to get the IDs from 'idx_EOBS_post2009'
                                    # dates2keep <- monthly_series[[m]][which(!(monthly_series[[m]]$obs == 0 & monthly_series[[m]]$E_OBS == 0)),"Date"]
                            
                                    # Use those to identify the IDs 'sub_EOBS'
                                    # sub_EOBS_dates2keep <- which(sub_EOBS$Date %in% dates2keep)
                                    # duplicated_rows <- duplicated(sub_EOBS_dates2keep) ; unique(duplicated_rows)

                                    # Combine: all pre-2009 + selective post-2009
                                    # correct_indices <- c(idx_EOBS_pre2009, sub_EOBS_dates2keep)

                                    ### 21/05/25: Correcting the following issue: 
                                    ### I applied each month-specific quantile model to ALL relevant dates (not filtered by month)
                                    ### -> Ended up with 12 monthly corrections but for each single 'Date'
                                    ### -> Ensure each month’s model is only applied to dates in that month (format(Date,"%m") == m)

                                    # Limit everything to month `m`
                                    # Format the current month index (m) as a two-digit string (e.g., "01", "02", ..., "12")
                                    month_str <- sprintf("%02d", m)
                                    
                                    # Identify the row indices in sub_EOBS corresponding to the current month
                                    sub_EOBS_month_idx <- which(format(sub_EOBS$Date, "%m") == month_str)
                                    
                                    # Identify indices in sub_EOBS that are:
                                    #   (1) before 2009-01-01 (pre-2009), &
                                    #   (2) within the current month
                                    idx_EOBS_pre2009_m <- intersect(which(sub_EOBS$Date < as.Date("2009-01-01")), sub_EOBS_month_idx)
                                    
                                    # For post-2009 data, identify the dates where neither obs nor E-OBS is zero
                                    # These are the "wet days" that should be corrected
                                    dates2keep <- monthly_series[[m]][which(!(monthly_series[[m]]$obs == 0 & monthly_series[[m]]$E_OBS == 0)),"Date"]
                                    
                                    # From sub_EOBS, identify indices of the dates that: are in dates2keep & belong to the current month 
                                    sub_EOBS_dates2keep_m <- intersect(which(sub_EOBS$Date %in% dates2keep), sub_EOBS_month_idx)

                                    # Combine the pre-2009 dates (from the same month) and the selected post-2009 wet days
                                    # These are the rows on which quantile mapping will be applied for month `m`
                                    correct_indices <- c(idx_EOBS_pre2009_m, sub_EOBS_dates2keep_m)
                            
                                    # Predict full TS based on 'qm_model' object using 'correct_indices'
                                    pred <- doQmapRQUANT(
                                                x = sub_EOBS[correct_indices,3],
                                                fobj = qm_model,
                                                type = "linear"
                                    ) # eo - doQmapRQUANT
                                    # sub_EOBS2 <- sub_EOBS
                                    # sub_EOBS2[correct_indices,3] <- pred

                                    # data.frame above replaced by this:
                                    ddf <- data.frame(
                                                EP = p,
                                                Date = sub_EOBS[correct_indices, "Date"],
                                                raw_value = sub_EOBS[correct_indices,3],
                                                corrected_value = pred
                                    ) # eo ddf
                                    
                                    ## Check corrected data against initial obs data quickly
                                    # summary( monthly_series[[m]][,3] )
                                    # summary(ddf[ddf$Date >= as.Date("2009-01-01"),"corrected_value"])
                                    # sum(duplicated(ddf[,c("EP","Date")]))

                                } else {

                                    # Fit QM model on TS of corresponding month m
                                    qm_model <- fitQmapRQUANT(
                                                obs = monthly_series[[m]][,3],
                                                mod = monthly_series[[m]][,4],
                                                qstep = 0.01,
                                                wet.day = FALSE
                                    ) # eo fitQmapRQUANT

                                    # Predict on monthly_EOBS_for_pred[[m]]
                                    pred <- doQmapRQUANT(
                                                x = monthly_EOBS_for_pred[[m]][,3],
                                                fobj = qm_model,
                                                type = "linear"
                                    ) # eo doQmapRQUANT

                                    # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                                    ddf <- data.frame(
                                                EP = p,
                                                Date = monthly_EOBS_for_pred[[m]][,"Date"],
                                                raw_value = monthly_EOBS_for_pred[[m]][,3],
                                                corrected_value = pred
                                    ) # eo ddf

                                } # eo if else loop - var = precipitation 
                        
                                # Clean & return ddf
                                rm(pred,qm_model)
                                return(ddf)

                            } # eo FUN 

                        ) # eo lapply - m

                        # Combine 'res_monthly_preds' back into one data frame & restore chronological order
                        tible <- dplyr::bind_rows(res_monthly_preds)
                        tible <- tible[order(tible$Date),]

                        # Check for duplicates?
                        # duplicated_rows <- duplicated(tible) ; unique(duplicated_rows)
                        # sum(duplicated(tible[,c("EP","Date")])) # should return 0

                        rm(monthly_series,monthly_EOBS_for_pred,sub_EOBS,sub_merged_df)
                        gc()

                        # Return
                        return(tible)

                }, mc.cores = 25

            ) # eo mclapply - plots

            # Rbind and save in proper dir
            table <- dplyr::bind_rows(pred_QM_all_plots)
            # dim(table) ; summary(table) ; table[12221:12281,]
            rm(pred_QM_all_plots); gc()

            # Add stat, vars, region before saving 
            table <- add_column(table, region = region, .after = "EP")
            table <- add_column(table, var = var, .after = "region")
            table <- add_column(table, stat = stat, .after = "var")

            ### 22/05/25: Check for rows with only NAs
            # only_na_rows <- apply(table, 1, function(x) all(is.na(x)))
            # any(only_na_rows) # OK, should be good
                
            # Save
            setwd("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/monthly") # dir()
            save(x = table, file = paste("table_pred_QM_forest_",method,"_",stat,"_",var,"_",region,".Rdata", sep = "") )
            # Remove some stuff
            rm(table); gc()

        } else if ( method == "mw" ) {
                
            # Define window size for the mowing window
            window_size <- w_size

            pred_QM_all_plots <- mclapply(plots, function(p) {

                        # p <- plots[13]
                        message(paste("Performing quantile mapping for ",p," based on the ",method," approach", sep = ""))
                        sub_merged_df <- merged_df[merged_df$EP == p,]
                        
                        # Add vector of unique days of year (1–366) to subset
                        sub_merged_df$doy <- as.integer(format(sub_merged_df$Date, "%j"))
                        doys <- sort(unique(sub_merged_df$doy))
                        
                        # Subset 'sub_EOBS_daily_stat' to have the full TS to priject the QM on below
                        sub_EOBS <- sub_EOBS_daily_stat[sub_EOBS_daily_stat$EP == p,]

                        # Add vector of unique days of year (1–366) to EOBS data (full TS needed for prediction)
                        sub_EOBS$doy <- as.integer(format(sub_EOBS$Date, "%j"))

                        # Use lapply to process each day-of-year (doy)
                        ### WARNING: moving window option takes a lot longer !
                        res_mw_preds <- lapply(doys, function(d) {
                                
                                    # Moving window: handle wrap-around at start/end of year
                                    window_days <- ((d - window_size):(d + window_size)) %% 366
                                    
                                    # if there is a 0 in the moving window, replace it 0 by 366 (loop around the year)
                                    if( 0 %in% window_days ) {
                                        window_days[window_days == 0] <- 366
                                    } # eo if loop - 0
  
                                    # Select those rows of 'sub_merged_df' that fall within the moving window
                                    in_window <- sub_merged_df$doy %in% window_days
  
                                    # Fit quantile mapping model
                                    # if loop: if var == "precipitation" use wet.day = 0, otherwise wet.day = FALSE
                                    if( var == "precipitation" ) {
                                        
                                        qm_model <- fitQmapRQUANT(
                                            obs = sub_merged_df$obs[in_window],
                                            mod = sub_merged_df$E_OBS[in_window],
                                            qstep = 0.01,
                                            wet.day = 0
                                        ) # eo - fitQmapRQUANT

                                    } else {

                                        qm_model <- fitQmapRQUANT(
                                            obs = sub_merged_df$obs[in_window],
                                            mod = sub_merged_df$E_OBS[in_window],
                                            qstep = 0.01,
                                            wet.day = FALSE
                                        ) # eo - fitQmapRQUANT

                                    } # eo if else loop 
                                    
                                    # Apply model to the full DOY values from 'sub_EOBS'
                                    # First, identify rows of full EP's TS that correspond to DOY
                                    apply_idx <- which(sub_EOBS$doy == d)
                                    
                                    # Second, predict based on the QM model for these
                                    pred <- doQmapRQUANT(
                                            x = sub_EOBS[apply_idx,3],
                                            fobj = qm_model,
                                            type = "linear"
                                    ) # eo - doQmapRQUANT

                                    # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                                    ddf <- data.frame(
                                            EP = p,
                                            Date = sub_EOBS[apply_idx,2],
                                            DOY = d,
                                            raw_value = sub_EOBS[apply_idx,3],
                                            corrected_value = pred
                                    ) # eo ddf
                        
                                    # Clean & return ddf
                                    rm(pred,qm_model,in_window,window_days,apply_idx)
                                    return(ddf)
    
                                } # eo FUN

                        ) # eo lapply - day of the year (doy)

                        # Combine 'res_monthly_preds' back into one data frame & restore chronological order
                        table_plot <- dplyr::bind_rows(res_mw_preds)
                        table_plot <- table_plot[order(table_plot$Date),]
                        # table_plot[2721:2921,]
                        
                        # Clean and return
                        rm(res_mw_preds,sub_merged_df)
                        gc()

                        # Return
                        return(table_plot)

                }, mc.cores = 25
                
            ) # eo mclapply - plots
                
            # Rbind and save in proper dir
            table <- dplyr::bind_rows(pred_QM_all_plots)
            rm(pred_QM_all_plots); gc()

            # Add stat, vars, region before saving 
            table <- add_column(table, region = region, .after = "EP")
            table <- add_column(table, var = var, .after = "region")
            table <- add_column(table, stat = stat, .after = "var")
            # head(table); summary(table)

            # Save
            setwd("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/mw") # dir()
            save(x = table, file = paste("table_pred_QM_forest_",method,"_",stat,"_",var,"_",region,".Rdata", sep = "") )
            # Remove some stuff
            rm(table); gc()

        } else if ( method == "anoms" ) {
                
            # Set training period
            train_start <- as.Date("2009-01-01")
            train_end <- max(merged_df$Date)

            # Run anomaly-based QM in parallel with mclapply()
            pred_QM_all_plots <- mclapply(plots, function(p) {
  
                        # p <- plots[sample(c(1:50),1,1)] # sample random EP for testing script below
                        message(paste("Performing quantile mapping for ",p," based on the ",method," approach", sep = ""))
                        # EP-level ddf with Date/obs & model data 
                        sub_merged_df <- merged_df[merged_df$EP == p,]
                        # summary(sub_merged_df)

                        # Subset 'sub_EOBS_daily_stat' to have the full TS to priject the QM on below
                        sub_EOBS <- sub_EOBS_daily_stat[sub_EOBS_daily_stat$EP == p,]
                        
                        # Compute monthly means from 'sub_EOBS' to get the full TS of daily anomalies to the monthly mean
                        # Add month column for grouping
                        sub_EOBS$Month <- format(sub_EOBS$Date, "%m")
                        monthly_means_mod <- tapply(sub_EOBS[,3], sub_EOBS$Month, mean, na.rm = TRUE)
                        # monthly_means_mod

                        # Compute anomalies to the month for the FULL model TS (for QM prediction later)
                        sub_EOBS$anom <- sub_EOBS[,3] - monthly_means_mod[sub_EOBS$Month]
                        # summary(sub_EOBS)
                        ### NOTE: sub_EOBS is not to be used to train the QM model, just for the prediction

                        # Remove 2008 because unreliable year for observations - not to be used to train QM model
                        # sub_merged_df2 <- sub_merged_df[which(sub_merged_df$Date >= "2009-05-01"),]
  
                        # Add month column for grouping
                        sub_merged_df$Month <- format(sub_merged_df$Date, "%m")
  
                        # Get training data (overlap)
                        sub_merged_df_train <- sub_merged_df[sub_merged_df$Date >= train_start & sub_merged_df$Date <= train_end,]
  
                        # Compute monthly means from OBS in training period
                        monthly_means_obs <- tapply(sub_merged_df_train$obs, sub_merged_df_train$Month, mean, na.rm = TRUE)
                        # monthly_means_obs ; monthly_means_mod
  
                        # Compute anomalies for the post 2009 period
                        sub_merged_df_train$mod_anom <- sub_merged_df_train$E_OBS - monthly_means_mod[sub_merged_df_train$Month]
                        sub_merged_df_train$obs_anom <- sub_merged_df_train$obs - monthly_means_obs[sub_merged_df_train$Month]
  
                        # Loop over months (can later parallelize if needed)
                        # m <- "10"
                        pred_QM_mon_anoms <- lapply(sprintf("%02d",1:12), function(m) {

                                # Subset obs anomalies for month 'm'
                                obs_anom_m <- sub_merged_df_train$obs_anom[sub_merged_df_train$Month == m]
                                
                                # Subset obs anomalies for month 'm'
                                mod_anom_m <- sub_merged_df_train$mod_anom[sub_merged_df_train$Month == m]
                                # summary(obs_anom_m) ; summary(mod_anom_m)
                                
                                # Dates to be corrected
                                dates <- sub_merged_df_train$Date[sub_merged_df_train$Month == m]
                                
                                # IDs of the point corresponding to month 'm'
                                idx_month <- which(sub_merged_df$Month == m)
                                
                                # Fit QM model        
                                qmodel <- fitQmapRQUANT(
                                                obs = obs_anom_m,
                                                mod = mod_anom_m,
                                                qstep = 0.01,
                                                wet.day = FALSE
                                ) # eo fitQmapRQUANT
                                
                                # Correct model-based anomalies with QM model
                                pred <- doQmapRQUANT(
                                                x = sub_EOBS[sub_EOBS$Month == m,"anom"],
                                                fobj = qmodel,
                                                type = "linear"
                                ) # eo doQmapRQUANT
                                # summary(mod_anom_m); summary(pred)
                                
                                # Add anomaly back to monthly mean from obs or mod
                                # Here, we want to match observations as closely as possible so we add 'monthly_means_obs[m]'
                                corrected_value_means_obs <- pred + monthly_means_obs[m]
                                # Or monthly_means_mod[m] if wan to preserve model's climatology while correcting variability
                                corrected_value_means_mod <- pred + monthly_means_mod[m]	

                                # Make a ddf containing: EP, date, full uncorrected TS + correct TS based on the global QM
                                ddf <- data.frame(
                                                EP = p,
                                                Date = sub_EOBS[sub_EOBS$Month == m,"Date"],
                                                raw_value = sub_EOBS[sub_EOBS$Month == m,3],
                                                corrected_value_means_obs = corrected_value_means_obs,
                                                corrected_value_means_mod = corrected_value_means_mod
                                ) # eo ddf
                                # summary(ddf)

                                ### 25/04/25: In the case of total precipitation with the 'anoms' approach, 
                                ### you may ignore the 'wet.day' argument (leave it to FALSE) since you are 
                                ### working with anomalies.
                                ### HOWEVER, after adding the 'monthly_means_obs', you may never find values == 0 
                                ### (i.e., dry days) or even negative precip. values! 
                                ### -> Need to correct that here by concerting the smallest values > 0 (i.e., 0.02766) to 0
                                
                                ## Find said min value & replace by 0
                                if( var == "precipitation" ) {
                                    
                                    min_val <- min(ddf$corrected_value_means_obs[ddf$corrected_value_means_obs > 0]) # min_val
                                    
                                    # Convert them to 0 in 'ddf'
                                    ddf[ddf$corrected_value_means_obs == min_val,"corrected_value_means_obs"] <- 0

                                    ### 29/04/25: Same for v2 (model-based monthly means used to get final corrected values)
                                    min_val <- min(ddf$corrected_value_means_mod[ddf$corrected_value_means_mod > 0]) # min_val
                                    ddf[ddf$corrected_value_means_mod == min_val,"corrected_value_means_mod"] <- 0
                                
                                    ### If there are any negative values...convert to 0 too
                                    ddf[ddf$corrected_value_means_obs < 0,"corrected_value_means_obs"] <- 0
                                    ddf[ddf$corrected_value_means_mod < 0,"corrected_value_means_mod"] <- 0

                                } # eo var == "precipitation"
                                
                                # Return & clean
                                return(ddf)
                                rm(min_val,pred,qmodel,idx_month,obs_anom_m,mod_anom_m,
                                    corrected_value_means_obs,corrected_value_means_mod)

                            } # eo FUN
                            
                        ) # eo lapply - m

                        # Combine 'res_monthly_preds' back into one data frame & restore chronological order
                        table <- dplyr::bind_rows(pred_QM_mon_anoms)
                        table <- table[order(table$Date),]
                        # dim(table); str(table); summary(table)

                        rm(sub_merged_df2,sub_EOBS,sub_merged_df)
                        gc()

                        # Return
                        return(table)
  
                }, mc.cores = 25
                    
            ) # eo mclapply - p in plots 

            # Rbind and save in proper dir
            table <- dplyr::bind_rows(pred_QM_all_plots)
            # dim(table) ; summary(table) ; table[136050:136100,]
            rm(pred_QM_all_plots); gc()

            # Add stat, vars, region before saving 
            table <- add_column(table, region = region, .after = "EP")
            table <- add_column(table, var = var, .after = "region")
            table <- add_column(table, stat = stat, .after = "var")

            # Save
            setwd("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/anomalies") # dir()
            save(x = table, file = paste("table_pred_QM_forest_",method,"_",stat,"_",var,"_",region,".Rdata", sep = "") )
            # Remove some stuff
            rm(table); gc()

        } # eo if else loop

} # eo master FUN - quantile_forest_mapper


for(r in c("HND","SCH","SWA")) {
    for(m in c("global","monthly","mw","anoms")) {
        quantile_forest_mapper(var = "precipitation", stat = "total", region = r, method = m, w_size = 15)
    } # eo for loop - m
} # eo for loop - r


### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------