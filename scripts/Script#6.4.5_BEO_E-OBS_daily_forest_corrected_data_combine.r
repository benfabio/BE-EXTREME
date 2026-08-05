### ------------------------------------------------------------------------------------------------------------

### 19/01/26 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to combine final (i.e., corrected) the EP-level forest climate values with the local measurements
### to generate the final long-term reconstructions table, like you did for grasslands (R Script#6.4.4).
### Will use these 'final' TS for ECE detection and quantification as well as SPEI computation once the final 
### precipitation downscaling and corrections through quantile mapping are done  

### Recyling script 6.4.4 to write a master FUN that will:
### - Load the corrected daily TS data from R Script#4.5.3 ('correcte.r')
### - Load the corresponding local measurements at EP level
### - Combine them into one final TS (all < 2009-05 -> corrected outputs; all > 2009-05 -> local obs)

### WARNING: Flag and remove those EPs with uneven age managment -> unreliable age estimates -> no stand age correction -> remove

### Last update: 16/02/26 (Recycle code from R Script#6.4.3 (plot_daily_stat_comparison()) to evaluate final forest predictions)

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
library("metrica") # for RIA

# Dir to save the combined TS in:
files.dir <- "/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/forests"
# Grasslands' final reconstructed TS are in: "/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/grasslands"

### ------------------------------------------------------------------------------------------------------------

### 1°) Load one of the grasslands' final files and examine their structure again so your forets files follow the same one
# setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/grasslands") #; dir()
# test <- get(load("table_combined_obs+corr_global_max_Ta_200_grasslands_SWA.Rdata"))
# dim(test); head(test)
# str(test)
#'data.frame':   1369700 obs. of  5 variables:
# $ EP         : chr  "AEG01" "AEG01" "AEG01" "AEG01" ...
# $ Date       : Date, format: "1950-01-01" "1950-01-02" ...
# $ obs        : num  NA NA NA NA NA NA NA NA NA NA ...
# $ corr_model : num  -1.588 4.732 1.221 0.744 2.78 ...
# $ final_value: num  -1.588 4.732 1.221 0.744 2.78 ...

### Too simple for the forests though. Lacking some key metadata like tree type, main tree species, age and so on...
### Go to 'files.dir' and examine which metadata to keep in the final files
# setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables") # dir()
# test2 <- readRDS("table_predict_GAMM_corrected_max_Ta_200_15_01_26.rds")
# dim(test2); head(test2)
# str(test2)
# Columns to add as metadata in the final full TS files
# StandAge
# TreeType
# ForestType (aka 'Type')
# MainTreeSpecies ('MTS')

### -> Then, split per 'Region' ('SWA', 'HND' and 'SCH') based on 'EP' and save separately under this nomenclature:
### "table_combined_obs+corr_max_Ta_200_forests_SWA.Rdata"

### ----------------------------------------------------

### 2°) Write master FUN to combine fial corrected values with original gap-filled and corrected observed values.
### Inspired from combine_quantile_maps() in R Script#6.4.4 but simpler since no quantile mapping directly involed here. 

# For testing combine_corrected_forest_values() while you write it: 
# stat = "max"
# var = "Ta_200"

combine_corrected_forest_values <- function(stat, var) {

        #' This function takes four arguments and returns a formatted data.frame:
        #' @param var the climate variable to process (character) - one of the following: 
        #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20", "Ts_50", "precipitation" or "SM_10"
        #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character): 'mean', 'max' or 'min'
        #' @return A formatted data.frame combining the daily statistics

        message(paste("Loading the modelled and observational TS of ",paste(stat,var, sep = " "), sep = ""))

        ### 1. Load the final corrected reconstructed data
        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables") # dir()
        corr <- readRDS(paste("table_predict_GAMM_filled_",stat,"_",var,"_10_02_26.rds", sep = ""))

        ### 2. Load the initial gap-filled corrected observational data to combine to the modelled historic data 
        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data")
        obs <- readRDS(paste("table_aggr_mean_offsets+metadata_",stat,"_",var,"_28_10_25.rds", sep = ""))
        colnames(obs)[5] <- "Value" # str(obs)
        # Add 'Date' to 'obs' based on Year + DOY
        obs <- obs %>% mutate(Date = as.Date(DOY - 1, origin = paste0(Year,"-01-01")))
        # summary(obs$Date); str(obs$Date)

        ### 3. Combine them and make sure they follow the same format as the graslands' tables + metadata
        # str(obs); str(corr)
        final <- data.frame(
            EP = corr$EP,
            Date = corr$Date,
            model_value = corr$Mean_Forest_value, # will need to put 'obs' before this one 
            final_value = corr$Mean_Forest_value,
            Stand_age = corr$StandAge,
            Tree_type = corr$TreeType,
            Main_tree_species = corr$MTS,
            Forest_type = corr$Type
        ) # eo ddf
        # dim(final); summary(final)

        # Merge 'obs'
        final <- final %>% full_join(obs[,c("EP","Date","Value")], by = c("EP","Date")) # str(final)
        # Change 'Value' to 'obs' and relocate before 'corr_model'
        colnames(final)[9] <- "obs"
        # Relocate
        final <- final %>% relocate("obs", .before = "model_value")
        # Integrate 'obs' values into 'final_value'
        final[!is.na(final$obs),"final_value"] <- final[!is.na(final$obs),"obs"]
        # summary(final) ; head(final)

        ### 4. Split 'final' into 3 regional tables. Will save all 3 in 3 different tables, like for the grasslands data
        # Add 'Region' based on 'EP'
        final$Region <- NA
        for(i in unique(final$EP)) {
            message(i)
            # if else loop
            if( grepl("HEW",i) ) {
                final[final$EP == i,"Region"] <- "HND"
            } else if( grepl("AEW",i) ) {
                final[final$EP == i,"Region"] <- "SWA"
            } else if( grepl("SEW",i) ) {
                final[final$EP == i,"Region"] <- "SCH"
            } # eo if els eloop
        } # eo for loop
        # unique(final$Region); summary(factor(final$Region))
        # Relocate
        final <- final %>% relocate("Region", .before = "EP")
        
        ### 5. Sanity checks and exclude forets EPs with uneven age management (cf. conversations with Forest Core Team)
        
        # Delete rows with only NAs = keep rows where at least one column is not NA
        final <- final %>% filter(if_any(everything(), ~ !is.na(.)))
        # (summary(factor(final$Forest_type))/nrow(final))*100
        
        # Remove obs if any Stand_age is < 0 
        if( nrow(final[final$Stand_age < 1,]) > 1 ) {
            message(paste("!!! WARNING: Some EPs show Stange ages < 0 - Correcting"))
            final <- final[final$Stand_age >= 1,]
        } # eo if loop

        # Remove the 26 EPs of the HND that do not have an age class ('AC') management system
        # Other words: 24 EPs of the HND that show an '_AC_' Forest_type
        types2keep <- unique(final$Forest_type)[grepl("_AC_",unique(final$Forest_type))]
        ep2keep <- unique( final[(final$Forest_type %in% types2keep) & final$Region == "HND","EP"] )
        # Check EPs that do not have 'AC' management system
        # unique( final[!(final$Forest_type %in% types2keep),"EP"] )
        # unique( final[!(final$Forest_type %in% types2keep),"Forest_type"] )
        ### Save only those HND EPs that belong to 'ep2keep'

        ### 6. Save and clean
        # Go to files.dir
        message(paste("Saving the final regional reconstruction of ",paste(stat,var, sep = " "),"\n", sep = ""))
        setwd(files.dir)
        saveRDS(final[final$Region == "SWA" & !is.na(final$Region),], file = paste("table_combined_obs+corr_",stat,"_",var,"_forests_SWA.rds", sep = ""))
        saveRDS(final[final$Region == "HND" & !is.na(final$Region) & final$EP %in% ep2keep,], file = paste("table_combined_obs+corr_",stat,"_",var,"_forests_HND.rds", sep = ""))
        saveRDS(final[final$Region == "SCH" & !is.na(final$Region),], file = paste("table_combined_obs+corr_",stat,"_",var,"_forests_SCH.rds", sep = ""))

        # Clean and move to next item
        rm(final,obs,corr,types2keep,ep2keep); gc()

} # eo FUN - combine_corrected_forest_values

# Test
# combine_corrected_forest_values(stat = "max", var = "Ta_200")

### Apply combine_corrected_forest_values() in a double for loop, as usual
for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")) {
    for(s in c("max","min")) {
        combine_corrected_forest_values(var = v, stat = s)
    } # eo for loop - s
} # eo for loop - v

### ----------------------------------------------------

### 3°) Checking outputs saved above. Plot TS of the reconstructed forest data to check for long-term trends etc.

### 3.A.1) Compute final evaluation metrics to report for each regin and variable

### In order to be consistent with the QM outputs for the grasslands, compute the following evaluation metrics for the forests:
# - mean bias error (MBE)
# - mean absolute error (MAE)
# - RMSE
# - corr coeff
# - refined index of agreement (RIA)

### Add a switch to compute metrics on anomalies ot the monthly mean (anoms == TRUE)

### Evaluation metrics in:
eval.metric.dir <- "/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables/evaluation_tables"

# To test evaluate_quantile_maps while you're writing it: 
# var = "Ta_10"
# stat = "max"
# anoms = TRUE

evaluate_corrected_forests <- function(var, stat, anoms = FALSE) {

        #' This function takes four arguments and returns a formatted data.frame:
        #' @param var the climate variable to process (character) - one of the following: 
        #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20", "Ts_50", "precipitation" or "SM_10"
        #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character): 'mean', 'max' or 'min'
        #' @param anoms Whether evaluation metrics should be computed on anomalies to the monthly mean (BOOLEAN). Default = FALSE
        #' @return A formatted data.frame combining the daily statistics

        ## Go to directory and extract TS of interest
        setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/forests") # dir()
        # Read the file containing your observed daily stats, for all 3 regions
        files <- dir()[grepl(paste("corr_",stat,"_",var, sep = ""),dir())]
        res <- lapply(files, function(f) {
                dat <- readRDS(f); return(dat)
            } # eo FUN
        ) # eo 
        # Rbind
        table <- dplyr::bind_rows(res)
        rm(res); gc()
        # str(table); dim(table); head(table); summary(table)

        # Sanity check
        if( exists("table") == FALSE ) {
            stop(
              paste("!!! ERROR: Could not load observed ",paste(var,stat, sep = "_")," file", sep = "")
            )
        } # eo if loop - sanity check

        # Remove NAs for evaluation (usually data < 2008 of course)
        table <- na.omit(table)
        # str(table)
        # Convert EP and Region to factors
        table$EP <- as.factor(table$EP)
        table$Region <- as.factor(table$Region)

        # Vector of EPs
        plots <- unique(table$EP)

        ### Calculate evaluation metrics: mean bias, corr coef, RMSE, mean absolute bias (MAE) for each EP separately
        if( anoms ) {

          # Compute monthly means for each data source
          table_anoms <- table %>%
              mutate(month = format(Date, "%Y-%m")) %>%
              group_by(EP,month) %>%
              mutate(
                  mon_mean_obs = mean(obs, na.rm = TRUE),
                  mon_mean_mod = mean(model_value, na.rm = TRUE)
              ) %>% 
              ungroup()
          # summary(table_anoms)

          # Compute monthly anomalies
          table_anoms <- table_anoms %>%
              mutate(
                  mon_anom_obs = obs - mon_mean_obs,
                  mon_anom_mod = model_value - mon_mean_mod
              )
          # summary(table_anoms)
            
          # Subset merged_df per EP and return evaluation metrics
          eval_metrics <- lapply(plots, function(p) {
                  
                  # p <- plots[29]
                  message(paste("Computing evaluation metrics for EP: ",p," based on monthly anomalies", sep = ""))
                  ## Subset
                  sub_table_anoms <- table_anoms[table_anoms$EP == p,]
                  ## Only keep Dates > '2009-05-01' because of obs data
                  sub_table_anoms <- sub_table_anoms[sub_table_anoms$Date >= '2009-05-01',]
                  ## Mean bias error (MBE)
                  mbe <- mean( (sub_table_anoms$mon_anom_mod) - (sub_table_anoms$mon_anom_obs) )
                  ## Mean absolute error (MAE)
                  mae <- mean(abs(sub_table_anoms$mon_anom_mod - sub_table_anoms$mon_anom_obs))
                  ## Root mean square error (RMSE)
                  rmse <- sqrt(mean((sub_table_anoms$mon_anom_mod - sub_table_anoms$mon_anom_obs)^2, na.rm = TRUE))
                  ## Corr coeff
                  corr_coeff <- cor(sub_table_anoms$mon_anom_obs, sub_table_anoms$mon_anom_mod, use = "complete.obs")
                  ## RIA
                  ria <- metrica::d1r(obs = sub_table_anoms$mon_anom_obs, pred = sub_table_anoms$mon_anom_mod, tidy = TRUE)
                      
                  # Return evaluation metrics in a data.frame
                  eval_metrics <- data.frame(
                        Region = unique(sub_table_anoms$Region),
                        EP = p,
                        mbe = mbe,
                        mae = mae, 
                        rmse = rmse, 
                        corr = corr_coeff,
                        ria = ria
                  ) # eo ddf
                  # Discard subset
                  rm(sub_table_anoms); gc()
                  # Return
                  return(eval_metrics)

            } # eo FUN - plots

          ) # eo lapply - eval_metrics
          # Rbind
          table_eval_metrics <- dplyr::bind_rows(eval_metrics)
          rm(eval_metrics); gc()
          # summary(table_eval_metrics) # quick check

          # Return and clean
          table_eval_metrics$variable <- var
          table_eval_metrics$stat <- stat
            
          # Save in dir as .Rdata
          message(paste("Saving evaluation metrics table for ",stat," ",var, sep = ""))
          setwd(eval.metric.dir)
          save(x = table_eval_metrics, file = paste("table_evaluation_metrics_anoms_",stat,"_",var,"_forests_10.02.26.Rdata", sep = ""))
          rm(table_anoms,table_eval_metrics); gc()

        } else {

          # Subset merged_df per EP and return evaluation metrics
          eval_metrics <- lapply(plots, function(p) {

                  # p <- sample(plots,1,1); p # for testing
                  message(paste("Computing evaluation metrics for ",p, sep = ""))
                  
                  ## subset
                  sub_table <- table[table$EP == p,]
                  ## Only keep Dates > '2009-05-01' because of obs data
                  sub_table <- sub_table[sub_table$Date >= '2009-05-01',]
                  ## Mean bias error (MBE)
                  mbe <- mean(sub_table$model_value - sub_table$obs)
                  ## Mean absolute error (MAE)
                  mae <- mean(abs(sub_table$model_value - sub_table$obs))
                  ## Root mean square error (RMSE)
                  rmse <- sqrt(mean((sub_table$model_value - sub_table$obs)^2, na.rm = TRUE))
                  ## Corr coeff
                  corr_coeff <- cor(sub_table$obs, sub_table$model_value, use = "complete.obs")
                  ## RIA
                  ria <- metrica::d1r(obs = sub_table$obs, pred = sub_table$model_value, tidy = TRUE)
                      
                  # Return evaluation metrics in a data.frame
                  eval_metrics <- data.frame(
                        Region = unique(sub_table$Region),
                        EP = p,
                        mbe = mbe,
                        mae = mae,
                        rmse = rmse,
                        corr = corr_coeff,
                        ria = ria
                  ) # eo ddf
                  # Discard subset
                  rm(sub_table); gc()
                  # Return
                  return(eval_metrics)

            } # eo FUN - plots

          ) # eo lapply - eval_metrics

          # Rbind
          table_eval_metrics <- dplyr::bind_rows(eval_metrics)
          rm(eval_metrics); gc()
          # summary(table_eval_metrics) # quick check

          # Return and clean
          table_eval_metrics$variable <- var
          table_eval_metrics$stat <- stat
            
          message(paste("\n","Saving evaluation metrics table for ",stat," ",var, sep = ""))
          setwd(eval.metric.dir)
          save(x = table_eval_metrics, file = paste("table_evaluation_metrics_",stat,"_",var,"_forests_10.02.26.Rdata", sep = ""))
          rm(table,table_eval_metrics); gc()

      } # eo if else loop - anomalies

} # eo master FUN - evaluate_corrected_forests


### Apply evaluate_quantile_maps() for all data (variables, stats, QM method, regions...)
### 10-11/02/26: Re-do with newest corrected predictions
for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")) {
    for(s in c("max","min")) {
        evaluate_corrected_forests(var = v, stat = s, anoms = FALSE)
        evaluate_corrected_forests(var = v, stat = s, anoms = TRUE)
    } # eo for loop - s
} # eo for loop - v


### 20/01/26: Examine distribution of of evalutaion metrics per region and variable

### 10-11/02/26: Re-do with newest corrected predictions
setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables/evaluation_tables/")


### 3.A.1) Examine metrics based on normal values
files <- dir()[!grepl("_anoms_",dir())]
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
table <- bind_rows(res)
rm(res); gc()
summary(table)

# Plots
table$var <- factor(paste(table$stat,table$variable, sep = "_")) # unique(table$var)
# Reorder: 
table$var <- factor(table$var, levels = c("max_Ta_10","min_Ta_10","max_Ta_200","min_Ta_200",
        "max_Ts_05","min_Ts_05","max_Ts_10","min_Ts_10","max_Ts_20","min_Ts_20","max_SM_10","min_SM_10"))

# summary(table$var)
p <- ggplot(data = table, aes(x = Region, y = d1r, fill = factor(Region))) + geom_violin() + 
    geom_boxplot(width = .2, colour = "black", fill = "white") +
    scale_fill_brewer(name = "Region", palette = "Paired") +
    xlab("") + ylab("Refined index of agreement") +
    theme_bw() + facet_wrap(.~var, ncol = 4, nrow = 3, scales = "free_y")
ggsave(plot = p, filename = "boxplots_RIA_corr_forests_10.02.26.jpg", dpi = 300, width = 6, height = 6)


p <- ggplot(data = table, aes(x = Region, y = rmse, fill = factor(Region))) + geom_violin() + 
    geom_boxplot(width = .2, colour = "black", fill = "white") +
    scale_fill_brewer(name = "Region", palette = "Paired") +
    xlab("") + ylab("RMSE") + theme_bw() + facet_wrap(.~var, ncol = 4, nrow = 3, scales = "free_y")
ggsave(plot = p, filename = "boxplots_RMSE_corr_forests_10.02.26.jpg", dpi = 300, width = 6, height = 6)


p <- ggplot(data = table, aes(x = Region, y = corr, fill = factor(Region))) + geom_violin() + 
    geom_boxplot(width = .2, colour = "black", fill = "white") +
    scale_fill_brewer(name = "Region", palette = "Paired") +
    xlab("") + ylab("Correlation coefficient") + theme_bw() + facet_wrap(.~var, ncol = 4, nrow = 3, scales = "free_y")
ggsave(plot = p, filename = "boxplots_CORR_corr_forests_10.02.26.jpg", dpi = 300, width = 6, height = 6)


p <- ggplot(data = table, aes(x = Region, y = mbe, fill = factor(Region))) + geom_violin() + 
    geom_boxplot(width = .2, colour = "black", fill = "white") +
    scale_fill_brewer(name = "Region", palette = "Paired") +
    xlab("") + ylab("Mean bias error") + theme_bw() + facet_wrap(.~var, ncol = 4, nrow = 3, scales = "free_y")
ggsave(plot = p, filename = "boxplots_MBE_corr_forests_10.02.26.jpg", dpi = 300, width = 6, height = 6)

### All very good fit except SM_10 which show poorer performance although still OK. 
### Note on RIA: RIA = 0 means that the magnitude of model errors equals the magnitude 
### of deviations from the observed mean -> no better prediction than using the mean
### In practice, RIA = 0.8 means the model accounts for 80% of the observed variation
### in a statistically meaningful way; RIA < 0.5 signals significant issues.

### 10-11/02/26: Summarize range of evaluation metrics based on raw values
# - max Ta_10: mean bias error ranging between -0.5 and +0.5; all correlation coefficients > 0.965; RMSE between 1.6 and 2.4 across all EPs; all RIA > 0.86
# - min Ta_10: mean bias error ranging between -1 and +2; all correlation coefficients > 0.95; RMSE between 1 and 2.5 across all EPs; all RIA > 0.84 (max = 0.92)
# - max Ta_200: mean bias error ranging between -0.5 and +1; all correlation coefficients > 0.990; RMSE between 0.8 and 1.4 across all EPs; all RIA > 0.92
# - min Ta_200: mean bias error ranging between -1 and +2; all correlation coefficients > 0.98; RMSE between 0.5 and 1.5 across all EPs; all RIA > 0.89 (max = 0.95)
# - max Ts_05: mean bias error ranging between -0.5 and +1; all correlation coefficients > 0.975; RMSE between 0.5 and 1.7 across all EPs; all RIA > 0.85
# - min Ts_05: mean bias error ranging between -0.5 and +1; all correlation coefficients > 0.98; RMSE between 0.5 and 1.5 across all EPs; all RIA > 0.85
# - max Ts_10: mean bias error ranging between -0.5 and +1; all correlation coefficients > 0.98; RMSE between 0.5 and 1.7 across all EPs; all RIA > 0.85
# - min Ts_10: mean bias error ranging between -0.5 and +1; all correlation coefficients > 0.98; RMSE between 0.5 and 1.7 across all EPs; all RIA > 0.85
# - max Ts_20: mean bias error ranging between -0.5 and +1; all correlation coefficients > 0.99; RMSE between 0.5 and 1.25 across all EPs; all RIA > 0.85
# - min Ts_20: mean bias error ranging between -0.5 and +0.8; all correlation coefficients > 0.99; RMSE between 0.5 and 1.25 across all EPs; all RIA > 0.84
# - max SM_10 and min SM_10: mean bias error ranging between -7.5% and +7%; correlation coefficients ranging between 0.5 and 0.8; RMSE between 5 and 10; 
#   RIA > 0 (between 0.25 and 0.75 for 2 regions: HND and SWA), BUT RIA mostly < 0 for the SCH

### Interpretation
### For temperatures: Excellent performance; No evidence of structural artefacts; Capping StandAge solved the long-term issue
### For SM10: Moderate skill overall; Region-dependent performance -> 1 region (SCH) seems problematic
### This suggests that your model structure is appropriate for temperature.
### Soil moisture likely needs additional structure (maybe precipitation lags, cumulative rainfall, or soil-specific modifiers)...


### ----------------------------------------------------


### 3.A.2) Examine metrics based on anomalies ot the monthy mean
files <- dir()[grepl("_anoms_",dir())]
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
table <- bind_rows(res)
rm(res); gc()
summary(table)

# Plots
table$var <- factor(paste(table$stat,table$variable, sep = "_")) # unique(table$var)
# Reorder: 
table$var <- factor(table$var, levels = c("max_Ta_10","min_Ta_10","max_Ta_200","min_Ta_200",
        "max_Ts_05","min_Ts_05","max_Ts_10","min_Ts_10","max_Ts_20","min_Ts_20","max_SM_10","min_SM_10"))

p <- ggplot(data = table, aes(x = Region, y = d1r, fill = factor(Region))) + geom_violin() + 
    geom_boxplot(width = .2, colour = "black", fill = "white") +
    scale_fill_brewer(name = "Region", palette = "Paired") +
    xlab("") + ylab("Refined index of agreement") +
    theme_bw() + facet_wrap(.~var, ncol = 4, nrow = 3, scales = "free_y")
ggsave(plot = p, filename = "boxplots_RIA_anoms_corr_forests_20.01.26.jpg", dpi = 300, width = 6, height = 6)


p <- ggplot(data = table, aes(x = Region, y = rmse, fill = factor(Region))) + geom_violin() + 
    geom_boxplot(width = .2, colour = "black", fill = "white") +
    scale_fill_brewer(name = "Region", palette = "Paired") +
    xlab("") + ylab("RMSE") + theme_bw() + facet_wrap(.~var, ncol = 4, nrow = 3, scales = "free_y")
ggsave(plot = p, filename = "boxplots_RMSE_anoms_corr_forests_20.01.26.jpg", dpi = 300, width = 6, height = 6)


p <- ggplot(data = table, aes(x = Region, y = corr, fill = factor(Region))) + geom_violin() + 
    geom_boxplot(width = .2, colour = "black", fill = "white") +
    scale_fill_brewer(name = "Region", palette = "Paired") +
    xlab("") + ylab("Correlation coefficient") + theme_bw() + facet_wrap(.~var, ncol = 4, nrow = 3, scales = "free_y")
ggsave(plot = p, filename = "boxplots_CORR_anoms_corr_forests_20.01.26.jpg", dpi = 300, width = 6, height = 6)


p <- ggplot(data = table, aes(x = Region, y = mbe, fill = factor(Region))) + geom_violin() + 
    geom_boxplot(width = .2, colour = "black", fill = "white") +
    scale_fill_brewer(name = "Region", palette = "Paired") +
    xlab("") + ylab("Mean bias error") + theme_bw() + facet_wrap(.~var, ncol = 4, nrow = 3, scales = "free_y")
ggsave(plot = p, filename = "boxplots_MBE_anoms_corr_forests_20.01.26.jpg", dpi = 300, width = 6, height = 6)


### Observations:
### Evaluation metrics based on anomalies to the monthly mean still quite good
### (scores usually between 0.8 and 0.9, even for RIA). max and min SM_10 show much lower scores again
### but RIA still around 0.6 which means decent performance.
### -> I would avoid using SM_10 reconstructions for extremes. Use only for exploration maybe?


### 10-11/02/26: Summarize range of evaluation metrics based on anomalies to the monthly means
# - max Ta_10:  mean bias error all very close to 0; most correlation coefficients > 0.92 (max = 0.96); RMSE between 1.2 and 1.4; most RIA > 0.81
# - min Ta_10:  mean bias error all very close to 0; most correlation coefficients > 0.85;  RMSE between 0.9 and 1.5; most RIA > 0.725 (max = 0.85)
# - max Ta_200: mean bias error all very close to 0; most correlation coefficients > 0.98;  RMSE between 0.65 and 0.9; most RIA > 0.91
# - min Ta_200: mean bias error all very close to 0; most correlation coefficients > 0.925; RMSE between 0.5 and 1.5; most RIA > 0.8 (max = 0.95)
# - max Ts_05:  mean bias error all very close to 0; most correlation coefficients > 0.875; RMSE between 0.5 and 0.8; most RIA > 0.74 (max = 0.825)
# - min Ts_05:  mean bias error all very close to 0; most correlation coefficients > 0.925; RMSE between 0.35 and 0.55; most RIA > 0.80 (max = 0.87)
# - max Ts_10:  mean bias error all very close to 0; most correlation coefficients > 0.90;  RMSE between 0.4 and 0.55; most RIA > 0.78 (max = 0.86)
# - min Ts_10:  mean bias error all very close to 0; most correlation coefficients > 0.92;  RMSE between 0.35 and 0.45; most RIA > 0.80 (max = 0.875)
# - max Ts_20:  mean bias error all very close to 0; most correlation coefficients > 0.93;  RMSE between 0.27 and 0.4; most RIA > 0.825
# - min Ts_20:  mean bias error all very close to 0; most correlation coefficients > 0.95;  RMSE between 0.25 and 0.4; most RIA > 0.825
# - max SM_10 and min SM_10: mean bias error all very close to 0; most correlation coefficients ranging between 0.2 and 0.7; 
#   RMSE between 1 and 3; most RIA > 0 (between 0.5 and 0.65 for 2 regions: HND and SWA), BUT RIA mostly < 0.5 for the SCH, sometimes even < 0

### Interpretation
### The model captures true temporal variability in temperature.
### It does not rely on seasonal smooths to inflate skill.
### It reproduces event-scale deviations realistically.
### The age cap did not degrade performance.
### The GAMM framework seems structurally viable and dynamically skillful for thermal offset modelling.

### Bottom-line conclusions. The capped-age GAMM:
### -> Produces stable long-term hindcasts
### -> Accurately reproduces seasonal structure
### -> Successfully captures temperature anomaly dynamics
### -> Shows physically consistent depth-dependent soil temperature behaviour
### -> Has moderate to low region-specific skill for soil moisture variability -> I wouldn't trust the SM_10 reconstructions so much for the forests

### Methods summary: 
### "To assess the model’s ability to reproduce temporal variability independently of seasonal mean biases, we evaluated performance using 
### monthly anomalies computed separately for observations and model predictions. For both datasets, anomalies were defined relative to 
### their respective monthly climatologies over the training period (2009–2024). This way we isolate interannual and synoptic-scale
### variability by removing differences in mean seasonal structure between observations and predictions. Model skill was quantified using
### the correlation coefficient (r), root mean square error (RMSE), mean bias error (MBE), and the refined index of agreement (RIA).
### For air temperature (Ta_10, Ta_200), anomaly correlations were consistently high (generally r > 0.85 and frequently > 0.92,
### reaching up to 0.98 for Ta_200), with low RMSE and negligible mean bias. RIA values were predominantly above 0.8, indicating strong agreement
### in both amplitude and timing of temperature deviations. Soil temperature anomalies (Ts_05–Ts_20) showed similarly robust performance,
### with correlations typically between 0.88 and 0.95 and decreasing RMSE with depth, reflecting physically consistent thermal buffering in deeper soil layers.
### In contrast, soil moisture anomalies (SM_10) exhibited moderate to low and region-dependent skill (r ≈ 0.2–0.7).
### While mean bias remained negligible, correlations and RIA values were lower, particularly in the SCH region, suggesting that short-term
### hydrological variability is only partially captured by the current model structure.

### Overall, the anomaly-based evaluation demonstrates that the GAMM framework successfully reproduces temperature variability beyond the seasonal cycle,
### confirming that predictive skill is not solely driven by climatological structure but reflects realistic dynamic responses to climate fluctuations."

### -------------------------------------------------

### 16/02/26: Preparing Group Meeting of the next GA in Wernigerode :) 

### 3.A.3) Recycle code from R Script#6.4.3 (plot_daily_stat_comparison()) to evaluate final forest predictions,
### for each region x var x stat, in the same style of plot than for the grasslands' QM predictions

# To test while you're writing it:
r <- "HND" 
v <- "Ta_10"
s <- "max"
mon.anoms <- FALSE
boxp <- TRUE
histo <- TRUE
scatt <- TRUE
time_series <- TRUE
biases <- TRUE 
std.anoms <- TRUE


plot_daily_stat_comparison <- function(region, var, stat, mon.anoms, boxp, histo, scatt, time_series) {

    #' This function takes 9 arguments and returns a formatted data.frame:
    #' @param var the climate variable to process (character) - one of the following: 
    #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20", "Ts_50", "precipitation" or "SM_10"
    #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character): 'mean', 'max' or 'min' 
    #' @param region Region name (character): 'SCH' or 'HND' or 'SWA'
    #' @param mon.anoms Switch (BOOLEAN) - Whether plots should be made based on anomalies to the monly mean 
    #' @param boxp Switch (BOOLEAN) - Whether boxplots should be made
    #' @param histo Switch (BOOLEAN) - Whether gghistograms should be made
    #' @param scatt Switch (BOOLEAN) - Whether scaterr plots should be made
    #' @param time_series Switch (BOOLEAN) - Whether time series plots should be made
    #' @return the plots selected through the arguments above

      ## Go to local obs directory and extract TS of interest
      setwd(paste("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily/forests/",region, sep = ""))

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

      # Read the file containing your observed daily stats
      message(paste("Loading local observations of ",stat," ",var, sep = ""))
      obs_daily_stat <- read.csv(files, h = T, sep = ",", dec = ".")
      colnames(obs_daily_stat)[3] <- paste(var,stat, sep = "_")

      # Sanity check
      if( exists("obs_daily_stat") == FALSE ) {
          stop(
            paste("!!! ERROR: Could not load observed ",paste(var,stat, sep = "_")," file for the ",region, sep = "")
          )
      } # eo if loop - sanity check

      ## Go to quantile mapping outputs directory and load corresponding dataset
      message(paste("Loading predictions of ",stat," ",var, sep = ""))
      setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables")
      file <- dir()[grepl(paste("filled",s,v, sep = "_"), dir())]
      corr_data <- readRDS(file) # dim(corr_data) ; str(corr_data) ; colnames(corr_data)

      # To make sure Date format is homogeneous across all 3 tables
      if( class(obs_daily_stat$Date) != "Date" ) {
          obs_daily_stat$Date <- as.Date(obs_daily_stat$Date)
      } # eo if loop

      if( class(corr_data$Date) != "Date" ) {
          corr_data$Date <- as.Date(corr_data$Date)
      } # eo if loop

      # Define names to merge by
      names <- c("EP","Date",paste(var,stat, sep = "_"))

      # Adjust colname in 'corr_data'
      colnames(corr_data)[21] <- paste(var,stat, sep = "_")

      # Subset region of interest: 
      corr_data <- corr_data[corr_data$Region == r,]

      # Merge
      message(paste("Merging both datasets...", sep = ""))
      merged_df <- reduce(
                    list(obs_daily_stat[,names],
                         corr_data[,names]
                    ), full_join, by = c("EP","Date")
      ) # eo reduce
      # summary(merged_df) # OK
      rm(obs_daily_stat,corr_data) ; gc()

      # Adjust colnames
      colnames(merged_df)[c(3,4)] <- c("obs","model")
        
      # Remove NAs for evaluation (usually data < 2008 of course)
      merged_df <- na.omit(merged_df)

      # Vector of EP IDs - use to subset 'merged_df'
      plots <- unique(merged_df$EP)

      ## 26/05/25: Restrict to Dates > '2009-05-01' because of faulty obs 2008 data
      merged_df <- merged_df[merged_df$Date >= '2009-05-01',]

      ### SWITCHES:
      ### If anoms == TRUE --> make plots based on anomalies to the monthly mean
      ### else --> make plots on normal data
      message(paste("Drawing plots\n", sep = ""))
      if( mon.anoms ) {

          # Compute monthly means for each data source
          merged_df <- merged_df %>%
                mutate(month = format(Date, "%Y-%m")) %>%
                group_by(EP,month) %>%
                mutate(
                    mon_mean_obs = mean(obs, na.rm = TRUE),
                    mon_mean_mod = mean(model, na.rm = TRUE)
                ) %>% 
                ungroup()

          # Compute monthly anomalies
          merged_df <- merged_df %>%
                mutate(
                    mon_anom_obs = obs - mon_mean_obs,
                    mon_anom_mod = model - mon_mean_mod
                )
            
          ### Save plots in plot.dir
          colors_sources <- c("obs" = "#00AFBB", "model" = "#E7B800")

          # 1st: boxplots 
          if( boxp ) {
                # Boxplot of stat distribution - needs melting
                m_merged_df <- reshape2::melt(merged_df[,c("EP","Date","mon_anom_obs","mon_anom_mod")], id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Re-name source levels so they match the 'colors_sources' palette above; unique(m_merged_df$source)
                m_merged_df$source <- as.character(m_merged_df$source)
                m_merged_df[m_merged_df$source == "mon_anom_obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "mon_anom_mod","source"] <- "model"
                # Make the plot
                plot <- ggplot(data = m_merged_df, aes(x = factor(source), y = value, fill = factor(source))) +
                  geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = .2) + 
                  geom_hline(yintercept = 0, linetype = "dashed") + 
                  scale_fill_manual(values = colors_sources) + labs(y = paste(stat,var,sep = " "), x = "") + 
                  theme_bw() + theme(legend.position = "none") +
                  ggtitle(paste("Distribution of ",stat," ",var," monthly anomalies in the ",region, sep = ""))
                # Save plot 
                ggsave(plot = plot, filename = paste("boxplot_evaluation_mon_anomalies_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 5, height = 5)  
          } # eo if loop - boxp

          # 2nd: histograms - also need melting
          if( histo ) {
                # Needs melting too
                m_merged_df <- reshape2::melt(merged_df[,c("EP","Date","mon_anom_obs","mon_anom_mod")], id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Re-name source levels so they match the 'colors_sources' palette above; unique(m_merged_df$source)
                m_merged_df$source <- as.character(m_merged_df$source)
                m_merged_df[m_merged_df$source == "mon_anom_obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "mon_anom_mod","source"] <- "model"
                # Make the plot
                plot <- gghistogram(m_merged_df, x = "value", add = "mean", rug = FALSE,
                          fill = "source", palette = colors_sources,
                          ylab = "Count", xlab = paste(stat,var,sep = " "),
                          title = paste("Distribution of ",stat," ",var," monthly anomalies in the ",region, sep = "")
                )
                # Save plot in dir  
                ggsave(plot = plot, filename = paste("gghistogram_evaluation_mon_anomalies_",stat,"_",var,"_",region,".jpg",sep = ""), dpi = 300, width = 7, height = 5)  
          } # eo if loop - histo

          # 3rd: scatt - make 2: obs vs corrected model
          if( scatt ) {
                # Make plot 1
                plot <- ggplot(merged_df, aes(x = mon_anom_obs, y = mon_anom_mod)) + geom_point(alpha = 0.2, colour = "grey50") +
                        geom_abline(intercept = 0, slope = 1, color = "#d53e4f", linetype = "dashed") + 
                        geom_smooth(method = "lm", se = TRUE, colour = "#00AFBB") + 
                        geom_hline(yintercept = 0, linetype = "dashed") + 
                        geom_vline(xintercept = 0, linetype = "dashed") + 
                        labs(title = paste("Local vs model daily ",stat," ",var,"\nmonthly anomalies", sep = ""),
                            x = "Local measurements anomalies", y = "model anomalies") +
                        theme_bw()
                # Save
                ggsave(plot = plot, filename = paste("scatter_evaluation_mon_anomalies_",stat,"_",var,"_",region,".jpg",sep = ""), dpi = 300, width = 5, height = 5)  
          } # eo if loop - scatt

          # 4th: time_series
          if( time_series ) {
                p <- ggplot(merged_df, aes(x = Date)) + 
                  geom_line(aes(y = mon_anom_obs), linewidth = 1, colour = "#00AFBB", alpha = .25) +
                  geom_line(aes(y = mon_anom_mod), linewidth = 1, colour = "#E7B800", alpha = .25) +
                  geom_hline(yintercept = 0, linetype = "dashed") + 
                  labs(title = paste("Time series of ",stat," ",var," monthly anomalies in the ",region," (all EPs)", sep = ""),
                      x = "Date", y = paste(stat,var,sep = " ")) +
                  theme_minimal()
                # Save
                 ggsave(plot = p, filename = paste("time_series_evaluation_mon_anomalies_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 10, height = 4)  
          } # eo if loop - time_series

      } else {

          ### Save plots depending on the switches activated as function's arguments: boxp, histo, scatt, time_series
          colors_sources <- c("obs" = "#00AFBB", "model" = "#E7B800")

          # 1st: boxplots 
          if( boxp ) {
                # Boxplot of stat distribution - needs melting
                m_merged_df <- reshape2::melt(merged_df, id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Re-name source levels so they match the 'colors_sources' palette above; unique(m_merged_df$source)
                m_merged_df$source <- as.character(m_merged_df$source)
                m_merged_df[m_merged_df$source == "obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "model","source"] <- "model"
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
                m_merged_df[m_merged_df$source == "obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "model","source"] <- "model"
                # Make the plot
                plot <- gghistogram(m_merged_df, x = "value", add = "mean", rug = FALSE,
                          fill = "source", palette = colors_sources,
                          ylab = "Count", xlab = paste(stat,var,sep = " "),
                          title = paste("Distribution of ",stat," ",var," in the ",region, sep = "")
                )
                # Save plot in dir  
                ggsave(plot = plot, filename = paste("gghistogram_evaluation_",stat,"_",var,"_",region,".jpg",sep = ""), dpi = 300, width = 7, height = 5)  
          } # eo if loop - histo

          # 3rd: scatt - make 2: obs vs corrected model outputs
          if( scatt ) {
                # Make plot 1
                plot <- ggplot(merged_df, aes(x = obs, y = model)) + geom_point(alpha = 0.2, colour = "grey50") +
                        geom_abline(intercept = 0, slope = 1, color = "#d53e4f", linetype = "dashed") + 
                        geom_smooth(method = "lm", se = TRUE, colour = "#00AFBB") + 
                        labs(title = paste("Local vs model daily ",stat," ",var, sep = ""),
                            x = "Local measurements", y = "Model (corrected)") +
                        theme_bw()
                # Save
                ggsave(plot = plot, filename = paste("scatter_evaluation_",stat,"_",var,"_",region,".jpg",sep = ""), dpi = 300, width = 5, height = 5)  
          } # eo if loop - histo

          # 4th: time_series
          if( time_series ) {
                p <- ggplot(merged_df, aes(x = Date)) + 
                  geom_line(aes(y = obs), linewidth = 1, colour = "#00AFBB", alpha = .25) +
                  geom_line(aes(y = model), linewidth = 1, colour = "#E7B800", alpha = .25) +
                  labs(title = paste("Time series of ",stat," ",var," in the ",region," (all EPs)", sep = ""),
                      x = "Date", y = paste(stat,var,sep = " ")) +
                  theme_minimal()
                # Save
                 ggsave(plot = p, filename = paste("time_series_evaluation_",stat,"_",var,"_",region,".jpg",sep = ""), dpi = 300, width = 10, height = 4)  
            } # eo if loop - time_series

          } # eo if else loop - anoms == T

} # eo plotting FUN - plot_daily_stat_comparison


for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")) {
    for(s in c("max","min")) {
        for(r in c("HND","SCH","SWA")) {
                
                plot_daily_stat_comparison(
                      region = r,
                      var = v,
                      stat = s,
                      mon.anoms = FALSE,
                      boxp = TRUE,
                      histo = TRUE,
                      scatt = TRUE,
                      time_series = TRUE
                )

        } # eo for loop - r
    } # eo for loop - s
} # eo for loop - v


### -------------------------------------------------

### 3.A.4) Plotting predictions' quality (RMSE, corr. coeff, FIA? MBE) across variables - for the GA 2026

## Go to dir 
setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables/evaluation_tables")

## For evaluation meteics based on raw values (main plot)
files <- dir()[!grepl("anoms",dir())]
res <- lapply(files, function(f) { d <- get(load(f)) ; return(d) } )
table <- bind_rows(res)
# dim(table); str(table)
## Add variable name 
table$var <- factor(paste(table$stat,table$variable, sep = " "))
# Need to reverse after computing order
table$var <- with(table, reorder(var, d1r, FUN = median))
table$var <- factor(table$var, levels = rev(levels(table$var)))
p <- ggplot(table, aes(x = d1r, y = var)) + geom_boxplot(fill = "#66bd63") + 
    xlab("RIA") + ylab("Variables") + theme_minimal()
ggsave(plot = p, filename = "boxplots_rankings_variables_forests_raw_RIA_16.02.26.jpg", dpi = 300, width = 5.5, height = 5)


## For evaluation meteics based on anomalies to the monthly mean (suppl. plot)
files <- dir()[grepl("anoms",dir())]
res <- lapply(files, function(f) { d <- get(load(f)) ; return(d) } )
table <- bind_rows(res)
## Add variable name 
table$var <- factor(paste(table$stat,table$variable, sep = " "))
# Need to reverse after computing order
table$var <- with(table, reorder(var, d1r, FUN = median))
table$var <- factor(table$var, levels = rev(levels(table$var)))
p <- ggplot(table,aes(x = d1r, y = var)) + geom_boxplot(fill = "#66bd63") + 
    xlab("RIA") + ylab("Variables") + theme_minimal()
ggsave(plot = p, filename = "boxplots_rankings_variables_forests_mon_anoms_RIA_16.02.26.jpg", dpi = 300, width = 5.5, height = 5)


### 16/02/26: BONUS - Same plots as above but for the grasslands
setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/evaluation_metrics")

## For evaluation meteics based on raw values (main plot)
files <- dir()[!grepl("anoms",dir())]
files <- files[grepl("mw",files)] # pick those based on the chosen 'mw' QM strategy
# files
res <- lapply(files, function(f) { d <- get(load(f)) ; return(d) } )
table <- bind_rows(res)
# dim(table); str(table)
## Add variable name 
table$var <- factor(paste(table$stat,table$variable, sep = " ")) # unique(table$var)
# Need to reverse after computing order
table$var <- with(table, reorder(var, d1r, FUN = median))
table$var <- factor(table$var, levels = rev(levels(table$var)))
p <- ggplot(table, aes(x = d1r, y = var)) + geom_boxplot(fill = "#fdae61") + 
    xlab("RIA") + ylab("Variables") + theme_minimal()
ggsave(plot = p, filename = "boxplots_rankings_variables_grasslands_raw_RIA_16.02.26.jpg", dpi = 300, width = 5.5, height = 5)


## For evaluation meteics based on anomalies to the monthly mean
files <- dir()[grepl("anoms",dir())]
files <- files[grepl("mw",files)] # pick those based on the chosen 'mw' QM strategy
# files
res <- lapply(files, function(f) { d <- get(load(f)) ; return(d) } )
table <- bind_rows(res)
# dim(table); str(table)
## Add variable name 
table$var <- factor(paste(table$stat,table$variable, sep = " ")) # unique(table$var)
# Need to reverse after computing order
table$var <- with(table, reorder(var, d1r, FUN = median))
table$var <- factor(table$var, levels = rev(levels(table$var)))
p <- ggplot(table, aes(x = d1r, y = var)) + geom_boxplot(fill = "#fdae61") + 
    xlab("RIA") + ylab("Variables") + theme_minimal()
ggsave(plot = p, filename = "boxplots_rankings_variables_grasslands_anoms_RIA_16.02.26.jpg", dpi = 300, width = 5.5, height = 5)



### -------------------------------------------------

### 3.B) Plot TS of the final corrected daily values for each region and variable
setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/forests") # dir()

# For testing
# v = "precipitation"
# s = "max"

for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10","precipitation")) {

    for(s in c("max","min")) {

        if( v == "precipitation" ) {
            files <- dir()[grepl(paste("mw_total_precipitation", sep = ""),dir())]
            # Load files
            res <- lapply(files, function(f) { dat <- get(load(f)); return(dat) } ) # eo lapply
            # Rbind
            table <- dplyr::bind_rows(res)
            rm(res); gc()
            # str(table); summary(table)
            # Add Region to 'table'
            table$Region <- NA
            for(i in unique(table$EP)) {
                # if else loop
                message(i)
                if( grepl("HEW",i) ) {
                    table[table$EP == i,"Region"] <- "HND"
                } else if( grepl("AEW",i) ) {
                    table[table$EP == i,"Region"] <- "SWA"
                } else if( grepl("SEW",i) ) {
                    table[table$EP == i,"Region"] <- "SCH"
                } # eo if els eloop
            } # eo for loop 
            # unique(table$Region)
        } else {
            files <- dir()[grepl(paste("corr_",s,"_",v, sep = ""),dir())]
            # Load files
            res <- lapply(files, function(f) { dat <- readRDS(f); return(dat) } ) # eo lapply
            # Rbind
            table <- dplyr::bind_rows(res)
            rm(res); gc()
            # str(table); summary(table)
        } # eo if else loop - v
        
        # Plot TS per Region after computing mean value + sd
        avg <- table %>% group_by(Date,Region) %>% summarise(mean = mean(final_value, na.rm = TRUE), sd = sd(final_value, na.rm = TRUE))
        # summary(avg)
        avg <- avg %>% group_by(Region) %>% mutate(MEAN = mean(mean)) %>% ungroup() # head(avg); str(avg)
        
        p <- ggplot(data = avg) + geom_ribbon(aes(x = Date, ymin = mean-sd, ymax = mean+sd), fill = "grey") + 
            geom_line(aes(x = Date, y = mean), color = "black", linewidth = .2) + 
            geom_smooth(aes(x = Date, y = mean), color = "#d53e4f", method = "gam") + 
            geom_hline(aes(yintercept = MEAN), colour = "#3288bd", linetype = "dashed") + 
            xlab("") + theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6)) +
            scale_x_date(date_breaks = "1 year", labels = date_format("%Y"), expand = c(0,0)) +
            facet_wrap(.~Region, nrow = 3, ncol = 1)
        
        # Save
        if(v == "precipitation") {
            ggsave(plot = p, filename = paste("plot_TS_final_forest_value_total_precipitation_11.02.26.jpg", sep = ""), dpi = 300, width = 10, height = 6.5)
        } else {
            ggsave(plot = p, filename = paste("plot_TS_final_forest_value_",s,"_",v,"_11.02.26.jpg", sep = ""), dpi = 300, width = 10, height = 6.5)
        } # eo if else loop - 

    } # eo for loop - s
    
} # eo for loop - v

### -> HUGE linear increase in soil temperatures since 1950!

### 10/02/26: Already plotted for newest predictions - no need to run


### -------------------------------------------------

### 26/01/26: 3.C) Plot the UNCORRECTED daily TS too and compare. May be needed to check why we have this enormous warming trends in Ts...

# v = "Ts_10"
# s = "max"

#for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")) {
#    for(s in c("max","min")) {
#
#        file <- dir()[grepl(paste("filled_",s,"_",v, sep = ""),dir())] # file
#        table <- readRDS(file)
#        # str(table); summary(table)
#
#        # Plot TS per Region after computing mean value + sd
#        avg <- table %>% group_by(Date,Region) %>% summarise(mean = mean(Mean_Forest_value, na.rm = TRUE), sd = sd(Mean_Forest_value, na.rm = TRUE))
#        avg <- avg %>% group_by(Region) %>% mutate(MEAN = mean(mean)) %>% ungroup()
#        # head(avg); str(avg); summary(avg)
#        
#        p <- ggplot(data = avg) + geom_ribbon(aes(x = Date, ymin = mean-sd, ymax = mean+sd), fill = "grey") + 
#            geom_line(aes(x = Date, y = mean), color = "black", linewidth = .2) + 
#            geom_smooth(aes(x = Date, y = mean), color = "#d53e4f", method = "gam") + 
#            geom_hline(aes(yintercept = MEAN), colour = "#3288bd", linetype = "dashed") + 
#            xlab("") + theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6)) +
#            scale_x_date(date_breaks = "1 year", labels = date_format("%Y"), expand = c(0,0)) +
#            facet_wrap(.~Region, nrow = 3, ncol = 1)
#        
#        # Save
#        ggsave(plot = p, filename = paste("plot_TS_uncorrected_forest_value_",s,"_",v,"_26.01.26.jpg", sep = ""), dpi = 300, width = 10, height = 6.5)
#
#    } # eo for loop - s
#} # eo for loop - v


### Observations: Uncorrected Ts data still seem to warm up too strongly

### Compute warming rates per decade for max Ta_10, max Ta_200 and Ts_10 and compare against literature 
library("tidyr")
library("broom")

setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables")
files <- dir()[grepl(paste("corrected_","max","_","Ts_20", sep = ""),dir())]
# files
table <- readRDS(files)
# str(table); summary(table)
ts_annual <- table %>%
  group_by(EP,Year) %>%
  summarise(
    MeanTemp = mean(Mean_Forest_corr, na.rm = TRUE),
    .groups = "drop"
)
# summary(ts_annual)

### Estimate warming rates (linear trend, °C per decade)

## NOTE: Sen’s slope estimates the trend as the median of all pairwise slopes between observations in a time series.
## That means:
## - No single extreme year can dominate the trend,
## - The estimate reflects the central tendency of change over time.
## This is fundamentally different from OLS, which minimises squared residuals and is highly sensitive to outliers.

## Sen’s slope:
## - Is non-parametric -> No assumption of normality or homoscedasticity
## - Is robust to outliers -> One or two extreme years cannot create spurious multi-degree trends
## - Handles monotonic but noisy trends well -> Exactly what climate warming signals look like
## - Is widely used in climate & soil literature -> Especially when paired with the Mann–Kendall test

# "Trends were estimated using Sen’s slope, a non-parametric and outlier-robust estimator widely applied in
# climatological and soil temperature studies, as it provides reliable trend estimates for noisy,
# autocorrelated environmental time series without assuming normally distributed residuals."

# install.packages("trend")
library("trend")
ts_trends_sen <- data.frame(
    ts_annual %>%
    group_by(EP) %>%
    summarise(
        trend_C_per_decade = sens.slope(MeanTemp,Year)$estimates * 10,
        p_value = sens.slope(MeanTemp, Year)$p.value,
        .groups = "drop"
    )
) # eo ddf
# ts_trends_sen
summary(ts_trends_sen)

### Per variable:
## max Ta_10
#      EP            trend_C_per_decade    p_value         
# Length:150         Min.   :0.05907    Min.   :0.000e+00  
# Class :character   1st Qu.:0.33437    1st Qu.:0.000e+00  
# Mode  :character   Median :1.05307    Median :0.000e+00  
#                    Mean   :0.82949    Mean   :1.398e-02  
#                    3rd Qu.:1.11274    3rd Qu.:3.300e-07  
#                    Max.   :1.31332    Max.   :6.100e-01

## max Ta_200
#      EP            trend_C_per_decade    p_value         
# Length:150         Min.   :0.1156     Min.   :0.000e+00  
# Class :character   1st Qu.:0.3486     1st Qu.:0.000e+00  
# Mode  :character   Median :0.6596     Median :0.000e+00  
#                    Mean   :0.5471     Mean   :1.698e-02  
#                    3rd Qu.:0.6835     3rd Qu.:1.560e-06  
#                    Max.   :0.7326     Max.   :5.047e-01 

## max Ts_05
#      EP            trend_C_per_decade    p_value         
# Length:150         Min.   :0.2810     Min.   :0.000e+00  
# Class :character   1st Qu.:0.7747     1st Qu.:0.000e+00  
# Mode  :character   Median :1.0253     Median :0.000e+00  
#                    Mean   :0.9615     Mean   :3.019e-06  
#                    3rd Qu.:1.1359     3rd Qu.:0.000e+00  
#                    Max.   :1.3494     Max.   :1.751e-04  

## max Ts_10
#      EP            trend_C_per_decade    p_value        
# Length:150         Min.   :0.1088     Min.   :0.000000  
# Class :character   1st Qu.:0.7103     1st Qu.:0.000000  
# Mode  :character   Median :1.2582     Median :0.000000  
#                    Mean   :1.0894     Mean   :0.001895  
#                    3rd Qu.:1.3892     3rd Qu.:0.000000  
#                    Max.   :1.6306     Max.   :0.113111  

## max Ts_20
#      EP            trend_C_per_decade    p_value         
# Length:150         Min.   :0.2933     Min.   :0.000e+00  
# Class :character   1st Qu.:0.3601     1st Qu.:0.000e+00  
# Mode  :character   Median :0.3849     Median :0.000e+00  
#                    Mean   :0.3942     Mean   :1.455e-09  
#                    3rd Qu.:0.4207     3rd Qu.:0.000e+00  
#                    Max.   :0.5761     Max.   :7.511e-08 

### -------------------------------------------------

### 21/01/26: 3.D) Make the same plots but for the grasslands and compare 
setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/grasslands"); dir()

# For testing
# v = "Ta_200"
# s = "max"

for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")) {
    for(s in c("max","min")) {

        files <- dir()[grepl(paste("mw_",s,"_",v, sep = ""),dir())]

        res <- lapply(files, function(f) {
                dat <- get(load(f)); return(dat)
            } # eo FUN
        ) # eo lapply
        # Rbind
        table <- dplyr::bind_rows(res)
        rm(res); gc()
        # str(table); summary(table)

        # Add 'Region' based on EP
        table$Region <- NA
        for(i in unique(table$EP)) {
            message(i)
            # if else loop
            if( grepl("HEG",i) ) {
                table[table$EP == i,"Region"] <- "HND"
            } else if( grepl("AEG",i) ) {
                table[table$EP == i,"Region"] <- "SWA"
            } else if( grepl("SEG",i) ) {
                table[table$EP == i,"Region"] <- "SCH"
            } # eo if else loop
        } # eo for loop
        # unique(table$Region); summary(factor(table$Region))
        # Relocate
        table <- table %>% relocate("Region", .before = "EP")

        # Plot TS per Region after computing mean value + sd
        avg <- table %>% group_by(Date,Region) %>% summarise(mean = mean(final_value, na.rm = TRUE), sd = sd(final_value, na.rm = TRUE))
        avg <- avg %>% group_by(Region) %>% mutate(MEAN = mean(mean)) %>% ungroup() 
        # summary(avg); head(avg); str(avg)
        
        p <- ggplot(data = avg) + geom_ribbon(aes(x = Date, ymin = mean-sd, ymax = mean+sd), fill = "grey") + 
            geom_line(aes(x = Date, y = mean), color = "black", linewidth = .2) + 
            geom_smooth(aes(x = Date, y = mean), color = "#d53e4f", method = "gam") + 
            geom_hline(aes(yintercept = MEAN), colour = "#3288bd", linetype = "dashed") + 
            xlab("") + theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6)) +
            scale_x_date(date_breaks = "1 year", labels = date_format("%Y"), expand = c(0,0)) +
            facet_wrap(.~Region, nrow = 3, ncol = 1)
        
        # Save
        ggsave(plot = p, filename = paste("plot_TS_final_grassland_value_",s,"_",v,"_21.01.26.jpg", sep = ""), dpi = 300, width = 10, height = 6.5)

    } # eo for loop - s
} # eo for loop - v


### 26/01/26: Same above: compute warming trends for temperature variables in the grasslands
setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/grasslands")
files <- dir()[grepl(paste("mw_","min","_","SM_10", sep = ""),dir())]
# files
res <- lapply(files, function(f) { dat <- get(load(f)); return(dat) } ) # eo lapply
table <- dplyr::bind_rows(res)
rm(res); gc()
# str(table); summary(table)
table$EP <- as.factor(table$EP)

# Compute annual values
ts_annual <- table %>%
  mutate(Year = lubridate::year(Date)) %>%
  group_by(EP,Year) %>%
  summarise(
    MeanTemp = mean(final_value, na.rm = TRUE),
    .groups = "drop"
)
# summary(ts_annual)
# Compute trends
library("trend")
ts_trends_sen <- data.frame(
    ts_annual %>%
    group_by(EP) %>%
    summarise(
        trend_C_per_decade = sens.slope(MeanTemp,Year)$estimates * 10,
        p_value = sens.slope(MeanTemp, Year)$p.value,
        .groups = "drop"
    )
) # eo ddf
# ts_trends_sen
summary(ts_trends_sen)

## max Ta_10
#       EP      trend_C_per_decade    p_value         
# AEG01  :  1   Min.   :0.3195     Min.   :3.277e-10  
# AEG02  :  1   1st Qu.:0.3831     1st Qu.:3.105e-09  
# AEG03  :  1   Median :0.3963     Median :1.449e-08  
# AEG04  :  1   Mean   :0.3937     Mean   :1.871e-08  
# AEG05  :  1   3rd Qu.:0.4096     3rd Qu.:2.769e-08  
# AEG06  :  1   Max.   :0.4395     Max.   :7.864e-08 

## max Ta_200
#       EP      trend_C_per_decade    p_value         
# AEG01  :  1   Min.   :0.3054     Min.   :7.414e-10  
# AEG02  :  1   1st Qu.:0.3293     1st Qu.:6.623e-09  
# AEG03  :  1   Median :0.3391     Median :1.267e-08  
# AEG04  :  1   Mean   :0.3408     Mean   :6.411e-08  
# AEG05  :  1   3rd Qu.:0.3504     3rd Qu.:2.956e-08  
# AEG06  :  1   Max.   :0.3930     Max.   :3.837e-06  

## max Ts_05
#       EP      trend_C_per_decade    p_value         
# AEG01  :  1   Min.   :0.1786     Min.   :1.765e-11  
# AEG02  :  1   1st Qu.:0.2270     1st Qu.:2.587e-10  
# AEG03  :  1   Median :0.2428     Median :1.605e-09  
# AEG04  :  1   Mean   :0.2496     Mean   :6.989e-09  
# AEG05  :  1   3rd Qu.:0.2656     3rd Qu.:3.986e-09  
# AEG06  :  1   Max.   :0.3554     Max.   :2.134e-07 

## max Ts_10
#       EP      trend_C_per_decade    p_value         
# AEG01  :  1   Min.   :0.1493     Min.   :1.372e-11  
# AEG02  :  1   1st Qu.:0.2009     1st Qu.:3.968e-10  
# AEG03  :  1   Median :0.2156     Median :1.902e-09  
# AEG04  :  1   Mean   :0.2216     Mean   :1.549e-08  
# AEG05  :  1   3rd Qu.:0.2370     3rd Qu.:5.859e-09  
# AEG06  :  1   Max.   :0.3100     Max.   :4.634e-07 

## max Ts_20
#       EP      trend_C_per_decade    p_value         
# AEG01  :  1   Min.   :0.1205     Min.   :3.295e-11  
# AEG02  :  1   1st Qu.:0.1795     1st Qu.:3.528e-10  
# AEG03  :  1   Median :0.1955     Median :7.854e-10  
# AEG04  :  1   Mean   :0.1996     Mean   :4.540e-09  
# AEG05  :  1   3rd Qu.:0.2223     3rd Qu.:1.955e-09  
# AEG06  :  1   Max.   :0.2891     Max.   :1.668e-07  

## max SM_10
#       EP      trend_C_per_decade    p_value        
# AEG01  :  1   Min.   :-1.0259    Min.   :0.002113  
# AEG02  :  1   1st Qu.:-0.3678    1st Qu.:0.045408  
# AEG03  :  1   Median :-0.2585    Median :0.126559  
# AEG04  :  1   Mean   :-0.2890    Mean   :0.212497  
# AEG05  :  1   3rd Qu.:-0.1547    3rd Qu.:0.307712  
# AEG06  :  1   Max.   : 0.0374    Max.   :0.963515  
### -> mostly N.S.
### Those that are signif. are probably negatove trends.
ts_trends_sen[ts_trends_sen$p_value < 0.05,] # yes
nrow(ts_trends_sen[ts_trends_sen$p_value < 0.05,]) # 39 EPs showing signif. decreasing trends in max SM_10

## min SM_10 - more interesting
#       EP      trend_C_per_decade    p_value        
# AEG01  :  1   Min.   :-1.17712   Min.   :0.001649  
# AEG02  :  1   1st Qu.:-0.36539   1st Qu.:0.039549  
# AEG03  :  1   Median :-0.24734   Median :0.117737  
# AEG04  :  1   Mean   :-0.29807   Mean   :0.179623  
# AEG05  :  1   3rd Qu.:-0.16778   3rd Qu.:0.255663  
# AEG06  :  1   Max.   : 0.08644   Max.   :0.963515 
ts_trends_sen[ts_trends_sen$p_value < 0.05,] # yes
nrow(ts_trends_sen[ts_trends_sen$p_value < 0.05,]) # 46 EPs showing signif. decreasing trends in min SM_10


## total precip.
#       EP      trend_C_per_decade     p_value        
# AEG01  :  1   Min.   :-0.047758   Min.   :0.001022  
# AEG02  :  1   1st Qu.:-0.021924   1st Qu.:0.182402  
# AEG03  :  1   Median :-0.010514   Median :0.372405  
# AEG04  :  1   Mean   :-0.008369   Mean   :0.447678  
# AEG05  :  1   3rd Qu.: 0.005123   3rd Qu.:0.721245  
# AEG06  :  1   Max.   : 0.065312   Max.   :0.985402  
### -> mostly N.S.
ts_trends_sen[ts_trends_sen$p_value < 0.05,]
#      EP trend_C_per_decade     p_value
#10 AEG10         0.06531168 0.001022208
#54 HEG04        -0.04775781 0.026831705
#83 HEG33        -0.04569929 0.030845050
### -> 3 EPs show signif. trends


### 26/01/26: Let’s compare median max trends:
# Variable	    Forest	    Grassland
# Ta_10	        ~1.05	    ~0.40
# Ta_200	    ~0.66	    ~0.34
# Ts_05	        ~1.03	    ~0.24
# Ts_10	        ~1.26	    ~0.22
# Ts_20	        ~0.38	    ~0.20

### -> For grasslands: very much in line with current knowledge - great! :) 
### -> For forests: way too high rate of warming per decade! Should be closer or even lower than grasslands! 
### Literature implies that as Europe warms, forests both warm (in absolute terms) and strengthen their relative cooling in summer,
### particularly where soil moisture is still sufficient to support evapotranspiration.
### -> We should expect higher absolute warming rates in open grasslands than in forests at 2m and at 5–10cm soil depth,
### especially for maximum temperatures and extremes.

### References: 
# https://www.nature.com/articles/s41467-025-63556-2 


### 27/01/26: Potential technical explanations for this major issue! 

### Forest climate reconstruction is: Temp_f = Temp_g + Offset
### Where Offset = f(stange, tree type*DOY, region, DOY, anomaly of Temp_g to DOY)

### GAMM is trained on a short modern period (≈2009–2024) -> applied to historical climate conditions that lie outside the joint predictor space
### We trained a model where the temperature sensitivity of forest offsets is learned from modern anomaly distributions,
### then we apply it to past periods where the anomaly distribution might be fundamentally different. This may cause a systematic
### cold bias in the hindcast, which then inflates warming trends.

### In other words, the GAMM learns relationships like: “when grassland anomaly is +3 °C in summer, forest offset behaves like X”
### “when anomaly is −2 °C in spring, forest offset behaves like Y”. The GAMMs learn this under a warm-biased anomaly regime (modern days). 

### In the past, grassland anomalies relative to modern DOY climatology should be systematically more negative, compressed toward the lower tail.
### Therefore the model sees: anomaly values it interprets as “cooler-than-usual days”, even though, climatologically, they were normal at the time
### -> The GAMM responds by predicting too-strong cooling offsets -> GAMM responds by predicting too-strong cooling offsets -> artificially steep warming trends

### TEST 1
### -> CHECK THIS BY PLOTTING GRASSLAND ANOMALY TO DOY IN TIME FIRST, ALSO PLOT PREDICTED OFFSET IN THE PAST
### -> COMPARE THEM IN THE TRAINING PERIOD TO HINDCATS PERIOD   

### TEST 2
### -> PLOT THE SMOOTH TERM: 's(Tgrass_anomaly_to_DOY, by = DOY)'
### Does the offset become more negative for negative anomalies? Is this effect asymmetric?

### -----------------------------------------------

### 11/02/26: Re-compute and analyze trends (=/- per decade) for corrected forests predictions

setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/forests/") # dir()

# For testing
# v = "precipitation"

for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10","precipitation")) {
    
    if( v == "precipitation" ) {
            files <- dir()[grepl(paste("mw_total_precipitation", sep = ""),dir())]
            # Load files
            res <- lapply(files, function(f) { dat <- get(load(f)); return(dat) } ) # eo lapply
            # Rbind
            table <- dplyr::bind_rows(res)
            rm(res); gc()
            # Add Region to 'table'
            table$Region <- NA
            for(i in unique(table$EP)) {
                # if else loop
                message(i)
                if( grepl("HEW",i) ) {
                    table[table$EP == i,"Region"] <- "HND"
                } else if( grepl("AEW",i) ) {
                    table[table$EP == i,"Region"] <- "SWA"
                } else if( grepl("SEW",i) ) {
                    table[table$EP == i,"Region"] <- "SCH"
                } # eo if els eloop
            } # eo for loop 
    } else {
            files <- dir()[grepl(paste("corr_max_",v, sep = ""),dir())]
            # Load files
            res <- lapply(files, function(f) { dat <- readRDS(f); return(dat) } ) # eo lapply
            # Rbind
            table <- dplyr::bind_rows(res)
            rm(res); gc()
    } # eo if else loop - v

    # Compute annual values
    ts_annual <- table %>%
        mutate(Year = lubridate::year(Date)) %>%
        group_by(EP,Year) %>%
        summarise(
            MeanTemp = mean(final_value, na.rm = TRUE),
            .groups = "drop"
        )

    # Compute trends
    require("trend")
    ts_trends_sen <- data.frame(
        ts_annual %>%
        group_by(EP) %>%
        summarise(
            trend_C_per_decade = sens.slope(MeanTemp,Year)$estimates * 10,
            p_value = sens.slope(MeanTemp, Year)$p.value,
            .groups = "drop"
        )
    ) # eo ddf

    message("\n")
    message(paste("Trends for: ",v, sep = ""))
    print(summary(ts_trends_sen))
    message("\n")

} # eo for loop

### Outputs: 

#Trends for: Ta_10
#       EP      trend_C_per_decade    p_value         
#AEW01  :  1   Min.   :-0.01506   Min.   :0.000e+00  
#AEW02  :  1   1st Qu.: 0.27504   1st Qu.:0.000e+00  
#AEW03  :  1   Median : 0.35277   Median :2.000e-08  
#AEW04  :  1   Mean   : 0.34287   Mean   :2.404e-02  
#AEW05  :  1   3rd Qu.: 0.41640   3rd Qu.:1.143e-05  
#AEW06  :  1   Max.   : 0.56022   Max.   :1.000e+00  
#(Other):118                                         


#Trends for: Ta_200
#       EP      trend_C_per_decade    p_value         
# AEW01  :  1   Min.   :-0.4012    Min.   :0.000e+00  
# AEW02  :  1   1st Qu.: 0.2168    1st Qu.:1.000e-08  
# AEW03  :  1   Median : 0.2693    Median :1.500e-06  
# AEW04  :  1   Mean   : 0.2583    Mean   :4.812e-02  
# AEW05  :  1   3rd Qu.: 0.3341    3rd Qu.:2.263e-03  
# AEW06  :  1   Max.   : 0.4213    Max.   :7.693e-01  
# (Other):118                                         


#Trends for: Ts_05
#       EP      trend_C_per_decade      p_value         
# AEW01  :  1   Min.   :-0.0006489   Min.   :0.000e+00  
# AEW02  :  1   1st Qu.: 0.1926703   1st Qu.:0.000e+00  
# AEW03  :  1   Median : 0.2686842   Median :0.000e+00  
# AEW04  :  1   Mean   : 0.4221233   Mean   :8.186e-03  
# AEW05  :  1   3rd Qu.: 0.6075955   3rd Qu.:2.000e-08  
# AEW06  :  1   Max.   : 1.3771278   Max.   :1.000e+00  
# (Other):118                                           


#Trends for: Ts_10
#       EP      trend_C_per_decade    p_value         
# AEW01  :  1   Min.   :-0.3697    Min.   :0.000e+00  
# AEW02  :  1   1st Qu.: 0.1737    1st Qu.:0.000e+00  
# AEW03  :  1   Median : 0.2469    Median :0.000e+00  
# AEW04  :  1   Mean   : 0.3551    Mean   :5.394e-03  
# AEW05  :  1   3rd Qu.: 0.5228    3rd Qu.:5.100e-07  
# AEW06  :  1   Max.   : 1.0995    Max.   :5.831e-01  
# (Other):118                                         


#Trends for: Ts_20
#       EP      trend_C_per_decade    p_value         
# AEW01  :  1   Min.   :0.06845    Min.   :0.000e+00  
# AEW02  :  1   1st Qu.:0.13643    1st Qu.:0.000e+00  
# AEW03  :  1   Median :0.19975    Median :0.000e+00  
# AEW04  :  1   Mean   :0.28512    Mean   :2.603e-04  
# AEW05  :  1   3rd Qu.:0.35351    3rd Qu.:3.710e-08  
# AEW06  :  1   Max.   :0.89981    Max.   :1.738e-02  
# (Other):118                                         


#Trends for: SM_10
#       EP      trend_C_per_decade    p_value         
# AEW01  :  1   Min.   :-5.09401   Min.   :0.000e+00  
# AEW02  :  1   1st Qu.:-0.69611   1st Qu.:4.670e-06  
# AEW03  :  1   Median : 0.03766   Median :3.272e-03  
# AEW04  :  1   Mean   :-0.18058   Mean   :9.483e-02  
# AEW05  :  1   3rd Qu.: 0.53677   3rd Qu.:3.016e-02  
# AEW06  :  1   Max.   : 2.58827   Max.   :9.773e-01  
# (Other):118         
### -> Significant trends but oscillate between losses and gains in moisture -> EP VARIABILITY DOMINATES

#Trends for: precipitation
#      EP            trend_C_per_decade     p_value        
# Length:150         Min.   :-0.048552   Min.   :0.000602  
# Class :character   1st Qu.:-0.020836   1st Qu.:0.126559  
# Mode  :character   Median :-0.008634   Median :0.357915  
#                    Mean   :-0.006756   Mean   :0.402148  
#                    3rd Qu.: 0.006889   3rd Qu.:0.667208  
#                    Max.   : 0.065782   Max.   :0.992701  
### -> N.S. in general

### 26/01/26: Let’s compare median trends in max values:
# Variable	    Forests	  Grasslands
# max Ta_10	    ~0.352	    ~0.40
# max Ta_200	~0.269	    ~0.34
# max Ts_05	    ~0.269	    ~0.24
# max Ts_10	    ~0.247	    ~0.22
# max Ts_20	    ~0.199	    ~0.20

### -> Comparable now! =) Forests even warm up less fast than grasslands when looking at surface and air temperatures
### -> Skin temperatures wrm up faster than air temperature and soil temperatures
### -> In forests, air temperatures warm at a similar pace as soil temperatures

### -> Looks good! Let's run the SPEI/ECEs computations based on those final forests values!

### ------------------------------------------------------------------------------------------------------------

### 16/02/26: For the GA 2026 in Wernigerode: plot distribution of long)term trends for each variable

### 1°) Grasslands

setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/grasslands/"); dir()

# For testing
v = "SM_10"

for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","precipitation","SM_10")) {
    
    if( v == "precipitation" ) {
            files <- dir()[grepl(paste("mw_total_precipitation", sep = ""),dir())]
            # Load files
            res <- lapply(files, function(f) { dat <- get(load(f)); return(dat) } ) # eo lapply
            # Rbind
            table <- dplyr::bind_rows(res)
            rm(res); gc()
            # Add Region to 'table'
            table$Region <- NA
            for(i in unique(table$EP)) {
                # if else loop
                message(i)
                if( grepl("HEW",i) ) {
                    table[table$EP == i,"Region"] <- "HND"
                } else if( grepl("AEW",i) ) {
                    table[table$EP == i,"Region"] <- "SWA"
                } else if( grepl("SEW",i) ) {
                    table[table$EP == i,"Region"] <- "SCH"
                } # eo if els eloop
            } # eo for loop 
    } else {
            files <- dir()[grepl(paste("_mw_max_",v, sep = ""),dir())]
            # Load files
            res <- lapply(files, function(f) { dat <- get(load(f)); return(dat) } ) # eo lapply
            # Rbind
            table <- dplyr::bind_rows(res)
            rm(res); gc()
    } # eo if else loop - v

    # Compute annual values
    ts_annual <- table %>%
        mutate(Year = lubridate::year(Date)) %>%
        group_by(EP,Year) %>%
        summarise(
            Mean = mean(final_value, na.rm = TRUE),
            .groups = "drop"
        )

    # Compute trends
    require("trend")
    ts_trends_sen <- data.frame(
        ts_annual %>%
        group_by(EP) %>%
        summarise(
            trend_C_per_decade = sens.slope(Mean,Year)$estimates * 10,
            p_value = sens.slope(Mean,Year)$p.value,
            .groups = "drop"
        )
    ) # eo ddf

    message(paste("Plotting trends distribution for ",v, sep = ""))
    if( v == "precipitation" ) { 
        plot <- ggplot(ts_trends_sen, aes(x = trend_C_per_decade)) +
            geom_density(fill = "#fdae61") +
            geom_rug(alpha = 0.2) + geom_vline(xintercept = 0, linetype = "dashed") +
            geom_vline(xintercept = median(ts_trends_sen$trend_C_per_decade), linewidth = 1) +
            xlab("Trend (mm per decade)") + ylab("Density") + theme_classic()
    } else { 
        plot <- ggplot(ts_trends_sen, aes(x = trend_C_per_decade)) +
            geom_density(fill = "#fdae61") +
            geom_rug(alpha = 0.2) + geom_vline(xintercept = 0, linetype = "dashed") +
            geom_vline(xintercept = median(ts_trends_sen$trend_C_per_decade), linewidth = 1) +
            xlab("Trend (°C per decade)") + ylab("Density") + theme_classic()
    }    

    # Saving plot
    ggsave(plot = plot, filename = paste("plots_density_trends_decades_grasslands_max_",v,"_16.02.26.jpg", sep = ""), dpi = 300, height = 3, width = 4)
    message("\n")

} # eo for loop


### 2°) Forests
setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/forests/"); dir()
# For testing
# v = "Ta_10"

for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20")) {
    
    files <- dir()[grepl(paste("max_",v, sep = ""),dir())]
    # Load files
    res <- lapply(files, function(f) { dat <- readRDS(f); return(dat) } ) # eo lapply
    # Rbind
    table <- dplyr::bind_rows(res)
    rm(res); gc()

    # Compute annual values
    ts_annual <- table %>%
        mutate(Year = lubridate::year(Date)) %>%
        group_by(EP,Year) %>%
        summarise(
            Mean = mean(final_value, na.rm = TRUE),
            .groups = "drop"
        )

    # Compute trends
    require("trend")
    ts_trends_sen <- data.frame(
        ts_annual %>%
        group_by(EP) %>%
        summarise(
            trend_C_per_decade = sens.slope(Mean,Year)$estimates * 10,
            p_value = sens.slope(Mean,Year)$p.value,
            .groups = "drop"
        )
    ) # eo ddf

    message(paste("Plotting trends distribution for ",v, sep = ""))
    plot <- ggplot(ts_trends_sen, aes(x = trend_C_per_decade)) +
        geom_density(fill = "#66bd63", alpha = 0.5) +
        geom_rug(alpha = 0.2) + geom_vline(xintercept = 0, linetype = "dashed") +
        geom_vline(xintercept = median(ts_trends_sen$trend_C_per_decade), linewidth = 1) +
        xlab("Trend (°C per decade)") + ylab("Density") + theme_classic()
  
    # Saving plot
    ggsave(plot = plot, filename = paste("plots_density_trends_decades_forests_max_",v,"_16.02.26.jpg", sep = ""), dpi = 300, height = 3, width = 4)
    message("\n")

} # eo for loop


### ------------------------------------------------------------------------------------------------------------

### 12/02/26: Trouble shoot issue of extremely long (Duration > 500) ECE of max Ts_05 and max Ts_10 at EP 'AEW11'

# Go to directory and extract TS of interest
# setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/forests") # dir()
# Read the file containing your observed daily stats, for all 3 regions
# files <- dir()[grepl(paste("corr_","max","_","Ts_05", sep = ""),dir())]; files
# table <- readRDS("table_combined_obs+corr_max_Ts_05_forests_SWA.rds")
# table <- table[table$EP == "AEW11",]
# str(table); head(table); summary(table)

#setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables")
#table <- readRDS("table_predict_GAMM_max_Ts_05_AEW11_06_02_26.rds")
#dim(table); str(table)
#summary(table)
#head(table)
#head( table[table$Date %in% seq(as.Date("2015-10-01"), as.Date("2018-10-01"), by = "days"),c("Date","Anom_flag_extrap")] )

### -> AEW11 was 2 yo in 1950 -> too young to be reliable in the chosen fixed baseline period (1950-1980)
### -> Unreliable stats for defining ECEs probably -> Add an argument or a threshold in ECE definition that forbids 
### the inclusion of forest EPs younger than 10-15 yo.

### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------