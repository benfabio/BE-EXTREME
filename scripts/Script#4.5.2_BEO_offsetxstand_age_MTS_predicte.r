### ------------------------------------------------------------------------------------------------------------

### 08/01/26 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to re-load the prediction tables from Script#4.5.1 and the final GAMM objects (Script#4.4.3) to predict 
### daily offsets of max/min Ta_10, Ta_200, SM_10, Ts_05, Ts_10, Ts_20 back to 1950. 

### In addition, we try to develop a function to evaluate the uncertainty of these back tracking predictions (e.g., use 95% CI width). 
### To quantify predictive uncertainty, try to apply a simulation-based approach that propagates multiple sources of variability:
# - (i) estimation uncertainty of fixed-effect coefficients and smooth terms
# - (ii) variability in EP-level random intercepts, and 
# - (iii) residual observation-level error

### Model objects are stored in:
# /home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data/models 

### Prediction tables to use are stored in: 
# /home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables

### Last update: 10/02/26 (Re-run evaluate_forest_predict() on latest predictions)

### ------------------------------------------------------------------------------------------------------------

library("purrr")
library("tidyr")
library("dplyr")
library("data.table")
library("reshape2")
library("lubridate")
library("ggplot2")
library("ggpubr")
library("viridis")
library("mgcv")
library("parallel")
library("MASS") # for mvrnorm
library("Metrics") # for rmse function in prediction evaluation step

### ------------------------------------------------------------------------------------------------------------

### Principle:
## GAMMs predict daily Offset as a function of: StandAge, DOY (with smooths), Anom/Mean_Grassland_Value (smooth),
## TreeType, Region (fixed effects), EP (random intercept), and AR1 temporal correlation (correlation within EP)
## Because the GAMMs include EP as a random effect, the predictions depend on EP, meaning: 
## predictions for known EPs can include conditional predictions (with the random effect for that EP), and
## predictions for new/unseen EPs will either: Use marginal predictions (random effect = 0) or simulate
## a new random effect (as we did in the uncertainty function)

## Note: it’s better to combine all EP-specific prediction tables into a single newdata first because:
# - GAMM’s predict function is vectorized (=it predicts for all rows at once)
# - we do nt need to loop over 150 EPs individually
# - it ensures that factor levels for TreeType, Region, and EP are consistent
# - makes it easier to use the simulation-based uncertainty function later

## Note: in the 'newdata', some EPs have been used to train the GAMM, some not...therefore, we should use marginal predictions.
## Those ignores the EP-specific random effect (= sets it to zero), but they work for all EPs and ensures consistency.
## In addition, for back-in-time projections, even “known” EPs are being extrapolated beyond the training data.
## Marginal predictions are more robust for extrapolation.


### Write a function that will fill in your prediction tables with the GAMM-based predictions of Offset
# To test FUN while writing it
# var = "SM_10"
# stat = "max"
# para = TRUE
# cores = 25


fill_predict_table <- function(var, stat, para, cores) {
        
        #' This function takes four arguments and returns a formatted data.frame:
        #' @param var the climate variable to process (character) - one of the following: 
        #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20" and "SM_10"
        #' @param stat the daily statistic (character): 'max' or 'min'
        #' @param para switch, whether to sue parallel computing when loading the EP_specific prediction tables (boolean)
        #' @param cores number of cores to use for parallel computing (integer)
        #' @return A formatted data.frame that is the full prediction table including the GAMM predictions of dailt Offset
    
        ### First, load the GAMM to be used for predictions
        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data/models")
        message(paste("\nLoading the ",stat," ",var," GAMM", sep = ""))
        mod <- readRDS(paste("model_GAMM_offsets_",stat,"_",var,"_subset_15yrs_03.02.26.rds", sep = "")) 
        # summary(mod$gam)

        ## Second, load all EP-specific prediction tables and rbind them into one ddf
        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables")
        files2load <- dir()[grepl(paste("table_predict_GAMM_",stat,"_",var, sep = ""),dir())]
        
        ### 11/01/26: Ignore '_filled_' file if it exists
        if( length( files2load[grepl("_filled_",files2load)] ) > 0 ) {
            files2load <- files2load[!grepl(files2load[grepl("_filled_",files2load)],files2load)]
        } # eo if loop

        ## You may use parallel computing depending on 'para'
        if( para ) {
            message(paste("Loading the EP-specific prediction tables for ",stat," ",var, sep = ""))
            res <- mclapply(files2load, function(f) {
                    d <- readRDS(f); return(d)
                }, mc.cores = cores
            ) # eo mclapply
            # Rbind
            table <- dplyr::bind_rows(res)
            rm(res); gc(files2load)
        } else {
            message(paste("Loading the EP-specific prediction tables for ",stat," ",var, sep = ""))
            res <- lapply(files2load, function(f) {
                    d <- readRDS(f); return(d)
                }
            ) # eo lapply
            # Rbind
            table <- dplyr::bind_rows(res)
            rm(res); gc(files2load)
        } # eo if else loop - para

        # Make sure factors are factors in 'table' ; str(table)
        table$Region <- as.factor(table$Region)
        table$EP <- as.factor(table$EP)
        table$MTS_type <- as.factor(table$MTS_type)
        
        # Change 'MTS_type' to 'TreeType' as in the GAMM
        colnames(table)[16] <- "TreeType"
        table <- table %>% dplyr::select(-OffsetPred)

        ### 02/02/26: Anoms are encoded as s(Anom) in the mod$gam. But need to use 'Anom_adj' for prediction!!
        ### Drop 'Anom' from 'table' and rename 'Anom_adj' to 'Anom'
        table <- table %>% dplyr::select(-Anom)
        colnames(table)[11] <- "Anom"

        ### Marginal prediction for all EPs
        message(paste("Predicting the daily offsets of ",stat," ",var," based on the chosen GAMM", sep = ""))
        table$Offset_mod <- predict(mod$gam, newdata = table, type = "response", exclude = "s(EP)")
        # summary(table) # table[is.na(table$Mean_Grassland_Value),]
        # summary( table[table$EP == "SEW47",] ) # some NAs sometimes but not so much, < 2%

        ### 12/01/26: Add forest climate prediction simply based on offset prediction and grassland predictions as follows:
        ###             Offset_mod = Grassland_mod - Forest_mod
        ### Meaning:    Forest_mod = Grassland_mod - Offset_mod
        table$Mean_Forest_value <- (table$Mean_Grassland_Value) - (table$Offset_mod)
        # summary(table$Mean_Forest_value) ; head(table)
        # For now, values make sense for max Ta_10, min Ta_200
        # For now, values make less sense for max Ta_200 & max Ts_10 (because modelled offset is mainly negative instead of positive?)
        # summary(table[,c("Mean_Grassland_Value","Mean_Forest_value")])

        ## 12/01/26: Post-processing correction of negative values of SM_10 (not realistic but expected)
        ## Negatives are an artefact of: additive modelling of offsets / anomalies, bias correction steps, 
        ## extrapolation outside the calibration range.
        ## A few options are possible:
        ## 1) For soil moisture, 0 is physically meaningful -> Replace negative values by 0 with 'SM_corr <- pmax(SM_forest_model, 0)'
        ## 2) If you dislike hard 0, use last-observation-carried-forward (LOCF) within EP with na.locf() # ?na.locf
        ## Option 2) is defensible in my opinion but it genertes new NAs (not that many but still)
        ## 3) Replace with the first non-negative value in the same EP found in the same month to keep seasonal realism
        ## (e.g., winter SM stays low, summer SM higher)& avoids unrealistic jumps at the start of the series

        if( var == "SM_10" ) {
            
            # Warning message
            message("Correcting potential predicted SM_10 values < 0")
            require("zoo")
            
            table <- table %>%
                group_by(EP,var,Month) %>% # group by EP and Month
                arrange(Date, .by_group = TRUE) %>%
                mutate(
                    # Flag negative SM_10 values only
                    SM_neg = var == "SM_10" & !is.na(Mean_Forest_value) & Mean_Forest_value < 0,
                    # Temporary series: only negatives become NA
                    SM_tmp = if_else(SM_neg, NA_real_, Mean_Forest_value),
                    # LOCF within EP and Month
                    SM_tmp_filled = zoo::na.locf(SM_tmp, na.rm = FALSE),
                    # Identify first non-negative value in this EP and Month
                    first_valid = first(Mean_Forest_value[Mean_Forest_value >= 0 & !is.na(Mean_Forest_value)]),
                    # Replace negative values
                    Mean_Forest_value = case_when(
                        SM_neg & !is.na(SM_tmp_filled) ~ SM_tmp_filled, # LOCF worked
                        SM_neg &  is.na(SM_tmp_filled) ~ first_valid, # use first valid in same EP+Month
                        TRUE ~ Mean_Forest_value
                    )
                ) %>%
                dplyr::select(-SM_neg, -SM_tmp, -SM_tmp_filled, -first_valid) %>%
                ungroup()
            # summary(table)
            # data.frame( table[!is.na(table$Mean_Grassland_Value) & table$Mean_Forest_value < 10,] )[1:100,]
            # dim(table[table$Mean_Forest_value < 1,])

        } # eo if loop - correcting SM_10 values < 0

        ### Final check: delete rows with only NAs = keep rows where at least one column is not NA
        table_clean <- table %>% filter(if_any(everything(), ~ !is.na(.)))

        ### Save predictions in a different .rds file but in the same directory
        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables")
        saveRDS(table_clean, file = paste0("table_predict_GAMM_filled_",stat,"_",var,"_10_02_26.rds") )
        # Clean and stop
        rm(table,table_clean,mod) ; gc()

} # eo fun - fill_predict_table


### Apply fill_predict_table to all variables
for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")) {
        for(s in c("max","min")) {
            fill_predict_table(var = v, stat = s, para = TRUE, cores = 25)
        } # eo 2nd for loop
} # eo 1st for loop


### 10/02/26: Running fill_predict_table() for c("Ta_10","Ta_200","Ts_05","Ts_10") - check outputs
# dir()[grepl("_filled_",dir())]
# [1] "table_predict_GAMM_filled_max_SM_10_10_02_26.rds" 
# [2] "table_predict_GAMM_filled_max_Ta_10_10_02_26.rds" 
# [3] "table_predict_GAMM_filled_max_Ta_200_10_02_26.rds"
# [4] "table_predict_GAMM_filled_max_Ts_05_10_02_26.rds" 
# [5] "table_predict_GAMM_filled_max_Ts_10_10_02_26.rds" 
# [6] "table_predict_GAMM_filled_max_Ts_20_10_02_26.rds" 
# [7] "table_predict_GAMM_filled_min_SM_10_10_02_26.rds" 
# [8] "table_predict_GAMM_filled_min_Ta_10_10_02_26.rds" 
# [9] "table_predict_GAMM_filled_min_Ta_200_10_02_26.rds"
# [10] "table_predict_GAMM_filled_min_Ts_05_10_02_26.rds" 
# [11] "table_predict_GAMM_filled_min_Ts_10_10_02_26.rds" 
# [12] "table_predict_GAMM_filled_min_Ts_20_10_02_26.rds" 
t <- readRDS(dir()[grepl("_filled_",dir())][3])
dim(t)
str(t)
head(data.frame(t))
summary(t)
t[is.na(t$Offset_mod),][1:100,]
#unique(t[is.na(t$OffsetPred),"EP"]) # 8 EPs: SEW07 SEW10 SEW11 SEW12 SEW19 SEW20 SEW46 SEW47
#summary(t[is.na(t$Offset_mod),]) # 22 to 172 yo; 2011-2024, Date
#unique(t[is.na(t$OffsetPred),"DOY"]) # all days
#unique(t[is.na(t$OffsetPred),"Year"]) # 2011 2012 2014 2022 2023 2024 only
#unique(t[is.na(t$OffsetPred),"Month"]) # all 12 months 
### Remove & clean
#rm(t); gc()

### 02/02/26 -> R Script#6.4.6 to diagnose potential issues in predicted Offset and Anoms

### 10/02/26 -> R Script#6.4.6 to diagnose potential issues in predicted Offset and Anoms

### ------------------------------------------------------------------------------------------------------------

### 12/01/26: After finishing the primary reconstruction of 'Mean_Forest_value' based on 'Mean_Grassland_value' and 
### predicted mean offset, let's run some diagnostics against OBSERVED forets values to assess how good we've done and 
### whether we still need to perform a trend-preserving quantile mapping on those modelled forest values. 

### So far, the back-transformation (Forest_mod = Grassland_mod - Offset_mod) is internally consistent, 
### preserves the coherence between grassland and forest temperatures, and ensures that forest values inherit
### long-term trends and synoptic variability from the grasslands. 

### Now, applying quantile mapping on 'Mean_Forest_value' (Tforest_final = QM(Tforest_model, Tforest_obs)) is justified if: 
### - Small distributional mismatches are left by offset modelling, scale differences between grassland-derived predictors 
###   and forest microclimate
### - There are biases in variance and tails of the distrbutions (e.g. underrepresented cold extremes under canopy) or 
###   if residual temporal effects were not fully captured by the GAMM coefficeints (e.g., DOY, DOYxTreeType)
### -> May still be relevant especially for extremes (tails of the distrbution)

### But, before committing to this extra QM, we'd need to check for the 2009-2024 period: 
# - QQ-plots: Tforest_model vs Tforest_obs
# - Seasonal variance ratios
# - Basic spell metrics
### If we see: compressed tails, seasonally compressed variances, very different spell metrics 
### -> QM on Tforest_model is justified
### If not, you may skip it and (proudly) stat that “no statistical post-processing was required” ;)

### Let's make a new FUN that does this for us -> evaluate_forest_predict() that will: 
# To test FUN whilre we're writing it:
# var = "Ta_200"
# stat = "max"
# monthly = TRUE
# plots = TRUE


evaluate_forest_predict <- function(stat, var, monthly = TRUE, plots = TRUE) {
        
        #' This function takes four arguments and returns a formatted data.frame:
        #' @param var the climate variable to process (character) - one of the following: 
        #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20" and "SM_10"
        #' @param stat the daily statistic (character): 'max' or 'min'
        #' @param monthly switch: whether to compute diagnostic indices per month as well (BOOLEAN)
        #' @param plots switch: whether to draw and save plots as .jpg in directory (BOOLEAN)
        #' @return A formatted data.frame containing the evaluation metrics - also returns some plots 
        
        ### Useless message
        message(paste("Evaluating modelled daily forest values of",stat,var,"\n", sep = " "))

        ### 1. Load filled predictin table, load initial forest observations and join into one table (2009-2024 period only) 
        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables")
        table_pred <- readRDS( dir()[grepl(paste("filled_",stat,"_",var, sep = ""),dir())] )

        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data")
        table_obs <- readRDS( dir()[grepl(paste("metadata_",stat,"_",var, sep = ""),dir())] )
        # dim(table_obs); dim(table_pred); head(table_obs)

        ## Join "table_pred$Mean_Forest_value" to 'table_obs'. Need to reate Date in Tforest_obs though:
        table_obs <- table_obs %>% mutate(Date = as.Date(DOY - 1, origin = paste0(Year,"-01-01"))) # summary(table_obs$Date) # guet
        # Join Mean_Forest_value from Tforest_model by EP & Date
        table_obs <- table_obs %>% left_join(table_pred %>% dplyr::select(EP, Date, Mean_Forest_value), by = c("EP","Date"))
        # dim(table_obs); str(table_obs); summary(table_obs)
        rm(table_pred); gc()

        ### 2 Compute & summarise diagnostics
        # Quickly need to rename the observatinal variable column
        colnames(table_obs)[5] <- "Obs"

        # Compute overall (across all Dates) diagnostics per EP
        stats <- data.frame(
            table_obs %>%
            group_by(EP) %>%
            summarise(
                n = sum(!is.na(Obs) & !is.na(Mean_Forest_value)),
                Obs_mean = mean(Obs, na.rm = TRUE),
                Obs_sd   = sd(Obs, na.rm = TRUE),
                Obs_min  = min(Obs, na.rm = TRUE),
                Obs_max  = max(Obs, na.rm = TRUE),
                Model_mean = mean(Mean_Forest_value, na.rm = TRUE),
                Model_sd   = sd(Mean_Forest_value, na.rm = TRUE),
                Model_min  = min(Mean_Forest_value, na.rm = TRUE),
                Model_max  = max(Mean_Forest_value, na.rm = TRUE),
                RMSE = sqrt(mean((Obs - Mean_Forest_value)^2, na.rm = TRUE)),
                MeanBias = mean(Mean_Forest_value - Obs, na.rm = TRUE),
                Correlation = cor(Obs, Mean_Forest_value, use = "complete.obs")
            ) # eo summarise
        ) # eo ddf
        # summary(stats)

        ## if monthly == TRUE, also compute diagnostics per month - same but with group_by(EP,Month) 
        if( monthly ) {
            
            table_obs$Month <- lubridate::month(table_obs$Date)

            stats_mon <- data.frame(
                table_obs %>%
                group_by(EP,Month) %>%
                summarise(
                    n = sum(!is.na(Obs) & !is.na(Mean_Forest_value)),
                    Obs_mean = mean(Obs, na.rm = TRUE),
                    Obs_sd   = sd(Obs, na.rm = TRUE),
                    Obs_min  = min(Obs, na.rm = TRUE),
                    Obs_max  = max(Obs, na.rm = TRUE),
                    Model_mean = mean(Mean_Forest_value, na.rm = TRUE),
                    Model_sd   = sd(Mean_Forest_value, na.rm = TRUE),
                    Model_min  = min(Mean_Forest_value, na.rm = TRUE),
                    Model_max  = max(Mean_Forest_value, na.rm = TRUE),
                    RMSE = sqrt(mean((Obs - Mean_Forest_value)^2, na.rm = TRUE)),
                    MeanBias = mean(Mean_Forest_value - Obs, na.rm = TRUE),
                    Correlation = cor(Obs, Mean_Forest_value, use = "complete.obs")
                ) # eo summarise
            ) # eo ddf
            # summary(stats_mon)
             
        } # eo if loop - monthly
        
        ### 3. Draw diagnostocs plots, if plots == TRUE
        if( plots ) {
            
            setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables")

            ## Scatterplot
            p1 <- ggplot(table_obs, aes(x = Obs, y = Mean_Forest_value)) + geom_point(alpha = 0.1) +
                geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#F21A00FF") +
                labs(
                    title = paste("Scatterplot of Obs vs. Modelled",stat,var, sep = " "),
                    x = paste("Observed",stat,var, sep = " "),
                    y = paste("Modelled",stat,var, sep = " ")) +
                theme_minimal()
  
            ## gghistogram - need to melt
            m_table_obs <- reshape2::melt(table_obs[,c("Obs","Mean_Forest_value")], value.name = "Value")
            # Adjust some names for fanciness; unique(m_table_obs[,"variable"])
            m_table_obs$variable <- as.character(m_table_obs$variable)
            m_table_obs[m_table_obs$variable == "Mean_Forest_value" & !is.na(m_table_obs$variable),"variable"] <- "Modelled"
            m_table_obs[m_table_obs$variable == "Obs" & !is.na(m_table_obs$variable),"variable"] <- "Observed"
            p2 <- gghistogram(m_table_obs, x = "Value", add = "mean", rug = TRUE, color = "variable", fill = "variable", palette = c("#EBCC2AFF","#78B7C5FF"))
            rm(m_table_obs); gc()

            ## Print plots as .jpg
            ggsave(plot = p1, filename = paste("plot_scatter_eval_forest_predict_",stat,"_",var,"_overall_10.02.26.jpg", sep = ""), dpi = 300, width = 5.5, height = 4)
            ggsave(plot = p2, filename = paste("gghist_eval_forest_predict_",stat,"_",var,"_overall_10.02.26.jpg", sep = ""), dpi = 300, width = 5.5, height = 4)

            ## Facet per month if monthly == TRUE
            if( monthly ) {
            
                table_obs$Month <- lubridate::month(table_obs$Date)
            
                # scatterplot
                p1_mon <- ggplot(table_obs, aes(x = Obs, y = Mean_Forest_value)) + geom_point(alpha = 0.3) +
                    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#F21A00FF") +
                    labs(
                        title = paste("Scatterplot of Obs vs. Modelled",stat,var, sep = " "),
                        x = paste("Observed",stat,var, sep = " "),
                        y = paste("Modelled",stat,var, sep = " ")) +
                    theme_minimal() + facet_wrap(.~factor(Month), nrow = 4, ncol = 3)
  
                # gghistogram - need to melt
                m_table_obs <- reshape2::melt(table_obs[,c("Obs","Mean_Forest_value","Month")], value.name = "Value", id.vars = "Month")
                # Adjust some names for fanciness; unique(m_table_obs[,"variable"])
                m_table_obs$variable <- as.character(m_table_obs$variable)
                m_table_obs[m_table_obs$variable == "Mean_Forest_value" & !is.na(m_table_obs$variable),"variable"] <- "Modelled"
                m_table_obs[m_table_obs$variable == "Obs" & !is.na(m_table_obs$variable),"variable"] <- "Observed"
            
                p2_mon <- gghistogram(
                        m_table_obs, x = "Value", add = "mean", rug = TRUE, color = "variable", fill = "variable", palette = c("#EBCC2AFF","#78B7C5FF")
                ) + facet_wrap(.~factor(Month), nrow = 3, ncol = 4)
    
                rm(m_table_obs); gc()

                # Save
                ggsave(plot = p1_mon, filename = paste("plot_scatter_eval_forest_predict_",stat,"_",var,"_mon_10.02.26.jpg", sep = ""), dpi = 300, width = 6, height = 8)
                ggsave(plot = p2_mon, filename = paste("gghist_eval_forest_predict_",stat,"_",var,"_mon_10.02.26.jpg", sep = ""), dpi = 300, width = 10, height = 7.5)

                ### 14/01/25: Examine distrbution of RMSE and average bias per month (boxplot)
                p1 <- ggplot(data = stats_mon, aes(x = factor(Month), y = MeanBias)) + geom_violin(fill = "grey", colour = "black") + 
                    geom_boxplot(fill = "white", colour = "black", width = .2) + xlab("Month") + ylab("Mean bias (modelled - observed)") +
                    ggtitle(paste("Evaluating",stat,var, sep = " ")) + theme_bw()

                p2 <- ggplot(data = stats_mon, aes(x = factor(Month), y = RMSE)) + geom_violin(fill = "grey", colour = "black") +
                    geom_boxplot(fill = "white", colour = "black", width = .2) + xlab("Month") + ylab("RMSE") + theme_bw()

                panel <- ggarrange(p1,p2, align = "hv", nrow = 2, ncol = 1)
                ggsave(plot = panel, filename = paste("boxplots_eval",stat,var,"RMSE+MeanBias_10.02.26.jpg", sep = "_"), dpi = 300, width = 8, height = 6)

            } # eo if loop - monthly

        } # eo if loop - plots

        ### 4. Return summaries for programmatic use
        if( monthly ) { 
            list(overall = stats, monthly = stats_mon)
        } else {
            list(overall = stats)
        } # eo - monthly

} # eo FUN - evaluate_forest_predict


### Apply evaluate_forest_predict()
list_stats <- evaluate_forest_predict(stat = "max", var = "Ta_10", monthly = TRUE, plots = TRUE)
# str(list_stats); list_stats[[2]]
# mon_stats <- list_stats[[2]]
# summary(mon_stats)

### 10/02/26: Generate plots in a for loop for the updated forets predictions
for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")) {
    for(s in c("max","min")) {
        evaluate_forest_predict(stat = s, var = v, monthly = TRUE, plots = TRUE)
    } # eo 2nd for loop
} # eo 1st for loop


### 10/02/26: Observations based on the newest predictions

### In case we observe a systematic bias (regardless of the month) and a very high correlation coeff. (e.g., systematic bias of +2–3 °C with r ≈ 0.98)
### -> We got the physics right but the baseline slightly wrong -> simple bias correction per month
### Good because: additive, month-specific, does not depend strongly on quantile
### Allows to preserve: dynamics, trends, extremes and physical meaning
### Additional QM may overwrite that...

### Additional final check to assess whether QM is really needed or not: plot, per month, the bias vs. observed values
### If bias flat -> QM unjustified
### If bias NOT flat -> QM may be justified

# For testing
# stat = "min"
# var = "Ta_200"

plot_forest_bias <- function(stat, var, EP = FALSE) {
        
        #' This function takes four arguments and returns a formatted data.frame:
        #' @param var the climate variable to process (character) - one of the following: 
        #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20" and "SM_10"
        #' @param stat the daily statistic (character): 'max' or 'min'
        #' @param EP switch - Whether to draw the plot per EP instead of per month (BOOLEAN). Default == FALSE
        #' @return A formatted data.frame containing the evaluation metrics - also returns some plots 
        
        ### Useless message
        message(paste("Evaluating modelled daily forest values of",stat,var,"\n", sep = " "))

        ### 1. Load filled predictin table, load initial forest observations and join into one table (2009-2024 period only) 
        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables")
        table_pred <- readRDS( dir()[grepl(paste("filled_",stat,"_",var, sep = ""),dir())] )

        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data")
        table_obs <- readRDS( dir()[grepl(paste("metadata_",stat,"_",var, sep = ""),dir())] )
        # dim(table_obs); dim(table_pred); head(table_obs)

        ## Join "table_pred$Mean_Forest_value" to 'table_obs'. Need to reate Date in Tforest_obs though:
        table_obs <- table_obs %>% mutate(Date = as.Date(DOY - 1, origin = paste0(Year,"-01-01"))) # summary(table_obs$Date) # guet
        # Join Mean_Forest_value from Tforest_model by EP & Date
        table_obs <- table_obs %>% left_join(table_pred %>% dplyr::select(EP, Date, Mean_Forest_value), by = c("EP","Date"))
        # dim(table_obs); str(table_obs); summary(table_obs)
        rm(table_pred); gc()

        # Add month
        table_obs$Month <- lubridate::month(table_obs$Date)

        # Adjust colname
        colnames(table_obs)[5] <- "Obs"

        ### 2. Create 'bias' variable per month ina. new ddf and plot it per month
        diag <- table_obs %>% mutate(Bias = Mean_Forest_value - Obs)
        # head(diag); summary(diag)

        ## Switch: if EP == TRUE: Plot not per month but per EP (to assess whether MonthxEP corrections are further needed)
        if( EP ) {
            
            plot <- ggplot(diag, aes(x = Obs, y = Bias)) + geom_point(alpha = 0.1) +
                geom_smooth(method = "lm", se = TRUE, color = "#d53e4f") +
                facet_wrap(~ EP, ncol = 15, nrow = 10) +
                geom_hline(yintercept = 0, linetype = "dashed") +
                labs(
                    x = paste("Observed",stat,var, sep = " "),
                    y = "Bias (modelled - observed)",
                    title = paste("EP-level bias vs. observed for",stat,var, sep = " ")
                ) + theme_bw()

            # Save
            setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables")
            ggsave(plot = plot, filename = paste("plot_eval_EP_biasxobs",stat,var,"02.02.26.jpg", sep = "_"), dpi = 300, width = 15, height = 12.5)        
        
            # Clean before exiting loop
            rm(plot,diag,table_obs,table_pred); gc()

        } else {
            
            plot <- ggplot(diag, aes(x = Obs, y = Bias)) + geom_point(alpha = 0.1) +
                geom_smooth(method = "lm", se = TRUE, color = "#d53e4f") +
                facet_wrap(~ Month, ncol = 4) +
                geom_hline(yintercept = 0, linetype = "dashed") +
                labs(
                    x = paste("Observed",stat,var, sep = " "),
                    y = "Bias (modelled - observed)",
                    title = paste("Monthly bias vs. observed for",stat,var, sep = " ")
                ) + theme_bw()

            # Save
            setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables")
            ggsave(plot = plot, filename = paste("plot_eval_mon_biasxobs",stat,var,"02.02.26.jpg", sep = "_"), dpi = 300, width = 7, height = 5)        
        
            # Clean before exiting loop
            rm(plot,diag,table_obs,table_pred); gc()

        } # eo if else loop - EP
        
} # eo FUN - evaluate_forest_bias


### Apply evaluate_forest_bias() 
### 02/02/26: Re-run with new adjusted ANom data
for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")) {
    for(s in c("max","min")) {
        plot_forest_bias(stat = s, var = v, EP = FALSE)
    } # eo 2nd for loop
} # eo 1st for loop


### Observations:

## Ta_10: 
# Mean bias (Modelled max/min Ta_10) and RMSE seem to be strongly conserved across months for both min and max Ta_10
# Modelled max/min Ta_10 is consistently higher than observations. RMSE of max Ta_10 is around 5. RMSE of min Ta_10 is around 2.75. 
# Plots of Biases as a function of observed max/min Ta_10 (per month) show strong linear decrease in bias with higher observed values, across all 12 months. 
# Slopes seem to change across months though. Slopes are stronger for biases in min Ta_10 than in max Ta_10.
# Yet, plots of Biases as a function of observed max/min Ta_10 per forest (EP; 150 values) show broader variability in slopes but most of them also show 
# a linear decrease in biases with higher max/min Ta_10 values. Some arer EPs show linear increase in Biases with observed values. 
# Some EPs show weak linear trend (nearly flat lines) but high absolute biase relative to observations. 

## Ta_200: 
# Mean bias (Modelled max/min Ta_200) and RMSE are strongly conserved across months for both min and max Ta_200
# Modelled max/min Ta_200 is consistently higher than observations, but less than for Ta_10.
# RMSE of max Ta_200 is around 2.8-3. RMSE of min Ta_200 is around 1.75-2.
# Plots of Biases as a function of observed max/min Ta_200 (per month) show linear decrease in bias with higher observed values, across all 12 months. 
# Slopes seem to change across months. Slopes are stronger for biases in min Ta_10 than in max Ta_10.
# Some monthly slopes are relatively flat for max Ta_200.
# Yet, plots of Biases as a function of observed max/min Ta_10 per forest (EP; 150 values) show broader variability in slopes than per month
# (12 values only though), but most of them also show a linear decrease in biases with higher max/min Ta_10 values.
# Some arer EPs show linear increase in Biases with observed values. 
# Some EPs show weak linear trend (nearly flat lines) but high absolute biase relative to observations. 

## Ts_05: 
# Mean bias (Modelled max/min Ts_05) and RMSE are strongly conserved across months for both min and max Ts_05.
# Modelled max/min Ts_05 is consistently higher than observations.
# RMSE of max Ts_05 is around 7-7.5. RMSE of min Ts_05 is around 6. 
# Distribution plots of monthly modelled min/max Ts_05 show some kind of bimodality (2 peaks) that is not present in the observed Ts_05 values.
# Plots of Biases as a function of observed max/min Ts_05 (per month) show strong linear decrease in bias with higher observed values, across all 12 months. 
# Slopes do not seem to change across months. Slopes are similar for biases in both min Ta_10 than in max Ta_10.
# Yet, plots of Biases as a function of observed max/min Ts_05 per forest (EP; 150 values) show very broad variability in slopes than per month,
# but most of them also show a linear decrease in biases with higher max/min Ts_05 values. Some arer EPs show linear increase in Biases with observed values. 
# Some EPs show weak linear trend (nearly flat lines) but high absolute biase relative to observations (meaning: a flat cloud of point but clearly located above the bias == 0 line). 
# Variability in slopes per EP is larger than for air temperature variables (Ta_10 & Ta_200), which is surprising because the GAMMs of 
# soil temperature variables (Ts_05, Ts_10 & Ts_20) showed higher predictive power! EP-level bias plots are very simialr for min Ts_05 & max Ts_05.

## Ts_10: 
# Mean bias (Modelled max/min Ts_10) and RMSE are strongly conserved across months for both min and max Ts_10
# Modelled max/min Ts_10 is consistently higher than observations. 
# RMSE of max Ts_10 is around 7-7.5. RMSE of min Ts_10 is also around 7-7.5. 
# Distribution plots of monthly modelled min/max Ts_05 show some kind of bimodality (2 peaks) that is not present in the observed Ts_05 values.
# Rets of the observations are very simialr to what we observed for min/max Ts_05!

## Ts_20:
# Mean bias (Modelled max/min Ts_20) and RMSE are strongly conserved across months for both min and max Ts_20.
# Modelled max/min Ts_20 is consistently higher than observations, but the bias is not as bad as for Ts_05 and Ts_10.
# RMSE of max Ts_20 is around 2-3. RMSE of min Ts_20 is around 2-2.5. 
# Plots of Biases as a function of observed max/min Ts_20 (per month) show strong linear decrease in bias with higher observed values, across all 12 months. 
# Slopes seem to change across months a bit. Slopes are similar for both min Ta_10 & max Ta_10.
# Same observations as above regarding the EP-level plots of biases against observed values. 

## SM_10:
# Mean bias (Modelled max/min SM_10) and RMSE are strongly conserved across months for both min and max SM_10.
# Modelled max/min SM_10 is consistently higher than observations, but the bias is not as bad as for SM_10 and SM_10. 
# RMSE of max SM_10 is around 10-13. RMSE of min SM_10 is around 10-11. 
# Plots of Biases as a function of observed max/min SM_10 (per month) show linear decrease in bias with higher observed values, across all 12 months. 
# The slopes of the bias seems conserved across all 12 months for both min SM_10 and max SM_10.
# Same observations as above for the EP-levem plots of bias against obsefvations: Higher variability in slopes at the EP level, but most are decreasing linearly
# with higher observed values of SM_10. Some rare cases of linear increases with higher observed values. Some flat clouds of points above the bias == 0 line. 
# Very simialr for both min SM_10 and max SM_10.

### Synthesis:
# Systematic positive bias (model > obs)
# Bias decreases with increasing observed values
# Bias–Obs relationship is approximately linear
# Monthly structure exists but is secondary
# EP-level variability exists but is noisy and often very variables
# RMSE stable across months
# This is not random error. It is a structural scaling problem.

## Formally: 
# cold / dry / low-T states -> overestimated
# warm / wet / high-T states -> closer to observations


### 02/02/26: Same structural problems remain after changing anomalies to DOY though! Overestimation of forests' values by several degrees
### -> Likely caused by a too strong warming trend due to decreasing offsets in time (s(StandAge) probably the driver)
### -> Let's correct them with EP-specific slopes again and re-plot final corrected vales to assess warming bias etc.


### How to test whether to perform Month-specific corrections or EP-specific corrections? 
### -> Out Of Sample RMSE comparison
# Hold out some years of data (e.g. 8 years).
# Fit the correction model on all other years.
# Predict the held-out year.
# Compute RMSE between predictions and observations.
# This simulates what you actually care about:
# - applying a correction outside the calibration window
# - especially relevant because you reconstruct 1950–2008
## = Checks whether EP- or month-specific corrections actually improve prediction of unseen years,
## which is exactly the situation you face when reconstructing climate back to 1950.

## Why RMSE (not p-values)? p-values tell you whether something exists, RMSE tells you whether it helps
## So, let's compute RMSE for each model:
#       Obs ~ Model * EP
#       Obs ~ Model * Month
# Then compare:
# rmse_EP  < rmse_mon   -> EP structure helps more
# rmse_mon < rmse_EP    -> Month structure helps more
# -> The smaller RMSE wins.

### Small FUN to compute cross validation RMSE
#data = table_obs
#model_formula = Obs ~ Mean_Forest_value * EP
#
#cv_rmse <- function(model_formula, data) {
#    years <- unique(data$Year)
#    rmse <- sapply(years, function(y) {
#        train <- filter(data, Year != y)
#        test  <- filter(data, Year == y)
#        fit <- lm(model_formula, data = train)
#        pred <- predict(fit, newdata = test)
#        rmse(test$Obs, pred)
#        }
#    )
#  mean(rmse, na.rm = TRUE)
#}

# For testing: 
# stat = "max"
# var = "Ta_200"

test_modl_slopes <- function(stat,var) {

        ### Useless message
        message(paste("Performing RMSE comparison for",stat,var, sep = " "))

        ### 1. Load filled predictin table, load initial forest observations and join into one table (2009-2024 period only) 
        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables")
        table_pred <- readRDS( dir()[grepl(paste("filled_",stat,"_",var, sep = ""),dir())] )

        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data")
        table_obs <- readRDS( dir()[grepl(paste("metadata_",stat,"_",var, sep = ""),dir())] )
        # dim(table_obs); dim(table_pred); head(table_obs)

        ## Join "table_pred$Mean_Forest_value" to 'table_obs'. Need to reate Date in Tforest_obs though:
        table_obs <- table_obs %>% mutate(Date = as.Date(DOY - 1, origin = paste0(Year,"-01-01"))) # summary(table_obs$Date) # guet
        # Join Mean_Forest_value from Tforest_model by EP & Date
        table_obs <- table_obs %>% left_join(table_pred %>% dplyr::select(EP, Date, Mean_Forest_value), by = c("EP","Date"))
        # dim(table_obs); str(table_obs); summary(table_obs)
        rm(table_pred); gc()

        # Add month
        table_obs$Month <- lubridate::month(table_obs$Date)
        table_obs$Year <- lubridate::year(table_obs$Date)

        # Adjust colname
        colnames(table_obs)[5] <- "Obs"

        ### 2. Perform CV RMSE comparison based on test/train blocks 
        train <- table_obs %>% filter(Year <= 2018)
        test  <- table_obs %>% filter(Year > 2018)
        fit_EP  <- lm(Obs ~ Mean_Forest_value * EP, data = train)
        fit_mon <- lm(Obs ~ Mean_Forest_value * Month, data = train)
        pred_EP  <- predict(fit_EP,  newdata = test)
        pred_mon <- predict(fit_mon, newdata = test)
        # summary(pred_EP); summary(pred_mon)
        rmse_EP  <- sqrt(mean((test$Obs - pred_EP)^2, na.rm = TRUE))
        rmse_mon <- sqrt(mean((test$Obs - pred_mon)^2, na.rm = TRUE))

        message(paste("RMSE of EP model = ",rmse_EP, sep = ""))
        message(paste("RMSE of Month model = ",rmse_mon, sep = ""))
        message("\n")
        
} # eo FUN - evaluate_forest_bias


### Apply test_modl_slopes()
for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")) {
    for(s in c("max","min")) {
        test_modl_slopes(stat = s, var = v)
    } # eo 2nd for loop
} # eo 1st for loop

# Performing RMSE comparison for max Ta_10
# RMSE of EP model = 1.730
# RMSE of Month model = 3.69

# Performing RMSE comparison for min Ta_10
# RMSE of EP model = 1.32
# RMSE of Month model = 1.94

# Performing RMSE comparison for max Ta_200
# RMSE of EP model = 0.96
# RMSE of Month model = 2.02

# Performing RMSE comparison for min Ta_200
# RMSE of EP model = 1.05
# RMSE of Month model = 1.49

# Performing RMSE comparison for max Ts_05
# RMSE of EP model = 1.05
# RMSE of Month model = 3.38

# Performing RMSE comparison for min Ts_05
# RMSE of EP model = 0.85
# RMSE of Month model = 2.98

# Performing RMSE comparison for max Ts_10
# RMSE of EP model = 1.06
# RMSE of Month model = 3.44

# Performing RMSE comparison for min Ts_10
# RMSE of EP model = 0.89
# RMSE of Month model = 3.16

# Performing RMSE comparison for max Ts_20
# RMSE of EP model = 0.66
# RMSE of Month model = 1.52

# Performing RMSE comparison for min Ts_20
# RMSE of EP model = 0.62
# RMSE of Month model = 1.54

# Performing RMSE comparison for max SM_10
# RMSE of EP model = 5.94
# RMSE of Month model = 9.16

# Performing RMSE comparison for min SM_10
# RMSE of EP model = 5.53
# RMSE of Month model = 8.75


### 15/01/26: Interpretation
### Across all variables and both min/max, we have:
# - EP model RMSE ≪ Month model RMSE, often by a factor of 2–3
# - Always in the same direction
### -> Forest-specific scaling differences are real, persistent, and predictive.
### -> EP ≠ noise !!

### Monthly slopes look similar on average but EP-level plots showed: wide slope dispersion, some reversed trends
### flat-but-shifted clouds. Meanwhile, monthly models: force one slope per month, average over incompatible EP-specific physics, 
### and cannot correct site-level scaling errors

### How to perform EP-level corrections safely: compare models of EP-specific intercepts vs. EP-specific intercepts
### In other words, compare RMSE of: 
### 'Obs ~ Mean_Forest_value + EP'    vs.     'Obs ~ Mean_Forest_value * EP'    vs.     Obs ~ Mean_Forest_value * EP + Mean_Forest_value:Month

### Simply write a second verison of the test_modl_slopes() function above to compare these 3 models based on RMSE: 
var = "Ta_10"
stat = "max"

test_modl_slopes2 <- function(stat, var) {

        #' This function takes four arguments and returns a formatted data.frame:
        #' @param var the climate variable to process (character) - one of the following: 
        #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20" and "SM_10"
        #' @param stat the daily statistic (character): 'max' or 'min'
        #' @return Print the RMSE and AIC values of each of the 3 linear models
    
        message(paste("Performing RMSE comparison for", stat, var))

        ### 1. Load prediction table
        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables")
        table_pred <- readRDS(dir()[grepl(paste("filled_", stat, "_", var, sep = ""), dir())])

        ### 2. Load observation table
        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data")
        table_obs <- readRDS(dir()[grepl(paste("metadata_", stat, "_", var, sep = ""), dir())])

        ### 3. Create Date and join predictions
        table_obs <- table_obs %>%
            mutate(Date = as.Date(DOY - 1, origin = paste0(Year, "-01-01"))) %>%
            left_join(
                table_pred %>% dplyr::select(EP, Date, Mean_Forest_value),
                by = c("EP", "Date")
            )

        rm(table_pred); gc()

        ### 4. Add time variables
        table_obs <- table_obs %>%
            mutate(
                Month = lubridate::month(Date),
                Year  = lubridate::year(Date)
            )

        ### 5. Rename observed column
        colnames(table_obs)[5] <- "Obs"

        ### 6. Train / test split (block-based)
        train <- table_obs %>% filter(Year <= 2018)
        test  <- table_obs %>% filter(Year > 2018)

        ### 7. Fit models
        fit_T1 <- lm(Obs ~ Mean_Forest_value + EP, data = train)
        fit_T2 <- lm(Obs ~ Mean_Forest_value * EP, data = train)
        fit_T3 <- lm(Obs ~ Mean_Forest_value * EP + Mean_Forest_value:Month, data = train)

        ### 8. Predict
        pred_T1 <- predict(fit_T1, newdata = test)
        pred_T2 <- predict(fit_T2, newdata = test)
        pred_T3 <- predict(fit_T3, newdata = test)

        ### 9. Compute RMSE & AIC
        rmse <- function(obs, pred) { sqrt(mean((obs - pred)^2, na.rm = TRUE)) }

        rmse_T1 <- rmse(test$Obs, pred_T1)
        rmse_T2 <- rmse(test$Obs, pred_T2)
        rmse_T3 <- rmse(test$Obs, pred_T3)

        AIC_T1 <- AIC(fit_T1)
        AIC_T2 <- AIC(fit_T2)
        AIC_T3 <- AIC(fit_T3)

        ### 10. Report
        message(sprintf("RMSE of Tier 1 (EP intercepts)          = %.3f", rmse_T1))
        message(sprintf("RMSE of Tier 2 (EP slopes)              = %.3f", rmse_T2))
        message(sprintf("RMSE of Tier 3 (EP slopes + Month mod.) = %.3f", rmse_T3))
        message(sprintf("AIC Tier of 1 (EP intercepts)          = %.3f", AIC_T1))
        message(sprintf("AIC Tier of 2 (EP slopes)              = %.3f", AIC_T2))
        message(sprintf("AIC Tier of 3 (EP slopes + Month mod.) = %.3f", AIC_T3))

        message("\n")
        rm(pred_T3,pred_T2,pred_T1,fit_T1,fit_T2,fit_T3,test,train,table_obs,table_pred); gc()

} # eo FUN - test_modl_slopes2

### Apply in a double for loop as usual
for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")) {
    for(s in c("max","min")) {
        test_modl_slopes2(stat = s, var = v)
    } # eo 2nd for loop
} # eo 1st for loop

### Observations
# Performing RMSE comparison for max Ta_10
# RMSE of Tier 1 (EP intercepts)          = 1.768
# RMSE of Tier 2 (EP slopes)              = 1.730
# RMSE of Tier 3 (EP slopes + Month mod.) = 1.730
# AIC Tier of 1 (EP intercepts)          = 2011018.641
# AIC Tier of 2 (EP slopes)              = 1949034.890
# AIC Tier of 3 (EP slopes + Month mod.) = 1948919.049

# Performing RMSE comparison for min Ta_10
# RMSE of Tier 1 (EP intercepts)          = 1.345
# RMSE of Tier 2 (EP slopes)              = 1.325
# RMSE of Tier 3 (EP slopes + Month mod.) = 1.326
# AIC Tier of 1 (EP intercepts)          = 1827432.060
# AIC Tier of 2 (EP slopes)              = 1814893.156
# AIC Tier of 3 (EP slopes + Month mod.) = 1814003.265

# Performing RMSE comparison for max Ta_200
# RMSE of Tier 1 (EP intercepts)          = 0.973
# RMSE of Tier 2 (EP slopes)              = 0.960
# RMSE of Tier 3 (EP slopes + Month mod.) = 0.959
# AIC Tier of 1 (EP intercepts)          = 1357329.262
# AIC Tier of 2 (EP slopes)              = 1332152.320
# AIC Tier of 3 (EP slopes + Month mod.) = 1332095.409

# Performing RMSE comparison for min Ta_200
# RMSE of Tier 1 (EP intercepts)          = 1.052
# RMSE of Tier 2 (EP slopes)              = 1.049
# RMSE of Tier 3 (EP slopes + Month mod.) = 1.049
# AIC Tier of 1 (EP intercepts)          = 1522616.618
# AIC Tier of 2 (EP slopes)              = 1519941.313
# AIC Tier of 3 (EP slopes + Month mod.) = 1518827.808

# Performing RMSE comparison for max Ts_05
# RMSE of Tier 1 (EP intercepts)          = 1.111
# RMSE of Tier 2 (EP slopes)              = 1.055
# RMSE of Tier 3 (EP slopes + Month mod.) = 1.058
# AIC Tier of 1 (EP intercepts)          = 1624039.117
# AIC Tier of 2 (EP slopes)              = 1553767.663
# AIC Tier of 3 (EP slopes + Month mod.) = 1544564.758

# Performing RMSE comparison for min Ts_05
# RMSE of Tier 1 (EP intercepts)          = 0.903
# RMSE of Tier 2 (EP slopes)              = 0.854
# RMSE of Tier 3 (EP slopes + Month mod.) = 0.856
# AIC Tier of 1 (EP intercepts)          = 1389059.106
# AIC Tier of 2 (EP slopes)              = 1330856.586
# AIC Tier of 3 (EP slopes + Month mod.) = 1329163.865

# Performing RMSE comparison for max Ts_10
# RMSE of Tier 1 (EP intercepts)          = 1.113
# RMSE of Tier 2 (EP slopes)              = 1.063
# RMSE of Tier 3 (EP slopes + Month mod.) = 1.069
# AIC Tier of 1 (EP intercepts)          = 1512790.911
# AIC Tier of 2 (EP slopes)              = 1443927.342
# AIC Tier of 3 (EP slopes + Month mod.) = 1440646.104

# Performing RMSE comparison for min Ts_10
# RMSE of Tier 1 (EP intercepts)          = 0.945
# RMSE of Tier 2 (EP slopes)              = 0.889
# RMSE of Tier 3 (EP slopes + Month mod.) = 0.893
# AIC Tier of 1 (EP intercepts)          = 1379196.628
# AIC Tier of 2 (EP slopes)              = 1312513.434
# AIC Tier of 3 (EP slopes + Month mod.) = 1311033.145

# Performing RMSE comparison for max Ts_20
# RMSE of Tier 1 (EP intercepts)          = 0.735
# RMSE of Tier 2 (EP slopes)              = 0.658
# RMSE of Tier 3 (EP slopes + Month mod.) = 0.660
# AIC Tier of 1 (EP intercepts)          = 1119286.266
# AIC Tier of 2 (EP slopes)              = 960362.494
# AIC Tier of 3 (EP slopes + Month mod.) = 951346.269

# Performing RMSE comparison for min Ts_20
# RMSE of Tier 1 (EP intercepts)          = 0.695
# RMSE of Tier 2 (EP slopes)              = 0.620
# RMSE of Tier 3 (EP slopes + Month mod.) = 0.622
# AIC Tier of 1 (EP intercepts)          = 1057036.562
# AIC Tier of 2 (EP slopes)              = 898911.063
# AIC Tier of 3 (EP slopes + Month mod.) = 888082.682

# Performing RMSE comparison for max SM_10
# RMSE of Tier 1 (EP intercepts)          = 6.112
# RMSE of Tier 2 (EP slopes)              = 5.942
# RMSE of Tier 3 (EP slopes + Month mod.) = 5.951
# AIC Tier of 1 (EP intercepts)          = 3148239.907
# AIC Tier of 2 (EP slopes)              = 3085746.755
# AIC Tier of 3 (EP slopes + Month mod.) = 3085271.658

# Performing RMSE comparison for min SM_10
# RMSE of Tier 1 (EP intercepts)          = 5.712
# RMSE of Tier 2 (EP slopes)              = 5.530
# RMSE of Tier 3 (EP slopes + Month mod.) = 5.540
# AIC Tier of 1 (EP intercepts)          = 3083034.158
# AIC Tier of 2 (EP slopes)              = 3011277.655
# AIC Tier of 3 (EP slopes + Month mod.) = 3010679.464

### Interpretation: 
## Across all 12 variables, Tier-2 models always improve RMSE relative to Tier-1, and Tier-2 models always show much lower AIC. 
## Meanwhile, Tier-3 models never really improve RMSE nor AIC.
## -> Tier-2 (EP-specific slope corrections) are the optimal correction strategy.
## -> bias structure is site-specific, mostly linear and stable across months

### Write a function that will help you evaluate the uncertainty of your temporal predictions of offset

## The GAMMs included smooth functions of StandAge, DOY by tree type, and anomaly metrics, as well as fixed effects
## for tree type and region, and random intercepts for ecological plots (EPs).

## To quantify predictive uncertainty, we test a simulation-based approach that propagates multiple sources of variability:
## (i) estimation uncertainty of fixed-effect coefficients and smooth terms,
## (ii) variability in EP-level random intercepts, and
### (iii) residual observation-level error.
## For each observation in the new dataset, we first constructed the linear predictor matrix relative to the fitted GAMM and
## generated 1,000 realizations of the fixed effects by sampling from a multivariate normal distribution defined by the 
## estimated coefficients and their covariance. Random intercepts were drawn from a normal distribution with variance equal to
## the estimated EP-level variance, and residual noise was added from a normal distribution with variance equal to the residual
## variance of the model. This yielded an ensemble of simulated predictions, from which we computed the predicted mean,
## 95% prediction intervals, absolute uncertainty (interval width), standard deviation of simulations, and relative uncertainty (CV).

## This approach provides observation-level estimates of predictive uncertainty that account for both model and residual variability.
## It is suitable for predictions at unseen EPs and for temporal extrapolation, including back-in-time projections,
## where uncertainty naturally increases with extrapolation beyond the training data.

## Let's automatically computes a prediction uncertainty metric for each observation. We’ll include:
# Pred = predicted mean
# Lower / Upper = 95% prediction interval
# PredUncertainty = width of the prediction interval (Upper − Lower)
# PredSD = standard deviation of simulated predictions
# PredCV = coefficient of variation (optional, relative uncertainty)

# For each row/predcited value in newdata we can have:
# newdata$Lower   # 2.5% quantile
# newdata$Upper   # 97.5% quantile
# newdata$Pred    # mean prediction
# Therefore, a simple uncertainty metric is the prediction interval width:
#           newdata$PredUncertainty <- newdata$Upper - newdata$Lower
# Wide interval = high uncertainty
# Narrow interval = low uncertainty

### Optional: normalize by prediction magnitude?
### Sometimes, when predicting back in time, the mean Offset may change a lot. To compare uncertainty relative
### to signal, you can compute coefficient of variation (CV):
#           newdata$PredCV <- newdata$PredUncertainty / abs(newdata$Pred)

### Write function for evaluating predictions' uncertainty with chatGPT5

# First, load the GAMM to be used for predictions
stat = "max"
var = "Ta_10"
setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data/models")
mod <- readRDS(paste("model_GAMM_offsets_",stat,"_",var,"_subset_14yrs_25.11.25.rds", sep = "")) 
# Second, load all EP-specific prediction tables and rbind them into one ddf
setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables")
files2load <- dir()[grepl(paste(stat,var, sep = "_"),dir())]
res <- mclapply(files2load, function(f) { d <- readRDS(f); return(d) }, mc.cores = 25 ) # eo mclapply
table <- dplyr::bind_rows(res)
rm(res); gc(files2load)
# Make sure factors are factors in 'table' ; str(table)
table$Region <- as.factor(table$Region)
table$EP <- as.factor(table$EP)
table$MTS_type <- as.factor(table$MTS_type)
colnames(table)[12] <- "TreeType"

# For testing as you write it
gamm_model = mod
newdata = table
n_sims = 999
include_random = TRUE
include_resid = TRUE
compute_cv = TRUE

evaluate_gamm_uncertainty <- function(gamm_model, newdata, n_sims = 999, include_random = FALSE, include_resid = TRUE, compute_cv = TRUE) {

        #' Predict offsets from a GAMM with uncertainty
        #'
        #' This function generates predicted offsets and associated uncertainty metrics
        #' from a generalized additive mixed model (GAMM) fitted with `mgcv::gamm`. 
        #' It propagates uncertainty from fixed-effect and smooth estimates, optional
        #' EP-level random effects, and residual observation-level variance.

        #' @param gamm_model A fitted GAMM object returned by `mgcv::gamm`.
        #' @param newdata A data.frame containing the new observations for prediction.
        #'   Must include all variables used in the model (predictors and factors). 
        #'   Factor levels should match the original training data.
        #' @param n_sims Integer. Number of simulations to perform for uncertainty estimation.
        #'   Default is 999
        #' @param include_random Logical. If TRUE, incorporates EP-level random effect 
        #'   uncertainty into the simulations. Default is TRUE.
        #' @param include_resid Logical. If TRUE, incorporates residual (observation-level) 
        #'   uncertainty into the simulations. Default is FALSE
        #' @param compute_cv Logical. If TRUE, computes the coefficient of variation 
        #'   (PredSD / Pred) as a relative uncertainty metric. Default is TRUE.
        #'
        #' @return A data.frame identical to `newdata` with additional columns:
        #'   \itemize{
        #'     \item \code{Pred}: predicted mean offset
        #'     \item \code{Lower}: 2.5th percentile of simulated predictions (95% CI lower)
        #'     \item \code{Upper}: 97.5th percentile of simulated predictions (95% CI upper)
        #'     \item \code{PredUncertainty}: width of the 95% prediction interval (Upper - Lower)
        #'     \item \code{PredSD}: standard deviation of simulated predictions
        #'     \item \code{PredCV}: relative uncertainty (PredSD / Pred), if \code{compute_cv = TRUE}
        #'   }
        #'
        #' @examples
        #' # Generate predictions with uncertainty for newdata
        #' predictions <- evaluate_gamm_uncertainty(
        #'     gamm_model = gamm_model,
        #'     newdata = newdata_all,
        #'     n_sims = 1000,
        #'     include_random = FALSE,  # marginal predictions
        #'     include_resid = TRUE,
        #'     compute_cv = TRUE
        #' )
        #'
        #' @export

        # --- 0. Check for unseen EPs ---------------------------------------------
        if ( "EP" %in% names(newdata) ) {
            ep_train <- levels(gamm_model$lme$data$EP)
            ep_new <- unique(as.character(newdata$EP))
            ep_unseen <- setdiff(ep_new, ep_train)
            if ( length(ep_unseen) > 0 ) {
                warning(
                    paste0(
                        "Predictions include ", length(ep_unseen),
                        " EP(s) not present in the training data.\n",
                        "For these EPs, predictions should be interpreted as population-level ",
                        "effects (random intercept = 0).\n",
                        "Consider using include_random = FALSE for consistency."
                    ),
                    call. = FALSE
                )
            }
        } # eo if loop

        # --- 1. Fixed effects / smooth uncertainty ---
        X <- predict(gamm_model$gam, newdata = newdata, type = "lpmatrix") # linear predictor; dim(X); str(X)
        beta <- coef(gamm_model$gam)
        Vb <- vcov(gamm_model$gam) # covariance of fixed effects
  
        set.seed(123) # for reproducible simulations
        beta_sim <- MASS::mvrnorm(n_sims, mu = beta, Sigma = Vb)
        pred_sims <- X %*% t(beta_sim)  # nrow(newdata) x n_sims matrix
        # dim(pred_sims); head(pred_sims)
  
        # --- 2. Random effects ---
        if (include_random) {
            sigma_EP <- as.numeric(VarCorr(gamm_model$lme)$EP[1])
            re_sim <- matrix(rnorm(n_sims * nrow(newdata), 0, sigma_EP), nrow = nrow(newdata), ncol = n_sims)
            pred_sims <- pred_sims + re_sim
        } # eo if loop
  
        # --- 3. Residual / observation-level uncertainty ---
        if (include_resid) {
            sigma_resid <- summary(gamm_model$lme)$sigma
            resid_sim <- matrix(rnorm(n_sims * nrow(newdata), 0, sigma_resid), nrow = nrow(newdata), ncol = n_sims)
            pred_sims <- pred_sims + resid_sim
        } # eo if loop
  
        # --- 4. Compute summary metrics ---
        newdata$Pred <- rowMeans(pred_sims)
        newdata$Lower <- apply(pred_sims, 1, quantile, 0.025)
        newdata$Upper <- apply(pred_sims, 1, quantile, 0.975)
  
        # Uncertainty metrics
        newdata$PredUncertainty <- newdata$Upper - newdata$Lower
        newdata$PredSD <- apply(pred_sims, 1, sd)
  
### Helps explain why Biases decrease with increasing observed value: 
### GAMMs seem to slightly over-amplify warm / SM_10 conditions and that amplification differs by site
### But it does not depend strongly on month once DOY effects are already in the GAMM
### The fact that EP-level slopes are > Monthly slopes supports the view that microclimate attenuation/amplification 
### is structural, liejly because it is driven by canopy, soil, exposure, and because it is stable across seasons

### Concrete next steps;
### - Perform Tier-2 type corrections on all variables ('Obs ~ Mean_Forest_value * EP')
### - Apply corrected values to full TS (1950-2024) - 'Tforest_corrected = a_EP + b_EP * Tforest_model'
### with a_EP = EP-specific intercept, and b_EP = EP-specific slope
### DO NOT re-apply QM

### -> R Script#4.5.3_correcte.R 

### BONUS: How to write this in the Methods of the paper: 
### "We evaluated three bias-correction strategies of increasing complexity: (i) EP-specific intercept corrections,
### (ii) EP-specific linear corrections, and (iii) EP-specific linear corrections with month-specific slope modulation.
### Model skill was evaluated using block-based out-of-sample RMSE, and model parsimony using AIC.
### Across all variables, EP-specific linear corrections substantially reduced RMSE relative to intercept-only models,
### while month-specific slope modulation provided negligible additional improvement in predictive performance despite
### slightly lower AIC values. We therefore adopted EP-specific linear corrections as the optimal bias-correction strategy,
### balancing predictive skill, robustness, and parsimony."


### 02/02/26: Same observations with corrected anomalies to DOY...

### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------