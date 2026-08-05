### ------------------------------------------------------------------------------------------------------------

### 28/03/25 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to evaluate the daily climate data corrected through various quantile mapping (QM) strategies
### against the locally-measured climate data.   
### Will help determine which QM strategy (global vs. monthly vs. mowing window vs. anomalies) provides
### the most accurate correction of E-OBS/ERA5-Land biases --> Which TS to use for ECEs detection and quantification?

### Recyling script 6.3.2 to write a master FUN that will:
###  - Load the corrected daily TS data from R script 6.4.2 (quantile_mapper outputs)
###  - Load the corresponding local measurements at EP level
###  - Compare the QM outputs against the observed daily data (see Script#6.3.1)
###  - Store the evaluation metrics for each EP and variable in a list (or some other type of object)
###  - Make plots to determine which QM strategy provided the best bias corrections

### Last update: 16/04/26 (Summarizing the mean ± sd values of all 5 evaluation metrics for Table 4 of the ESSD data paper)

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
# install.packages("metrica") # for RIA
library("metrica")

# Directories to save functions' outputs: 
plot.dir <- "/home/fbenedetti/plots/quantile_mapping_evaluation"
eval.metric.dir <- "/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/evaluation_metrics"

### ------------------------------------------------------------------------------------------------------------

### Write master FUN to return table containing the evaluation metrics. What argument should the FUN comprise?
### - variable
### - stat (min/max/mean)
### - region (SCH/HND/SWA)
### - QM method (global, monthly, mw, anoms)
### - anomalies (TRUE/FALSE) - should evaluation metrics be computed on anomalies to the monthly mean

### Master FUN - evaluate_quantile_maps()
# To test evaluate_quantile_maps while you're writing it:
var <- "precipitation"
stat <- "total"
region <- "SCH"
method <- "anomalies"
anoms <- TRUE

evaluate_quantile_maps <- function(var, stat, region, method, anoms) {

        #' This function takes four arguments and returns a formatted data.frame:
        #' @param var the climate variable to process (character) - one of the following: 
        #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20", "Ts_50", "precipitation" or "SM_10"
        #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character): 'mean', 'max' or 'min'
        #' @param region Region name (character): 'SCH' or 'HND' or 'SWA'
        #' @param method Which quantile mapping approach to run (character): 'global' or 'monthly' or 'mw' or 'anomalies'
        #' @param anoms Whether evaluation metrics should be computed on anomalies to the monthly mean (BOOLEAN)
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
        file <- dir()[grepl(paste(stat,var,region, sep = "_"),dir())]
        # Load 'file'
        corr_data <- get(load(file)) # dim(corr_data) ; str(corr_data) ; colnames(corr_data)
        
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
        # Clean ram
        rm(obs_daily_stat,corr_data)
        gc()

        # Adjust colnames
        colnames(merged_df)[c(3,4)] <- c("obs","corr_model") 
        
        # Remove NAs for evaluation (usually data < 2008 of course)
        merged_df <- na.omit(merged_df)

        # Vector of EP IDs - use to subset 'merged_df'
        plots <- unique(merged_df$EP)

        ### Calculate evaluation metrics: mean bias, corr coef, RMSE, mean absolute bias (MAE) for each EP separately!
        ### 04/06/25: ADDING REFINEX IDEX OF AGREEMENT (RIA) TO THE LIST OF EVALUATION METRICS

        if( anoms ) {

          # Compute monthly means for each data source
          merged_df <- merged_df %>%
              mutate(month = format(Date, "%Y-%m")) %>%
              group_by(EP,month) %>%
              mutate(
                  mon_mean_obs = mean(obs, na.rm = TRUE),
                  mon_mean_mod = mean(corr_model, na.rm = TRUE)
              ) %>% 
              ungroup()
          # summary(merged_df)

          # Compute monthly anomalies
          merged_df <- merged_df %>%
              mutate(
                  mon_anom_obs = obs - mon_mean_obs,
                  mon_anom_mod = corr_model - mon_mean_mod
              )
          # summary(merged_df)
            
          # Subset merged_df per EP and return evaluation metrics
          eval_metrics <- lapply(plots, function(p) {
                  
                  # p <- plots[13]
                  message(paste("Computing evaluation metrics for EP: ",p," based on monthly anomalies", sep = ""))
                  ## Subset
                  sub_merged_df <- merged_df[merged_df$EP == p,]
                  ## Only keep Dates > '2009-05-01' because of obs data
                  sub_merged_df <- sub_merged_df[sub_merged_df$Date >= '2009-05-01',]
                  ## Mean bias error (MBE)
                  mbe <- mean(sub_merged_df$mon_anom_mod - sub_merged_df$mon_anom_obs)
                  ## Mean absolute error (MAE)
                  mae <- mean(abs(sub_merged_df$mon_anom_mod - sub_merged_df$mon_anom_obs))
                  ## Root mean square error (RMSE)
                  rmse <- sqrt(mean((sub_merged_df$mon_anom_mod - sub_merged_df$mon_anom_obs)^2, na.rm = TRUE))
                  ## Corr coeff
                  corr_coeff <- cor(sub_merged_df$mon_anom_obs, sub_merged_df$mon_anom_mod, use = "complete.obs")
                  ## RIA
                  ria <- d1r(obs = sub_merged_df$mon_anom_obs, pred = sub_merged_df$mon_anom_mod, tidy = TRUE)
                      
                  # Return evaluation metrics in a data.frame
                  eval_metrics <- data.frame(EP = p,
                        mbe = mbe,
                        mae = mae, 
                        rmse = rmse, 
                        corr = corr_coeff,
                        ria = ria
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
          table_eval_metrics$variable <- var
          table_eval_metrics$stat <- stat
          table_eval_metrics$qm_method <- method
            
          # Save in dir as .Rdata
          message(paste("Saving evaluation metrics table for QM outputs (anomalies-based) based on ",method," ",
                  stat," ",var," in the ",region, sep = ""))
          setwd(eval.metric.dir)
          save(x = table_eval_metrics, file = paste("table_evaluation_metrics_anoms_",method,"_",stat,"_",var,"_grasslands_",region,".Rdata", sep = ""))
          rm(merged_df,table_eval_metrics); gc()

        } else {

          # Subset merged_df per EP and return evaluation metrics
          eval_metrics <- lapply(plots, function(p) {

                  # p <- sample(plots,1,1); p # for testing
                  message(paste("Computing evaluation metrics for ",p, sep = ""))
                  
                  ## subset
                  sub_merged_df <- merged_df[merged_df$EP == p,]
                  ## Only keep Dates > '2009-05-01' because of obs data
                  sub_merged_df <- sub_merged_df[sub_merged_df$Date >= '2009-05-01',]
                  ## Mean bias error (MBE)
                  mbe <- mean(sub_merged_df$corr_model - sub_merged_df$obs)
                  ## Mean absolute error (MAE)
                  mae <- mean(abs(sub_merged_df$corr_model - sub_merged_df$obs))
                  ## Root mean square error (RMSE)
                  rmse <- sqrt(mean((sub_merged_df$corr_model - sub_merged_df$obs)^2, na.rm = TRUE))
                  ## Corr coeff
                  corr_coeff <- cor(sub_merged_df$obs, sub_merged_df$corr_model, use = "complete.obs")
                  ## RIA
                  ria <- d1r(obs = sub_merged_df$obs, pred = sub_merged_df$corr_model, tidy = TRUE)
                      
                  # Return evaluation metrics in a data.frame
                  eval_metrics <- data.frame(EP = p,
                        mbe = mbe,
                        mae = mae,
                        rmse = rmse,
                        corr = corr_coeff,
                        ria = ria
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
          table_eval_metrics$variable <- var
          table_eval_metrics$stat <- stat
          table_eval_metrics$qm_method <- method
            
          message(paste("\n","Saving evaluation metrics table for QM outputs based on ",method," ",
                  stat," ",var," in the ",region, sep = ""))
          setwd(eval.metric.dir)
          save(x = table_eval_metrics, file = paste("table_evaluation_metrics_",method,"_",stat,"_",var,"_grasslands_",region,".Rdata", sep = ""))
          rm(merged_df,table_eval_metrics); gc()

      } # eo if else loop - anomalies

} # eo master FUN - evaluate_quantile_maps

### Apply evaluate_quantile_maps() for all data (variables, stats, QM method, regions...)

for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","precipitation","SM_10")) {
    for(s in c("max","min")) {
        for(r in c("HND","SCH","SWA")) {
            for(m in c('global','monthly','mw','anomalies')) {

                evaluate_quantile_maps(var = v, stat = s, region = r, method = m, anoms = FALSE)

                evaluate_quantile_maps(var = v, stat = s, region = r, method = m, anoms = TRUE)

            } # eo for loop - m
        } # eo for loop - r
    } # eo for loop - s
} # eo for loop - v


### ------------------------------------------------------------------------------------------------------------

### 26/05/2025: RE-load all evaluation metrics table generated above and make plots (geom_boxplot probably)
### to compare the outputs from the various QM strategies and help pick the best one (i.e., the one that leads
### to the lowest RMSE, MAB and highest corr coeff.) - after correcting

setwd(eval.metric.dir)
### NOTE:
###  '_anomalies_' -> metrics based on the anomalies QM stategy
###  '_anoms_' -> metrics based on anomalies to the mean (whatever the QM strategy)

### A) Evaluation metrics based on raw data

files <- dir()[!grepl("_anoms_",dir())]
files <- files[!grepl("_forests_",files)]
# files # gut
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
# Rbind
table <- bind_rows(res); rm(res); gc()
# dim(table) ; head(table) ; summary(table)

### The goal is twofold here:
### - 1) get a sense of the quality of the QM corrections overall per variables
### - 2) determine which QM strategy is the best overall
### Note: it could be that one strategy is not the best across ALL variables 

## Need to make two types of plots:
## (i) Distribution of RMSE/corr coeff/MAE per variable (+ facet per regions)
## (ii) Distribution of RMSE/corr coeff/MAE per variable and QM strategy (+ facet per regions)

## Add a 'variable' factor (stat x var combination)
table$variable2 <- as.factor(paste(table$stat, table$variable, sep = " "))
# summary(table$variable2) # guet


p1 <- ggplot(data = table, aes(x = variable2, y = mae)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("Mean Absolute Error (MAE)") +
    theme_bw() + theme(axis.text = element_text(size = 6))

p2 <- ggplot(data = table, aes(x = variable2, y = mae)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("Mean Absolute Error (MAE)") +
    theme_bw() + theme(axis.text = element_text(size = 6)) + 
    facet_wrap(~factor(region), ncol = 1, nrow = 3)

# Save them in plot.dir
setwd(plot.dir)
ggsave(plot = p1, filename = "boxplot_MAE_QM_all_strats_26_05_25.jpg", dpi = 300, width = 10, height = 3)
ggsave(plot = p2, filename = "boxplot_MAE_QM_all_strats_regions_26_05_25.jpg", dpi = 300, width = 10, height = 8)


## Same for the other 2 metrics: RMSE and corr coefficients
p1 <- ggplot(data = table, aes(x = variable2, y = rmse)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("RMSE") +
    theme_bw() + theme(axis.text = element_text(size = 6))

p2 <- ggplot(data = table, aes(x = variable2, y = rmse)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("RMSE") +
    theme_bw() + theme(axis.text = element_text(size = 6)) + 
    facet_wrap(~factor(region), ncol = 1, nrow = 3)

# Save them in plot.dir
ggsave(plot = p1, filename = "boxplot_RMSE_QM_all_strats_26_05_25.jpg", dpi = 300, width = 10, height = 3)
ggsave(plot = p2, filename = "boxplot_RMSE_QM_all_strats_regions_26_05_25.jpg", dpi = 300, width = 10, height = 8)

### --> MAE and RMSE are basically the same. Keep RMSE only.

p1 <- ggplot(data = table, aes(x = variable2, y = corr)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("Correlation coeff.") +
    theme_bw() + theme(axis.text = element_text(size = 6))

p2 <- ggplot(data = table, aes(x = variable2, y = corr)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("Correlation coeff.") +
    theme_bw() + theme(axis.text = element_text(size = 6)) + 
    facet_wrap(~factor(region), ncol = 1, nrow = 3)

# Save them in plot.dir
ggsave(plot = p1, filename = "boxplot_corr_QM_all_strats_26_05_25.jpg", dpi = 300, width = 10, height = 3)
ggsave(plot = p2, filename = "boxplot_corr_QM_all_strats_regions_26_05_25.jpg", dpi = 300, width = 10, height = 8)


### CONCLUSIONS:
### -> temperature variables show extremely good fits! Very encouraging
### -> total precip. more prone to errors but still ok (mean corr = 0.77)
### -> SM_10 even more prone to errors (mean corr = 0.70) - still OK or not? 


### Now plot distributions per QM strategy
p1 <- ggplot(data = table, aes(x = factor(qm_method), y = rmse)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("RMSE") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(factor(variable2) ~ factor(region), scales = "free_y")

p2 <- ggplot(data = table, aes(x = factor(qm_method), y = corr)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("Correlation coeff.") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(factor(variable2) ~ factor(region), scales = "free_y")

ggsave(plot = p1, filename = "boxplot_RMSE_QM_all_strats_regionsxstrats_26_05_25.jpg", dpi = 300, width = 12, height = 10)
ggsave(plot = p2, filename = "boxplot_corr_QM_all_strats_regionsxstrats_26_05_25.jpg", dpi = 300, width = 12, height = 10)

### CONCLUSIONS:
# max Ta_10: Global approach worse than the others (not by that much though), especially in the HND
# min Ta_10: No real differences...

# max Ta_200: Global slightly worse, but not huge differences
# min Ta_200: Global slightly worse, but not huge differences

# max Ts_05: Global approach worse than the others - in the HND, the 'mw' seems slightly better
# min Ts_05: Global approach worse than the others again (though still not by that much)

# max Ts_10: Global approach worse than the others - in the HND, the 'mw' seems slightly better?
# min Ts_10: Same

# max Ts_20: Global approach worse than the others - in the HND, the 'mw' seems slightly better?
# min Ts_20: Same

# max SM_10: No clear differences...Let's check metrics based on anomalies (see below)
# min SM_10: No clear differences...

# total precipitation: 'monthly' performs worst, quite clearly, in all 3 regions

### monthly approach and anomalies approach give VERY similar results


### PCA plot and analyze PC scores per strategies
library("FactoMineR")
# str(table); colnames(table)
## Input columns 2:5 in PCA
pca <- PCA(X = table[,c(2:5)], scale.unit = TRUE, graph = FALSE)
summary(pca)
#                        Dim.1   Dim.2
# Variance               2.903   0.930
# % of var.             70.981  24.889
# Cumulative % of var.  70.981  95.870

#         Dim.1    ctr   cos2    Dim.2    ctr   cos2 
# mbe  |  0.323  3.584  0.104 |  0.946 96.281  0.896 |
# mae  |  0.970 32.415  0.941 | -0.090  0.869  0.008 |
# rmse |  0.989 33.680  0.978 | -0.095  0.972  0.009 |
# corr | -0.938 30.321  0.880 |  0.132  1.878  0.017 |

## The higher the corr coeff, the lower the RMSE/MAE (as expected)
## --> either MAXIMIZE or MINIMIZE this, furthest away from 0
## PC2 = mbe, you will want to MINIMIZE that (= lowest mean errors)

## Plot PC1/2 per strategy/variables
table[,c("PC1","PC2")] <- pca$ind$coord[,c(1,2)]

p1 <- ggplot(data = table, aes(x = factor(qm_method), y = PC1)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    #geom_hline(yintercept = 0, linetype = "dashed") + 
    xlab("Variables") + ylab("PC1") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(.~factor(variable2), scales = "free_y")

p2 <- ggplot(data = table, aes(x = factor(qm_method), y = PC2)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    #geom_hline(yintercept = 0, linetype = "dashed") + 
    xlab("Variables") + ylab("PC2") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(.~factor(variable2), scales = "free_y")

setwd(plot.dir)
ggsave(plot = p1, filename = "boxplot_PC1_QM_all_stratsxvariables_26_05_25.jpg", dpi = 300, width = 7.5, height = 7.5)
ggsave(plot = p2, filename = "boxplot_PC2_QM_all_stratsxvariables_26_05_25.jpg", dpi = 300, width = 7.5, height = 7.5)


### CONCLUSIONS:
## -> Same as above obviously
## -> Avoid the 'global' strategy in general, especially for the soil temperature (Ts) variables
## -> looks like 'monthly' is slightly better for precipitation
## -> PC2/MBE shows novel aspects! the 'mw' strategy tends to be deliver more consistent 
##    negative (good) and lower MBE. Other strategies may deliver lower MBE in some plots
##    but also much higher MBE in other plots. Meanwhile, the 'mv' is less variable in terms
##    of MBE

### For now, I would say that the 'mv' approach tends to be better, especially compared to the 
### 'global' approach. Let's see for metrics based on anomalies though

### ---------------------------------------------------

# Go back to dir where eval metrics are stored
setwd(eval.metric.dir)

### B) Evaluation metrics based on the anomalies to the mean
files <- dir()[grepl("_anoms_",dir())]
# files
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
# Rbind
table <- bind_rows(res); rm(res); gc()

## Add a 'variable' factor (stat x var combination)
table$variable2 <- as.factor(paste(table$stat, table$variable, sep = " "))

setwd(plot.dir)

p1 <- ggplot(data = table, aes(x = variable2, y = rmse)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("RMSE") +
    theme_bw() + theme(axis.text = element_text(size = 6))

p2 <- ggplot(data = table, aes(x = variable2, y = rmse)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("RMSE") +
    theme_bw() + theme(axis.text = element_text(size = 6)) + 
    facet_wrap(~factor(region), ncol = 1, nrow = 3)

# Save them in plot.dir
ggsave(plot = p1, filename = "boxplot_RMSE_anomalies_QM_all_strats_26_05_25.jpg", dpi = 300, width = 10, height = 3)
ggsave(plot = p2, filename = "boxplot_RMSE_anomalies_QM_all_strats_regions_26_05_25.jpg", dpi = 300, width = 10, height = 8)

p1 <- ggplot(data = table, aes(x = variable2, y = corr)) +
    geom_violin(fill = "#00AFBB", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("Correlation coeff.") +
    theme_bw() + theme(axis.text = element_text(size = 6))

p2 <- ggplot(data = table, aes(x = variable2, y = corr)) +
    geom_violin(fill = "#00AFBB", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("Correlation coeff.") +
    theme_bw() + theme(axis.text = element_text(size = 6)) + 
    facet_wrap(~factor(region), ncol = 1, nrow = 3)

# Save them in plot.dir
ggsave(plot = p1, filename = "boxplot_corr_anomalies_QM_all_strats_26_05_25.jpg", dpi = 300, width = 10, height = 3)
ggsave(plot = p2, filename = "boxplot_corr_anomalies_QM_all_strats_regions_26_05_25.jpg", dpi = 300, width = 10, height = 8)



### CONCLUSIONS:
### -> Very similar to the eval metrics that were not based on anomalies
### -> SM_10 even worse than before as expected

### Now plot distributions per QM strategy
p1 <- ggplot(data = table, aes(x = factor(qm_method), y = rmse)) +
    geom_violin(fill = "#00AFBB", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("RMSE") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(factor(variable2) ~ factor(region), scales = "free_y")

p2 <- ggplot(data = table, aes(x = factor(qm_method), y = corr)) +
    geom_violin(fill = "#00AFBB", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("Correlation coeff.") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(factor(variable2) ~ factor(region), scales = "free_y")

ggsave(plot = p1, filename = "boxplot_RMSE_anomalies_QM_all_strats_regionsxstrats_26_05_25.jpg", dpi = 300, width = 12, height = 10)
ggsave(plot = p2, filename = "boxplot_corr_anomalies_QM_all_strats_regionsxstrats_26_05_25.jpg", dpi = 300, width = 12, height = 10)


### CONCLUSIONS:
# max Ta_10: Global approach worse than the others again (not by that much though - again)
# min Ta_10: No real differences...again

# max Ta_200: No real differences - all perform very high
# min Ta_200: Same - no notable differences

# max Ts_05: 'mw' seems best in the HND and SCH but NOT in SWA
# min Ts_05: 'mw' seems best overall, but not by much

# max Ts_10: 'mw' seems best overall
# min Ts_10: 'mw' seems best overall - very clearly for once

# max Ts_20: 'mw' seems best overall
# min Ts_20: 'mw' seems best overall - 'global' the worse

# max SM_10: 'global' seems best in the HND and SWA! surprising
# min SM_10: 'global' seems best in the HND and SWA too

# total precipitation: 'mw' seems best but not by much; 'monthly' is worst


### PCA plot and analyze PC scores per strategies
pca <- PCA(X = table[,c(2:5)], scale.unit = TRUE, graph = FALSE, ncp = 5)
summary(pca)
#                        Dim.1
# Variance               2.753
# % of var.             91.836
### -> Look at PC1 only since 92% of variance already

#         Dim.1    ctr   cos2  
# mbe  |  0.000  0.000  0.001 |
# mae  |  0.974 34.431  0.948 |
# rmse |  0.976 34.632  0.953 |
# corr | -0.923 30.937  0.852 |

## The higher the corr coeff, the lower the RMSE/MAE (as expected)
## --> either MAXIMIZE or MINIMIZE this, furthest away from 0

## Plot PC1/2 per strategy/variables
table[,c("PC1","PC2")] <- pca$ind$coord[,c(1,2)]

p1 <- ggplot(data = table, aes(x = factor(qm_method), y = PC1)) +
    geom_violin(fill = "#00AFBB", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    #geom_hline(yintercept = 0, linetype = "dashed") + 
    xlab("Variables") + ylab("PC1") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(.~factor(variable2), scales = "free_y")

setwd(plot.dir)
ggsave(plot = p1, filename = "boxplot_PC1_anomalies_QM_all_stratsxvariables_26_05_25.jpg", dpi = 300, width = 7.5, height = 7.5)


### MBE is not well captured by that PCA somehow. Examine it's distribution directly.
p2 <- ggplot(data = table, aes(x = factor(qm_method), y = mbe)) +
    geom_violin(fill = "#00AFBB", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    #geom_hline(yintercept = 0, linetype = "dashed") + 
    xlab("Variables") + ylab("Mean bias") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(.~factor(variable2), scales = "free_y")

ggsave(plot = p2, filename = "boxplot_MBE_anomalies_QM_all_stratsxvariables_26_05_25.jpg", dpi = 300, width = 7.5, height = 7.5)


### CONCLUSIONS:
## -> 'global' is worst, 'mw' is slightly better than the rest except for total precipitation 
##    for which 'monthly' is better
## -> MBE gives no additional information about which QM strategy provides the most 'constrained' variability in errors

### -------------------------------------------------

### 04/06/2025: Examining distribution of Refined Index of Agreement (RIA) across QM strategies and variables
setwd(eval.metric.dir) # dir()

### A) RIA distrbution based on normal data
files <- dir()[!grepl("_anoms_",dir())] # files # gut
res <- lapply(files, function(f) { d <- get(load(f)); return(d) } )
# Rbind
table <- bind_rows(res); rm(res); gc()
# dim(table) ; head(table) ; summary(table)
## Add a 'variable' factor (stat x var combination)
table$variable2 <- as.factor(paste(table$stat, table$variable, sep = " "))

p <- ggplot(data = table, aes(x = factor(qm_method), y = d1r)) +
    geom_violin(fill = "#e31a1c", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("Refined Index of Agreement") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(factor(variable2) ~ factor(region), scales = "free_y")

setwd(plot.dir) # dir()
ggsave(plot = p, filename = "boxplot_RIA_QM_all_strats_regionsxstrats_04_06_25.jpg", dpi = 300, width = 12, height = 10)


### B) RIA distribution based on anomalies to the monthly means
setwd(eval.metric.dir)
files <- dir()[grepl("_anoms_",dir())]
res <- lapply(files, function(f) { d <- get(load(f)); return(d) } )
table <- bind_rows(res); rm(res); gc()
# dim(table) ; head(table) ; summary(table)
## Add a 'variable' factor (stat x var combination)
table$variable2 <- as.factor(paste(table$stat, table$variable, sep = " "))

p <- ggplot(data = table, aes(x = factor(qm_method), y = d1r)) +
    geom_violin(fill = "#fc4e2a", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Variables") + ylab("Refined Index of Agreement") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(factor(variable2) ~ factor(region), scales = "free_y")

setwd(plot.dir)
ggsave(plot = p, filename = "boxplot_RIA_QM_all_strats_regionsxstrats_anomalies_04_06_25.jpg", dpi = 300, width = 12, height = 10)


### CONCLUSIONS: 

## -> 'monthly' is worst for total precipitation - discard x
## -> 'global' is worst for many temperature-related variables, but it looks OK for SM_10 or total precip. - discard x
## -> use either 'mw' or 'anomalies'

### ------------------------------------------------------------------------------------------------------------

### 29/04/25: Write another FUN to make plots such as histograms, time series or scatter plots
### The FUN will understand its arguments as 'switches' to display certain type of plots 
### It will use the same plots arguments to make them based on biases (model - obs) and
### monthly anomalies too thanks to the 'mon.anoms' and 'biases' arguments


###  Associated plotting FUN
# To test while you're writing it:
# region <- "HND" 
# var <- "Ta_10"
# stat <- "max"
# method <- "global"
# mon.anoms <- FALSE
# boxp <- TRUE
# histo <- TRUE
# scatt <- TRUE
# time_series <- TRUE
# biases <- TRUE 
# std.anoms <- TRUE


plot_daily_stat_comparison <- function(region, var, stat, method, mon.anoms, boxp, histo, scatt, time_series) {

    #' This function takes 9 arguments and returns a formatted data.frame:
    #' @param var the climate variable to process (character) - one of the following: 
    #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20", "Ts_50", "precipitation" or "SM_10"
    #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character): 'mean', 'max' or 'min' 
    #' @param region Region name (character): 'SCH' or 'HND' or 'SWA'
    #' @param method Which quantile mapping approach to run (character): 'global' or 'monthly' or 'mw' or 'anoms' 
    #' @param mon.anoms Switch (BOOLEAN) - Whether plots should be made based on anomalies to the monly mean 
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

      # Read the file containing your observed daily stats
      obs_daily_stat <- read.csv(files, h = T, sep = ",", dec = ".")
      colnames(obs_daily_stat)[3] <- paste(var,stat, sep = "_")

      # Sanity check
      if( exists("obs_daily_stat") == FALSE ) {
          stop(
            paste("!!! ERROR: Could not load observed ",paste(var,stat, sep = "_")," file for the ",region, sep = "")
          )
      } # eo if loop - sanity check

      ## Go to quantile mapping outputs directory and load corresponding dataset
      setwd(paste("/home/fbenedetti/E-OBS/Explos/quantile_mapping_outputs/",method, sep = ""))
      file <- dir()[grepl(paste(stat,var,region, sep = "_"), dir())]
      corr_data <- get(load(file)) # dim(corr_data) ; str(corr_data) ; colnames(corr_data)
        
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

      # Merge
      merged_df <- reduce(
                    list(obs_daily_stat[,names],
                         corr_data[,names]
                    ), full_join, by = c("EP","Date")
      ) # eo reduce
      # summary(merged_df) # OK
      rm(obs_daily_stat,corr_data) ; gc()

      # Adjust colnames
      colnames(merged_df)[c(3,4)] <- c("obs","corr_model") 
        
      # Remove NAs for evaluation (usually data < 2008 of course)
      merged_df <- na.omit(merged_df)

      # Vector of EP IDs - use to subset 'merged_df'
      plots <- unique(merged_df$EP)

      ## 26/05/25: Restrict to Dates > '2009-05-01' because of faulty obs 2008 data
      merged_df <- merged_df[merged_df$Date >= '2009-05-01',]

      ### SWITCHES:
      ### If anoms == TRUE --> make plots based on anomalies to the monthly mean
      ### else --> make plots on normal data

      if( mon.anoms ) {

          # Compute monthly means for each data source
          merged_df <- merged_df %>%
                mutate(month = format(Date, "%Y-%m")) %>%
                group_by(EP,month) %>%
                mutate(
                    mon_mean_obs = mean(obs, na.rm = TRUE),
                    mon_mean_mod = mean(corr_model, na.rm = TRUE)
                ) %>% 
                ungroup()

          # Compute monthly anomalies
          merged_df <- merged_df %>%
                mutate(
                    mon_anom_obs = obs - mon_mean_obs,
                    mon_anom_mod = corr_model - mon_mean_mod
                )
            
          ### Save plots in plot.dir
          setwd(plot.dir)
          colors_sources <- c("obs" = "#00AFBB", "model" = "#E7B800")

          # 1st: boxplots 
          if( boxp ) {
                # Boxplot of stat distribution - needs melting
                m_merged_df <- reshape2::melt(merged_df[,c("EP","Date","mon_anom_obs","mon_anom_mod")],
                                        id.vars = c("EP","Date"))
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
                ggsave(plot = plot, filename = paste("boxplot_evaluation_mon_anomalies_",method,"_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 5, height = 5)  
          } # eo if loop - boxp

          # 2nd: histograms - also need melting
          if( histo ) {
                # Needs melting too
                m_merged_df <- reshape2::melt(merged_df[,c("EP","Date","mon_anom_obs","mon_anom_mod")],
                                        id.vars = c("EP","Date"))
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
                ggsave(plot = plot, filename = paste("gghistogram_evaluation_mon_anomalies_",method,"_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 7, height = 5)  
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
                ggsave(plot = plot, filename = paste("scatter_evaluation_mon_anomalies_",method,"_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 5, height = 5)  
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
                 ggsave(plot = p, filename = paste("time_series_evaluation_mon_anomalies_",method,"_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 10, height = 4)  
          } # eo if loop - time_series

      } else {

          ### Save plots depending on the switches activated as function's arguments: boxp, histo, scatt, time_series
          setwd(plot.dir)
          colors_sources <- c("obs" = "#00AFBB", "model" = "#E7B800")

          # 1st: boxplots 
          if( boxp ) {
                # Boxplot of stat distribution - needs melting
                m_merged_df <- reshape2::melt(merged_df, id.vars = c("EP","Date"))
                colnames(m_merged_df)[c(3)] <- c("source")
                # Re-name source levels so they match the 'colors_sources' palette above; unique(m_merged_df$source)
                m_merged_df$source <- as.character(m_merged_df$source)
                m_merged_df[m_merged_df$source == "obs","source"] <- "obs"
                m_merged_df[m_merged_df$source == "corr_model","source"] <- "model"
                # Make the plot
                plot <- ggplot(data = m_merged_df, aes(x = factor(source), y = value, fill = factor(source))) +
                  geom_violin(colour = "black") + geom_boxplot(colour = "black", fill = "white", width = .2) + 
                  scale_fill_manual(values = colors_sources) + labs(y = paste(stat,var,sep = " "), x = "") + 
                  theme_bw() + theme(legend.position = "none") +
                  ggtitle(paste("Distribution of ",stat," ",var," in the ",region, sep = ""))
                # Save plot in dir  
                ggsave(plot = plot, filename = paste("boxplot_evaluation_",method,"_",stat,"_",var,"_",region,".jpg",sep = ""),
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
                m_merged_df[m_merged_df$source == "corr_model","source"] <- "model"
                # Make the plot
                plot <- gghistogram(m_merged_df, x = "value", add = "mean", rug = FALSE,
                          fill = "source", palette = colors_sources,
                          ylab = "Count", xlab = paste(stat,var,sep = " "),
                          title = paste("Distribution of ",stat," ",var," in the ",region, sep = "")
                )
                # Save plot in dir  
                ggsave(plot = plot, filename = paste("gghistogram_evaluation_",method,"_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 7, height = 5)  
          } # eo if loop - histo

          # 3rd: scatt - make 2: obs vs corrected model outputs
          if( scatt ) {
                # Make plot 1
                plot <- ggplot(merged_df, aes(x = obs, y = corr_model)) + geom_point(alpha = 0.2, colour = "grey50") +
                        geom_abline(intercept = 0, slope = 1, color = "#d53e4f", linetype = "dashed") + 
                        geom_smooth(method = "lm", se = TRUE, colour = "#00AFBB") + 
                        labs(title = paste("Local vs model daily ",stat," ",var, sep = ""),
                            x = "Local measurements", y = "Model (corrected)") +
                        theme_bw()
                # Save
                ggsave(plot = plot, filename = paste("scatter_evaluation_",method,"_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 5, height = 5)  
          } # eo if loop - histo

          # 4th: time_series
          if( time_series ) {
                p <- ggplot(merged_df, aes(x = Date)) + 
                  geom_line(aes(y = obs), linewidth = 1, colour = "#00AFBB", alpha = .25) +
                  geom_line(aes(y = corr_model), linewidth = 1, colour = "#E7B800", alpha = .25) +
                  labs(title = paste("Time series of ",stat," ",var," in the ",region," (all EPs)", sep = ""),
                      x = "Date", y = paste(stat,var,sep = " ")) +
                  theme_minimal()
                # Save
                 ggsave(plot = p, filename = paste("time_series_evaluation_",method,"_",stat,"_",var,"_",region,".jpg",sep = ""),
                      dpi = 300, width = 10, height = 4)  
            } # eo if loop - time_series

          } # eo if else loop - anoms == T

} # eo plotting FUN - plot_daily_stat_comparison


for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","precipitation","SM_10")) {
    for(s in c("max","min")) {
        for(r in c("HND","SCH","SWA")) {
            for(m in c('global','monthly','mw','anomalies')) {

                plot_daily_stat_comparison(
                      region = r,
                      var = v,
                      stat = s,
                      method = m,
                      mon.anoms = TRUE, ### 26/05/25: Need to re-run with mon.anoms = TRUE
                      boxp = TRUE,
                      histo = TRUE,
                      scatt = TRUE,
                      time_series = TRUE
                ) 

            } # eo for loop - m
        } # eo for loop - r
    } # eo for loop - s
} # eo for loop - v


### ------------------------------------------------------------------------------------------------------------

### 21/01/26: Write evaluate_forest_quantile_maps() for forest precipitation reconstruction 
# To test evaluate_quantile_maps while you're writing it:
# var <- "precipitation"
# stat <- "total"
# region <- "HND"
# method <- "mw"
# anoms <- FALSE

evaluate_forest_quantile_maps <- function(var = "precipitation", stat = "total", region, method, anoms) {

        #' This function takes four arguments and returns a formatted data.frame:
        #' @param var the climate variable to process (character) - Default = "precipitation"
        #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character) - Default = "total"
        #' @param region Region name (character): 'SCH' or 'HND' or 'SWA'
        #' @param method Which quantile mapping approach to run (character): 'global' or 'monthly' or 'mw' or 'anomalies'
        #' @param anoms Whether evaluation metrics should be computed on anomalies to the monthly mean (BOOLEAN)
        #' @return A formatted data.frame combining the daily statistics

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
        file <- dir()[grepl(paste("forest", sep = ""),dir())]
        file <- file[grepl(region,file)]
        # Load 'file'
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
        # Clean ram
        rm(obs_daily_stat,corr_data)
        gc()

        # Adjust colnames
        colnames(merged_df)[c(3,4)] <- c("obs","corr_model") 
        
        # Remove NAs for evaluation (usually data < 2008 of course)
        merged_df <- na.omit(merged_df)

        # Vector of EP IDs - use to subset 'merged_df'
        plots <- unique(merged_df$EP)

        ### Calculate evaluation metrics: mean bias, corr coef, RMSE, mean absolute bias (MAE) for each EP separately!

        if( anoms ) {

          # Compute monthly means for each data source
          merged_df <- merged_df %>%
              mutate(month = format(Date, "%Y-%m")) %>%
              group_by(EP,month) %>%
              mutate(
                  mon_mean_obs = mean(obs, na.rm = TRUE),
                  mon_mean_mod = mean(corr_model, na.rm = TRUE)
              ) %>% 
              ungroup()
          # summary(merged_df)

          # Compute monthly anomalies
          merged_df <- merged_df %>%
              mutate(
                  mon_anom_obs = obs - mon_mean_obs,
                  mon_anom_mod = corr_model - mon_mean_mod
              )
          # summary(merged_df)
            
          # Subset merged_df per EP and return evaluation metrics
          eval_metrics <- lapply(plots, function(p) {
                  
                  # p <- plots[13]
                  message(paste("Computing evaluation metrics for EP: ",p," based on monthly anomalies", sep = ""))
                  ## Subset
                  sub_merged_df <- merged_df[merged_df$EP == p,]
                  ## Only keep Dates > '2009-05-01' because of obs data
                  sub_merged_df <- sub_merged_df[sub_merged_df$Date >= '2009-05-01',]
                  ## Mean bias error (MBE)
                  mbe <- mean(sub_merged_df$mon_anom_mod - sub_merged_df$mon_anom_obs)
                  ## Mean absolute error (MAE)
                  mae <- mean(abs(sub_merged_df$mon_anom_mod - sub_merged_df$mon_anom_obs))
                  ## Root mean square error (RMSE)
                  rmse <- sqrt(mean((sub_merged_df$mon_anom_mod - sub_merged_df$mon_anom_obs)^2, na.rm = TRUE))
                  ## Corr coeff
                  corr_coeff <- cor(sub_merged_df$mon_anom_obs, sub_merged_df$mon_anom_mod, use = "complete.obs")
                  ## RIA
                  ria <- d1r(obs = sub_merged_df$mon_anom_obs, pred = sub_merged_df$mon_anom_mod, tidy = TRUE)
                      
                  # Return evaluation metrics in a data.frame
                  eval_metrics <- data.frame(EP = p,
                        mbe = mbe,
                        mae = mae, 
                        rmse = rmse, 
                        corr = corr_coeff,
                        ria = ria
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
          table_eval_metrics$variable <- var
          table_eval_metrics$stat <- stat
          table_eval_metrics$qm_method <- method
            
          # Save in dir as .Rdata
          message(paste("Saving evaluation metrics table for QM outputs (anomalies-based) based on ",method," ", stat," ",var," in the forests of the ",region, sep = ""))
          setwd(eval.metric.dir)
          save(x = table_eval_metrics, file = paste("table_evaluation_metrics_anoms_",method,"_",stat,"_",var,"_forests_",region,".Rdata", sep = ""))
          rm(merged_df,table_eval_metrics); gc()

        } else {

          # Subset merged_df per EP and return evaluation metrics
          eval_metrics <- lapply(plots, function(p) {

                  # p <- sample(plots,1,1); p # for testing
                  message(paste("Computing evaluation metrics for ",p, sep = ""))
                  
                  ## subset
                  sub_merged_df <- merged_df[merged_df$EP == p,]
                  ## Only keep Dates > '2009-05-01' because of obs data
                  sub_merged_df <- sub_merged_df[sub_merged_df$Date >= '2009-05-01',]
                  ## Mean bias error (MBE)
                  mbe <- mean(sub_merged_df$corr_model - sub_merged_df$obs)
                  ## Mean absolute error (MAE)
                  mae <- mean(abs(sub_merged_df$corr_model - sub_merged_df$obs))
                  ## Root mean square error (RMSE)
                  rmse <- sqrt(mean((sub_merged_df$corr_model - sub_merged_df$obs)^2, na.rm = TRUE))
                  ## Corr coeff
                  corr_coeff <- cor(sub_merged_df$obs, sub_merged_df$corr_model, use = "complete.obs")
                  ## RIA
                  ria <- d1r(obs = sub_merged_df$obs, pred = sub_merged_df$corr_model, tidy = TRUE)
                      
                  # Return evaluation metrics in a data.frame
                  eval_metrics <- data.frame(EP = p,
                        mbe = mbe,
                        mae = mae,
                        rmse = rmse,
                        corr = corr_coeff,
                        ria = ria
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
          table_eval_metrics$variable <- var
          table_eval_metrics$stat <- stat
          table_eval_metrics$qm_method <- method
            
          message(paste("\n","Saving evaluation metrics table for QM outputs based on ",method," ",stat," ",var," in the ",region, sep = ""))
          setwd(eval.metric.dir)
          save(x = table_eval_metrics, file = paste("table_evaluation_metrics_",method,"_",stat,"_",var,"_forests_",region,".Rdata", sep = ""))
          rm(merged_df,table_eval_metrics); gc()

      } # eo if else loop - anomalies

} # eo master FUN - evaluate_quantile_maps

### Apply evaluate_quantile_maps() for all data (variables, stats, QM method, regions...)

for(r in c("HND","SCH","SWA")) {
    for(m in c("global","monthly","mw","anomalies")) {
        evaluate_forest_quantile_maps(var = "precipitation", stat = "total", region = r, method = m, anoms = FALSE)
        evaluate_forest_quantile_maps(var = "precipitation", stat = "total", region = r, method = m, anoms = TRUE)
    } # eo for loop - m
} # eo for loop - r


### --------------------------------------------------

### 21/01/26: Examine distrbution of evaluation metrics for the forest precipitation data 
setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/evaluation_metrics")
dir()[grepl("forest",dir())] # 24 tables because 3 regions x 4 QM methods x 2 evaluation scales (anomalies to the monthly mean or not)


### A) Evaluation metrics based on raw data
files <- dir()[grepl("forest",dir())]
files <- files[!grepl("_anoms_",files)] # files # should be length 12 (3x4)
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
table <- bind_rows(res); rm(res); gc()
# dim(table) ; head(table) ; summary(table)

### The goal is twofold here:
### - 1) get a sense of the quality of the QM corrections overall
### - 2) determine which QM strategy is the best overall
### Note: it could be that one strategy is not the best across ALL variables 

## Need to make two types of plots:
## (i)  Distribution of RMSE/corr coeff/MAE (+ facet per regions)
## (ii) Distribution of RMSE/corr coeff/MAE and QM strategy (+ facet per regions)

p1 <- ggplot(data = table, aes(x = factor(region), y = mae)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Region") + ylab("Mean Absolute Error (MAE)") +
    theme_bw() + theme(axis.text = element_text(size = 6))

# Save them in plot.dir
setwd(plot.dir)
ggsave(plot = p1, filename = "boxplot_MAE_QM_forest_precip_all_strats_regions_21_01_26.jpg", dpi = 300, width = 4, height = 4)


## Same for the other 2 metrics: RMSE and corr coefficients
p1 <- ggplot(data = table, aes(x = factor(region), y = rmse)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Region") + ylab("RMSE") +
    theme_bw() + theme(axis.text = element_text(size = 6)) 

# Save them in plot.dir
ggsave(plot = p1, filename = "boxplot_RMSE_QM_forest_precip_all_strats_regions_21_01_26.jpg", dpi = 300, width = 4, height = 4)

### --> MAE and RMSE are basically the same. Keep RMSE only.

p1 <- ggplot(data = table, aes(x = factor(region), y = corr)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Region") + ylab("Correlation coeff.") +
    theme_bw() + theme(axis.text = element_text(size = 6)) 

# Save them in plot.dir
ggsave(plot = p1, filename = "boxplot_CORR_QM_forest_precip_all_strats_regions_21_01_26.jpg", dpi = 300, width = 4, height = 4)

### MBE
p1 <- ggplot(data = table, aes(x = factor(region), y = mbe)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Region") + ylab("Mean bias error.") +
    theme_bw() + theme(axis.text = element_text(size = 6)) 

# Save them in plot.dir
ggsave(plot = p1, filename = "boxplot_MBE_QM_forest_precip_all_strats_regions_21_01_26.jpg", dpi = 300, width = 4, height = 4)

### RIA
p1 <- ggplot(data = table, aes(x = factor(region), y = d1r)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("Region") + ylab("RIA") +
    theme_bw() + theme(axis.text = element_text(size = 6)) 

# Save them in plot.dir
ggsave(plot = p1, filename = "boxplot_RIA_QM_forest_precip_all_strats_regions_21_01_26.jpg", dpi = 300, width = 4, height = 4)

### CONCLUSIONS:
### -> SWA shows the highest errors, and HND the lowest
### -> HND shows the lowest mean bias errors (very close to 0)
### -> SWA shows the best corr coeff overall though (~0.8) and best RIA so not that bad
### -> corr coeff overall range between 0.7 and 0.85 so good to reasonable :)
### -> same with RIA actually. Good distrbution as mostly > 0.75

### Now plot distributions per QM strategy
p1 <- ggplot(data = table, aes(x = factor(qm_method), y = rmse)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("QM method") + ylab("RMSE") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(.~ factor(region), scales = "free_y")

p2 <- ggplot(data = table, aes(x = factor(qm_method), y = corr)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("QM method") + ylab("Correlation coeff.") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(.~ factor(region), scales = "free_y")

p3 <- ggplot(data = table, aes(x = factor(qm_method), y = mbe)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("QM method") + ylab("Mean bias error") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(.~ factor(region), scales = "free_y")

p4 <- ggplot(data = table, aes(x = factor(qm_method), y = d1r)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("QM method") + ylab("RIA") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(.~ factor(region), scales = "free_y")

ggsave(plot = p1, filename = "boxplot_RMSE_QM_forest_precip_all_strats_regionsxstrats_21_01_26.jpg", dpi = 300, width = 6, height = 3)
ggsave(plot = p2, filename = "boxplot_corr_QM_forest_precip_all_strats_regionsxstrats_21_01_26.jpg", dpi = 300, width = 6, height = 3)
ggsave(plot = p3, filename = "boxplot_MBE_QM_forest_precip_all_strats_regionsxstrats_21_01_26.jpg", dpi = 300, width = 6, height = 3)
ggsave(plot = p4, filename = "boxplot_RIA_QM_forest_precip_all_strats_regionsxstrats_21_01_26.jpg", dpi = 300, width = 6, height = 3)


### CONCLUSIONS:
### -> 'monthly' QM strat is the worst performing, as for the grasslands obviously
### -> 'anomalies' and 'global' QM strat are the top performers
### -> But, in the SWA, 'mw' shows the lowest MBE
### -> Keep 'mw' as the standard to be consistent



### B) Evaluation metrics based on ANOMALIES TO THE MONTHLY MEAN
setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/evaluation_metrics")
files <- dir()[grepl("forest",dir())]
files <- files[grepl("_anoms_",files)] # files # should be length 12 (3x4)
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
table <- bind_rows(res); rm(res); gc()
# dim(table) ; head(table) ; summary(table)

# Save plots in plot.dir
setwd(plot.dir)

### Now plot distributions per QM strategy
p1 <- ggplot(data = table, aes(x = factor(qm_method), y = rmse)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("QM method") + ylab("RMSE") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(.~ factor(region), scales = "free_y")

p2 <- ggplot(data = table, aes(x = factor(qm_method), y = corr)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("QM method") + ylab("Correlation coeff.") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(.~ factor(region), scales = "free_y")

p3 <- ggplot(data = table, aes(x = factor(qm_method), y = mbe)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("QM method") + ylab("Mean bias error") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(.~ factor(region), scales = "free_y")

p4 <- ggplot(data = table, aes(x = factor(qm_method), y = d1r)) +
    geom_violin(fill = "#E7B800", colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    xlab("QM method") + ylab("RIA") +
    theme_bw() + theme(axis.text = element_text(size = 5)) +
    facet_wrap(.~ factor(region), scales = "free_y")

ggsave(plot = p1, filename = "boxplot_RMSE_QM_forest_precip_anoms_all_strats_regionsxstrats_21_01_26.jpg", dpi = 300, width = 6, height = 3)
ggsave(plot = p2, filename = "boxplot_corr_QM_forest_precip_anoms_all_strats_regionsxstrats_21_01_26.jpg", dpi = 300, width = 6, height = 3)
ggsave(plot = p3, filename = "boxplot_MBE_QM_forest_precip_anoms_all_strats_regionsxstrats_21_01_26.jpg", dpi = 300, width = 6, height = 3)
ggsave(plot = p4, filename = "boxplot_RIA_QM_forest_precip_anoms_all_strats_regionsxstrats_21_01_26.jpg", dpi = 300, width = 6, height = 3)

### Same observations and conclusions as above with normal values

### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------