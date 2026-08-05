### ------------------------------------------------------------------------------------------------------------

### 31/07/25 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to load the datasets created by Script#4.4.1 and apply the GAMM approach developed in
### Script#4.4.1 to variables other than daily max Ta_200 (those variables listed at the beginning 
### of Script#4.4.2: "Ta_10","Ts_05","Ts_10","Ts_20","SM_10" 

### See Script#4.4.1 and Script#4.4.2 for exhaustive details about the approach and why we carry it out.  
### Based on the analyses carried out in Script#4.4.2, the empirical regressive models should follow this
### structure: 
#mod <- gamm(q
#    Offset ~ s(StandAge, k = 20) + s(DOY, bs = "cc", by = TreeType) + TreeType + Region,
#    random = list(EP = ~1),
#    correlation = corAR1(form = ~ 1|EP),
#    data = df
#) # takes a ~30min, even without correlation terms

### Last update: 30/06/26 (Summarizing GAMM skill metrics in Table S3 for ESSD draft)

### ------------------------------------------------------------------------------------------------------------

library("purrr")
library("tidyr")
library("dplyr")
library("data.table")
library("reshape2")
library("lubridate")
library("ggplot2")
library("ggpubr")
library("cowplot")
library("viridis")
library("nlme")
library("lme4")
library("mgcv")
library("MuMIn")
library("performance")
library("broom.mixed")
# sessionInfo() # to see the versions of the package used

### ------------------------------------------------------------------------------------------------------------

### Add daily total precipitation as a side data.frame for Ts
#res <- lapply(c("HND","SCH","SWA"), function(r) {
#
#      # Useless message
#      message(paste("Extracting precipitation and relative humidity data for the grasslands and forests of the ",r, sep = ""))
#
#      ## Grasslands
#      # Load precipitation data
#      setwd(paste("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily/grasslands/",r,"/", sep = ""))
#      precip_grass <- read.csv(dir()[grepl("interp_daily_precipitation",dir())], h = TRUE, sep = ",", dec = ".")
#      precip_grass$Date <- as.Date(precip_grass$Date) # str(precip); precip[is.na(precip$precipitation_radolan),][1:100,] # Keep NaN for now
#
#      ## Forests
#      setwd(paste("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily/forests/",r,"/", sep = ""))
#      # Load precipitation data
#      precip_forest <- read.csv(dir()[grepl("interp_daily_precipitation",dir())], h = TRUE, sep = ",", dec = ".")
#      precip_forest$Date <- as.Date(precip_forest$Date)
#
#      ## Combine
#      df <- rbind(precip_grass,precip_forest)
#      # str(df) ; unique(df$EP)
#      colnames(df)[3] <- "precipitation"
#
#      # Return
#      rm(precip_grass,precip_forest); gc()
#      return(df)
#
#  } # eo FUN - r
#
#) # eo lapply - region
#
## Rbind
#table_precip <- bind_rows(res)
#rm(res); gc()
### Add precipitation and rH data based on EP and Date
#colnames(table_precip)[3] <- "precipitation"

### ------------------------------------------------------------------------------------------------------------

### Write a master FUN that will use the structure of the GAMM developed at the end of Script#4.4.2 to perform a
### similar approach on all variables. The FUN should save the model objects issued by gamm()

### Files should be in:
# setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data"); dir()

### 29/01/26: IDENTIFIED ISSUE FOR THE SM GAMMs: INTRODUCING CIRCULARITY IN THE LONG TERM SHIFTS BY ADDING MEAN GRASSLAND VALUES 
### -> GRASSLAND VALUES WILL BE USED FOR THE FINAL RECONSTRUCTION SO THE LONG TERM DRIFT IS IMPRINTED TWICE
### -> NEED TO CORRECT AND RE-RUN THE FUNCTION BY REMOVING THE SM_10 IF ELSE LOOP

# For testing master FUN below while you write it
# var = "Ta_10"
# stat = "min"
# age = 75
# years = 14

gamm_offset <- function(var, stat, age, years, cap = 40) {

        #' This function takes three arguments and returns a model object of class 'gamm':
        #' @param var The variable to model (character): "Ta_200","Ta_10","Ts_05","Ts_10","Ts_20","SM_10" 
        #' @param stat The daily stat of the associated variable (character): 'max' or 'min'
        #' @param age Maximum stand age to account for in the GAMM approach (integer) 
        #' @param years Number of years of daily data to train the GAMM (integer): 6, 8, 10, 12, 14 or 16
        #' 'age' Should be between 50 and 70 based on R Script#4.4.2
        #' @param cap Age at which canopy closure should happen (s(StangeAge) smooth is fixed beyond this age for GAMM
        #' training and prediction) (integer): Default is 40. Should be a value between 30 or 50. 
        #' @return A formatted data.frame combining the inputs.
      
        # Useless message
        message(paste("Loading the aggregated data for ",stat," ",var, sep = ""))

        # Read in the data after identifying the corresponding
        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data") #; dir()

        file <- dir()[grepl(paste("metadata",stat,var,sep = "_"),dir())] # file
        df <- readRDS(file)
        # dim(df); str(df) # should be 858'600 daily measurements

        # Add month and date from DOY and Year
        df$Date <- as.Date(df$DOY - 1, origin = paste0(df$Year, "-01-01"))
        df$Month <- lubridate::month(df$Date)

        # Add sine and cosine of DOY
        df$sin_DOY <- sin(2 * pi * df$DOY / 365)
        df$cos_DOY <- cos(2 * pi * df$DOY / 365)

        # Rarefy to have a bablanced dataset to mode offsets as a function of stand age and tree types
        decid_df <- df[df$TreeType == "Deciduous" & df$StandAge < age,]
        coni_df <- df[df$TreeType == "Coniferous" & df$StandAge < age,]

        # Rbind 
        combin_df <- rbind(decid_df,coni_df)
        rm(coni_df,decid_df,file,df,file); gc()

        # To model a different seasonal pattern (smooth over DOY) for each TreeType, use the by = argument inside the smooth function
        combin_df$TreeType <- factor(combin_df$TreeType)
        combin_df$Region <- factor(combin_df$Region)
        colnames(combin_df)[4] <- "Offset"

        ### 14/10/25: Remove rows with NaN
        combin_df <- combin_df %>% drop_na(StandAge)
        combin_df <- combin_df %>% drop_na(Offset)

        ### 19/11/25: For temperature-related variables: Add anomalies to DOY to better capture extremes and add information about finer scale variability
        ### Is it important to include a first-order autoregressive structure in residuals within each EP ("residuals at day t are correlated with day t–1, with correlation ρ") ? 
        ### The AR(1) structure affects: standard errors and inference (p-values, edf) and sometimes the level of smoothing chosen by the model.
        ### But it rarely changes the shape of the DOY seasonal cycle very much and thus the resulting anomalies
        
        ### 19/11/25: Let's test this for a 6-year subset of the data (so it goes faster) based on a simple acf() test
        ### Fit both models. Compare ACF of residuals from Model A. If strong autocorrelation remains -> use Model B. If not -> Model A is fine and simpler.
        #test_df <- combin_df[combin_df$Year %in% c(2019:2024),]
        #m1 <- gamm(Mean_Grassland_Value ~ s(DOY, bs = "cc"), random = list(EP = ~1), data = test_df)
        #m2 <- gamm(Mean_Grassland_Value ~ s(DOY, bs = "cc"), random = list(EP = ~1), correlation = corAR1(form = ~ 1|EP), data = test_df)
        ## summary(m1$gam) ; summary(m2$gam)
        #acf_m1 <- acf(residuals(m1$lme), main = "ACF – DOY Model without AR(1)")
        #acf_m2 <- acf(residuals(m2$lme), main = "ACF – DOY Model with AR(1)")
        #pdf(file = paste(getwd(),"/",paste("plot_acf_DOY_anom_noAR1_19.11.25.pdf", sep = ""),sep = ""), width = 5, height = 3.5)
        #acf_m1
        #dev.off()
        #pdf(file = paste(getwd(),"/",paste("plot_acf_DOY_anom_withAR1_19.11.25.pdf", sep = ""),sep = ""), width = 5, height = 3.5)
        #acf_m2
        #dev.off()
        ### -> ACF shows the same results: AR(1) structure did not change the residuals
        ### -> DOY anomalies (observed – smooth prediction) will be the same to within noise
        ### -> Most likely the cyclic smooth s(DOY, bs="cc") captures nearly all temporal structure

        ### 29/01/26: Computing anomalies to DOY based on a simple GAM FOR ALL VARIABLES - NOT JUST TEMPERATURE
        ### WHY s(Year) SHOULD NOT BE ADDED: 
        ### - Anoms become residuals around a detrended (unreal) climate
        ### - Anoms will show weaker variance and extremes because they will learn "damped/buffered" anomalies
        ### ADVICE: Never “clean” the training data to fix a prediction-time problem.
        ###         Instead, constrain what the model is allowed to see when extrapolating.
        message(paste("Computing anomalies to DOY based on a simple GAM", sep = ""))
        tmp <- gamm(Mean_Grassland_Value ~ s(DOY, bs = "cc"), random = list(EP = ~1), data = combin_df) 
        # combin_df$Anom <- with(combin_df, Mean_Grassland_Value - predict(tmp$gam, newdata = combin_df))
        # combin_df$Anom <- NA
        # combin_df[!is.na(combin_df$Mean_Grassland_Value),"Anom"] <- resid(tmp$gam)
        combin_df$Anom <- resid(tmp$gam)
        # summary(tmp$gam) ; summary(combin_df$Anom)
        
        ### 15/10/25-22/10: Train model on 6 representative years of data instead of 16 -> reduce risk of OOM errors due to the LME in the GAMM
        # combin_df <- combin_df[combin_df$Year %in% c(2010,2013,2016,2019,2022,2024),] # 6 years version, 70'114 rows
        # combin_df <- combin_df[combin_df$Year %in% c(2009,2010,2011,2012,2014,2016,2018,2020,2022,2024),] # 10 years version, 116'839
        # combin_df <- combin_df[combin_df$Year %in% c(2009:2023),] # dim(combin_df)
        # dim(combin_df) ; summary(combin_df) ; unique(combin_df$EP)

        ### 20/11/25: Set number of years of daily obervations to train the GAMM - main factor influencing computational cost
        ### Possible values are: 6, 8, 10, 12, 14 or 16 (16 = full time series)
        y <- 2009:2024

        if( years == 6 ) {
            message(paste("Subsetting ",years," of daily data", sep = ""))
            combin_df <- combin_df[combin_df$Year %in% c(2010,2013,2016,2019,2022,2024),]
        } else if( years == 8 ) {
            message(paste("Subsetting ",years," of daily data", sep = ""))
            combin_df <- combin_df[combin_df$Year %in% y[seq(1, length(y), length.out = 8)],]
        } else if( years == 10 ) {
            message(paste("Subsetting ",years," of daily data", sep = ""))
            combin_df <- combin_df[combin_df$Year %in% c(2009,2010,2011,2012,2014,2016,2018,2020,2022,2024),]
        } else if( years == 12 ) {
            message(paste("Subsetting ",years," of daily data", sep = ""))
            combin_df <- combin_df[combin_df$Year %in% y[seq(1, length(y), length.out = 12)],]
        } else if( years == 14 ) {
            message(paste("Subsetting ",years," of daily data", sep = ""))
            combin_df <- combin_df[combin_df$Year %in% y[seq(1, length(y), length.out = 14)],]
        } else if( years == 15 ) {
            message(paste("Subsetting ",years," of daily data", sep = ""))
            combin_df <- combin_df[combin_df$Year %in% y[seq(1, length(y), length.out = 15)],]
        } else {
            message(paste("Using all 16 years of daily data", sep = ""))
        } # eo if else loop to subset years

        # Make some plots to examine distrbution of Offsets and stand age (the latter should not change across dataset)
        pdf(NULL) # https://stackoverflow.com/questions/6535927/how-do-i-prevent-rplots-pdf-from-being-generated 

        p1 <- gghistogram(
            combin_df, x = "StandAge", y = "..density..",
            add = "mean", rug = FALSE, fill = "TreeType",
            palette = c("#238443","#d9f0a3"),
            add_density = TRUE, bins = 50
        ) + facet_wrap(.~factor(Region))

        p2 <- gghistogram(
            combin_df, x = "Offset", y = "..density..",
            add = "mean", rug = FALSE, fill = "TreeType",
            palette = c("#238443","#d9f0a3"),
            add_density = TRUE, bins = 50
        ) + facet_wrap(.~factor(Region))

        panel <- ggarrange(p1,p2, align = "hv", ncol = 1, nrow = 2)
        ggsave(plot = panel, filename = paste("panel_distribution_subset_",stat,"_",var,"_regions_03.02.26.jpg", sep = ""), dpi = 300, height = 5, width = 5.5)

        rm(panel,p1,p2)
        gc()

        # Perform GAMM & save
        message(paste("Running the GAMM for ",stat," ",var, sep = ""))

        ### "To ensure computational feasibility of the GAMM while retaining temporal representativeness,
        ### the model was fitted on a stratified subset of six years (2010, 2013, 2016, 2019, 2022, and 2024).
        ### These years were selected to evenly span the study period and encompass the full range of climatic
        ### variability observed in the dataset. This subsampling approach preserved the seasonal signal and
        ### interannual variation necessary for estimating the temporal correlation structure, while reducing data
        #### volume to a level compatible with the AR(1) correlation model implemented in gamm()"

        ### gamm() is a wrapper around lme() from nlme, and that combination is not designed for datasets with > hundreds of thousands of rows.
        ### The crash after a few minutes is very likely an Out of Memory (OOM) kill (not a coding error).

        ### 19/11/25: Modify this part to include Anom or Mean_Grassland_Value as a predictor in the GAMM depending on variable modelled

        ### 15/10/25: This was a controlled test to confirm that your crashes are indeed caused by the corAR1 correlation structure
        ### (and not by the random effect or data size itself).
        ## Take a manageable subset of the combin_df data (e.g., 10,000–20,000 rows, still with multiple EPs).
        ## Ensures you still have several EPs and multiple time points per EP for correlation estimation.
        # test_df <- combin_df |>
        #    dplyr::group_by(EP) |>
        #    dplyr::slice_sample(n = 300, replace = FALSE) |>
        #    dplyr::ungroup()
        # str(test_df) ; dim(test_df)

        ## Fit Model A - with correlation structure
        #mod_A <- gamm(
        #    Offset ~ s(StandAge, k = 8, fx = TRUE) +
        #    s(DOY, bs = "cc", by = TreeType, k = 12) +
        #    TreeType + Region,
        #    random = list(EP = ~1),
        #    correlation = corAR1(form = ~ 1 | EP),
        #    data = test_df,
        #    control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
        #)

        ## Fit Model B - without correlation structure
        #mod_B <- gamm(
        #    Offset ~ s(StandAge, k = 8, fx = TRUE) +
        #    s(DOY, bs = "cc", by = TreeType, k = 12) +
        #    TreeType + Region,
        #    random = list(EP = ~1),
        #    data = test_df,
        #    control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100)
        #)

        ## Once both models converge, compare their structure and performance:
        # summary(mod_A$lme)
        # summary(mod_B$lme)
        # summary(mod_A$gam) # R2 = 0.38
        # summary(mod_B$gam) # R2 = 0.383
        # r.squaredGLMM(mod_A$lme)
        #   R2m  R2c
        # 0.065 0.665
        # r.squaredGLMM(mod_B$lme)
        #   R2m  R2c
        # 0.062 0.631
        # AIC(mod_A$lme, mod_B$lme) # A has lower AIC by 700 points

        ## You can quantify how much the corAR1 actually helps using the autocorrelation of residuals:
        # pdf(file = "/home/fbenedetti/acf_modeB.pdf", width = 5, height = 5)
        # acf(residuals(mod_B$lme, type = "normalized"))
        # dev.off()
        # If the residual autocorrelation is small (lags drop quickly), you can safely drop the AR1 term without much loss in model validity
        # Strong autocorrelation at lag 1 and persistent correlation beyond lag 1 -> residuals are not independent
        # -> model violates one of the key assumptions of standard GAMs/GAMMs -> Model B is underfitting the temporal dependence
        # Large autocorrelation only at lag 1 for Model A, but then it drops below significance for higher lags
        # -> AR1 successfully captured the main serial correlation structure
        # -> Model A correctly models temporal dependence; the AR(1) structure is statistically justified and materially improves the fit.
        ## Therefore, we have a theoretically necessary structure that is also computationally too heavy in its current form...

        ### Train the GAMM
        message(paste("Training the GAMM", sep = ""))
        # Drop unused vars and factor levels
        combin_df <- droplevels(combin_df)
        vars_used <- c("Region","EP","Offset","Mean_Grassland_Value","Anom","StandAge","DOY","TreeType")
        combin_df <- combin_df %>% select(all_of(vars_used))
        combin_df <- combin_df %>% drop_na(all_of(vars_used))
        # summary(combin_df) ; dim(combin_df)

        ### 03/02/26: Because forest structural effects on microclimate are expected to stabilize after canopy closure,
        ### and because GAMMs were trained on stands aged 12–75 years only, we capped StandAge at 40 years for offsets modelling
        ### and prediction. This prevents extrapolation beyond the training support and enforces ecologically
        ### realistic saturation of stand-age effects -> Should allow to avoid any long)term warming drift

        combin_df$StandAge_cap <- pmin(combin_df$StandAge, cap) # literature-based subjective choice, could be 30 or 45 yo

        gamm_model <- gamm(
                Offset ~ s(StandAge_cap, k = 4, bs = "ts") +  # bs = "ts" means “Let age matter if and where the data support it, otherwise fade it out.”
                s(DOY, bs = "cc", by = TreeType, k = 8, fx = TRUE) +
                s(Anom, k = 5, fx = TRUE, bs = "tp") +
                TreeType + Region,
                random = list(EP = ~1),
                correlation = corAR1(form = ~ 1|EP),
                data = combin_df
        ) # eo gamm

        ### 03/02/26: Test the effect of capping StandAge
        # summary(gamm_model$gam)
        # require("gratia")
        # p <- draw(gamm_model$gam, select = "s(StandAge_cap)", rug = TRUE)
        # ggsave(plot = p, filename = "test_smooth_term_StandAge_cap_min_Ta_10_gratia.jpg", dpi = 300, width = 4, height = 4)
        
        ### Plot potential prediction until 100 yo? 
        # age_seq <- seq(0,100, by = 1)
        # newdat <- tibble(
        #    StandAge_cap = age_seq,
        #    Anom = 0,
        #    DOY = 180,
        #    TreeType = levels(combin_df$TreeType)[1],
        #    Region = levels(combin_df$Region)[1]
        # )
        # pred <- predict(gamm_model$gam, newdata = newdat, type = "terms", terms = "s(StandAge_cap)", se.fit = TRUE)
        # plot_df <- newdat %>%
        #    mutate(
        #        fit = as.numeric(pred$fit),
        #        se = as.numeric(pred$se.fit),
        #        upper = fit + 2 * se,
        #        lower = fit - 2 * se
        #    ) # eo mutate
        # plot <- ggplot(plot_df, aes(x = StandAge_cap, y = fit)) +
        #    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.4) +
        #    geom_line(linewidth = 1) + geom_vline(xintercept = 40, linetype = "dashed", colour = "red", linewidth = 0.8) +
        #    labs(x = "Stand age (years)", y = "Partial effect on Offset", title = "Effect of stand age on offsets",
        #        subtitle = "Effect increases during forest development and stabilises after canopy closure") +
        #    theme_minimal(base_size = 13)
        # Save
        #ggsave(plot = plot, filename = "test_smooth_term_StandAge_cap_min_Ta_10.jpg", dpi = 300, width = 6, height = 4.5)

        ### Save GAMM object as .RDS
        message(paste("Saving the GAMM for ",stat," ",var, sep = ""))
        saveRDS(gamm_model, file = paste("model_GAMM_offsets_",stat,"_",var,"_subset_",years,"yrs_03.02.26.rds", sep = ""))
        
        # Clean and go next
        rm(combin_df,gamm_model); gc()

} # eo FUN - model_offset_gamm

### Apply model_offset_gamm() in for loops: c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")
for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")) {
    for(s in c("max","min")) {
        gamm_offset(var = v, stat = s, age = 75, years = 15, cap = 40)
    } # eo 2nd for loop
} # eo 1st for loop

### 20/10/25: Was re-run on the 15/10/25 with models of reduced complexity and years 

### 19-20/11/25: Was re-run on the 19/11/25 with models including Anomalies to DOY or Mean_Grassland_Value for SM_10
### based on the results of the tests performed in the 4.4.2 R scripts - Try with 14 years of data

### 29/01/26: Re-running the SM_10 models with Anoms instead of mean grassland values
# gamm_offset(var = "SM_10", stat = "max", age = 75, years = 14)
# gamm_offset(var = "SM_10", stat = "min", age = 75, years = 14)

### 03/02/36: Re-running with capped s(StandAge) for all variables -> Avoid long-term drift due to extrapolation beyond training range


### ------------------------------------------------------------------------------------------------------------

### 20/10/25: Analyzing outputs from new runs of gamm_offset() ran on the 15/10/25 before the group retreat in Kiental

### 21/10/25: Analyzing outputs from new runs of gamm_offset() ran on the 20/10/25 (same subset of 6 years but with fx = T and without precip.)

### 22/10/25: Analyzing outputs from new runs of gamm_offset() ran on the 21-22/10/25 (subset of 10 years but with fx = T and without precip.)

### 23/10/25: Analyzing outputs from new runs of gamm_offset() ran on the 21-22/10/25 (subset of 10 years but with fx = T and without precip.)

### 19-20/11/25: Need to analyze the outputs from new runs of gamm_offset() ran on the 19/11/25
### (15 years and refined formula to include grass land values and anomalies)

### 29/01/26: Analyzing outputs of the corrected SM_10 GAMMs

### 06/02/26: Analyzing outputs of the new GAMMs based on new s(StandAge, k = 4) with clipping ages > 40 yo


## Go to dir on local machine where GAMM outputs are stored
# setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/BEO-UniBern/Data/Climate data/offsets")
# dir on the climcal server
setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data/models/"); dir()
## When gamm_offset() is done running, there should be 12 model objects to analyze (6 variables x 2 stats)

### Re-use content of R Script#4.4.2 to analyze outputs (plots, stats etc.)

### 1) Ta_10
## 1.A) max
mod <- readRDS("model_GAMM_offsets_max_Ta_10_subset_15yrs_03.02.26.rds") # class(mod); str(mod)
summary(mod$gam) # R-sq.(adj) = 0.594  ### fixed effects & smooths are capturing a lot of structure - same as GAMM runs from August!
# 'RegionSWA' non signif. 
# summary(mod$lme)
### 21/10: With fx = FALSE: R-sq.(adj) = 0.594 (same as before)
### 22/10: With fx = TRUE and 10 years: 0.585

### 03/02/26: With s(StandAge, k = 4, fx = F): 0.65 ! Best model so far! All terms signif.

# gam.check
gam.check(mod$gam) # gut - k should be too low for s(StandAge)
# 'TreeTypeConiferous' & 'TreeTypeDeciduous' too low in terms of 'k'
### BEWARE outputs of gam.check may vary: gam.check() is not a pure “read‐only” diagnostic,
### it actually refits parts of your model internally. It re‐evaluates your smooth terms with
### k‐index tests to check if your chosen basis size k was sufficient.
### To do that, it resamples residuals (via simulation) and sometimes uses randomization for the p‐values in those checks.
### It also recomputes EDF, GCV scores, and other fit summaries during this check.

### 03/02/26: With s(StandAge, k = 4, fx = F): s(StandAge_cap) looks good in terms of k choice :) 

## concurvity() checks nonlinear functional dependence between smooth terms
# 0.0 → no concurvity (completely independent from other terms)
# 1.0 → perfect concurvity (term is an exact function of other terms — complete redundancy)
# < 0.3 → not a problem
# 0.3–0.6 → moderate; watch for inflated SEs, unstable shapes
# > 0.8 → serious concurvity; one term can be almost entirely predicted from others, making estimates unstable
concurvity(mod$gam, full = TRUE)
## The fixed-effect part ('TreeType' & 'Region') is 91% predictable from other terms in the model.
## For 's(StandAge)', concurvity is 0.34-0.39. This is moderate, not terrible!
### 03/02/26: With s(StandAge, k = 4, fx = F): No concurvity issues with 'StandAge_cap'

# Adjusted R2
r.squaredGLMM(mod$lme)
#  R2m   R2c
# 0.61  0.66
### 20/10/25 -> Much higher R2m as before! :) Gut

### 21/10: With fx = FALSE 
#  R2m       R2c
# 0.072     0.650
### -> fx = T really is the cause for changing the R2m

### 22/10: With fx = TRUE and 10 years: 
#  R2m       R2c
# 0.596     0.632

### 03/02/26: With s(StandAge, k = 4, fx = F):
#  R2m      R2c
# 0.578    0.752
### R2c increased even furtgher, good.

# When & why is the R2m low?
# Because most of your “signal” might be in:
# - Nonlinear smooth terms (GAM part) - that is already 60% of the deviance explained
# - Or random effects (LME part)
# - Or the autocorrelation structure

### -> Here, our model is learning a lot from the smooth terms (GAM), from the grouping structure (R2c = 0.66) 
### (random effects in the LME), AND from the fixed effects (R2m = 0.61)

### 08/25: 
### Combination of low R2m and high R2c is very common in GAMMs: the smooths do the heavy lifting for structured variation,
### and the random intercepts soak up group‐level shifts. A low R2m reflects that the LME’s fixed effects aren’t where the main fit power is.

summary(mod$lme)

### Plots
library("gratia")
pdf(NULL) 
smooths.plot <- draw(mod$gam)
resids.plot <- appraise(mod$gam, line_col = "#238443", point_alpha = .01)
# class(smooths.plot) ; class(resids.plot)
ggsave(plot = smooths.plot, filename = "plot_smooths_GAMM_max_Ta_10_06.02.26.jpg", dpi = 300, height = 7.5, width = 7.5)
ggsave(plot = resids.plot, filename = "plot_residuals_GAMM_max_Ta_10_06.02.26.jpg", dpi = 300, height = 7, width = 7.5)

## To plot fixed effects estimates (± standard errors)
coefs <- tidy(mod$lme, effects = "fixed")
# data.frame(coefs)
coeffs.plot <- ggplot(coefs, aes(x = estimate, y = reorder(term, estimate))) + geom_point() +
    geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
    labs(x = "Estimate", y = "Fixed Effect", title = "Fixed Effects (GAMM - lme)") +
    theme_minimal()
ggsave(plot = coeffs.plot, filename = "plot_coeffs_GAMM_max_Ta_10_06.02.26.jpg", dpi = 300, height = 5, width = 5)

### Smooths check out. Similar to daily max Ta_200. This model seems good :)

rm(coefs,coeffs.plot,resids.plot,smooths.plot,mod); gc()

### -------------------------------------------------

## 1.B) min Ta_10
mod <- readRDS("model_GAMM_offsets_min_Ta_10_subset_15yrs_03.02.26.rds")
summary(mod$gam) # R-sq.(adj) = 0.057   ### not great based on this...even slightly worse than before
# 'RegionSCH' & 'RegionSWA' non signif.
# summary(mod$lme)

### 21/10: With fx = FALSE -> 0.057 (same as before)
### 22/10: With fx = TRUE and 10 years: 0.084 (best so far, but still very weak model)
### 03/02/26: With s(StandAge, k = 4, fx = F): 0.145, sill quite weak but best so far! 

# gam.check
gam.check(mod$gam)
# looks good

# To report Phi
summary(mod$lme)

# Adjusted R2
r.squaredGLMM(mod$lme)
#  R2m   R2c
# 0.092  0.20 (used to be 0.69)
### 20/10 -> Great reduction in 'R2c'...
### GAMM from 07/25 was better!

### 21/10 With fx = FALSE:
#  R2m   R2c
# 0.04  0.311
# R2c increased again but not as much as inital GAMM from 07/25

### 22/10: With fx = TRUE and 10 years:
#  R2m   R2c
# 0.127 0.226
### GAMM from 07/25 was better with a R2c of 0.69! 
# Check again
old_mod <- get(load("model_GAMM_offsets_min_Ta_10_31.07.25.RData"))
r.squaredGLMM(old_mod$lme)
#   R2m  R2c
# 0.042 0.690
rm(old_mod,mod_gamm); gc()

### 03/02/26: With s(StandAge, k = 4, fx = F):
#  R2m     R2c
# 0.286   0.426

### CCL: R2c decreased with new capped StandAge smooth but R2m (part explained by fixed effects) improved massively, good.

### Again, most variance explanation in the LME comes from random intercepts
### (differences between EPs/Regions) - not the linear fixed effects. 
### This model is NOT learning a lot from the smooth terms (GAM) but it is from the grouping structure
### (random effects in the LME). Again, the linear fixed effects in the LME aren’t explaining much by themselves.
### -> This model finds most of its power in the random effects (EP) and autocorr structure 

## Plots
library("gratia")
pdf(NULL) 
smooths.plot <- draw(mod$gam)
resids.plot <- appraise(mod$gam, line_col = "#238443", point_alpha = .01)
# class(smooths.plot) ; class(resids.plot)
ggsave(plot = smooths.plot, filename = "plot_smooths_GAMM_min_Ta_10_06.06.26.jpg", dpi = 300, height = 7.5, width = 7.5)
ggsave(plot = resids.plot, filename = "plot_residuals_GAMM_min_Ta_10_06.06.26.jpg", dpi = 300, height = 7, width = 7.5)

### 03/02/26: With s(StandAge, k = 4, fx = F): Resdiual plots still look OK though :)

## To plot fixed effects estimates (± standard errors)
coefs <- tidy(mod$lme, effects = "fixed")
# data.frame(coefs)
coeffs.plot <- ggplot(coefs, aes(x = estimate, y = reorder(term, estimate))) + geom_point() +
    geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
    labs(x = "Estimate", y = "Fixed Effect", title = "Fixed Effects (GAMM - lme)") +
    theme_minimal()
ggsave(plot = coeffs.plot, filename = "plot_coeffs_GAMM_min_Ta_10_20.10.25.jpg", dpi = 300, height = 5, width = 5)

### Stand age smooth is somehow reversed compared to max Ta_10? Why?

### Potential explanations:
### - Canopy insulation reduces nocturnal radiative cooling: At night, canopies block longwave sky radiation
###   and reduce wind mixing, so older forests often have warmer minima than open grassland
###   → ΔT_min = Tgrass_min − Tforest_min will be smaller or even negative with age
###   (i.e. min temp. in forests warmer relative to grasslands). So the sign/shape may be opposite.

### - Soil moisture & evapotranspiration: In daytime, evapotranspiration in mature forests may further cool surfaces;
###   at night, moisture effects are less important for minima vs radiative trapping.

### -> "Stand age increases daytime buffering (larger cooling of maxima, see max Ta_10 above),
###    but increases nighttime warming (reduces negative offset of minima), consistent with
###    canopy shading reducing daytime heating while canopy-mediated longwave trapping and
###    reduced mixing raise night-time surface temperatures"

### HOWEVER, could also be due to concurvity (collinearity) in the model - not an ecological reason then: 
## -> check concurvity and variable dependence

## concurvity() checks nonlinear functional dependence between smooth terms
# 0.0 → no concurvity (completely independent from other terms)
# 1.0 → perfect concurvity (term is an exact function of other terms — complete redundancy)
# < 0.3 → not a problem
# 0.3–0.6 → moderate; watch for inflated SEs, unstable shapes
# > 0.8 → serious concurvity; one term can be almost entirely predicted from others, making estimates unstable
concurvity(mod$gam, full = TRUE)

## The fixed-effect part ('TreeType' & 'region') is 91% predictable from other terms in the model.
## For 's(StandAge)', concurvity is 0.30. This is moderate, not terrible! 

### Given that your parametric part has 0.91 concurvity, we can suspect TreeType is strongly tied to StandAge in the data
### (as evidenced before). That means our model is trying to separate effects that, in reality, barely vary independently...

### But remember GAM explains almost no variance anyway..
### The smooth terms might be flexible but still not capturing strong trends -> the relationship is weak, noisy, or dominated by random variation
### Most of the predictive power actually comes from random effects, not from the smooths or fixed effects.
### -> This is consistent with ecological or environmental data where site-level variation dominates.

rm(coefs,coeffs.plot,resids.plot,smooths.plot,mod); gc()

### Same smooths as before! So just worse models for min Ta_10 in 10/25 compared to 07/10 with all years as input

### -------------------------------------------------

### 2) Ta_200
## 2.A) max
mod <- readRDS("model_GAMM_offsets_max_Ta_200_subset_15yrs_03.02.26.rds")
summary(mod$gam) # R-sq.(adj) = 0.39 - good! Same with new GAMM
### all terms signif.
# summary(mod$lme)

### 21/10: With fx = FALSE -> R-sq.(adj) = 0.39 (same)
### 22/10: With fx = TRUE and 10 years: 0.385 (same)
### 03/02/26: With s(StandAge, k = 4, fx = F): 0.418 (best so far!)

# gam.check
gam.check(mod$gam) # All good here

# concurvity check
concurvity(mod$gam, full = TRUE)
## same observations as the models above

# To report Phi in Table S3
summary(mod$lme) 

# Adjusted R2
r.squaredGLMM(mod$lme)
#       R2m             R2c
# 0.40 (was 0.32)  0.45 (was 0.83)
# Increase in fixed effects but decrease in local factors
# Still good :)

## 21/10: With fx = FALSE
#  R2m   R2c
# 0.12  0.683
### R2m even smaller than before, but R2c increased back slightly (so fx = TRUE decreases R2c the most, but increases R2m)

### 22/10: With fx = TRUE and 10 years:
#  R2m   R2c
# 0.396 0.444
### Similar to version from 15/10

### 03/02/26: With s(StandAge, k = 4, fx = F):
#  R2m   R2c
# 0.104 0.877
### Lowest R2m ever but highest, most variance has been explained by random effects when switching the StandAge smooth

## Plots
library("gratia")
pdf(NULL) 
smooths.plot <- draw(mod$gam)
resids.plot <- appraise(mod$gam, line_col = "#238443", point_alpha = .01)
# class(smooths.plot) ; class(resids.plot)
ggsave(plot = smooths.plot, filename = "plot_smooths_GAMM_max_Ta_200_06.02.26.jpg", dpi = 300, height = 7.5, width = 7.5)
ggsave(plot = resids.plot, filename = "plot_residuals_GAMM_max_Ta_200_06.02.26.jpg", dpi = 300, height = 7, width = 7.5)

## To plot fixed effects estimates (± standard errors)
#coefs <- tidy(mod$lme, effects = "fixed")
# data.frame(coefs)
#coeffs.plot <- ggplot(coefs, aes(x = estimate, y = reorder(term, estimate))) + geom_point() +
#    geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
#    labs(x = "Estimate", y = "Fixed Effect", title = "Fixed Effects (GAMM - lme)") +
#    theme_minimal()
#ggsave(plot = coeffs.plot, filename = "plot_coeffs_GAMM_max_Ta_200_20.10.25.jpg", dpi = 300, height = 5, width = 5)

### Decent model overall. Makes sense - similar to max Ta_10 too.

rm(coefs,resids.plot,smooths.plot,mod); gc()

### -------------------------------------------------

## 2.B) min Ta_200
mod <- readRDS("model_GAMM_offsets_min_Ta_200_subset_15yrs_03.02.26.rds")
summary(mod$gam) # R-sq.(adj) = 0.14 (was 0.136) - same as before
### RegionSCH non signif.
# summary(mod$lme)

## 21/10: With fx = FALSE -> 0.142 (slight increase with fx = F compared to fx = T)
### 22/10: With fx = TRUE and 10 years: 0.14 (as always)
### 03/02/26: With s(StandAge, k = 4, fx = F): 0.187 (best so far again!)

# gam.check
gam.check(mod$gam) # k too small overall

# concurvity
concurvity(mod$gam, TRUE)
# same as always: concurvity issue in the parametric part

# To report 'Phi'
summary(mod$lme)

# Adjusted R2
r.squaredGLMM(mod$lme)
#   R2m            R2c
#  0.12 (0.24)  0.26 (0.55)
### 20/10 -> Older GAMM was better - But smooth plots look the same honestly

## 21/10: With fx = FALSE:
#  R2m       R2c
# 0.104     0.305
### -> decrease in R2m but intermediary R2c (model from 07/25 was better)

### 22/10: With fx = TRUE and 10 years:
#   R2m     R2c
#  0.117   0.276
### OLDER GAMM FROM 07/10 WAS BETTER

### 03/02/26: With s(StandAge, k = 4, fx = F):
#  R2m      R2c
# 0.147    0.632
# Increase in both R2m and R2c, actually much higher R2c which makes sense for min Ta_200?

### In the min air temperature GAMM, R²m = 0.24 means the fixed effects explain ~24% of the variation (way more than the 8% earlier).
### Meanwhile, R²c = 0.55 means the full model with random effects explains ~55% — so site-level variation is still important,
### but now the smooths and factors are actually pulling real weight. It suggests that the relationship between StandAge
### and offsets is stronger (and more consistent) for air temperatures than for skin temperatures,
### making the sign reversal less likely to be noise

## Plots
library("gratia")
pdf(NULL) 
smooths.plot <- draw(mod$gam)
resids.plot <- appraise(mod$gam, line_col = "#238443", point_alpha = .01)
# class(smooths.plot) ; class(resids.plot)
ggsave(plot = smooths.plot, filename = "plot_smooths_GAMM_min_Ta_200_06.02.26.jpg", dpi = 300, height = 7.5, width = 7.5)
ggsave(plot = resids.plot, filename = "plot_residuals_GAMM_min_Ta_200_06.02.26.jpg", dpi = 300, height = 7, width = 7.5)

## To plot fixed effects estimates (± standard errors)
coefs <- tidy(mod$lme, effects = "fixed")
# data.frame(coefs)
coeffs.plot <- ggplot(coefs, aes(x = estimate, y = reorder(term, estimate))) + geom_point() +
    geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
    labs(x = "Estimate", y = "Fixed Effect", title = "Fixed Effects (GAMM - lme)") +
    theme_minimal()
ggsave(plot = coeffs.plot, filename = "plot_coeffs_GAMM_min_Ta_200_20.10.25.jpg", dpi = 300, height = 5, width = 5)

### Reversed effect of stand age on min temp. offsets again! Let's unravel this as it makes lots of sense :) 

## Offset < 0 → forest is warmer than grassland (forest retains/keeps more heat than grassland)
## Offset > 0 → forest is cooler than grassland (forest buffers daytime maxima, etc.)

## If your s(StandAge) smooth decreases (more negative) with age for nighttime minima, that means that 
## as stands get older they show larger negative offsets → i.e. older stands retain more heat at night relative to grasslands
## Matches physics: denser canopy / leaf-on conditions (especially for deciduous stands in summer) increase nighttime retention
## → more negative offsets! 

## For daytime maxima you saw the opposite shape — that is consistent too:
## Older stands shade more so forests are cooler in daytime (offset more positivewith time),
## while at night they retain more heat (offset more negative with time).

rm(coefs,coeffs.plot,resids.plot,smooths.plot,mod); gc()

### -------------------------------------------------

### 3) Ts_05
## 3.A) Ts_05 max
mod <- readRDS("model_GAMM_offsets_max_Ts_05_subset_15yrs_03.02.26.rds")
summary(mod$gam) # R-sq.(adj) = 0.764 (old GAMM: 0.75) ## Still an EXCELLENT fit
### 'TreeTypeDeciduous' & 'RegionSWA' non signif.
# summary(mod$lme)

### 21/10: With fx = FALSE -> 0.76 (same as usual, very good model)
### 22/10: With fx = TRUE and 10 years: 0.731
### 03/02/26: With s(StandAge, k = 4): 0.716 - slightly lower but still quite good

# gam.check
gam.check(mod$gam)
# s(StandAge) likely too low. But deviance explained is 75%. Let's see how the smooths look like later.

# concurvity
concurvity(mod$gam, TRUE)
# Same as usual. Concurvity issue in the parametric part

summary(mod$lme)

# Adjusted R2
r.squaredGLMM(mod$lme)
#  R2m           R2c
# 0.758 (0.01)  0.81 (0.953)
### 20/10 -> Great improvement in R2m! Good. Slight decrease in R2c but that is OK
### 20/10 -> new model is a better GAMM overall

### 21/10: With fx = FALSE:
#  R2m   R2c
# 0.044 0.738
### Huge decrease again in R2m and still decrease in R2c
### -> fx = FALSE discarded the variance attributed to the fixed effects - but it also somehow decreased total variance (R2c)

### 22/10: With fx = TRUE and 10 years: 
#  R2m   R2c
# 0.729 0.786
### This model and the one from 15/10 are the best :)

### 03/02/26: With s(StandAge, k = 4):
#  R2m      R2c
# 0.637    0.872
### Still very good

### -> Again, the GAMM is learning a lot from the smooth terms (GAM) and from the grouping structure
### (random effects in the LME), but the strictly linear fixed effects in the LME aren’t explaining much by themselves.
### That combination is very common in GAMMs: the smooths do the heavy lifting for structured variation,
### and the random intercepts soak up group‐level shifts.
### The low 𝑅𝑚2 just reflects that the LME’s fixed effects aren’t where the main fit power is.

## Plots
library("gratia")
pdf(NULL) 
smooths.plot <- draw(mod$gam)
resids.plot <- appraise(mod$gam, line_col = "#238443", point_alpha = .01)
# class(smooths.plot) ; class(resids.plot)
ggsave(plot = smooths.plot, filename = "plot_smooths_GAMM_max_Ts_05_06.02.26.jpg", dpi = 300, height = 7.5, width = 7.5)
ggsave(plot = resids.plot, filename = "plot_residuals_GAMM_max_Ts_05_06.02.26.jpg", dpi = 300, height = 7, width = 7.5)

## To plot fixed effects estimates (± standard errors)
#coefs <- tidy(mod$lme, effects = "fixed")
## data.frame(coefs)
#coeffs.plot <- ggplot(coefs, aes(x = estimate, y = reorder(term, estimate))) + geom_point() +
#    geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
#    labs(x = "Estimate", y = "Fixed Effect", title = "Fixed Effects (GAMM - lme)") +
#    theme_minimal()
#ggsave(plot = coeffs.plot, filename = "plot_coeffs_GAMM_max_Ts_05_20.10.25.jpg", dpi = 300, height = 5, width = 5)

### Interesting, the residuals look very good but most of the information must come from the DOY
### and DOY:TreeType terms (judging by their confidence intervals)!
### The s(StandAge) smooth has higher confidence intervals and is much more variable. 

### Considering the high % of deviance explained as well the VERY high R2c above, 
### it looks like offsets in max Ts_05 can mainly be predicted from local factors (EP)
### and DOY and its interactions with TreeType!

rm(resids.plot,smooths.plot,mod); gc()

### 20/10 -> Although new GAMM seems better overall, the s(precipitation) does not seem to make sense to me 
### as offsets should increase with higher precipitation
### Try to understand why R2m is increased because the increase in the GAM part is not that high (1/4% of deviance explained added by precip.)

### -------------------------------------------------

## 3.B) Ts_05 min
mod <- readRDS("model_GAMM_offsets_min_Ts_05_subset_15yrs_03.02.26.rds")
summary(mod$gam) # R-sq.(adj) = 0.69 (0.64), very good again
### 'TreeTypeDeciduous' & 'RegionSWA' non signif.
# summary(mod$lme)

### 21/10: With fx = FALSE -> 0.683 (same as usual)
### 22/10: With fx = TRUE and 10 years: 0.613, lower but still good
### 03/02/26: With s(StandAge, k = 4): 0.525, quite strong decrease (by 0.1 unit) in R2 of the GAMM. Due to simplification of the s(StandAge)

# gam.check
gam.check(mod$gam)
# s(StandAge) likely too low.

# concurvity check
concurvity(mod$gam, full = TRUE)
## same as usual

# To report 'Phi' for Table S3
summary(mod$lme)

# Adjusted R2
r.squaredGLMM(mod$lme)
#      R2m          R2c
# 0.681 (0.052) 0.81 (0.937)
### Like for max Ts_05 -> great improvement in R2m and sligh decrease in R2c

## 21/10: With fx = FALSE
#  R2m   R2c
# 0.054 0.741
### -> Same as max Ts_05 -> variance attributed to fixed effects returned to 0, decreased R2c
### -> Previous model from 15/10 actually seemed better

### 22/10: With fx = TRUE and 10 years:
#  R2m   R2c
# 0.635 0.776
### This model and the one from 15/07 are the best

### 06/02/26: With s(StandAge, k = 4):
#  R2m       R2c
# 0.416     0.862
### -> Decrease in R2m too but increase in R2c, EP soaked up more variance

## Plots
library("gratia")
pdf(NULL) 
smooths.plot <- draw(mod$gam)
resids.plot <- appraise(mod$gam, line_col = "#238443", point_alpha = .01)
# class(smooths.plot) ; class(resids.plot)
ggsave(plot = smooths.plot, filename = "plot_smooths_GAMM_min_Ts_05_06.02.26.jpg", dpi = 300, height = 7.5, width = 7.5)
ggsave(plot = resids.plot, filename = "plot_residuals_GAMM_min_Ts_05_06.02.26.jpg", dpi = 300, height = 7, width = 7.5)

## To plot fixed effects estimates (± standard errors)
#coefs <- tidy(mod$lme, effects = "fixed")
## data.frame(coefs)
#coeffs.plot <- ggplot(coefs, aes(x = estimate, y = reorder(term, estimate))) + geom_point() +
#    geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
#    labs(x = "Estimate", y = "Fixed Effect", title = "Fixed Effects (GAMM - lme)") +
#    theme_minimal()
#ggsave(plot = coeffs.plot, filename = "plot_coeffs_GAMM_min_Ts_05_20.10.25.jpg", dpi = 300, height = 5, width = 5)

### Very similar to max Ts_05. DOY shows a strong signal, but s(StandAge) also shows a strong, non linear,
### decreasing trend which makes sense based on what we observed before for air temperatures! 
### R2 is slightly lower and DOY have slightly bigger CI and the s(StandAge) pattern is clearer. 

### Maybe the effect of stand age on min Ts_05 is a bit stronger than for max Ts_05?

rm(resids.plot,smooths.plot,mod); gc()

### 20/10 -> s(precipitation) does not make much sense -> remove! 
### -> Try to understand origin of changes in R2m 0.681 (0.052) & R2c 0.81 (0.937)

### -------------------------------------------------

### 4) Ts_10
## 4.A) max Ts_10
mod <- readRDS("model_GAMM_offsets_max_Ts_10_subset_15yrs_03.02.26.rds")
summary(mod$gam) # R-sq.(adj) = 0.77 (0.751) - Excellent fit again!
### 'TreeTypeDeciduous' non signif. though
# summary(mod$lme)

### 21/10: With fx = FALSE -> 0.77
### 22/10: With fx = TRUE and 10 years: 0.732
### 06/02/26: With s(StandAge, k = 4): 0.712, slightly lower but still very good fit :)

# gam.check
gam.check(mod$gam)
# like for Ts_05 - k might be too low for s(StandAge)

# concurvity check
concurvity(mod$gam, full = TRUE)
# same as usual

# To report 'Phi'
summary(mod$lme)

# Adjusted R2
r.squaredGLMM(mod$lme)
#     R2m           R2c
# 0.76 (0.011)  0.82 (0.940)
# Same as max/min Ts_05
### 20/10 -> Same as other soil tempertaure variables - Strong increase in R2m & slight decrease in R2c
### Need to understand why.

## 21/10: With fx = FALSE
#  R2m   R2c
# 0.044 0.765
### -> Again, discarded back the variance attributed to fixed effects, and further decreased R2c (total variance)

### 22/10: With fx = TRUE and 10 years: 
#  R2m   R2c
# 0.729 0.798
### This model and the one from 15/10 are the bets again :)

### 06/02/26: With s(StandAge, k = 4):
#  R2m    R2c
# 0.528  0.886
### -> Lower R2m again, but increased R2c


## Plots
library("gratia")
pdf(NULL) 
smooths.plot <- draw(mod$gam)
pdf(NULL) 
resids.plot <- appraise(mod$gam, line_col = "#238443", point_alpha = .01)
# class(smooths.plot) ; class(resids.plot)
ggsave(plot = smooths.plot, filename = "plot_smooths_GAMM_max_Ts_10_06.02.26.jpg", dpi = 300, height = 7.5, width = 7.5)
ggsave(plot = resids.plot, filename = "plot_residuals_GAMM_max_Ts_10_06.02.26.jpg", dpi = 300, height = 7, width = 7.5)

## To plot fixed effects estimates (± standard errors)
#coefs <- tidy(mod$lme, effects = "fixed")
#coeffs.plot <- ggplot(coefs, aes(x = estimate, y = reorder(term, estimate))) + geom_point() +
#    geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
#    labs(x = "Estimate", y = "Fixed Effect", title = "Fixed Effects (GAMM - lme)") +
#    theme_minimal()
#ggsave(plot = coeffs.plot, filename = "plot_coeffs_GAMM_max_Ts_10_20.10.25.jpg", dpi = 300, height = 5, width = 5)

### Interesting again, very similar to max Ts_05 (as expected). DOY and EP have the strongest effects.
### s(StandAge) has a negative effect and not a positive saturating effect like Ta_10 & Ta_200. 
### This means older stands tend to retain heat better, relative to grasslands, no matter whether
### it is day or night time. 
### On average, offsets of daily max and min Ts_10 are positive: grasslands have warmer Ts than forests 
### at night and during the day. With time, it looks like as if older stands tend towards lower offsets 
### meaning that as the forests grows older, it's better able to retain the heat, probably because it 
### better retaind the heat at night. This is supported by the fact that the negative effect of stand age
### seems STRONGER for min Ts relative to max Ts - let's check again right below. 

rm(resids.plot,smooths.plot,mod); gc()

### -------------------------------------------------

## 4.B) min Ts_10
mod <- readRDS("model_GAMM_offsets_min_Ts_10_subset_15yrs_03.02.26.rds")
summary(mod$gam) # R-sq.(adj) = 0.731 (0.696) - very good again!
### 'TreeTypeDeciduous' non signif. & 'RegionSWA' barely signif.
# summary(mod$lme)

### 21/10: With fx = FALSE -> 0.73 (same as models from 07/25 - fx = TRUE decreases R2 of the GAM)
### 22/10: With fx = TRUE and 10 years: 0.665 (lowest R2 for now)
### 06/02/26: With s(StandAge, k = 4): 0.604 (R2 still decreasing as we simplify smooth terms)

# gam.check
gam.check(mod$gam)
# value of 'k' might be too low fo s(StandAge) again - like above

# concurvity check
concurvity(mod$gam, full = TRUE)
# same as usual

summary(mod$lme)

# Adjusted R2
r.squaredGLMM(mod$lme)
#     R2m          R2c
# 0.71 (0.02)  0.82 (0.913)
### 20/10 -> Same as other Ts variables

## 21/10: With fx = FALSE
#  R2m   R2c
# 0.062 0.758
### -> same pattern as usual for all Ts variables

### 22/10: With fx = TRUE and 10 years: 
#  R2m    R2c
# 0.670 0.788
### Same as Ts variables above

### 06/02/26: With s(StandAge, k = 4):
#  R2m   R2c
# 0.469 0.859
### -> strong decrease in R2m but variance is soaked by EP then

## Plots
library("gratia")
pdf(NULL) 
smooths.plot <- draw(mod$gam)
pdf(NULL) 
resids.plot <- appraise(mod$gam, line_col = "#238443", point_alpha = .01)
# class(smooths.plot) ; class(resids.plot)
ggsave(plot = smooths.plot, filename = "plot_smooths_GAMM_min_Ts_10_06.06.26.jpg", dpi = 300, height = 7.5, width = 7.5)
ggsave(plot = resids.plot, filename = "plot_residuals_GAMM_min_Ts_10_06.06.26.jpg", dpi = 300, height = 7, width = 7.5)

## To plot fixed effects estimates (± standard errors)
#coefs <- tidy(mod$lme, effects = "fixed")
#coeffs.plot <- ggplot(coefs, aes(x = estimate, y = reorder(term, estimate))) + geom_point() +
#    geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
#    labs(x = "Estimate", y = "Fixed Effect", title = "Fixed Effects (GAMM - lme)") +
#    theme_minimal()
#ggsave(plot = coeffs.plot, filename = "plot_coeffs_GAMM_min_Ts_10_20.10.25.jpg", dpi = 300, height = 5, width = 5)


### As expected, clear non linear weakening effect of stand age on offsets. Forests get better and better at retaining
### soil heat relative to grasslands. 

### NOTE: The distribution of offsets for soil temperatures all follow a similar profile and show consistent GAMM outputs.
### Their mean offset values are all > 0 (whatever the depth of the soil and the region and the dominant tree type)
### and these offsets (Ts grass - Ts forests) get lower with stand age. 
### DOY shows the strongest effect and this effect is of course seasonal: offsets get more positive
### (= warmer grasslands soils than forests soils) in summer. This is likely due to the shading effect of the canopy. 
### This easonal effect is slightly stronger for deciduous forests (also as expected). 

### As before, DOY and EP dominate the variance in daily soil temp offsets. Stand age plays a secondary (or even tertiary
### relative to tree type) role.

rm(resids.plot,smooths.plot,mod); gc()

### 20/10 -> Same as other Ts variables

### -------------------------------------------------

### 5) Ts_20
## 5.A) max Ts_20
mod <- readRDS("model_GAMM_offsets_max_Ts_20_subset_15yrs_03.02.26.rds")
summary(mod$gam) # R-sq.(adj) = 0.855 - very good fit again! 
### 'TreeTypeDeciduous' non signif. & 'RegionSWA' barely signif.
# summary(mod$lme)

### 21/10: With fx = FALSE -> 0.853 still very good as usual
### 22/10: With fx = TRUE and 10 years: 0.847
### 06/02/26: With s(StandAge, k = 4): 0.865 (excellent as usual)

# gam.check
gam.check(mod$gam)
# like above, k of s(StandAge) too low but we WANT it that way

# concurvity check
concurvity(mod$gam, full = TRUE)
# same as usual

# To report Phi
summary(mod$lme)

# Adjusted R2
r.squaredGLMM(mod$lme)
#  R2m    R2c
# 0.856  0.875
### 20/10 -> Same as other Ts variables

## 21/10: With fx = FALSE
#  R2m   R2c
# 0.055 0.788
### -> Same pattern as other Ts variables: fx = FALSE discard the variance attributed to fixed effects and even decreases R2c
### So worse than before I guess

### 22/10: With fx = TRUE and 10 years: 
#  R2m   R2c
# 0.847 0.871
### Excellent model

### 06/02/26: With s(StandAge, k = 4):
#  R2m   R2c
# 0.864 0.902
### EXCELLENT model again


## Plots
library("gratia")
pdf(NULL) 
smooths.plot <- draw(mod$gam)
pdf(NULL) 
resids.plot <- appraise(mod$gam, line_col = "#238443", point_alpha = .01)
# class(smooths.plot) ; class(resids.plot)
ggsave(plot = smooths.plot, filename = "plot_smooths_GAMM_max_Ts_20_06.02.26.jpg", dpi = 300, height = 7.5, width = 7.5)
ggsave(plot = resids.plot, filename = "plot_residuals_GAMM_max_Ts_20_06.02.26.jpg", dpi = 300, height = 7, width = 7.5)

### To plot fixed effects estimates (± standard errors)
#coefs <- tidy(mod$lme, effects = "fixed")
#coeffs.plot <- ggplot(coefs, aes(x = estimate, y = reorder(term, estimate))) + geom_point() +
#    geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
#    labs(x = "Estimate", y = "Fixed Effect", title = "Fixed Effects (GAMM - lme)") +
#    theme_minimal()
#ggsave(plot = coeffs.plot, filename = "plot_coeffs_GAMM_max_Ts_20_20.10.25.jpg", dpi = 300, height = 5, width = 5)

rm(resids.plot,smooths.plot,mod); gc()

### 20/10 -> Same as max Ts_05 & max Ts_10

### 06/02/26: s(StandAge) has a much weaker efefct of Ts_20 than upper soil levels. s(Anom) have much higher impact.

### -------------------------------------------------

## 5.B) min Ts_20
mod <- readRDS("model_GAMM_offsets_min_Ts_20_subset_15yrs_03.02.26.rds")
summary(mod$gam) # R-sq.(adj) = 0.858 - very good fit again! 
### 'TreeTypeDeciduous' non signif. & 'RegionSWA' non signif.
# summary(mod$lme)

### 21/10: With fx = FALSE -> 0.857
### 22/10: With fx = TRUE and 10 years: 0.843 - still very good
### 06/02/26: With s(StandAge, k = 4): 0.842 - still very good

# gam.check
gam.check(mod$gam)
# like above

# concurvity check
concurvity(mod$gam, full = TRUE)
# same as usual

# To report Phi
summary(mod$lme)

# Adjusted R2
r.squaredGLMM(mod$lme)
#  R2m    R2c
# 0.857  0.881
### 20/10 -> Same as other Ts variables

## 21/10: With fx = FALSE
#  R2m   R2c
# 0.084 0.78
### -> Same pattern as other Ts variables

### 22/10: With fx = TRUE and 10 years: 
#  R2m   R2c
# 0.842 0.874
### Very good model - like other Ts variables

### 06/02/26: With s(StandAge, k = 4):
#  R2m    R2c
# 0.834  0.896

## Plots
library("gratia")
pdf(NULL) 
smooths.plot <- draw(mod$gam)
pdf(NULL) 
resids.plot <- appraise(mod$gam, line_col = "#238443", point_alpha = .01)
# class(smooths.plot) ; class(resids.plot)
ggsave(plot = smooths.plot, filename = "plot_smooths_GAMM_min_Ts_20_06.02.26.jpg", dpi = 300, height = 7.5, width = 7.5)
ggsave(plot = resids.plot, filename = "plot_residuals_GAMM_min_Ts_20_06.02.26.jpg", dpi = 300, height = 7, width = 7.5)

## To plot fixed effects estimates (± standard errors)
#coefs <- tidy(mod$lme, effects = "fixed")
#coeffs.plot <- ggplot(coefs, aes(x = estimate, y = reorder(term, estimate))) + geom_point() +
#    geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
#    labs(x = "Estimate", y = "Fixed Effect", title = "Fixed Effects (GAMM - lme)") +
#    theme_minimal()
#ggsave(plot = coeffs.plot, filename = "plot_coeffs_GAMM_min_Ts_20_20.10.25.jpg", dpi = 300, height = 5, width = 5)

### 20/10 -> Same as Ts_05 & Ts_10

rm(resids.plot,smooths.plot,mod); gc()

### -------------------------------------------------

### On clicmal server
setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data/models/"); dir()

### 6) 29-30/01/26: RE-ANALYZE THE UPDATED max/min SM_10 AFTER RUNNING THEM ON 'Anoms'

## 6.A) max SM_10
mod <- readRDS("model_GAMM_offsets_max_SM_10_subset_15yrs_03.02.26.rds")
summary(mod$gam) # R-sq.(adj) = 0.37 (from 0.28) - not as good as soil temp. models but not bad at all. Reasonable fit.
### All smooth terms signif. 'TreeTypeDeciduous' & 'intercept' not signif.
### NOTE: a non-significant intercept means the predicted offset at the arbitrary baseline is not different from zero.
### It says nothing about the overall significance of the model or its smooth terms.
# summary(mod$lme)

### 21/10: With fx = FALSE -> 0.283 (same as previous 15/10 model) - fx = FALSE has no effect on this
### 22/10: With fx = TRUE and 10 years: 0.313, better
### 29/01/26: With fx = TRUE & 14 years & including s(Anoms, k = 2, fx = T): 0.37, good
### 30/01/26: With fx = TRUE & 14 years & including s(Anoms, k = 5, fx = T): 0.363, as good 
### 06/02/26: With s(StandAge, k = 4): 0.379, best model so far!

# gam.check
pdf(NULL)
gam.check(mod$gam)
# value of 'k' for s(StandAge may be too low again - but we want to avoid overfitting anyways)
# s(Anom) is good with k = 5

# concurvity check
concurvity(mod$gam, full = TRUE)
# s(Anom) doesnt show concurvity issues

# To report Phi
summary(mod$lme)

# Adjusted R2
r.squaredGLMM(mod$lme)
#     R2m           R2c
# 0.34 (0.083)   0.65 (0.899)
### 20/10 -> Again, signif. increase in R2m (good) but decrease in R2c. Why? 

## 21/10: With fx = FALSE 
#  R2m   R2c
# 0.259 0.724
### Similar to Ta variables, fx = FALSE decreased fixed effects' explained variance but restored some R2c
### R2c was highest for GAMM based on all 16 years + higher k values (0.90) but R2m was also lowest

### 22/10: With fx = TRUE and 10 years: 
#  R2m   R2c
# 0.324 0.622
### Best model so far, good R2c of the LME and better R2 of the GAM

### 29/01/26: With fx = TRUE and 14 years and s(Anom, k = 2, fx = T)
#  R2m       R2c
# 0.39      0.72
### Better fixed effects too! -> best model so far

### 30/01/26: With fx = TRUE and 14 years and s(Anom, k = 5, fx = T)
#  R2m       R2c
# 0.380     0.715

### 06/02/26: With s(StandAge, k = 4):
#  R2m     R2c
# 0.279  0.798
### Decrease in R2m in the LME though, but soaked up by EP, whoch is OK for SM_10

## Plots
library("gratia")
pdf(NULL) 
smooths.plot <- draw(mod$gam)
pdf(NULL) 
resids.plot <- appraise(mod$gam, line_col = "#238443", point_alpha = .01)
# class(smooths.plot) ; class(resids.plot)
ggsave(plot = smooths.plot, filename = "plot_smooths_GAMM_max_SM_10_06.02.26.jpg", dpi = 300, height = 7.5, width = 7.5)
ggsave(plot = resids.plot, filename = "plot_residuals_GAMM_max_SM_10_06.02.26.jpg", dpi = 300, height = 7, width = 7.5)

## To plot fixed effects estimates (± standard errors)
#coefs <- tidy(mod$lme, effects = "fixed")
#coeffs.plot <- ggplot(coefs, aes(x = estimate, y = reorder(term, estimate))) + geom_point() +
#    geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
#    labs(x = "Estimate", y = "Fixed Effect", title = "Fixed Effects (GAMM - lme)") +
#    theme_minimal()
#ggsave(plot = coeffs.plot, filename = "plot_coeffs_GAMM_max_SM_10_20.10.25.jpg", dpi = 300, height = 5, width = 5)


### To guide interpretation: SM offsets are > 0 on average: grasslands have higher soil moisture than nearby forests. 
### Why? 
### Forests are typically thirstier than grasslands. In many climates, especially those with limited
### or seasonal water availability, grasslands can show higher average soil moisture in the root zone
### than nearby forests. Here are some arguments why: 
### - Forests generally have higher leaf area index and deeper rooting systems than grasslands
### - This means that forests transpire more water through the canopy and extract more water from the soil profile
### - Grasslands, with shallower roots and less leaf area, often lose less water overall
### - Forest canopies intercept rainfall — water caught on leaves evaporates back to the atmosphere without entering the soil
### - In grasslands, rainfall reaches the soil more directly, which can maintain higher soil moisture after precipitation events
### - In dry seasons, forests often continue to transpire using deep water reserves, further lowering soil moisture compared to grasslands,
###   which may go dormant and reduce water use.
### - Grasslands can have dense root mats and high soil porosity in the top layer, promoting infiltration and water retention after rainfall.
### - In contrast, forest soils with high organic matter can also retain water, but deeper rooting means more extraction from all layers.

### So, still a very strong effect of EP (local variability/conditions). Residuals still look OK. 

### Stand age has (weak) negative effect on soil moisture offsets. This implies the following
## - Older forests’ soils are closer in moisture to grasslands
## - Since the gap is closing from above, the forest’s soil is getting relatively wetter compared to younger forests
## - This implies older forests are better at retaining soil moisture — potentially due to deeper rooting systems,
##   more developed organic layers, and improved canopy microclimate buffering against evaporation.

### 20/10 -> Smooths plots look the same but less complex (good) -> decrease in % deviance explained as expected
### But why such strong changes in the LME component? 

### 29/01/26: First test was with 14 years of data and including s(Anom, k = 2, fx = TRUE) but k = 2 is too restrictive! 
### -> leads to asymetric partial effects curve: increases when < 0 and then plateau -> irrealistic -> re-running with k = 5, fx = T
### -> for now, the model residual plots do not look aligned on the 1:1 plot, but look normally distrbuted with a long lower tail
### Waiting on the SM_10 models to re-run

rm(resids.plot,smooths.plot,mod); gc()

### 30/01/26: s(Anom, k = 5) still shows an asymetric response: strong slope in the negative part 
### and plateau when Anom > 0.
### When Anoms get more negative (= forests' soils get more moist than expected by the seasonal cycle)
### -> the offsets get more negative, which makes a lot of sense.  


### -------------------------------------------------

## 6.B) min SM_10
mod <- readRDS("model_GAMM_offsets_min_SM_10_subset_15yrs_03.02.26.rds")
summary(mod$gam) # R-sq.(adj) = 0.171 (weak GAM model)
# Only intercept is non significant

### 21/10: With fx = FALSE -> 0.171 (same)
### 22/10: With fx = TRUE and 10 years: 0.246 better!
### 29-30/01/26: With s(Anom, k = 5, fx = TRUE) and 14 years: 0.364
### 06/02/26: With s(StandAge, k = 4): 0.373; best model so far


# gam.check
gam.check(mod$gam)
# value of 'k' for s(StandAge may be too low again - but we want to avoid overfitting anyways)
# value of 'k' for s(Anom, k = 5) looks good

# concurvity check
concurvity(mod$gam, full = TRUE)
# Same as above for temperature variables
# no concurvity problem

summary(mod$lme)

# Adjusted R2
r.squaredGLMM(mod$lme)
#   R2m  R2c
# 0.368 0.678
### -> Most of the power of the GAMM comes from the LME, and especially the random effects

## 21/10: With fx = FALSE
#  R2m   R2c
# 0.235 0.739
### -> fx = FALSE decreases R2m (fixed effects power) but restores some R2c (random effects)

### 22/10: With fx = TRUE and 10 years: 
#  R2m   R2c
# 0.332 0.641
### Best model so far with the 15/10 model

### 30/01/26: With s(Anom, k = 5, fx = T) and 14 years of data 
#   R2m       R2c
#  0.417    0.4171
### Largest R2m increase with s(Anom, k = 5), but lower R2c.
### VERY WEIRD that they are equal!

### 06/02/26: With s(StandAge, k = 4):
#  R2m    R2c
# 0.405  0.737
### Best LME component so far, increased R2m and highest R2c

## Plots
library("gratia")
pdf(NULL) 
smooths.plot <- draw(mod$gam)
pdf(NULL) 
resids.plot <- appraise(mod$gam, line_col = "#238443", point_alpha = .01)
# class(smooths.plot) ; class(resids.plot)
ggsave(plot = smooths.plot, filename = "plot_smooths_GAMM_min_SM_10_06.02.25.jpg", dpi = 300, height = 7.5, width = 7.5)
ggsave(plot = resids.plot, filename = "plot_residuals_GAMM_min_SM_10_06.02.25.jpg", dpi = 300, height = 7, width = 7.5)
## To plot fixed effects estimates (± standard errors)
#coefs <- tidy(mod$lme, effects = "fixed")
#coeffs.plot <- ggplot(coefs, aes(x = estimate, y = reorder(term, estimate))) + geom_point() +
#    geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
#    labs(x = "Estimate", y = "Fixed Effect", title = "Fixed Effects (GAMM - lme)") +
#    theme_minimal()
#ggsave(plot = coeffs.plot, filename = "plot_coeffs_GAMM_min_SM_10_20.10.25.jpg", dpi = 300, height = 5, width = 5)

rm(coefs,resids.plot,smooths.plot,mod); gc()

### 30/01/26: Similar to max SM_10, but residuals look way better aligned along the expected 1:1 line
### s(Anom, k = 5, fx = T) shows the same smooth as max SM_10: strong decrease towards negative effects
### (offset gets lower as anomalies get more negative). Plateau for Anom > 0.

### 06/02/26: Very weak effetc of stand age, most of the power of the GAM comes from DOY and Anom 

### ------------------------------------------------------------------------------------------------------------

### 20/10/25: Above, I noticed a subtle (but important) shift in how variance is being allocated.
### (1) smaller / fixed spline bases (k reduced and fx = TRUE) and (2) subsampling years affect how much of the
### signal the fixed part can capture and how much is left to the random intercepts / residuals. 

## Remember:
# R2m (marginal) — proportion of total variance explained by fixed effects only
# R2c (conditional) — proportion explained by fixed + random effects

### Potential reasons for the increase in R2m and decrease in R2c:

## 1) Subsetting years changed the variance structure
## If EP-specific baseline differences or EP-specific time effects are smaller in this subset, the estimated
## random-intercept variance will shrink -> reduces the random component of explained variance (reducing R2c)
## while leaving fixed effects to explain relatively more of the remaining variation (increasing R2m)

## 2) Different estimated AR(1) / residual variance
## (I don't think so)

## 3) Fixing spline complexity (fx = TRUE) and reducing k
## fx = TRUE fixes the amount of wiggliness (i.e., you do not estimate the smoothing parameter λ).
## Reducing k and fixing smoothness reduces the GAM’s effective flexibility. You noted the GAM component R² decreased — that makes sense.
## A less flexible GAM can sometimes assign more of the structured variation to simple parametric fixed effects (e.g., TreeType, Region)
## because the splines are less able to soak up subtle, smooth variation. That increases R2m.
## If the GAM is less flexible, some remaining structure may be left in the residuals and picked up somewhat by random intercepts

## 4) Estimation/shrinkage effects from fewer data per EP
## With fewer years, each EP might have fewer observations contributing to its random-intercept estimate.
## That may increase shrinkage (random effects closer to zero), reducing random variance estimates and shifting
## apparent explanatory power to fixed effects


### Recommendations: 

## - If you care about population-level fixed effects (TreeType, Region, DOY shape): increase in R2m may be acceptable
##   (or even desirable!!). A simpler smooth may actually improve interpretability and reduce overfitting

## - If you care about capturing EP-specific variability (e.g., predictions per EP or variance partitioning): 
##   decrease in random variance is a concern! To recover EP variance you may use larger k or fx = FALSE
##   for at least some smooths, or fit on more years (if computationally feasible)

## - Keep AR(1) in the model; it seems justified by AIC and ACF. If memory is the bottleneck, keep the subset approach
##   but document it and run sensitivity analyses.

### -> Run the two-model experiment on the same subset (code above). That isolates subsampling vs smoothing choices.
### (re-run R code on same 6 years, and use fx = FALSE but keep same values of 'k')
### -> Remove s(precipitation)
### -> Inspect: VarCorr() differences, rho differences, sigma differences, EDF changes (summary(mod$gam)).

# To extract variance components
VarCorr(mod$lme) # or:
# vc <- as.data.frame(VarCorr(mod$lme))

# Extract residual standard deviation
summary(mod$lme)$sigma

# To compare rho (AR1) estimates
mod$lme$modelStruct$corStruct
# coef(mod$lme$modelStruct$corStruct, unconstrained = FALSE)

### Why comparing them tells you why the shift occurred: 
## If VarCorr (EP variance) decreased substantially, then the random intercepts simply have less room to explain variability
## -> this directly explains why R2m increased (more variance is being pushed into fixed effects)

## If at the same time sigma increased, then part of that “lost” random variance is now just unexplained noise,
## i.e. model simplicity is driving underfitting at that level

## If rho changed, it tells you whether the time correlation structure is absorbing variance that used to sit in
## the smooth of DOY or in random EP differences.
## Together with R2c, those changes tell you where the model reallocated variance after you:
# - reduced k
# - fixed smoothing (fx = TRUE)
# - subsampled years

## What you see in output               ## The interpretation
# VarCorr(EP) gets smaller	            Subset or smooth constraints (k, fx) reduced EP differences
# sigma gets bigger	                    The model is underfitting some signal (lost smooth complexity)
# rho gets bigger	                    Autocorrelation structure is pulling weight from GAM/random effects
# R2m increase and R2c decrease	        More of the structure is now attributed to fixed effects

### Do this here: 
setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data/models"); dir()

# For each variable, load the old and new model and run diagnostics mentioned above

vars <- c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")
stats <- c("max","min")
# For testing:
v <- "Ta_200"
s <- "max"

### WARNING: WILL FAIL AT min SM_10 SINCE NO OLD GAMM MODEL

### 21/10: Adding another test: true intraclass correlation (ICC).
### Goal: Is the GAMM learning ecological variation, or just “soaking up” variance via random intercepts?

###                   VarCorr(EP) / ( VarCorr(EP) + sigma^2 )

### If ICC > ~0.6 for a microclimate variable, it means that:
### - site explains more variation than stand age + DOY + tree type. Suggesting either:
### - missing covariates (e.g., local topography, LAI, groundwater access),
### - or heterogeneities you may model later

for(v in vars) {
    
    for(s in stats) {
        
        # Useless message
        message(paste("\n","Loading GAMMs for ",s," ",v,"\n", sep = ""))
        # mod_old <- get(load( dir()[grepl(paste("model_GAMM_offsets_",s,"_",v,"_31.07.25.RData", sep = ""),dir())] )) # For GAMMs ran in 07/25 (all 16 years, fx = FALSE, more complex k)
        # mod_old <- readRDS( dir()[grepl(paste("model_GAMM_offsets_",s,"_",v,"_subset_6yrs_15.10.25.rds", sep = ""),dir())] ) # GAMMs ran on the 15/10 (fx = TRUE, simpler k, 6 years)
        # mod_old <- readRDS( dir()[grepl(paste("model_GAMM_offsets_",s,"_",v,"_subset_10yrs_22.10.25.rds", sep = ""),dir())] ) # GAMMs ran on the 22/10 (fx = TRUE, simpler k, 10 years)
        mod_old <- readRDS( dir()[grepl(paste("model_GAMM_offsets_",s,"_",v,"_subset_15yrs_23.10.25.rds", sep = ""),dir())] ) # GAMMs ran on the 22/10 (fx = TRUE, simpler k, 15 years)
        mod_new <- readRDS( dir()[grepl(paste("model_GAMM_offsets_",s,"_",v,"_subset_14yrs_25.11.25.rds", sep = ""),dir())] ) # GAMMs ran on the 25-01/12 (fx = TRUE, simpler k, 14 years + mean grassland value)

        # VarCorr() = amount of variability explained by random intercepts
        vc_old <- as.data.frame(VarCorr(mod_old$lme)) # dim(vc_old); str(vc_old); class(vc_old)
        # message(paste("VarCorr (EP variance) of OLD model is = ",round(as.numeric(vc_old$x[39,"Variance"]),4), sep = ""))
        message(paste("VarCorr (EP variance) of OLD model is = ",round(as.numeric( VarCorr(mod_old$lme)[1,"Variance"] ),4), sep = ""))
        message(paste("VarCorr (EP variance) of NEW model is = ",round(as.numeric( VarCorr(mod_new$lme)[1,"Variance"] ),4),"\n", sep = ""))

        # Residual standard deviation = Sigma
        # Sigma increases when smooth complexity is lost
        # Aka; a smaller sigma means that at the level of the residuals, the model is performing better: there’s less variability left over that is not accounted for by the model
        # lower sigma -> higher conditional R² (usually) ; BUT a reduction in sigma does not necessarily mean the whole model is better
        message(paste("Sigma (residual standard deviation) of OLD model is = ",round(summary(mod_old$lme)$sigma,4), sep = ""))
        message(paste("Sigma (residual standard deviation) of NEW model is = ",round(summary(mod_new$lme)$sigma,4),"\n", sep = ""))
        
        # Importance of autocorr. structure (rho/Phi) -> increases if autocorr. structure gains more importance
        message(paste("Rho (weight of autocorr.) of OLD model is = ",round(mod_old$lme$modelStruct$corStruct,3), sep = ""))
        message(paste("Rho (weight of autocorr.) of NEW model is = ",round(mod_new$lme$modelStruct$corStruct,3),"\n", sep = ""))

        message(paste("R2m of OLD model is = ",round(r.squaredGLMM(mod_old$lme)[,1],3), sep = ""))
        message(paste("R2m of NEW model is = ",round(r.squaredGLMM(mod_new$lme)[,1],3),"\n", sep = ""))

        message(paste("R2c of OLD model is = ",round(r.squaredGLMM(mod_old$lme)[,2],3), sep = ""))
        message(paste("R2c of NEW model is = ",round(r.squaredGLMM(mod_new$lme)[,2],3),"\n", sep = ""))

        # And deviance explained of the GAM component
        message(paste("R squared of OLD model is = ",round(summary(mod_old$gam)$r.sq,3), sep = ""))
        message(paste("R squared of NEW model is = ",round(summary(mod_new$gam)$r.sq,3),"\n", sep = ""))

        # ICC for OLD model
        var_ep_old <- as.numeric(VarCorr(mod_old$lme)[1,"Variance"] )
        sigma_old <- summary(mod_old$lme)$sigma
        icc_old <- var_ep_old / (var_ep_old + sigma_old^2)
        
        # ICC for NEW model
        var_ep_new <- as.numeric( VarCorr(mod_new$lme)[1,"Variance"] )
        sigma_new <- summary(mod_new$lme)$sigma
        icc_new <- var_ep_new / (var_ep_new + sigma_new^2)
        message(paste("ICC of OLD model is = ", round(icc_old, 3), sep = ""))
        message(paste("ICC of NEW model is = ", round(icc_new, 3), "\n", sep = ""))

        # Remove and clean
        rm(mod_old,mod_new,vc_old,var_ep_old,sigma_old,icc_old,var_ep_new,sigma_new,icc_new) ; gc()
        message(" ------------------------------------------------------------------------------------ ")

    } # eo 2d for loop - s in stats

} # eo 1st for loop - v in vars 


### Analyzing outputs

## max Ta_200: 
# The GAM smooths are performing roughly the same; the changes in R2c/R2m are due to variance partitioning
# among fixed, random, and residual, not the GAM fit itself.
# Fixed effects explain more in the new GA%% (R²m increases), possibly due to differences in smoothing (fx = TRUE)
# or how subsetted data distributes variance. Conditional R2 drops because total variance in this subset/simpler model is smaller,
# so proportion explained by combined effects decreases.
# Here, fixed effects are stronger in the NEW model


## min Ta_200: 
# Decrease in R2m (0.241 -> 0.122) -> the fixed effects now explain much less variance,
# probably due to smaller k or fx = TRUE in the smooths. TreeType + Region + smooths are less flexible
# -> less variance explained by fixed structure
# We also see a large decrease in total explained variance (R2c). The combined effect of simpler smooths
# + smaller EP variance + subset of years -> less variance captured overall
# The GAM smooth itself is performing about the same or slightly better, even though total variance captured by
# fixed/random effects drops/
# -> Subsetting years + simpler/fixed smooths trades flexibility for computational stability, at the cost of less
# variance captured by fixed/random effects.


## max Ta_10:
# Substantial increase in EP variance (0.28 -> 0.50) -> the NEW model relies much more on EP-level random effects to explain variation.
# Could indicate that the subset of years emphasizes between-EP differences more strongly, or that fixed effects/smooths
# are now absorbing less of the variation at the population level.
# Large increase in fixed-effects explained variance.
# Fixed effects + smooths now capture much more structured variation,
# likely due to the choice of subset years and smooth simplification (fx = TRUE, smaller k).
# R2c essentially unchanged -> total variance explained by the combination of fixed + random effects is about the same
# GAM smooth itself is performing the same; again, the smooths are robust across OLD vs NEW models
# -> For max Ta_10, the new smoothing and subset choice improves fixed-effect capture without altering total fit!


## min Ta_10: 
# Increase in EP variance -> new model relies more on EP-level random effects to explain variation.
# This could reflect that the subset of years emphasizes between-EP differences more than the dataset with 16 years
# Slight increase in 'R2m', but very large decrease (0.69 -> 0.204) in 'R2c'
# -> subset of years / model simplification reduces overall variance captured by fixed + random effects combined
# NEW model shifts variance toward fixed effects and EP variance slightly.
# Total explained variance drops sharply (R²c ↓) → consequence of subset of years and simplified smooths.
# Residual variance remains dominant. GAM smooths lose some explanatory power here (dev.expl ↓).
## -> Contrary to max Ta_10, the new GAMM for min Ta_10 is not better. 


## max Ts_05:
# We observe a large increase in fixed-effects explained variance. The new model’s fixed effects
# (smooths + TreeType + Region) are now capturing most of the structured variation.
# Likely caused by fx = TRUE smooths and the subset of years, which allow the model to fit fixed structure more aggressively.
# Decrease in total explained variance (by 0.14 points only though). Even though fixed effects explain much more (R²m incr.),
# total variance explained drops because random effect contribution decreases and total variance in the subset is smaller
# -> max Ts_05 shows a shift from random effects toward fixed effects with the NEW model.
# The subset of years and fx/simpler smooths greatly amplify fixed-effect contribution.


## min Ts_05:
# min Ts_05 behaves very similarly to max Ts_05 -> new model shifts variance from random effects toward fixed effects, with smooths performing slightly better.

### Ts_10 and Ts_20: All the same as max Ts_05
### -> new models shift variance from random effects toward fixed effects, with smooths performing slightly better! 
### (better models overall)

## max SM_10:
# EP variance basically unchanged. Residual standard deviation also unchanged. Same for Phi.
# Large increase in R2m -> Substantial increase in variance explained by fixed effects.
# NEW model’s fixed effects and smooths capture more structured variation.
# BUT, decrease in R2c (total variance) because residual variance dominates and the dataset subset reduces overall variance
# Slight decrease in R2 -> confirms that, overall, the NEW model explains slightly less total variance.
# -> fixed effects capture more variance, but total explained variance (R2c) often decreases due to dataset subsetting


### Key takeaways are: 
# - Fixed effects in the NEW models generally explain more variance (R2m increases).
# - Random effects (EP variance) are mostly stable but vary slightly depending on the variable.
# - Residual variance (sigma) decreases slightly most cases.
# - Total variance explained (R2c) often decreases due to subsetting, even when fixed effects gain power
# - GAM smooth contributions (deviance explained) remain robust across model simplifications.
# - AR1 correlation (rho / Phi) shifts only slightly.

### Two things to do:
### X Re-run without precip. and turn 'fx' to 'FALSE' and re-assess differences (DONE)
### What fx = TRUE does: 
# Fixed regression spline: It tells the gam() function to treat the term as a fixed regression spline rather than
# a penalized regression spline, where smoothing parameters are estimated.
# Freezing smoothing: By fixing the smoothing parameter at zero, it effectively removes the penalty and turns
# the term into a simple regression spline with fixed degrees of freedom. 


### X Re-run without precip. keep 'fx' to 'TRUE', and try adding more years (10 instead of 8) and re-assess differences
### (DONE - see interpretation below line 1480)


### --------------------------------------------------

### 21/10/25: Re-ran without precip. and turn 'fx' to 'FALSE' and re-assess differences. 
### OLD model has fx = TRUE; NEW model has fx = FALSE.
### When we turned fx = FALSE, the smooths became fully penalised, which means the model is free to shrink them
### nearly to flat lines if it thinks the penalty demands it.
### This tests how much wiggliness was "forced" into the old models (fx = T), versus how much smoothing the data
### actually supports once smoothness penalties are active.

### Interpreting the main ouputs: 
## - Switching from fx = TRUE to fx = FALSE drastically weakens the fixed-effect smooths, causing:
## - Instead of capturing structured seasonal or stand-age effects in the GAM component, the new model is flattening the smooths
## - the signal is being reallocated to either the random intercepts or the temporal correlation
## - GAM smooths are explaining real structure — but they do so only when they are not penalised toward zero (fx = TRUE)
## -> fixed smooths should not be penalised, because their role is inferential / explanatory, not purely predictive

## NOTE: Why soil moisture behaves a bit differently?
# - R2c improves a bit (from 0.65 -> 0.725 for max SM), and variance shifts into RE/AR1, but the smooths still
#   collapse (R2m decreases), meaning:
# - the EP-level variation completely dominates SM
# - fixed smooths contribute very little regardless of penalisation
# -> This is consistent with SM being more site-driven than season-driven -> SM is actually telling a biologically meaningful story.

## NOTE: fx = FALSE is most useful when you’re building a purely predictive statistical model or trying to test whether the smooth exists

### CCL -> If the main purpose is explanation / interpretation -> automated smooth shrinking: x = TRUE version is the correct specification


### Interpreting the extra ICC tests: 
## - For air temperature variables, ICCs are quite low (0.09–0.16).
##   This indicates that most of the variability in daily temperature offsets occurs within EPs rather than between EPs.
##   Random intercepts capture some EP-level differences, but they are not dominant.

## - For soil/near surface temperatures, ICCs are moderate (0.22–0.40).
##   This suggests that EP-level differences explain a substantial portion of variation, particularly for minimum temperatures
##   Within-EP day-to-day variation still matters, but EP identity plays a more important role here than for air temperature

## - For soil moisture, ICCs are even higher (0.47–0.49).
##   Almost half of the variance is explained by differences between EPs. This is intuitive: soil moisture is strongly site-dependent,
##   influenced by local soil texture, water retention, and landscape position.
##   Temporal variation within EPs is smaller relative to between-EP differences.

### CCL: 
## -> Air temperature offsets -> mostly within-EP variability, low ICC.
## -> Soil/near surface temperature -> EP differences more important, moderate ICC.
## -> Soil moisture -> strong EP-level structuring, high ICC.
## -> Overall, the ICC values reflect how much the random intercept matters:
## -> For SM, EP identity is crucial. For Ta, the day-to-day signal dominates. For Ts, it’s intermediate.

### --------------------------------------------------

### 22/10/25: Re-ran without precip. and turn 'fx' to 'FALSE' and used 10 years of data instead of 6

### So, what happens when you use 10 years of data instead of 6? 

## - VarCorr and ICC drop -> EP effects become less dominant over time-scale variation
## -> Between-plot (EP) differences become relatively smaller 

## - sigma increases & autocorrelation becomes a bit stronger -> Within-plot temporal variability increases
## (i.e., more hot years, cold years, droughts) as one may have expected

## - The smooths do not break (GAM R² stays very stable) -> models ARE STABLE :) ->  ecological relationships remain
## structurally consistent

## - R2 doesn’t drop that much, despite more noise and variance being added -> our GAMMs may generalize,
## not overfitting to 6 years

### CCL: Extending from 6 to 10 years makes the models more realistic (more temporal heterogeneity),
### slightly reduces EP-driven structure (ICC decreases), but the functional ecological drivers remain stable.


### --------------------------------------------------

### 22/10/25: Computing model diagnostics like above but to compare the models from the 22/10/25
### (fx = TRUE, simpler k, 10 years of data) to the old legacy models from back in July (07/25) 
### (with: fx = FALSE, complex k (20 and 8) and 16 years of data). 

### By doing so, we ask to 2 things: 
### -> which variables gain most from 16y vs 10y?
### -> what does a change in ICC mean ecologically?
### (e.g., “offset is more site-structured (legacy effect)” vs “offset is more weather-driven”)

## - For temperature variables (Ta & Ts), the new 10-year fx = TRUE models show dramatic improvement in marginal R²
## meaning fixed effects explain much more variance.

## - For soil moisture (SM), fixed effects improve moderately, conditional variance slightly decreases.

## - ICC values slightly decrease in most cases, meaning less reliance on random effects.

## - Residual SD (Sigma) slightly higher in some new models, but generally comparable.

### -> Overall, the new models are more interpretable with respect to fixed effects,
### though conditional variance sometimes decreases slightly.

### Here’s a simplified traffic-light summary of your new 10-year fx = TRUE models compared to the legacy 16-year k = 20 models.
### I focused on overall improvement in interpretability and explanatory power, taking into account R2m, ICC, and GAM R2:

### Variable	Overall Assessment	    Comment
# max Ta_10	    🟢 Improved	            Fixed effects now explain most variation; ICC slightly lower
# min Ta_10	    🟢 Improved	            Marginal R² improved; more interpretable fixed effects

# max Ta_200	🟡 Similar	            Small changes; marginal R² slightly better, overall variance similar
# min Ta_200	🟡 Similar	            Minor changes; small decrease in R2m, ICC unchanged

# max Ts_05	    🟢 Improved	            Fixed effects dominate; model more interpretable
# min Ts_05	    🟢 Improved	            Huge improvement in R2m; better fixed-effect explanation
# max Ts_10	    🟢 Improved	            Fixed effects vastly improved; R2m much higher
# min Ts_10	    🟢 Improved	            R2m incqreased dramatically; better fixed-effect interpretability

# max SM_10	    🟡 Similar	            Moderate improvement in fixed effects; total variance slightly lower but still ok


### ------------------------------------------------------------------------------------------------------------

### 23/10/25: Computing model diagnostics like above but to compare the models from the 23/10/25
### (fx = TRUE, simpler k, 15 years of data) to the same models but based on 10 years of data and those
### based on 16 years (fx = TRUE, simpler k)

### Last quality checks. We will likely be using these models from the 23/10 to finish the forest microclimates
### reconstructions. 

### Main observations: 
## - Adding more years almost always improved ICC
## - The EP-level variance (VarCorr) got larger in most variables OR the residual variance shrank — both increase ICC.
## -> This means: between-site (EP) structure is more stable and better identified when training on a longer time series.

## - Marginal and conditional R² improved for most variables. Particularly for soil temperatures and soil moisture.
## -> This means the fixed effects are more robustly estimated with 15 y of data than 6 or 10

## - The biggest gains are in subsurface & slower processes, not surface air. Ts_05, Ts_10, Ts_20: ICC and R2 improve consistently, often meaningfully
## SM_10 also clearly improves in ICC and R2. Air temperature (Ta) shows only very small changes 
## -> expected, because air T has much higher short-term stochasticity and weaker site memory than soil processes

## - No sign of “overfitting to more time” -> the NEW models detect more stable site-level structure, residuals shrink slightly,
## autocorrelation (rho) stays stable -> no “spurious time trends”

### -> This means the longer time series adds signal rather than noise.

# Variable class	    Gains from more years?	        Why
# Air temperature	    🟢 small	                    Mostly short-memory processes
# Soil temperature	    🟢🟢 moderate	                Slow buffering, canopy/soil structure effects accumulate
# Soil moisture	        🟢🟢🟢 largest	                Strong persistence & hydrological identity of sites


### ------------------------------------------------------------------------------------------------------------

### 03/12/25: Computing model diagnostics like above for the final GAMMs (14 years of data, fx = TRUE, 
### lower 'k' values and including grassland values) - Models that were run between the 24/11/25 and the 01/12/25
### The model file were all named under the "subset_14yrs_25.11.25.rds" labeL.

### Reminder
## What you see in output               ## The interpretation
# VarCorr(EP) gets smaller	            Subset or smooth constraints (k, fx) reduced EP differences
# sigma gets bigger	                    The model is underfitting some signal (lost smooth complexity)
# rho gets bigger	                    Autocorrelation structure is pulling weight from GAM/random effects
# R2m increase and R2c decrease	        More of the structure is now attributed to fixed effects

### Overall summary

## Air temperature (Ta)
# - Substantial improvement across all metrics
# - Best gains in minimum temperatures (R² doubled)

## Shallow soil temperature (Ts)
# - Mixed results
# - R2 sometimes decreases (−5% to −20%)
# - ICC increases (30–120%)
# - Variance captured by EP increases a lot -> Suggests Mean_Grassland_Value is not a necessary predictor of Ts offsets.

## Deeper soil temperature (Ts_20)
# - Similar to Ts_05 and Ts_10: ICC and EP variance increase

## Soil moisture (SM_10)
# - Consistent improvement in predictive skill and errors.

### Higher ICC values imply that EP-level characteristics explain a larger proportion of the variability
### BUT is it that bad that ICC increases for soil-related variables? 
### Belowground, spatial variation may be abrupt and strongly tied to unmeasured local soil structure. 
### Low ICC would actually be suspicious, because it would imply that your model somehow explains soil heterogeneity
### without measuring it. When DOY × TreeType captures seasonal and vegetation-related patterns more effectively,
### the temporal variation is explained away by fixed effects. Once this is accounted for, EP-level variations remain 
### and “micro-site identity” becomes proportionally more important.

### Therefore: 
# -> High ICCs in soil-related variables are expected
# -> Increased ICC after adding strong temporal fixed effects may be conceptually correct
# -> This does not indicate model misspecification


### ------------------------------------------------------------------------------------------------------------

### 08/01/26: ISSUE WHILE WRITING R Script#4.5.1: We actually need to also save the DOY anomalies models (for non SM_10 variables)
### to build the prediction tables. You must compute 'Anom' for 'newdata' using the same anomaly-defining model that was used when 
### you created Anom in combin_df. That means:
# - Do NOT fit a new GAM or GAMM
# - Do NOT use a DOY-only model
# - Do NOT aggregate EPs
### Instead, save and re-load the original model created to derive the 'Anom' predictor.
### Anom used at prediction time must be generated by the same data-generating model as the Anom used at training time, meaning:
# - same fixed effects
# - same random‐effect structure
# - same smoothing basis / knots
# - same centering logic
### REMEMBER: Whenever a covariate is created by a model, that model becomes part of your data pipeline and must be saved!!

### Therefore, use the main FUN above to save the DOY GAMMs in: 
### /home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data/models

#var <- "SM_10"
#stat <- "max"
#age <- 75

gamm_DOY_anom <- function(var, stat, age) {

        #' This function takes three arguments and returns a model object of class 'gamm':
        #' @param var The variable to model (character): "Ta_200","Ta_10","Ts_05","Ts_10","Ts_20","SM_10" 
        #' @param stat The daily stat of the associated variable (character): 'max' or 'min'
        #' @param age Maximum stand age to account for in the GAMM approach (integer)
        #' @return A formatted data.frame combining the inputs.
      
        # Useless message
        message(paste("Loading the aggregated data for ",stat," ",var, sep = ""))

        # Read in the data after identifying the corresponding
        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data") #; dir()

        file <- dir()[grepl(paste("metadata",stat,var,sep = "_"),dir())] # file
        df <- readRDS(file)
        # dim(df); str(df) # should be 858'600 daily measurements

        # Add month and date from DOY and Year
        df$Date <- as.Date(df$DOY - 1, origin = paste0(df$Year, "-01-01"))
        df$Month <- lubridate::month(df$Date)

        # Add sine and cosine of DOY
        df$sin_DOY <- sin(2 * pi * df$DOY / 365)
        df$cos_DOY <- cos(2 * pi * df$DOY / 365)

        # Rarefy to have a bablanced dataset to mode offsets as a function of stand age and tree types
        decid_df <- df[df$TreeType == "Deciduous" & df$StandAge < age,]
        coni_df <- df[df$TreeType == "Coniferous" & df$StandAge < age,]

        # Rbind 
        combin_df <- rbind(decid_df,coni_df)
        rm(coni_df,decid_df,file,df,file); gc()

        # To model a different seasonal pattern (smooth over DOY) for each TreeType, use the by = argument inside the smooth function
        combin_df$TreeType <- factor(combin_df$TreeType)
        combin_df$Region <- factor(combin_df$Region)
        colnames(combin_df)[4] <- "Offset"

        ### 14/10/25: Remove rows with NaN
        combin_df <- combin_df %>% drop_na(StandAge)
        combin_df <- combin_df %>% drop_na(Offset)

        message(paste("Computing anomalies to DOY based on a simple GAM", sep = ""))
        anom_mod <- gamm(Mean_Grassland_Value ~ s(DOY, bs = "cc"), random = list(EP = ~1), data = combin_df)
        # summary(anom_mod$gam)

        ### Save GAMM object as .RDS
        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data/models")
        message(paste("Saving the DOY-anomaly GAMM for ",stat," ",var, sep = ""))
        saveRDS(anom_mod, file = paste("model_GAMM_DOY_anom_",stat,"_",var,"_29.01.26.rds", sep = "") )
        
        # Clean and go next
        rm(combin_df,anom_mod); gc()
    
} # eo FUN - model_offset_gamm

### Apply gamm_DOY_anom() in for loops: c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")
for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")) {
    for(s in c("max","min")) {
        gamm_DOY_anom(var = v, stat = s, age = 75)
    } # eo 2nd for loop
} # eo 1st for loop

### 29/01/26: Re-running DOY GAMs for SM_10 as well 
### -> Will be needed to standardize and detrend the 'Anoms' used for Offset predictions
gamm_DOY_anom(var = "SM_10", stat = "max", age = 75)
gamm_DOY_anom(var = "SM_10", stat = "min", age = 75)

### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------