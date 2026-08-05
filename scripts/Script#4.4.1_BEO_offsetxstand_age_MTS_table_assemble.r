### ------------------------------------------------------------------------------------------------------------

### 25/06/25 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to load the dataset created by R Script#1.4 (mean stand age and MTS of each forets EP) and combine it
### with estimates of mean temperature offsets between grasslands and forests (Tgrass - Tfor) based on the daily 
### temperatures of the top 5 closest grasslands (also given in the dataset from R Script#1.4). Compute the mean 
### offset for each day (do everything at the daily scale) based on the 5 neighbouring grasslands' data but weight
### that average based on the grassland plots' distance to the forest (closer = higher weight). 
### Once combined, you will use multivariate regression models (with random effects for plot ID nested within region)
### such as GAMM and MEM (lme() and gamm4()) to model the effect of stand and tree type on temperature offsets
### (-> R Script#4.4.2).
### Do this for air temperature, but also  skin and soil temperature. Including seasonality and random forest/site effects
### is critical to avoid spurious trends.

### Approach can be applied to near-ground (“skin”) temperature and shallow soil temperature.
### This because forest canopy strongly shades the ground, so skin temperature (10 cm above soil)
### will be cooler in forests, especially in summer, and this offset may depend on stand age
### because older stands = denser shade. Plus, soil at 10 cm is buffered by the soil itself, but is
### still influenced by shading and litter. Thus it makes sense to model these similar to air tempertaure at 2m.
### You may find that offsets for skin are large (e.g. grass heats much more in sun), while deep soil offsets are smaller.

### We may perform similarly analyzes on soil moisture differences (forest minus grassland) against stand age.
### Soil moisture influences evapotranspiration and can co-drive temperature buffering.
### For example, wetter soils can enhance cooling under trees.
### You would compute ΔSM = SM_grassland – SM_forest (or vice versa) and regress this on stand age and covariates,
### with the same mixed‐model structure. This makes sense as an auxiliary analysis, but keep in mind soil moisture
### is affected by precipitation, soil texture, and root uptake in complex ways.

# Main sources: empirical studies of forest–grassland microclimate differences and mixed‐model methods
# (e.g. Zellweger et al. 2019; Lindenmayer et al. 2022; Zhang et al. 2024)

### Last update: 30/10/25 (Re-running R Script#4.4.1 to add Precipitation -> might be useful for modelling SM_10 offsets)

### ------------------------------------------------------------------------------------------------------------

# Basic libraries
library("purrr")
library("tidyr")
library("dplyr")
library("reshape2")
library("lubridate")
library("ggplot2")
library("parallel")

### ------------------------------------------------------------------------------------------------------------

### 1°) Go to EP metadata directory and load dataset from R Script#1.4
setwd("/home/fbenedetti/Exploratories/EP metadata")
df <- get(load("table_stand_age_MTS_distances_elevation_forest_EPs_24_06_25.RData"))
# dim(df) ; str(df) ; summary(df)
# unique(df$Exploratory)

# Quickly change the 'Exploratory' levels so they match the usual labels I give them:
df[df$Exploratory == "ALB","Exploratory"] <- "SWA"
df[df$Exploratory == "HAI","Exploratory"] <- "HND"
# summary(factor(df$Exploratory))


### 2°) Create a function to be applied within a mclapply() that will help you create the final dataset
### that will be used for the multivariate regression models (MEM or GAMMs). 
### The final table should include:
## - Exploratory
## - Forest EP
## - Mean stand age
## - MTS
## - MTS_type
## - Date (D/M/Y)
## - Day of the year
## - Month
## - Year
## - Original value if the forest EP's target variable (min/max temperature)
## - Mean offset based on 5 daily offset values from the 5 closest grassland EPs

## NOTE:
## - Compute the offsets for EACH pairs of EPs and then average, weight by distance to the forest EP
## - Correct for differences in elevation because of adiabatic corrections (0.65°C/100m) in the case 
##   of air temperature - not applicable to skin & soil temperature


### Master FUN
# For testing FUN below: 
# EP <- "AEW04"
# var <- "precipitation"
# stat <- "total"

compute_table_offsets <- function(EP, var, stat) {

    #' This function takes four arguments and returns a formatted data.frame:
    #' @param EP the identity of the forest expeimental plot (EP), one of 'forest_ids'
    #' @param var the climate variable to process (character) - one of the following: 
    #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20" and "SM_10"
    #' @param stat the daily statistic (character): 'max' or 'min'
    #' @return A formatted data.frame combining the variables listed below section 2°)

    ## Useless message to start with
    setwd(paste("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/", sep = ""))

    if( var == "precipitation" ) { stat <- "total" }

    message(paste("\nLoading the ",stat," ",var," data for ",EP, sep = ""))

    ## Retrieve info of target EP from 'df'
    df_subset <- df[df$EP == EP,]
    # dim(df_subset) ; head(df_subset)

    # Extract region name out of EP name because directories are split per region
    if( grepl("A", EP) ) {
        region <- "SWA"
    } else if( grepl("S", EP) ) {
        region <- "SCH"
    } else if( grepl("H", EP) ) {
        region <- "HND"
    } # eo else if loop

    ## Go to dir and load the data of interest
    setwd(paste("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily/forests/",region, sep = ""))

    ### 30/10/25: Adding precipitation
    if( var == "precipitation" ) {
         file <- dir()[grepl(paste(var,sep = "_"), dir())]
        # Read the .csv file
        obs_daily_forest <- read.csv(file, h = T, sep = ",", dec = ".")
        # dim(obs_daily_forest) ; str(obs_daily_forest)
    } else {
        file <- dir()[grepl(paste(var,stat,sep = "_"), dir())]
        # Read the .csv file
        obs_daily_forest <- read.csv(file, h = T, sep = ",", dec = ".")
        # dim(obs_daily_forest) ; str(obs_daily_forest)
    } # eo if else loop

    ## Subset data from target EP
    obs_daily_forest <- obs_daily_forest[obs_daily_forest$EP == EP,]
    # dim(obs_daily_forest) ; summary(obs_daily_forest)

    ## Only keep data from May 2009 and beyond
    obs_daily_forest$Date <- as.Date(obs_daily_forest$Date)
    obs_daily_forest <- obs_daily_forest[which(obs_daily_forest$Date >= "2009-05-01"),]
    obs_daily_forest <- obs_daily_forest[which(obs_daily_forest$Date <= "2024-12-31"),]
    # summary(obs_daily_forest) ; obs_daily_forest[is.na(obs_daily_forest$),]

    # Need to adjust a colname first for next actions
    colnames(obs_daily_forest)[3] <- "variable"
    # obs_daily_forest[is.na(obs_daily_forest$variable),]

    # Sanity check
    #if( sum(is.na(obs_daily_forest$variable)) > 0 ) {
    #    stop(
    #        paste("!!! NAs still present in observed ",paste(var,stat, sep = "_")," file for ",EP," !!!", sep = "")
    #    )
    #} # eo if loop - sanity check

    ## Now, identify the top 5 closest grassland EP from 'df_subset'
    grass_ids <- unique(df_subset$EP2)

    ## And get their equivalent data for the "2009-05-01" to "2024-12-31" period
    setwd(paste("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily/grasslands/",region, sep = ""))
    ### 30/10/25: When adding precipitation:
    if( var == "precipitation" ) {
        file2 <- dir()[grepl(paste(var,sep = "_"), dir())] # file2
        # Read the .csv file
        obs_daily_grass <- read.csv(file2, h = T, sep = ",", dec = ".")
        # dim(obs_daily_grass) ; str(obs_daily_grass)
    } else {
        file2 <- dir()[grepl(paste(var,stat,sep = "_"), dir())] # file2
        # Read the .csv file
        obs_daily_grass <- read.csv(file2, h = T, sep = ",", dec = ".")
        # dim(obs_daily_grass) ; str(obs_daily_grass)
    } # eo if else loop
    
    ## Subset data from the 5 target grassland EPs
    obs_daily_grass <- obs_daily_grass[obs_daily_grass$EP %in% grass_ids,]
    # dim(obs_daily_grass) ; summary(obs_daily_grass)

    ## Only keep data from May 2009 and beyond
    obs_daily_grass$Date <- as.Date(obs_daily_grass$Date)
    obs_daily_grass <- obs_daily_grass[which(obs_daily_grass$Date >= "2009-05-01"),]
    obs_daily_grass <- obs_daily_grass[which(obs_daily_grass$Date <= "2024-12-31"),]
    
    ## Put the values as columns (dcast) per EP (5 different climate variable columns)
    # Need to adjust a colname first for next actions
    colnames(obs_daily_grass)[3] <- "variable"
    d_obs_daily_grass <- dcast(data = obs_daily_grass, formula = Date + year + month + day + dayOfYear ~ EP, value.var = "variable", fun.aggregate = mean)
    # dim(d_obs_daily_grass) ; summary(d_obs_daily_grass)

    ### IF VARIABLE OF INTEREST IS AIR TEMPERATURE (Ta-200) -> NEED TO PERFORM ADIABATIC LAPSE RATE CORRECTIONS!
    ### It’s the rate at which air temperature decreases with altitude due to decreasing atmospheric pressure (and expansion of rising air)
    ### Use 0.65°C/100m following Kim Weissig & al. (preprint from the Exploratories)

    ### To correct a temperature reading from a reference site (e.g., grassland) to match another site's elevation: 
    #               Tcorr = Tgrass − 0.0065×(elev.forest − elev.grassland)
    # If the forest is higher -> correction is positive (grassland temps are adjusted down to colder elevation)
    # If the grassland is higher -> correction is negative (adjusted up to warmer elevation)

    message(paste("\nComputing mean offsets of ",stat," ",var," for ",EP, sep = ""))

    if( var == "Ta_200" ) {

        # Sanity check: d_obs_daily_grass should have the same nrows as obs_daily_forest
        if( nrow(obs_daily_forest) != nrow(d_obs_daily_grass) ) {
            stop(
                paste("!!! Non matching dimensions between forest EP and grasslands EP for ",paste(var,stat, sep = "_")," @ ",EP," !!!", sep = "")
            )
        } # eo if loop - sanity check

        # Join forest Ta_200 and grasslands Ta_200
        d_obs_daily_grass$temp_forest <- obs_daily_forest$variable
        
        # Define the elevation fo the target forest EP as well as the elevation difference for each grassland EP
        elev_forest <- unique(df_subset$EP_elevation) # elev_forest
        elev_diff <- elev_forest - df_subset[1:5,"EP2_elevation"]
        names(elev_diff) <- df_subset[1:5,"EP2"]

        # Calculate Tcorr for each grassland EP
        d_obs_daily_grass[,paste("Tcorr", names(elev_diff)[1], sep = "_")] <- (d_obs_daily_grass[,names(elev_diff)[1]]) - 0.0065*(elev_diff[1])
        d_obs_daily_grass[,paste("Tcorr", names(elev_diff)[2], sep = "_")] <- (d_obs_daily_grass[,names(elev_diff)[2]]) - 0.0065*(elev_diff[2])
        d_obs_daily_grass[,paste("Tcorr", names(elev_diff)[3], sep = "_")] <- (d_obs_daily_grass[,names(elev_diff)[3]]) - 0.0065*(elev_diff[3])
        d_obs_daily_grass[,paste("Tcorr", names(elev_diff)[4], sep = "_")] <- (d_obs_daily_grass[,names(elev_diff)[4]]) - 0.0065*(elev_diff[4])
        d_obs_daily_grass[,paste("Tcorr", names(elev_diff)[5], sep = "_")] <- (d_obs_daily_grass[,names(elev_diff)[5]]) - 0.0065*(elev_diff[5])

        # Drop uncorrected columns (columns 6:10)
        d_obs_daily_grass <- d_obs_daily_grass[,c(1:5,11:length(d_obs_daily_grass))] # head(d2_obs_daily_grass) ; summary(d2_obs_daily_grass)
        # Adjust colnames
        colnames(d_obs_daily_grass)[c(7:11)] <- names(elev_diff)

        ### Compute offsets between each pairs of EPs: Tgrass - Tforest
        d_obs_daily_grass[,paste("Offset", names(elev_diff)[1], sep = "_")] <- d_obs_daily_grass[,names(elev_diff)[1]] - d_obs_daily_grass$temp_forest
        d_obs_daily_grass[,paste("Offset", names(elev_diff)[2], sep = "_")] <- d_obs_daily_grass[,names(elev_diff)[2]] - d_obs_daily_grass$temp_forest
        d_obs_daily_grass[,paste("Offset", names(elev_diff)[3], sep = "_")] <- d_obs_daily_grass[,names(elev_diff)[3]] - d_obs_daily_grass$temp_forest
        d_obs_daily_grass[,paste("Offset", names(elev_diff)[4], sep = "_")] <- d_obs_daily_grass[,names(elev_diff)[4]] - d_obs_daily_grass$temp_forest
        d_obs_daily_grass[,paste("Offset", names(elev_diff)[5], sep = "_")] <- d_obs_daily_grass[,names(elev_diff)[5]] - d_obs_daily_grass$temp_forest
        # summary(d_obs_daily_grass)

    } else {

        # Join forest Ta_200 and grasslands Ta_200
        d_obs_daily_grass$value_forest <- obs_daily_forest$variable
        # summary(d_obs_daily_grass)

        names <- as.character(df_subset[1:5,"EP2"])

        ### Compute offsets between each pairs of EPs: Tgrass - Tforest
        d_obs_daily_grass[,paste("Offset", names[1], sep = "_")] <- d_obs_daily_grass[,names[1]] - d_obs_daily_grass$value_forest
        d_obs_daily_grass[,paste("Offset", names[2], sep = "_")] <- d_obs_daily_grass[,names[2]] - d_obs_daily_grass$value_forest
        d_obs_daily_grass[,paste("Offset", names[3], sep = "_")] <- d_obs_daily_grass[,names[3]] - d_obs_daily_grass$value_forest
        d_obs_daily_grass[,paste("Offset", names[4], sep = "_")] <- d_obs_daily_grass[,names[4]] - d_obs_daily_grass$value_forest
        d_obs_daily_grass[,paste("Offset", names[5], sep = "_")] <- d_obs_daily_grass[,names[5]] - d_obs_daily_grass$value_forest
        # summary(d_obs_daily_grass)

    } # eo if else loop - Ta_200 or not
    
    ## Define the 5 weights to be used in the weighted average based on the distance to the forest EP in 'df_subset'
    weights <- 1 / df_subset[1:5,"Distance_to_EP2_km"] # weights

    ### Compute weighted mean offset
    d_obs_daily_grass <- data.frame(
        d_obs_daily_grass %>%
        rowwise() %>%
        mutate(
            Mean_Grassland_Value = weighted.mean(c_across(7:11), w = weights),
            Mean_Offset = weighted.mean(c_across(12:16), w = weights)
        ) %>%
        ungroup()
    ) # eo ddf
    # summary(d_obs_daily_grass)

    ### Looks good. Merge with data from 'df_subset'
    # First, put Region and EP as first columns
    d_obs_daily_grass$Region <- region
    d_obs_daily_grass$EP <- EP
    d_obs_daily_grass <- d_obs_daily_grass %>% relocate(Region)
    d_obs_daily_grass <- d_obs_daily_grass %>% relocate(EP, .after = Region)
    # head(d_obs_daily_grass) ; head(df_subset)

    if( var == "Ta_200" ) {
        
        ## Adjust some colnames too
        colnames(d_obs_daily_grass)[c(4:8)] <- c("Year","Month","Day","DOY",paste(stat,var, sep = "_"))
        colnames(df_subset)[1] <- "Region"

        ## Merge by EP
        t <- merge(df_subset[,c(1,2,5:7)], d_obs_daily_grass[,c(2:8,length(d_obs_daily_grass)-1,length(d_obs_daily_grass))], by = "EP")
        # str(t); dim(t)

    } else {

        ## Adjust some colnames too
        colnames(d_obs_daily_grass)[c(4:7,13)] <- c("Year","Month","Day","DOY",paste(stat,var, sep = "_"))
        colnames(df_subset)[1] <- "Region"

        ## Merge by EP
        t <- merge(df_subset[,c(1,2,5:7)], d_obs_daily_grass[,c(2:7,13,length(d_obs_daily_grass)-1,length(d_obs_daily_grass))], by = "EP")

    } # eo if else loop - var

    ## Re-provide mean stand age - got lost in th process
    Age_2012 <- unique(df_subset[df_subset$Year == 2012,"Age"]) # Age_2012
    t$Mean_stand_age <- Age_2012 + (t$Year - 2012)
    # summary(t) ; head(t)

    ## Useless message to finish
    message(paste("Saving mean offsets of ",stat," ",var," for ",EP,"\n", sep = ""))
    
    ## Return 't'
    setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling")
    saveRDS(t, file = paste0("table_mean_offsets+metadata_",stat,"_",var,"_",EP, "_28_10_25.rds"))
    
    # & clean
    rm(t,Age_2012,df_subset,d_obs_daily_grass,weights,obs_daily_grass,obs_daily_forest,file,file2)
    gc()

} # eo master FUN - compute_table_offsets

## Simply modify the compute_table_offsets() FUN to save individual files in the 'daily_offsets_for_microlimate_modelling' 
## directory through for loops - should not lead to swapping
forest_ids <- unique(df$EP)
# For testing:
# p <- "AEW03"
# v <- "Ts_20"
# s <- "min"

for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")) {
    for(s in c("max","min")) {
        
      mclapply(forest_ids, function(p) { 
                  compute_table_offsets(EP = p, var = "precipitation", stat = "total")
            }, mc.cores = 20
      ) # eo mclapply

    } # eo 2nd for loop
} # eo 1st for loop

### 28/10/25: Had to re-run script to save mean grassland values on top of the offset and forest values

### Deleting files in the dir that are from the 30/10
setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling"); dir()
# dir()[grep("30_06_25",dir())]
file.remove(dir()[grep("30_06_25",dir())])

### 30/10/25: Had to re-run script to add precipitation

### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------