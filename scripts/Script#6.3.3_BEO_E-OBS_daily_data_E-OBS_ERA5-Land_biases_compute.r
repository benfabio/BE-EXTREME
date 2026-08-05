### ------------------------------------------------------------------------------------------------------------

### 28/03/25 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to compute the biases in the daily climate data from E-OBS and ERA5 Land (model - obs)

### R script to create a function that:
###  - Loads daily data from E-OBS (Script#6.0) & ERA Land (Script#5.2) and the daily local Exploratories
###    climate data (Script#2.5)
###  - Computes biases in daily stat of each variable fo each EP, date and region
###  - Returns these biases in a table. May save them in: 
###    /home/fbenedetti/ERA5-Land-DEU-processed/daily/biases/ 
###    (mix all files in the same directory)

### Last update: 23/04/25 (Analyzing temporal patterns of SM_10 vs. volumetric_soil_water_layer_2 biases)

### ------------------------------------------------------------------------------------------------------------

# Libraries 
library("dplyr")
library("data.table")
library("ggplot2")
library("purrr")
library("reshape2")
library("lubridate")
library("parallel")

### ------------------------------------------------------------------------------------------------------------

### Master FUN - compute_daily_biases

compute_daily_biases <- function(region, system, var, stat) {

        #' This function takes four arguments and returns a formatted data.frame:
        #' @param region Region name (character): 'SCH' or 'HND' or 'SWA'
        #' @param system System within the region of interest (character): 'grasslands' or 'forests'
        #' @param var the climate variable to process (character) - one of the following: 
        #' "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20", "Ts_50", "precipitation" or "SM_10"
        #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character):
        #' 'mean', 'max' or 'min' or 'total' for precipitation
        #' @return A formatted data.frame combining the daily biases - saved as .Rdata file in 'biases' dir

        ## Go to local obs directory and extract TS of interest
        setwd(paste("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily/",system,"s/",region, sep = ""))

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

        # Message - depends on variable
        if( var %in% c("Ta_200","precipitation") ) {

            ## Go to E-OBS and ERA5-Land dirs and load their daily data too
            message(paste("Computing biases in daily ",stat," ",var," for the ",system," of the ",region," for E-OBS and ERA5-Land data", sep = ""))
            
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
                plots2subset <- unique(EOBS_daily_stat$EP)
                sub_EOBS_daily_stat <- EOBS_daily_stat[EOBS_daily_stat$EP %in% plots2subset,]
                rm(EOBS_daily_stat,plots2subset); gc()
            } else if( region == "HND" ) {
                plots2subset <- unique(EOBS_daily_stat$EP)
                sub_EOBS_daily_stat <- EOBS_daily_stat[EOBS_daily_stat$EP %in% plots2subset,]
                rm(EOBS_daily_stat,plots2subset); gc()
            } else if( region == "SWA" ) {
                plots2subset <- unique(EOBS_daily_stat$EP)
                sub_EOBS_daily_stat <- EOBS_daily_stat[EOBS_daily_stat$EP %in% plots2subset,]
                rm(EOBS_daily_stat,plots2subset); gc()
            } # eo if else loop - var
            
            # Re-name to avoid more of else loops later
            colnames(sub_EOBS_daily_stat) <- c("EP","Date","value")

            ## Go to ERA5-Land dir and load corresponding dataset
            setwd(paste("/home/fbenedetti/ERA5-Land-DEU-processed/daily/",system,"/",region, sep = ""))

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
              era5_var <- "volumetric_soil_water_layer_1"
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
            rm(era5_data,era5_files)
            gc()

            # Keep 'stat' of interest (discard the others)
            era5_ddf <- era5_ddf[,c("region","system","plot","date","day","month","year",stat)]

            # For each EP: match the daily data that overlap across all 3 sources and combine in a single ddf
            colnames(obs_daily_stat)[3] <- paste(var,stat, sep = "_")
            colnames(sub_EOBS_daily_stat)[3] <- paste(var,stat, sep = "_")
            colnames(era5_ddf)[c(3,4,length(era5_ddf))] <- c("EP","Date",paste(var,stat, sep = "_"))
            
            names <-  c("EP","Date",paste(var,stat, sep = "_")) # vector of colnames to join the 3 ddf by

            # To make sure date format is homogeneous across all tables
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
            rm(obs_daily_stat,sub_EOBS_daily_stat,era5_ddf)
            gc()

            # Adjust colnames
            colnames(merged_df)[c(3:5)] <- c("obs","E_OBS","ERA5_Land") 
            
            ### Calculate E-OBS & ERA5-Land biases in daily data for each EP 
            merged_df$biases_EOBS <- merged_df$E_OBS - merged_df$obs
            merged_df$biases_ERA5Land <- merged_df$ERA5_Land - merged_df$obs

            # Add system in 1st position
            merged_df <- merged_df %>% mutate(system = system) %>% relocate(system, .before = 1)

            # Save in dir, clean and return to wd
            message(paste("Saving biases table\n", sep = ""))
            setwd(paste("/home/fbenedetti/ERA5-Land-DEU-processed/daily/biases/", sep = ""))
            save(x = merged_df, file = paste("table_biases_daily_",stat,"_",var,"_",system,"_",region,"_",Sys.Date(),".Rdata", sep = ""))
            rm(merged_df)
            gc()

        } else {

            message(paste("Computing biases in daily ",stat," ",var," for the ",system," of the ",region," ERA5-Land data", sep = ""))

            ## Go to ERA5-Land dir and load corresponding dataset
            setwd(paste("/home/fbenedetti/ERA5-Land-DEU-processed/daily/",system,"/",region, sep = ""))
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
              era5_var <- "volumetric_soil_water_layer_2"  # 23/04/25: WARNING - modified volumetric_soil_water_layer_1 to volumetric_soil_water_layer_2 here to compare
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
            rm(era5_data,era5_files)
            gc()

            era5_ddf <- era5_ddf[,c("region","system","plot","date","day","month","year",stat)]

            # For each EP: match the daily data that overlap across both data sources and combine in a single ddf
            colnames(obs_daily_stat)[3] <- paste(var,stat, sep = "_")
            colnames(era5_ddf)[c(3,4,length(era5_ddf))] <- c("EP","Date",paste(var,stat, sep = "_"))
            names <-  c("EP","Date",paste(var,stat, sep = "_")) # vector of colnames to join the 3 ddf by

            # To make sure date format is homogeneous across both tables
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

            # Compute biases and return
            merged_df$biases_ERA5Land <- merged_df$ERA5_Land - merged_df$obs
            # > 0 --> ERA5-Land warmer/more humid than obs
            # < 0 --> ERA5-Land colder/less humid than obs

            # Add system in 1st position
            merged_df <- merged_df %>% mutate(system = system) %>% relocate(system, .before = 1)

            # Save in dir, clean and return to wd
            message(paste("Saving biases table\n", sep = ""))
            setwd(paste("/home/fbenedetti/ERA5-Land-DEU-processed/daily/biases/", sep = ""))
            save(x = merged_df, file = paste("table_biases_daily_",stat,"_",var,"_",system,"_",region,"_",Sys.Date(),".Rdata", sep = ""))
            rm(merged_df); gc()

        } # eo if else loop

} # eo master FUN - compute_daily_biases

### Running compute_daily_biases on the following variables: 

for(s in c("grassland","forest")) {
    for(r in c("HND","SCH","SWA")) {
        for(stat in c("max","min")) {
            compute_daily_biases(region = r, system = s, stat = stat, var = "SM_10")
        } # eo for loop - stat
    } # eo for loop - r
} # eo for loop - s

### ------------------------------------------------------------------------------------------------------------

### 07/04/25: Analyzing outputs of the FUN above (plots etc.) - FOR SM_10
### Check if there are temporal patterns in biases: interannual, seasonal, monthly etc.

### Focus on grassland plots only for now.
### Note: total precip. and Ta_200 have two biases: E-OBS & ERA5-Land

### 15/04/25: Re-plotting the biases for daily total precipitation but without ERA5-Land data 
### because those were so much higher one could not see the E-OBS biases well enough

### 16/04/25: Plotting the biases for Ts_20

### 23/04/25: Plotting the biases for SM_10 again, but against level_2 ERA5-Land data

### Create a simple FUN to read all grassland biases files of one variable
read_daily_biases <- function(var,stat) {
      # Identify the files 
      # Re-write 'stat' argument for precip. data
      if(var == "precipitation") {
          stat <- "total"
      } # eo if loop
      setwd(paste("/home/fbenedetti/ERA5-Land-DEU-processed/daily/biases/", sep = ""))
      files2keep <- dir()[grepl(paste(c(stat,var,"grassland"), collapse = "_"), dir())]
      # Load and rbind them
      res <- lapply(files2keep, function(f) {
                d <- get(load(f))
                # Add region as a factor
                if( grepl("HND",f) ) {
                    d$region <- "HND"
                } else if( grepl("SCH",f) ) {
                    d$region <- "SCH"
                } else if( grepl("SWA",f) ) {
                    d$region <- "SWA"
                } # eo if else loop - region's name
                # extract M/Y from 'Date' vector
                d$month <- lubridate::month(d$Date)
                d$year <- lubridate::year(d$Date)
                return(d)} 
      ) # eo lapply - files2keep
      # Rbind all files into one table
      data <- dplyr::bind_rows(res)
      rm(res); gc()
      # Return
      return(data)
} # eo FUN - read_daily_biases

### Per variable/stat examine biases distribution through time

### 1. Ta_200
##    1.A. max Ta_200
d <- na.omit( read_daily_biases(var = "Ta_200", stat = "max") )
## BEWARE: For Ta_200 & precip.: melt() to have data source (E-OBS vs. ERA5-Land) in a column
m.d <- melt(d, id.vars = c("EP","system","region","Date","month","year","obs","E_OBS","ERA5_Land"))
# head(m.d); str(m.d)
colnames(m.d)[c(10,11)] <- c("source","bias")
m.d$source <- as.character(m.d$source)
m.d[m.d$source == "biases_EOBS","source"] <- "E_OBS"
m.d[m.d$source == "biases_ERA5Land","source"] <- "ERA5_Land"

# Compute and plot yearly anomalies with boxplots (different colours per data source, facet_wrap())
yearly_biases <- m.d %>% group_by(source,EP,region,year) %>% summarize(mean_bias = mean(bias, na.rm = T))
# str(yearly_biases); summary(yearly_biases)
# Plot biases
colors_sources <- c("obs" = "#00AFBB", "E_OBS" = "#e31a1c", "ERA5_Land" = "#E7B800")
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias, fill = factor(source))) +
    geom_boxplot(colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    scale_fill_manual(values = colors_sources) + labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- m.d %>% group_by(source,EP,region,month) %>% summarize(mean_bias = mean(bias, na.rm = T))
# Plot biases
colors_sources <- c("obs" = "#00AFBB", "E_OBS" = "#e31a1c", "ERA5_Land" = "#E7B800")
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias, fill = factor(source))) +
    geom_boxplot(colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    scale_fill_manual(values = colors_sources) + labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 10, height = 3.5)  
# Good. looks like E-OBS is better in this case :) 

### -----------------------------------------------

##    1.B. min Ta_200
var <- "Ta_200"
stat <- "min"
d <- na.omit( read_daily_biases(var = var, stat = stat) )
# melt 
m.d <- melt(d, id.vars = c("EP","system","region","Date","month","year","obs","E_OBS","ERA5_Land"))
colnames(m.d)[c(10,11)] <- c("source","bias")
m.d$source <- as.character(m.d$source)
m.d[m.d$source == "biases_EOBS","source"] <- "E_OBS"
m.d[m.d$source == "biases_ERA5Land","source"] <- "ERA5_Land"
# Compute and plot yearly anomalies with boxplots (different colours per data source, facet_wrap())
yearly_biases <- m.d %>% group_by(source,EP,region,year) %>% summarize(mean_bias = mean(bias, na.rm = T))
# str(yearly_biases); summary(yearly_biases)
# Plot biases
colors_sources <- c("obs" = "#00AFBB", "E_OBS" = "#e31a1c", "ERA5_Land" = "#E7B800")
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias, fill = factor(source))) +
    geom_boxplot(colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    scale_fill_manual(values = colors_sources) + labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- m.d %>% group_by(source,EP,region,month) %>% summarize(mean_bias = mean(bias, na.rm = T))
# Plot biases
colors_sources <- c("obs" = "#00AFBB", "E_OBS" = "#e31a1c", "ERA5_Land" = "#E7B800")
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias, fill = factor(source))) +
    geom_boxplot(colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    scale_fill_manual(values = colors_sources) + labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 10, height = 3.5)  
## Weaker biases than with max Ta_200...However, considering the monthly biases, it looks like E-OBS again is the better one :) 

### -----------------------------------------------

### 2. Ta_10
##    2.A. max Ta_10
var <- "Ta_10"
stat <- "max"
d <- na.omit( read_daily_biases(var = var, stat = stat) )
# Compute and plot yearly anomalies with boxplots
yearly_biases <- d %>% group_by(EP,region,year) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# str(yearly_biases); summary(yearly_biases)
# Plot biases
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias)) +
    geom_boxplot(fill = "#E7B800", colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- d %>% group_by(EP,region,month) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias)) +
    geom_boxplot(colour = "black", fill = "#E7B800") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 10, height = 3.5)  
### ERA5-Land usually too cold (biases mostly negative); Weak interannulaity (2008 is off, as expected)
### Strong seasonality though, with stronger negative biases in summer (expected since ERA5-Land too cold)

### -----------------------------------------------

##    2.B. min Ta_10
var <- "Ta_10"
stat <- "min"
d <- na.omit( read_daily_biases(var = var, stat = stat) )
# Compute and plot yearly anomalies with boxplots
yearly_biases <- d %>% group_by(EP,region,year) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias)) +
    geom_boxplot(fill = "#E7B800", colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- d %>% group_by(EP,region,month) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias)) +
    geom_boxplot(colour = "black", fill = "#E7B800") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 10, height = 3.5)   
### Interesting: Quite somme interannuality in min Ta and this interannuality varies between regions! 
### Seasonality is evident in all regions but is again stronger in the SCH and the SWA relative to the HND.
### ERA5-Land too warm in the summer of the SCH (mean biases > 0) but too cold in winter in most regions

### -----------------------------------------------

### 3. total precipitation
var <- "precipitation"
stat <- "total"
d <- na.omit( read_daily_biases(var = var, stat = stat) )
## BEWARE: For Ta_200 & precip.: melt() to have data source (E-OBS vs. ERA5-Land) in a column
m.d <- melt(d, id.vars = c("EP","system","region","Date","month","year","obs","E_OBS","ERA5_Land"))
# head(m.d); str(m.d)
colnames(m.d)[c(10,11)] <- c("source","bias")
m.d$source <- as.character(m.d$source)
m.d[m.d$source == "biases_EOBS","source"] <- "E_OBS"
m.d[m.d$source == "biases_ERA5Land","source"] <- "ERA5_Land"

# Compute and plot yearly anomalies with boxplots (different colours per data source, facet_wrap())
yearly_biases <- m.d %>% group_by(source,EP,region,year) %>% summarize(mean_bias = mean(bias, na.rm = T))
# str(yearly_biases); summary(yearly_biases)
# Plot biases
colors_sources <- c("obs" = "#00AFBB", "E_OBS" = "#e31a1c", "ERA5_Land" = "#E7B800")
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias, fill = factor(source))) +
    geom_boxplot(colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    scale_fill_manual(values = colors_sources) + labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- m.d %>% group_by(source,EP,region,month) %>% summarize(mean_bias = mean(bias, na.rm = T))
# Plot biases
colors_sources <- c("obs" = "#00AFBB", "E_OBS" = "#e31a1c", "ERA5_Land" = "#E7B800")
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias, fill = factor(source))) +
    geom_boxplot(colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    scale_fill_manual(values = colors_sources) + labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 10, height = 3.5)  
### Simply: MUCH stronger biases in ERA5 precipitation compared to E-OBS
### (ERA5-Land overpredicts daily total precip.) 


### 15/04/25: Re-plotting the biases for daily total precipitation but without ERA5-Land data
yearly_biases <- m.d %>% group_by(source,EP,region,year) %>% summarize(mean_bias = mean(bias, na.rm = T))
# str(yearly_biases)
# Plot 
colors_sources <- c("obs" = "#00AFBB", "E_OBS" = "#e31a1c", "ERA5_Land" = "#E7B800")
plot <- ggplot(data = yearly_biases[yearly_biases$source == "E_OBS",], aes(x = factor(year), y = mean_bias, fill = factor(source))) +
    geom_boxplot(colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    scale_fill_manual(values = colors_sources) + labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,"_E-OBS_only.jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- m.d %>% group_by(source,EP,region,month) %>% summarize(mean_bias = mean(bias, na.rm = T))
# Plot biases
colors_sources <- c("obs" = "#00AFBB", "E_OBS" = "#e31a1c", "ERA5_Land" = "#E7B800")
plot <- ggplot(data = monthly_biases[monthly_biases$source == "E_OBS",], aes(x = factor(month), y = mean_bias, fill = factor(source))) +
    geom_boxplot(colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    scale_fill_manual(values = colors_sources) + labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,"_E-OBS_only.jpg",sep = ""), dpi = 300, width = 10, height = 3.5) 

### For E-OBS only data
## --> both monthly and yearly biases range beteen +1mm and -1mm, so very small biases 
## --> seasonal and annual biases in the same range, none really stronger than the other
## --> precipitation biases are weaker in the SCH compared to the HND and the SWA
##     they're actually very small in the SCH
## --> most prominent feature in the seasonal variation with E-OBS too dry in summer ()= negative biases)
##     and too wet in winter (= positive anomalies)


### -----------------------------------------------

### 4. Ts_05
##    4.A. max Ts_05
var <- "Ts_05"
stat <- "max"
d <- na.omit( read_daily_biases(var = var, stat = stat) )
# Compute and plot yearly anomalies with boxplots
yearly_biases <- d %>% group_by(EP,region,year) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias)) +
    geom_boxplot(fill = "#E7B800", colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- d %>% group_by(EP,region,month) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias)) +
    geom_boxplot(colour = "black", fill = "#E7B800") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 10, height = 3.5)  
### Some interannuality (2008 off in the SWA). Looks like mean biases increase a bit in the ast years for the SCH ad HND but not SWA
### Higher seasonality in biases in the SCH as always. ERA5-Land is too warm in summer and a bit too cold in winter. 
### Looks like seasoanlity > interannuality. Interannuality quite stable through time I'd say.  

### -----------------------------------------------

##    4.B. min Ts_05
var <- "Ts_05"
stat <- "min"
d <- na.omit( read_daily_biases(var = var, stat = stat) )
# Compute and plot yearly anomalies with boxplots
yearly_biases <- d %>% group_by(EP,region,year) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias)) +
    geom_boxplot(fill = "#E7B800", colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- d %>% group_by(EP,region,month) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias)) +
    geom_boxplot(colour = "black", fill = "#E7B800") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 10, height = 3.5)  
# Similar to max Ts_05: some small interannuality that is < seasonality.
# Seasonality weaker in the SCH for once. ERA5-Land seems too cold in general for min Ts. biases are mostly < 0.

### -----------------------------------------------

### 4. Ts_10
##    4.A. max Ts_10
var <- "Ts_10"
stat <- "max"
d <- na.omit( read_daily_biases(var = var, stat = stat) )
# Compute and plot yearly anomalies with boxplots
yearly_biases <- d %>% group_by(EP,region,year) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias)) +
    geom_boxplot(fill = "#E7B800", colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- d %>% group_by(EP,region,month) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias)) +
    geom_boxplot(colour = "black", fill = "#E7B800") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 10, height = 3.5)  
### VERY similar biases patterns as Ts_05 obviously (as expected).
### ERA5-Land too warm so most biases > 0. Interannuality > seasonality. 
### Stronger positve biases in summer in all regions. 

### -----------------------------------------------

##    4.B. min Ts_10
var <- "Ts_10"
stat <- "min"
d <- na.omit( read_daily_biases(var = var, stat = stat) )
# Compute and plot yearly anomalies with boxplots
yearly_biases <- d %>% group_by(EP,region,year) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias)) +
    geom_boxplot(fill = "#E7B800", colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- d %>% group_by(EP,region,month) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias)) +
    geom_boxplot(colour = "black", fill = "#E7B800") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 10, height = 3.5)  


### -----------------------------------------------

## 5.A. max SM_10
var <- "SM_10"
stat <- "max"
d <- na.omit( read_daily_biases(var = var, stat = stat) )
# Compute and plot yearly anomalies with boxplots
yearly_biases <- d %>% group_by(EP,region,year) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias)) +
    geom_boxplot(fill = "#E7B800", colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- d %>% group_by(EP,region,month) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias)) +
    geom_boxplot(colour = "black", fill = "#E7B800") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 10, height = 3.5)  


### -----------------------------------------------

## 5.B. min SM_10
var <- "SM_10"
stat <- "min"
d <- na.omit( read_daily_biases(var = var, stat = stat) )
# Compute and plot yearly anomalies with boxplots
yearly_biases <- d %>% group_by(EP,region,year) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias)) +
    geom_boxplot(fill = "#E7B800", colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- d %>% group_by(EP,region,month) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias)) +
    geom_boxplot(colour = "black", fill = "#E7B800") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 10, height = 3.5)  


### -----------------------------------------------

### 6. Ts_20
##    6.A. max Ts_20
var <- "Ts_20"
stat <- "max"
d <- na.omit( read_daily_biases(var = var, stat = stat) )
# Compute and plot yearly anomalies with boxplots
yearly_biases <- d %>% group_by(EP,region,year) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias)) +
    geom_boxplot(fill = "#E7B800", colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- d %>% group_by(EP,region,month) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias)) +
    geom_boxplot(colour = "black", fill = "#E7B800") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 10, height = 3.5)  
### 

### -----------------------------------------------

##    6.B. min Ts_20
var <- "Ts_20"
stat <- "min"
d <- na.omit( read_daily_biases(var = var, stat = stat) )
# Compute and plot yearly anomalies with boxplots
yearly_biases <- d %>% group_by(EP,region,year) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias)) +
    geom_boxplot(fill = "#E7B800", colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- d %>% group_by(EP,region,month) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias)) +
    geom_boxplot(colour = "black", fill = "#E7B800") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,".jpg",sep = ""), dpi = 300, width = 10, height = 3.5)  

### -----------------------------------------------

### 7. SM_10 v2: against level_2 ERA5 Land data

### First, modify the read_daily_biases() to only load files from the 23/04/25
# var <- "SM_10"
# stat <- "max"
read_daily_biases <- function(var,stat) {
      # Identify the files 
      # Re-write 'stat' argument for precip. data
      if(var == "precipitation") {
          stat <- "total"
      } # eo if loop
      setwd(paste("/home/fbenedetti/ERA5-Land-DEU-processed/daily/biases/", sep = ""))
      files2keep <- dir()[grepl(paste(c(stat,var,"grassland"), collapse = "_"), dir())]
      # Only keep those from the '2025-04-23'
      files2keep <- files2keep[grepl(paste("2025-04-23", collapse = "_"), files2keep)]
      # Load and rbind them
      res <- lapply(files2keep, function(f) {
                d <- get(load(f))
                # Add region as a factor
                if( grepl("HND",f) ) {
                    d$region <- "HND"
                } else if( grepl("SCH",f) ) {
                    d$region <- "SCH"
                } else if( grepl("SWA",f) ) {
                    d$region <- "SWA"
                } # eo if else loop - region's name
                # extract M/Y from 'Date' vector
                d$month <- lubridate::month(d$Date)
                d$year <- lubridate::year(d$Date)
                return(d)} 
      ) # eo lapply - files2keep
      # Rbind all files into one table
      data <- dplyr::bind_rows(res)
      rm(res); gc()
      # Return
      return(data)
} # eo FUN - read_daily_biases

## 7.A. max SM_10
var <- "SM_10"
stat <- "max"
d <- na.omit( read_daily_biases(var = var, stat = stat) )
# Compute and plot yearly anomalies with boxplots
yearly_biases <- d %>% group_by(EP,region,year) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias)) +
    geom_boxplot(fill = "#E7B800", colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,"_level2.jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- d %>% group_by(EP,region,month) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias)) +
    geom_boxplot(colour = "black", fill = "#E7B800") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,"_level2.jpg",sep = ""), dpi = 300, width = 10, height = 3.5)  


### -----------------------------------------------

## 7.B. min SM_10
var <- "SM_10"
stat <- "min"
d <- na.omit( read_daily_biases(var = var, stat = stat) )
# Compute and plot yearly anomalies with boxplots
yearly_biases <- d %>% group_by(EP,region,year) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = yearly_biases, aes(x = factor(year), y = mean_bias)) +
    geom_boxplot(fill = "#E7B800", colour = "black") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," annual biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_annual_biases_",stat,"_",var,"_level2.jpg",sep = ""), dpi = 300, width = 13, height = 3.5)  

# Compute and plot monthly anomalies with boxplots (different colours per data source, facet_wrap())
monthly_biases <- d %>% group_by(EP,region,month) %>% summarize(mean_bias = mean(biases_ERA5Land, na.rm = T))
# Plot biases
plot <- ggplot(data = monthly_biases, aes(x = factor(month), y = mean_bias)) +
    geom_boxplot(colour = "black", fill = "#E7B800") + geom_hline(yintercept = 0, linetype = "dashed") + 
    labs(y = paste("Biases in ",stat,var,sep = " "), x = "") + 
    theme_bw() + theme(
        legend.position = "bottom",
        axis.text.x = element_text(size = 6)) +
    facet_wrap(.~factor(region)) + 
    ggtitle(paste("Distribution of ",stat," ",var," monthly biases", sep = ""))
# Save
ggsave(plot = plot, filename = paste("boxplot_monthly_biases_",stat,"_",var,"_level2.jpg",sep = ""), dpi = 300, width = 10, height = 3.5)  


### ------------------------------------------------------------------------------------------------------------

### 31/03/2025: MAIN CONCLUSIONS from the plots made above: E-OBS > ERA5-Land. Quite some regional variability
### and interannuality; but seasonality stronger than interannuality overall. Especially in summer months.
### Should use EP-specific, monthly quantile mapping.

## max Ta_200: E-OBS has lower biases than ERA5-Land. Should keep E-OBS for air temperature.
##             ERA5-Land colder than obs. E-OBS warmer than obs. 

## min Ta_200: E-OBS has lower biases than ERA5-Land again.Keep E-OBS if possible.

## total precipitation: ERA5-Land shows MUCH higher (positive, too wet) biases than E-OBS. Keep E-OBS for sure. 

## max Ta_10: ERA5-Land colder than obs. but usually by less than 2°C.
## Strong seasonality: much higher (negative) biases in summer. Quite close to 0 in winter.
## ERA5-Land good for winter max Ta_10

## min Ta_10: Mostly biases < 0 except in the SCH! Quite some interannuality.
## Biases (< 0) stronger in winter for the HND and SWA. But higher (and > 0) in summer for the SCH.

## max Ts_05: Mostly positive biases (ERA5-Land warmer soil temperatures than obs. contrary to some other metrics like ait temp.)
## Lower biases in the SWA compared to HND and especially SCH. Stronger positive biases in summer. 

## min Ts_05: Mostly negative biases (but hard to interpret because of potential double negatives). 
## Seasonality > interannuality. Stronger biases in summer.

## max Ts_10: Mostly positive biases again. ERA5-Land warmer soil temperatures than obs.
## Same temporal patterns as Ts_05. Stronger seasonality in biases than interannuality. 
## Especially in summer. 

## min Ts_10: Same as min Ts_05.

## max Ts_20: Annual biases strcuture very similar to Ts_10, but weaker than Ts_10 in general.
## Ts_20 became too cold in the HND compared to obs while Ts_10 was too warm. 
## Monthly speaking, patterns also VERY similar to Ts_10 but < 2°C. Somehiw, Ts_20 seems easier to model than Ts_10.

## min Ts_20: closer to 0 than min Ts_05 & Ts_10 (deeper in the soil is easier to model by ERA5-Land?). 
## monthly patterns of biases are VERY different compared to those patterns of Ts_10 but also closer to 0 as above


### 07/04/25
## max SM_10: Mostly positive biases in SM_10 for the SCH and the SWA (ERA5-Land more humid than obs)
## and more variable for HND (SM_10 biases also lower in the HND than in the other 2 regions). 
## Biases surprisingly not too strong (-10%/+20%). Biases usually between +10% and -5% in the HND, +15%/-10% in the SCH
## and nearly zlways positive (= ERA5-Land too humid) in the SWA (between +20% and +5%). 
## Seasonality > interannuality in the HND (too humid in summer, too dry in winter). Similar temporal pettern in the SCH 
## and the SWA but just ALWAYS > 0 (always too humid).

## min SM_10: Extremely similar patterns as max SM_10 :-) Looks promising. ERA5-Land isn't as bad as I thought.


### 23/04/25: SM_10 vs. level_2 ERA5-Land data.

### Biases are VERY similar to level_1 data. level_2 data seem better in the HND but otherwise hard to tell.

### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------