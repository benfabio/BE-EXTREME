### ------------------------------------------------------------------------------------------------------------

### 12/05/25 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to extract the bias-corrected time series of daily temperature, moisture and precipitation (see Script#6.4.4) 
### and define daily climatologies and statistical thresholds for ECEs detection based on the chosen percentile (1st/5th/10th). 
### Then use the full time series of the E-OBS data to identify extremes and measure some of their features
### (i.e., duration, intensity, heterogeneity, frequency etc.). 

### NOTE: Run master FUN on all QM strategies (global, monthly, anoms and mw) and 3 percentiles to account for uncertainty
### NOTE: Make sure to enable the possibility to change the time period coverage of the reference period
###      (e.g., 1950-90, 1970-2000, etc.)
## NOTE: ECEs of min precipitation will have their own FUN

### Re-use R Script#6.2 to: 
###  - Make master FUN for detecting Xtremes and compute their features
###  - Run master FUN based on the bias-corrected daily TS and save outputs on clicmal

### Last update: 09/07/26 (Re-running find_extremes() for all forests and by removing those that were too young (< 10 yo) in the reference period)

### ------------------------------------------------------------------------------------------------------------

# Libraries 
library("dplyr")
library("data.table")
library("zoo")
library("reshape2")
require("scales")
library("lubridate")
library("parallel")

# Directory where to store ECEs tables
extremes_dir <- "/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/ECEs_tables"

### ------------------------------------------------------------------------------------------------------------

### Re-use the code from Script#6.2 to compute ECEs dynamics for every EP based on E-OBS data

# To test master FUN while you are writing it
# system = "forests"
# var = "precipitation"
# stat = "total"
# region = "SCH"
# method = "mw"
# perc = 0.99
# start = as.Date("1950-01-01")
# end = as.Date("1980-01-01")
# min_age = 15

### 12/02/26: Trouble shoot issue of extrmeley long (Duration > 500) ECE of max Ts_05 and max Ts_10 at AEW11 (see R Script#7.2)
### -> Add min_age argument to detect and remove those foretss that were too young/uncertain in the chosen baseline period

### 09/07/26: Re-running find_extremes() for all forests and by removing those that were too young (< 10 yo) in the reference period

find_extremes <- function(system, var, stat, region, method, perc, start, end, min_age = 10) {

        #' This function takes four arguments and returns a formatted data.frame:
        #' @param system the type of ecosystem (character): "grasslands" or "forests"
        #' @param var the climate variable to process (character): "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20", "Ts_50",
        #' "precipitation" or "SM_10"
        #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character): 'mean', 'max' or 'min'
        #' @param region Region name (character): 'SCH' or 'HND' or 'SWA'
        #' @param method Which quantile mapping approach to run (character): 'global' or 'monthly' or 'mw' or 'anomalies'
        #' @param perc Percentile to use from the refence period to define the statistical threshold (numeric): 
        #' 0 to 1 (0.90, 0.95, 0.99 for high ECEs such as heat events, or .01,.05,.1 for low events such as cold ones)
        #' @param start Date at which the reference basetine period should start (Date)
        #' @param end Date at which the reference basetine period should end (Date)
        #' @param min_age Age below which forets data should not be taken into account because StandAge is outside of GAMM
        #' traning range (10-40 years old) (integer) - Default == 10 because youngest stands fit into the GAMMs are 12 yo
        #' @return A formatted data.frame combining the daily statistics
    
        ## Message
        message(paste("\nDetecting extremes of ",paste(stat,var, sep = " ")," for the ",system," of the ",region,
            " based on the ",method," QM strategy and the ",perc," percentile\n", sep = ""))

        if( system == "grasslands" ) {
            setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/grasslands")
        } else (
            setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/forests")
        ) # eo if else loop
        
        ## Then, there will be a if else loop according to max/min below.
        ## So make sure precipitation has 'max' instead of 'min' or 'total'. 
        ## NOTE: ECEs of min precipitation will have their own FUN
        if( var == "precipitation" ) {
            stat <- "total"
        } # eo if else loop - precipitation
        
        ## Loading the dedicated file containing the full TS
        if( system == "grasslands" ) {
            d <- get(load(paste("table_combined_obs+corr_",method,"_",stat,"_",var,"_",system,"_",region,".Rdata", sep = "")))
            NULL
        } else (
            if( var == "precipitation" ) {
                d <- get(load(paste("table_combined_obs+corr_",method,"_",stat,"_",var,"_",system,"_",region,".Rdata", sep = "")))
                NULL
            } else {
                d <- readRDS(paste("table_combined_obs+corr_",stat,"_",var,"_",system,"_",region,".rds", sep = ""))
                NULL
            } # eo 2nd if loop
        ) # eo 1st if else loop

        ### WARNING: Only use 0.90, 0.95, 0.99 as perc for stat == 'max'
        ### Make sure to overwrite 'perc' if stat == "min" and vice-versa
        if( stat == "max" | var == "precipitation" ) {
            # Adjust perc
            if( perc == 0.1 ) {
                perc <- 0.9
            } 
            if( perc == 0.05 ) {
                perc <- 0.95
            }
            if( perc == 0.01 ) { 
                perc <- 0.99
            }
        } else if( stat == "min" ) {
            # Adjust perc in the opposite direction
            if( perc == 0.9 ) {
                perc <- 0.1
            }
            if( perc == 0.95 ) {
                perc <- 0.05
            }
            if( perc == 0.99 ) { 
                perc <- 0.01
            } # eo if else loop
        } # eo if else loop - perc

        ## In a mclappy, detect and quantify the ECEs at EP-level
        require("parallel")

        ### 22/01/26: For forests only, identify those EPs that do not span the whole 'start'/'end' period
        if( system == "forests" & var != "precipitation" ) {

            require("dplyr")
            message(paste("Removing EPs that do not have enough values for the ",start," to ",end," period", sep = ""))

            ### Remove data from stand ages < min_age
            d <- d[-which(d$Stand_age < min_age),]
            
            # Set the expected number of days based on 'start' and 'end'
            expected <- as.integer(end - start) + 1
            
            # Summarise coverage per EP
            coverage <- data.frame(
                d %>%
                filter(Date >= start, Date <= end) %>%
                group_by(EP) %>%
                summarise(
                        n_days_present = n_distinct(Date),
                        prop_coverage = n_days_present / expected,
                        first_date = min(Date),
                        last_date = max(Date),
                        .groups = "drop"
                    ) # eo summarise
                )
                # head(coverage) ; summary(coverage)

            # Check which EP exists at BOTH ends of the period
            ep2keep <- coverage %>% filter(first_date <= start, last_date >= end) %>% pull(EP)

        } else {
            
            ep2keep <- unique(d$EP)

        } # eo if loop - forests incomplete start/end period
        
        # For testing
        # Trouble shoot issue of extrmeley long (Duration > 500) ECE of max Ts_05 and max Ts_10 at AEW11
        # p <- ep2keep[1]; p
        # p = "AEW11"

        list_ECEs <- mclapply(X = ep2keep,
        
            FUN = function(p) {
                
                # Subset EP data
                message(paste("Computing ECEs for ",p,sep = ""))
                dd <- d[d$EP == p,]
                
                # Add D/M/Y
                dd$DOY <- lubridate::yday(dd$Date) # day of the month (1-31)
                dd$month <- lubridate::month(dd$Date)
                dd$year <- lubridate::year(dd$Date)
                # Subset reference period (start-end)
                baseline <- dd[which(dd$Date >= start & dd$Date <= end),]
                # summary(baseline)
                
                # Compute percentile thresholds per DOY
                clim <- baseline %>% group_by(DOY) %>% summarise(thresh = quantile(final_value, probs = perc, na.rm = TRUE) ) # summary(clim)
                # Join with full historical time series
                full_data <- dd %>% left_join(clim[,c("DOY","thresh")], by = "DOY") # summary(full_data)
                
                ### Adjust way of detecting ECEs and estimating their intensity depending on whether you are
                ### dealing with 'min' or 'max'
                if( stat == "max" | perc %in% c(0.9,0.95,0.99) ) {
                    
                    # Apply the threshold to identify heat extremes in the full time series (1950-2024)
                    full_data <- full_data %>% mutate(ECE = ifelse(final_value > thresh, 1, 0))
                    # Compute the intensity of ExHeEs
                    full_data <- full_data %>% mutate(Intensity = ifelse(ECE == 1, final_value - thresh, NA))

                } else if( stat == "min" | perc %in% c(0.1,0.05,0.01) ) {
                    
                    # Apply the threshold to identify cold extremes in the full time series (1950-2024)
                    full_data <- full_data %>% mutate(ECE = ifelse(final_value < thresh, 1, 0))
                    # Compute the intensity of ExHeEs
                    full_data <- full_data %>% mutate(Intensity = ifelse(ECE == 1, thresh - final_value, NA))

                } # eo if else loop - min or max

                # Compute duration of consecutive events
                duration <- full_data %>%
                    mutate(Event_ID = rleid(ECE)) %>%
                    group_by(Event_ID) %>%
                    mutate(Duration = ifelse(ECE == 1, n(), NA)) %>%
                    ungroup()
                # data.frame(duration[duration$ECE == 1,][100:110,])
                # summary(duration)

                # Summarize ECE features
                d_ECEs <- duration %>% filter(ECE == 1) %>%
                    group_by(Event_ID) %>%
                    summarise(
                        Start_Date = min(Date),
                        End_Date = max(Date),
                        Duration = n(),
                        Mean_Intensity = mean(Intensity, na.rm = TRUE),
                        Max_Intensity = max(Intensity, na.rm = TRUE),
                        Abruptness = (first(Max_Intensity) - 0) / Duration, # how fast the PEAK stress is reached
                        Heterogeneity = sd(Intensity, na.rm = TRUE) ) %>% 
                    arrange(Start_Date) %>%  # Ensure events are in chronological order
                    mutate(Recurrence_Interval = c(NA, diff(Start_Date)))
                # summary(d_ECEs); data.frame(d_ECEs)[1:200,]
                # head(data.frame(d_ECEs))
                
                ### Trouble shoot issue of extrmeley long (Duration > 500) ECE of max Ts_05 and max Ts_10 at AEW11
                # d_ECEs[d_ECEs$Duration > 500,]
                # full_data[full_data$Date %in% seq(as.Date("2015-10-01"), as.Date("2018-10-01"), by = "days"),c("Date","final_value","thresh","ECE")] 

                # Add PlotID
                d_ECEs$EP <- p

                # Clean and return
                rm(duration,full_data,dd,clim,baseline); gc()
                return(d_ECEs)

            } # eo FUN - EPs
            ,
            mc.cores = 25

        ) # eo mclapply
    
        ## Rbind
        table_ECEs <- dplyr::bind_rows(list_ECEs)
        table_ECEs <- table_ECEs %>% filter(!if_all(everything(), is.na))
        # dim(table_ECEs); head(data.frame(table_ECEs)); summary(table_ECEs)
        # data.frame(table_ECEs[7586:7686,])
        rm(list_ECEs); gc()
        
        ## Add metadata ; colnames(table_ECEs)
        table_ECEs$var <- var
        if( var == "precipitation" ) {
            table_ECEs$stat <- "total"
        } else {
            table_ECEs$stat <- stat
        } # eo if else loop - precipitation
        table_ECEs$region <- region
        table_ECEs$percentile <- perc
        table_ECEs$baseline_start <- start
        table_ECEs$baseline_end <- end

        ### 22/01/26: If system == 'grasslands' or var == 'precipitation' (for the forests), add QM strategy
        if( system == "grasslands" | var == "precipitation" ) {
            table_ECEs$qm_method <- method
        } else {
            table_ECEs$qm_method <- "No QM applied"
        } # eo if loop 

        ## And save in proper dir
        if( exists("table_ECEs") == FALSE ) {
            
            stop(
              paste("!!! ERROR: Could not find the final ECE table of ",paste(stat,var, sep = " "),
                " file for the ",system," of the ",region," based on the ",method," QM strategy and the ",perc," percentile\n", sep = "")
            )

        } else {

            setwd(extremes_dir)

            paste("Saving the final ECEs table of ",paste(stat,var, sep = " ")," file for the ",system," of the ",region,
                " based on the ",method," QM strategy and the ",perc," percentile\n", sep = "")
            
            if( system == "grasslands" ) {
                save(x = table_ECEs, file = paste("table_ECEs_",stat,"_",var,"_",system,"_",region,"_",method,"_q",perc,".Rdata", sep = "") )
            } else {
                if( var == "precipitation" ) {
                    save(x = table_ECEs, file = paste("table_ECEs_",stat,"_",var,"_",system,"_",region,"_",method,"_q",perc,".Rdata", sep = "") )
                } else {
                    save(x = table_ECEs, file = paste("table_ECEs_",stat,"_",var,"_",system,"_",region,"_q",perc,".Rdata", sep = "") )
                } # eo 2nd if else loop - var
            } # eo 1st if else loop - system

        } # eo if else loop - saving

} # eo FUN - find_extremes


### Apply find_extremes() to all possible combinations and examine outputs later (Script#7.2)
### When "forests" && != "precipitation", don't run all variations of 'qm_methods' (since no QM was directly applied)
### But make sure that "forests" && "precipitation" lead to normal settings (all thresh and QM methods)

### 09/07/26: Re-run find_extremes() for all forests and by removing those that were too young (< 10 yo) in the reference period
sys = "forests"

for(sys in c("forests","grasslands")) {

  for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10","precipitation")) {
    for(s in c("max","min")) {
      for(r in c("HND", "SCH", "SWA")) {
            m_vals <- if( sys == "forests" && v != "precipitation" ) { 
                "global"
        } else {
                c("global","monthly","mw","anomalies")
        }
        for (m in m_vals) {
          for (p in c(0.90,0.95,0.99)) {

            find_extremes(
                system = "forests",
                var = v,
                stat = s,
                region = r,
                method = m,
                perc = p,
                start = as.Date("1950-01-01"),
                end = as.Date("1980-01-01"),
                min_age = 10
            ) # eo find_extremes

          } # p
        } # m
      } # r
    } # s
  } # v

} # sys


### 09/07/26: Checking ECEs table generated
setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/ECEs_tables")
test <- get(load("table_ECEs_min_Ta_10_forests_SCH_q0.05.Rdata"))
dim(test)
head(test)
str(test)


### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
