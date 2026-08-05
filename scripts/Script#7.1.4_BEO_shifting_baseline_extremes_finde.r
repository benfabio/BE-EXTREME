### ------------------------------------------------------------------------------------------------------------

### 17/03/26 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to implement a shifting baseline approach to detect ECEs on top of the standard 
### fixed baseline approach. Fundamental to all ECEs analyses is a clearly defined climate background 
### – i.e., a temperature ‘baseline’ against which the heat event is defined.
### While a single approach to implementing a baseline may not be suitable for all extremes research applications,
### the choice of a baseline for analysing ECEs must be intentional as it affects research outcomes.
### In the case of the Biodiversity Exploratories, many researchers work on widely different processes and varities 
### of organisms. The temporal scales of these processes ad the life cycles of these organisms can span a few days
### (protists, bacteria) or seveal years (beetles, trees). They likely show heterogeneous responses to extremes and 
### likely display different adaptative responses. As a result, ECEs based on the fixed 1950-1980 basleine may not be 
### true extremes anymore to organisms and ecosystems living in 2020-2024 because they had time to adapt to the new
### climate variability. 
### Therefore, it is important that Explorers are free to switch the basleine definition based on their expertise, assumptions 
### target organisms/provcesses and study questions.

### Following Smith et al. 2025 (https://www.sciencedirect.com/science/article/pii/S0079661124002106), we will here briefly 
### summarize the shifting baseline approach and provide a R function (find_shifting_extremes()) for people in the project 
### to use. First, here are some basic definitions:
# - Baseline:	                    While reference period and baseline period are often used synonymously, we use baseline 
#                                   here to refer to the approach used for defining the reference period.

# - Fixed baseline:	                An approach based on an unchanging reference period (e.g. 1982–2011 or 1851–1900) whereby 
#                                   the climatology and threshold remain fixed over time (apart from possible seasonal variations, e.g. Hobday et al., 2016).

# - Detrended baseline:	            An approach in which temperature data are detrended prior to applying a fixed baseline. This is equivalent to
#                                   a ECE threshold that changes over time to remove the influence of slow changes in the climatological mean temperature 
#                                   (but not the variability). Note that this method has also been referred to as a shifting baseline.

# - Shifting baseline:	            An approach based on a frequently updated reference period prior to (or sometimes centred on) the analysis period 
#                                   (e.g. the climatology and threshold for a given year is based on the preceding 30-year reference period). 
#                                   Under a shifting baseline, the ECE threshold changes over time, removing the influence of slow changes in
#                                   climatological conditions (mean and variability). Detrending may also precede the use of a shifting baseline.

# - Adaptation-adjusted baseline:	An approach in which the ECE threshold changes over time in a prescribed manner, to account for the assumed adaptation 
#                                   rate of organisms (e.g. threshold increases linearly at 0.01 °C/decade).
#                                   Threshold evolution would typically be species specific (potentially informed by manipulative experiments).

### "Common approaches include a fixed baseline (Fig. 2a) where the reference period remains static for all time throughout the analysis period
### (Frölicher et al., 2018, Hobday et al., 2018, Hobday et al., 2016), a shifting baseline (Fig. 2b) where the reference period updates over time,
### keeping pace with the current analysis period (Cheung et al., 2021, Burger et al., 2022, Amaya et al., 2023), a detrended baseline (Fig. 2c)
### where long-term trends are removed from temperature data prior to calculating thresholds (Amaya et al., 2023; e.g. Jacox et al., 2020, Xu et al., 2022)
### and an adaptation-adjusted baseline (Fig. 2d) where the threshold is updated over time to reflect species adaptation potential (Logan et al., 2014, Li et al., 2023).
### Since the latter baseline approach would be tailored to understanding consequences of MHWs on a specific species or population,
### the details of the approach would be distinct in each case. As such, we focus on the other, more well-defined approaches." (Smith et al. 2025)

### Consideration of baseline type: Researchers should consider carefully what baseline type is appropriate for their application,
### including considerations of the science question at hand as well as how the information will be used.
### For example, using a:
# a. ‘Fixed baseline’ to understand changes in ecological risk under the assumption of limited adaptation or in attribution studies aimed at 
#     understanding the increased likelihood of extremes because of human-caused warming -> OUR STANDARD 
# b. ‘Shifting baseline’ to understand ecological risk in the case of rapidly adapting species, or to investigate local drivers of 
#     extremes relative to new normal conditions
# c. ‘Detrended baseline’ to separate the effect of long-term mean anthropogenic warming from variability changes in driving extremes
# d. ‘Adaptation adjusted baseline’ to account for empirical or assumed adaptation rates of specific organisms, WHERE AND WHEN SUCH DATA EXISTS
# e. ‘Periodically updated baseline’ for near term and operational assessment following adaptation of systems to contemporary conditions

### Re-use R Script#7.1 and maybe Script#7.4.1 to implement a shifting baseline approach that will cover the last 20-30 years preceding a target 
### date 'd' (to be selected by the user; couod be a sampling date for instance), instead of the fixed 1950-1980 baseline used in R Script#7.1.

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
library("data.table")

# Directory where to store ECEs tables
extremes_dir <- "/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/ECEs_tables/shifting"

### ------------------------------------------------------------------------------------------------------------

### Re-use the code from Script#6.2 to compute ECEs dynamics for every EP based on E-OBS data

### Rolling quantile helper function
roll_quantile <- function(x, p, n) {

        #' @param x = numeric vector (e.g., daily values for a given DOY across years)
        #' @param p = percentile to compute (e.g., 0.90, 0.95, 0.99)
        #' @param n = window length (number of years, e.g., 30)

        frollapply(
            x,
            n = n,
            align = "right", # ensures trailing window, DOES INCLUDE THE DOY BEYOND THE DATE OF INTEREST
            FUN = function(z) quantile(z, probs = p, type = 7, na.rm = TRUE),
            # function applied to each rolling window:
            # z = subset of x within the current window, returns the p-th quantile of that window
            # type = 7 -> standard R default quantile definition
            # na.rm = TRUE -> ignores missing values
            fill = NA # fills the first (n-1) positions with NA because there are not enough past values to form a full window
        ) # eo - frollapply

} # eo FUN - roll_quantile


## To test master FUN while you are writing it
# system = "grasslands"
# var = "precipitation"
# stat = "toal"
# region = "SWA"
# method = "mw"
# perc = 0.95
# min_age = 15
# analysis_start = base::as.Date("01-01-2004", format = "%d-%m-%Y")
# analysis_end = base::as.Date("31-12-2024", format = "%d-%m-%Y")
# window_length = 30


find_shifting_extremes <- function(system, var, stat, region, method, perc, analysis_start, analysis_end, window_length = 30, min_age = 10) { 

        #' This function takes four arguments and returns a formatted data.frame:
        #' @param system the type of ecosystem (character): "grasslands" or "forests"
        #' @param var the climate variable to process (character): "Ta_10", "Ta_200", "Ts_05", "Ts_10", "Ts_20", "Ts_50",
        #' "precipitation" or "SM_10"
        #' @param stat the daily statistic to evaluate against E-OBS and ERA5 Land (character): 'mean', 'max' or 'min'
        #' @param region Region name (character): 'SCH' or 'HND' or 'SWA'
        #' @param method Which quantile mapping approach to run (character): 'global' or 'monthly' or 'mw' or 'anomalies'
        #' @param perc Percentile to use from the refence period to define the statistical threshold (numeric): 
        #' 0 to 1 (0.90, 0.95, 0.99 for high ECEs such as heat events, or .01,.05,.1 for low events such as cold ones)
        #' @param analysis_start Date marking the start of the analysis period to detect ECE for (Date)
        #' @param analysis_end Date marking the end of the analysis period to detect ECE for (Date)
        #' @param window_length Length, in years, of the shifting baseline to derive the thresholds from (integer) - default = 30
        #' @param min_age Age below which forets data should not be taken into account because StandAge is outside of GAMM
        #' traning range (10-40 years old) (integer) - Default == 10 as youngest forest stands included in GAMM training were 12 yo
        #' @return A formatted data.frame combining the daily statistics

        ### Message
        message(
            paste("\nDetecting shifting extremes of ", stat, " ", var," for the ", system, " of the ", region,
            " using a shifting ", window_length, "-year trailing window and a ",perc, " percentile threshold\n", sep = "")
        )

        ### Load the data
        if(system == "grasslands") {
            setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/grasslands")
        } else {
            setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/combined_time_series_full/forests")
        } # eo if else loop

        if(var == "precipitation") { 
            stat <- "total"
        } # eo if loop - precip.

        if(system == "grasslands") {
            d <- get(load(paste("table_combined_obs+corr_",method,"_",stat,"_",var,"_",system,"_",region,".Rdata", sep = "")))
        } else {
            if(var == "precipitation") {
                d <- get(load(paste("table_combined_obs+corr_",method,"_",stat,"_",var,"_",system,"_",region,".Rdata", sep = "")))
            } else {
                d <- readRDS(paste("table_combined_obs+corr_",stat,"_",var,"_",system,"_",region,".rds", sep = ""))
            }
        } # eo if else loop - system

        ### Adjust percentiles when needed
        if(stat == "max" | var == "precipitation") {
            if(perc == 0.1) perc <- 0.9
            if(perc == 0.05) perc <- 0.95
            if(perc == 0.01) perc <- 0.99
        } else if(stat == "min") {
            if(perc == 0.9) perc <- 0.1
            if(perc == 0.95) perc <- 0.05
            if(perc == 0.99) perc <- 0.01
        } # eo if else loop

        ### Filtering for forest EPs
        if( system == "forests" & var != "precipitation" ) {
            d <- d[-which(d$Stand_age < min_age),]
            ep2keep <- unique(d$EP)
        } else {
            ep2keep <- unique(d$EP)
        } # eo if else loop

        ### Paralell looping over the EPs
        # ep <- ep2keep[43] ; ep # to test mclapply below
        list_ECEs <- mclapply(X = ep2keep, FUN = function(ep) {
            
                message(paste("Computing ECEs with a shifting baseline for ", ep))
                dd <- d[d$EP == ep,]

                ## Convert to data.table
                setDT(dd) # class(dd); str(dd)

                 ## Add time variables
                dd[, DOY := yday(Date)]
                dd[, year := year(Date)]

                ## Order data - Ensures rolling window is applied correctly within each DOY
                setorder(dd, DOY, year) # head(dd)

                ## Compute rolling thresholds using the roll_quantile() helper FUN from above
                ## Computes a time-varying threshold for each single DOY, based on the previous 30 years
                # roll_quantile(x = dd$final_value, p = perc, n = window_length) # test
                dd[, thresh := roll_quantile(x = final_value, p = perc, n = window_length), by = DOY]
                # head(dd); summary(dd); dim(dd)
                # dd[19354:19394,]
                ### There will be NA in the 'thresh' vector when:
                # - not enough past values for the rolling window
                # - missing values inside the window (secondary effect)
                # - incomplete TS (likely not the case here)
                # - irregular TS or ordering issues (likely not the case here)

                ## Remove years without full window
                dd <- dd[!is.na(thresh)]

                ## Restrict to analysis period
                dd <- dd[Date >= analysis_start & Date <= analysis_end]
                # dd[6235:6285,]

                ### Detect ECEs
                if( stat == "max" | stat == "total" | perc %in% c(0.9,0.95,0.99) ) {
                    dd[, ECE := ifelse(final_value > thresh, 1, 0)]
                    dd[, Intensity := ifelse(ECE == 1, final_value - thresh, NA)]
                    # as.data.frame(dd[dd$ECE == 1,][400:430,])
                } else {
                    dd[, ECE := ifelse(final_value < thresh, 1, 0)]
                    dd[, Intensity := ifelse(ECE == 1, thresh - final_value, NA)]
                } # eo if else loop

                ### IMPORTANT: reorder by Date BEFORE computing events
                setorder(dd, Date)

                ## Compute ECE duration (fix: enforce date continuity)
                dd[, Event_ID := rleid(ECE, as.numeric(Date - shift(Date)))]
                dd[, Duration := ifelse(ECE == 1, .N, NA), by = Event_ID]

                ## Summarize ECE data
                d_ECEs <- dd[ECE == 1, .(
                    Start_Date = min(Date),
                    End_Date = max(Date),
                    Duration = .N,
                    Mean_Intensity = mean(Intensity, na.rm = TRUE),
                    Max_Intensity = max(Intensity, na.rm = TRUE),
                    Abruptness = max(Intensity, na.rm = TRUE) / .N,
                    Heterogeneity = sd(Intensity, na.rm = TRUE)
                ), by = Event_ID]

                ## Sort ECEs chronologically
                setorder(d_ECEs, Start_Date)
                # head(d_ECEs); summary(d_ECEs)

                # Compute time between consecutive ECEs: the number of days between the start dates of two consecutive extreme events
                # = tells us how frequently extreme events occur over time
                # Pad the first value (since the first event has no previous one)
                d_ECEs[, Recurrence_Interval_days := c(NA, as.numeric(diff(Start_Date)))]
                
                # Add the EP ID to each ECE
                d_ECEs[, EP := ep]
                setcolorder(d_ECEs, c("EP", setdiff(names(d_ECEs), "EP")))

                ## Return
                rm(dd); gc()
                return(as.data.frame(d_ECEs))

            }, mc.cores = 25

        ) # eo mclapply - list_ECEs

        ### Combine results
        table_ECEs <- dplyr::bind_rows(list_ECEs)

        # Remove empty rows uf they exist
        table_ECEs <- table_ECEs %>% filter(!if_all(everything(), is.na))
        
        # Add metadata
        table_ECEs$var <- var
        table_ECEs$stat <- stat
        table_ECEs$region <- region
        table_ECEs$percentile <- perc
        table_ECEs$window_length <- window_length
        table_ECEs$analysis_start <- analysis_start
        table_ECEs$analysis_end <- analysis_end

        if( system == "grasslands" | var == "precipitation" ) {
            table_ECEs$qm_method <- method
        } else {
            table_ECEs$qm_method <- "No QM applied"
        } # eo if else loop

        ## Save in proper dir
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
                save(x = table_ECEs, file = paste("table_ECEs_shifting_",window_length,"yr_",stat,"_",var,"_",system,"_",region,"_",method,"_q",perc,"_from:",analysis_start,"_to:",analysis_end,".Rdata", sep = "") )
            } else {
                if( var == "precipitation" ) {
                    save(x = table_ECEs, file = paste("table_ECEs_shifting_",window_length,"yr_",stat,"_",var,"_",system,"_",region,"_",method,"_q",perc,"_from:",analysis_start,"_to:",analysis_end,".Rdata", sep = "") )
                } else {
                    save(x = table_ECEs, file = paste("table_ECEs_shifting_",window_length,"yr_",stat,"_",var,"_",system,"_",region,"_q",perc,"_from:",analysis_start,"_to:",analysis_end,".Rdata", sep = "") )
                } # eo 2nd if else loop - var
            } # eo 1st if else loop - system

        } # eo if else loop - saving

} # eo FUN - find_shifting_extremes


### 18/03/26: Test find_shifting_extremes() before re-running it for all combinations possible 
find_shifting_extremes(
    system = "grasslands",
    var = "Ta_10",
    stat = "max",
    region = "SCH",
    method = "mw",
    perc = 0.95,
    analysis_start = base::as.Date("01-01-2004", format = "%d-%m-%Y"),
    analysis_end = base::as.Date("31-12-2024", format = "%d-%m-%Y"), 
    window_length = 30,
    min_age = 15
)

find_shifting_extremes(
    system = "grasslands",
    var = "Ts_20",
    stat = "min",
    region = "HND",
    method = "mw",
    perc = 0.99,
    analysis_start = base::as.Date("01-01-2004", format = "%d-%m-%Y"),
    analysis_end = base::as.Date("31-12-2024", format = "%d-%m-%Y"), 
    window_length = 30,
    min_age = 15
)

find_shifting_extremes(
    system = "forests",
    var = "precipitation",
    stat = "total",
    region = "HND",
    method = "global",
    perc = 0.99,
    analysis_start = base::as.Date("01-01-2004", format = "%d-%m-%Y"),
    analysis_end = base::as.Date("31-12-2024", format = "%d-%m-%Y"), 
    window_length = 20, # check if it works too 
    min_age = 15
)

### Check outputs
# setwd(extremes_dir) # dir()
# t <- get(load("table_ECEs_shifting_20yr_total_precipitation_forests_HND_global_q0.99_from:2004-01-01_to:2024-12-31.Rdata"))
# class(t)
# str(t)
# head(t)
# summary(t)
### -> All test cases seem to be running smoothly for now...


### 18-19/03/26: Apply find_shifting_extremes() to all cases of interest in a multi-level for loop 
### 23/03/26: Same but for 1980-2024 period

### 09/07/26: Re-running find_extremes() for all forests and by removing those that were too young (< 10 yo) in the reference period
sys = "forests"

for(sys in c("grasslands","forests")) {

  for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10","precipitation")) {
    for(s in c("max","min")) {
      for(r in c("HND", "SCH", "SWA")) {
            m_vals <- if( sys == "forests" && v != "precipitation" ) { 
                "global"
        } else {
                c("global","monthly","mw","anomalies")
        }
        for(m in m_vals) {
          for(p in c(0.90,0.95,0.99)) {

            find_shifting_extremes(
                system = sys,
                var = v,
                stat = s,
                region = r,
                method = m,
                perc = p,
                analysis_start = base::as.Date("01-01-2004", format = "%d-%m-%Y"),
                analysis_end = base::as.Date("31-12-2024", format = "%d-%m-%Y"), 
                window_length = 30,
                min_age = 10
            ) # ep find_shifting_extremes

          } # p
        } # m
      } # r
    } # s
  } # v

} # sys


### 09/07/26: Checking ECEs table generated
setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/ECEs_tables/shifting/1980_2024")
test <- get(load("table_ECEs_shifting_30yr_max_Ta_200_forests_HND_q0.95_from:1980-01-01_to:2024-12-31.Rdata"))
dim(test)
head(test)
str(test)


### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------