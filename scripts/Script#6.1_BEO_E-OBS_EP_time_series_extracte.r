
### ------------------------------------------------------------------------------------------------------------

### 09/01/25 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to extract the daily temperature and precipitation data from E-OBS (see Script#6.0) at the EP's precise location
### and define daily climatologies and statistical thresholds of temperature and precipitation extremes (10th & 90th percentiles). 
### Then use the full time series of the E-OBS data to identify extremes and measure some of their basic features
### (duration, intensity)

### R script to: 
###  - Extract daily max/min temperature and daily precipitation data from E-OBS at the EP location
###  - Make gif or a video of the daily time series data (progressive time series shown through a line chart
###    or heatmaps)
###  - Compute daily climatology over the 1950-1980 period and set extreme events thresholds
###  - Identify extremes over the full time series. 
###  - Compute number/frequency/duration and intensity
###  - Make plots.

### Last update: 20/01/25 (Make master FUN for detecting Xtremes and compute their features (intensity, duration, abtruptness and heterogeneity))

### ------------------------------------------------------------------------------------------------------------

# Libraries 
library("geodata")
library("sf")
library("tidyverse")
library("data.table")
library("zoo")
library("reshape2")
library("viridis")
require("scales")
require("RColorBrewer")
library("lubridate")

### ------------------------------------------------------------------------------------------------------------

### 1°) Go to dedicated dir and read the EP metadat file (dataset #1000_10)
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/BEO-UniBern/Data/EP metadata/1000_10_Dataset")
plots <- read.csv("1000_10_data.csv", h = T)
# dim(plots); str(plots)

# Check EP IDs
# unique(plots$EP_Plot_ID)
# Keep those that are not 'na'
plots <- plots[plots$EP_Plot_ID != "na",]
# dim(plots); head(plots)
unique(plots$EP_Plot_ID) # OK, still need to add some 0 though

# FUN from Marc Beringer
BEplotZeros <- function(dat, column, plotnam = "PlotSTD"){
      dat <- as.data.frame(dat)
      funz <- function(x) ifelse((nchar(as.character(x))==4), gsub("(.)$", "0\\1", x), as.character(x)) # eo funz
      dat[,plotnam] <- sapply(dat[,column],funz)
      return(dat)
} # eo FUN
# Apply
plots <- BEplotZeros(plots, "EP_Plot_ID", plotnam = "EP")
# unique(plots$EP) # Gut


### 2°) Go to E-OBS dir and read in daily data
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/BEO-UniBern/Data/Climate data/E-OBS/") #; dir()

### 2.A) Daily max temp series
tiff_files <- dir()[grep("raster_max_temp",dir())]; tiff_files
temp_ras <- lapply(tiff_files, rast)
# Combine all rasters into one multi-layer raster
temp_ras_all <- do.call(c, temp_ras)
rm(temp_ras,tiff_files); gc()
# temp_ras_all; class(temp_ras_all)
# plot(temp_ras_all)

## Now, extract values of ma xtemp at the EP location 
## First, convert coordinates to a SpatVector object
coords <- plots[,c("EP","Longitude","Latitude")]
vect_plots <- vect(coords, geom = c("Longitude","Latitude"), crs = crs(temp_ras_all))
# class(vect_plots); plot(vect_plots)
# use extract() to extract values
max.temp.val <- terra::extract(x = temp_ras_all, y = vect_plots)
# head(max.temp.val[,1:10]) # need to drop column 1
# Cbind with coords
data <- cbind(coords[,"EP"], max.temp.val[,c(2:length(max.temp.val))])
colnames(data)[1] <- "EP"
# head(data[,1:10]) # OK
# dim(data) # OK
# summary(data)

# Melt and save as .Rdata
m.data <- melt(data, id.var = "EP")
# head(m.data) 
# Adjust colnames
colnames(m.data) <- c("EP","Date","maxT")
# str(m.data)
m.data$Date <- as.Date(m.data$Date)
# summary(m.data$Date)

# Save as .Rdata
save(x = m.data, file = "table_daily_max_temp_all_EPs_E-OBS_1950-2024_09.01.25.Rdata")
# Good. Clean and do the same for min temp and precipitation
rm(m.data,data,max.temp.val,vect_plots,coords,temp_ras_all); gc()


### ---------------------------------------------------

### 2.B) Daily min temp series
tiff_files <- dir()[grep("raster_min_temp",dir())]; tiff_files
temp_ras <- lapply(tiff_files, rast)
# Combine all rasters into one multi-layer raster
temp_ras_all <- do.call(c, temp_ras)
rm(temp_ras,tiff_files); gc()

## Now, extract values of ma xtemp at the EP location 
## First, convert coordinates to a SpatVector object
coords <- plots[,c("EP","Longitude","Latitude")]
vect_plots <- vect(coords, geom = c("Longitude","Latitude"), crs = crs(temp_ras_all))
# class(vect_plots); plot(vect_plots)
# use extract() to extract values
min.temp.val <- terra::extract(x = temp_ras_all, y = vect_plots)
# head(min.temp.val[,1:10]) # need to drop column 1
# Cbind with coords
data <- cbind(coords[,"EP"], min.temp.val[,c(2:length(min.temp.val))])
colnames(data)[1] <- "EP"
# head(data[,1:10]) # OK
# dim(data) # OK

# Melt and save as .Rdata
m.data <- melt(data, id.var = "EP")
# head(m.data) 
# Adjust colnames
colnames(m.data) <- c("EP","Date","minT")
# str(m.data)
m.data$Date <- as.Date(m.data$Date) 
# summary(m.data$Date)

# Save as .Rdata
save(x = m.data, file = "table_daily_min_temp_all_EPs_E-OBS_1950-2024_09.01.25.Rdata")
# Good. Clean and do the same for min temp and precipitation
rm(m.data,data,min.temp.val,vect_plots,coords,temp_ras_all); gc()

### ---------------------------------------------------

### 2.C) Daily precipitation
tiff_files <- dir()[grep("raster_precip_",dir())]; tiff_files
precip_ras <- lapply(tiff_files, rast)
# Combine all rasters into one multi-layer raster
precip_ras_all <- do.call(c, precip_ras)
rm(precip_ras,tiff_files); gc()

## Now, extract values of ma xtemp at the EP location 
## First, convert coordinates to a SpatVector object
coords <- plots[,c("EP","Longitude","Latitude")]
vect_plots <- vect(coords, geom = c("Longitude","Latitude"), crs = crs(precip_ras_all))
# class(vect_plots); plot(vect_plots)
# use extract() to extract values
precip.val <- terra::extract(x = precip_ras_all, y = vect_plots)
# head(precip.val[,1:40])
# Cbind with coords
data <- cbind(coords[,"EP"], precip.val[,c(2:length(precip.val))])
colnames(data)[1] <- "EP"
# head(data[,1:10]) # OK
# dim(data) # OK

# Melt and save as .Rdata
m.data <- melt(data, id.var = "EP")
# head(m.data) 
# Adjust colnames
colnames(m.data) <- c("EP","Date","precip")
# str(m.data)
m.data$Date <- as.Date(m.data$Date) 
# summary(m.data$Date); summary(m.data$precip)

# Save as .Rdata
save(x = m.data, file = "table_daily_precip_all_EPs_E-OBS_1950-2024_09.01.25.Rdata")
# Good. Clean and do the same for min temp and precipitation
rm(m.data,data,precip.val,vect_plots,coords,precip_ras_all); gc()

### Plot time series of all EPs with heatmap or geom_path
max.temp.ts <- get(load("table_daily_max_temp_all_EPs_E-OBS_1950-2024_09.01.25.Rdata"))
precip.ts <- get(load("table_daily_precip_all_EPs_E-OBS_1950-2024_09.01.25.Rdata"))
# dim(max.temp.ts); dim(precip.ts) # OK

ggplot(max.temp.ts, aes(x = Date, y = maxT)) + 
    geom_line() + labs(title = "Daily maximum temperature (°C)", x = "Date") +
    theme_minimal()

# Make an animated verison of it
library("gganimate")
# ?gganimate::animate
anim <- ggplot(max.temp.ts[max.temp.ts$EP == "AEG31",], aes(x = Date, y = maxT)) + 
    geom_line() + labs(title = "Daily maximum temperature at AEG31", x = "Date", y = "Max. temperature (°C)") +
    theme_minimal() + transition_reveal(Date)

# Save the animation in a 2 step process with animate() first
my.anim <- animate(anim, duration = 15, fps = 10, width = 1000, height = 400, renderer = gifski_renderer(), res = 100, type = "cairo")
anim_save("test.gif", animation = my.anim, path = getwd())


### ------------------------------------------------------------------------------------------------------------

### 3°) Use the TS data made above to compute daily climatologies and percentile thresholds (90th) for each EP 
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/BEO-UniBern/Data/Climate data/E-OBS/") #; dir()

# Read the daily data
max.temp.ts <- get(load("table_daily_max_temp_all_EPs_E-OBS_1950-2024_09.01.25.Rdata"))

# Add a "day of year" (day) column
max.temp.ts <- max.temp.ts %>% mutate(day = yday(Date))

# For talk at GA & other presentations
p <- ggplot(max.temp.ts[max.temp.ts$EP == "HEG45",], aes(x = Date, y = maxT)) +
  geom_line(colour = "grey50") + labs(title = "Daily maximum temperature (°C)\n(HEG45)",
    x = "Time", y = "Temperature (°C)") +
  theme_minimal()
ggsave(plot = p, filename = "time_series_daily_max.temp_single_EP_01.02.25.jpg", dpi = 300, height = 2.5, width = 6)


# Compute the daily climatology
daily_clim_max_temp <- max.temp.ts %>%
    group_by(EP,day) %>%
    summarize(
        clim = mean(maxT, na.rm = TRUE),
        .groups = "drop"
    )
# Check results
# dim(daily_clim_max_temp) 
# head(daily_clim_max_temp)
# summary(daily_clim_max_temp)
# Looks gut

# unique(daily_clim_max_temp$day) # 366? remove
daily_clim_max_temp <- daily_clim_max_temp[-which(daily_clim_max_temp$day == 366),]

# Test plot for the climatology
ggplot(daily_clim_max_temp, aes(x = day, y = clim)) +
  geom_line() + labs(title = "Daily Climatologies for Vegetation Plots",
    x = "Day of Year", y = "Mean Maximum Temperature (°C)") +
  theme_minimal()
# gut, that is for all EP

# For 1 EP only
p <- ggplot(daily_clim_max_temp[daily_clim_max_temp$EP == "HEG45",], aes(x = day, y = clim)) +
  geom_line() + labs(title = "Daily temperature thresholds for warm extremes\n(HEG45)",
    x = "Day of Year", y = "Average daily maximum temp. (°C)\n= threshold") +
  theme_minimal()

ggsave(plot = p, filename = "time_series_threhsolds_single_EP_01.02.25.jpg", dpi = 300, height = 3, width = 5)

### Now, to define the statitsical thresholds needed to classify extreme events (heatwaves in this case), one can use 
### several approaches:
# Approach	                      Pros	                                        Cons
# Fixed Day                       Percentiles	Simple to calculate.	            Sensitive to outliers and short-term fluctuations.
# Moving Window                   Percentiles	More robust to variability.	      Requires more computation.
# Longer Windows (e.g., 31 days)	More stable thresholds over time.	            May smooth out real extremes too much.

### Why use a moving window?
# Using a 15-day (or similar) moving window is a common practice in climate research. It’s recommended by the
# World Meteorological Organization (WMO) for calculating extreme temperature thresholds.
# For example: Heatwaves are often defined using a 90th percentile threshold based on a moving window around
# each calendar day. Cold spells use a similar approach with the 10th percentile threshold.
# Without a moving window: You might classify January 14 as an extreme cold day, but not January 15,
# due to slight temperature variations.
# With a moving window: The thresholds are smoother, making extreme events more consistent and realistic.
# By using a moving window, you avoid basing thresholds on small sample sizes.
# If you calculate percentiles based on data for a single day (e.g., Jan. 15), you're using fewer data points,
# which can make the threshold sensitive to outliers.

## All in all, If you calculate a percentile threshold for Jan 15 using only the data from January 15 across all years,
# you risk:
# - Missing nearby context (e.g., temperatures from January 14 or January 16).
# - Overfitting to a specific date, ignoring natural short-term fluctuations.
# A 15-day moving window smooths this variability by considering data from a window of days around each date,
# giving more stable and reliable thresholds.

### Let's use the historical (1950-1980) tempertaure data of one chosen EP to define those thresholds based on a
### 15-day rolling window from 'max.temp.ts' (and other data) then.
# str(max.temp.ts)

### Subset one random grassland plot
# unique(max.temp.ts$EP)
max.temp.ts <- max.temp.ts[max.temp.ts$EP == "HEG45",]
# dim(max.temp.ts); tail(max.temp.ts)

max.temp.ts$day <- yday(max.temp.ts$Date) # day of the year (1-366)
max.temp.ts$day2 <- lubridate::day(max.temp.ts$Date) # day of the month (1-31)
max.temp.ts$month <- lubridate::month(max.temp.ts$Date)
max.temp.ts$year <- lubridate::year(max.temp.ts$Date)
# Check
# summary(max.temp.ts)
# unique(max.temp.ts$year)

hist.max.temp.ts <- max.temp.ts[max.temp.ts$year < 1980,]
# summary(hist.max.temp.ts$Date)
# unique(hist.max.temp.ts$Date)

### Define rolling window FUN
#require("slider")
#moving_window_perc <- function(data, percentile) {
#    data %>%
#    group_by(day) %>%
#    mutate(
#      threshold = slider::slide_dbl(
#        .x = maxT,
#        .f = ~ quantile(.x, probs = percentile / 100, na.rm = TRUE),
#        .before = 7,  # 7 days before the current day
#        .after = 7    # 7 days after the current day
#      )
#    )
#} # eo FUN - moving_window_perc

## Calculate the 90th percentile threshold using a 15-day window
#p90_thres <- moving_window_perc(data = hist.max.temp.ts, percentile = 90)
# dim(p90_thres)
# head(p90_thres)
# There should be one threshold value per 'day'
#ggplot(aes(x = day, y = threshold), data = p90_thres) + geom_point() + theme_classic()
# summary(p90_thres)
# data.frame(p90_thres[3500:4000,]) # Looks good
#unique(p90_thres$threshold)


### Alternarive method with the 'zoo' package 
# If you want, apply a rolling window directly to historical data for each DOY
# Compute one 90th percentile threshold per DOY
clim <- hist.max.temp.ts %>% 
      group_by(day) %>%
      summarise(
        t90 = quantile(maxT, probs = 0.90, na.rm = TRUE)
      )
#dim(clim) # data.frame(clim)
# Quick plot to check
# ggplot(aes(x = day, y = t90), data = clim) + geom_point() + theme_void()
# Good enough? You may use those threhsolds for defining heat extreme events at EP of interets now!



### OPTIONAL ##################################################################################
                                                                                              #
# Apply a circular rolling window 15-day rolling window to smooth thresholds.                 #
# May be needed if you don't have enouhg data to define robust thresholds.                    #
# Here, we have 30 years of daily data so it should be fine.                                  #
                                                                                              #
# It needs to be CIRCULAR otherwise you loose the first and last seven days                   #
# of the year with the normal rolling 15 day window.                                          #
# To do so, repeat the first 7 and last 7 DOYs to ensure full smoothing                       # 
clim_wrapped <- bind_rows(                                                                    #
  clim %>% slice((n() - 7 + 1):n()),  # Add last 7 DOYs to the start                          #
  clim,                                                                                       #
  clim %>% slice(1:7)                # Add first 7 DOYs to the end                            #
)                                                                                             #
                                                                                              #
# Apply smoothing with 15-day rolling window                                                  #
s.clim <- clim_wrapped %>%                                                                    #
  mutate(                                                                                     #
    s_t90 = rollapply(                                                                        #
      c(t90, t90, t90),  # Wrap data around                                                   # 
      width = 15,                                                                             #
      FUN = mean,                                                                             #
      fill = NA,                                                                              #
      align = "center"                                                                        #  
    )[seq_along(t90)]  # Keep only the middle part (original data)                            #
  )                                                                                           #
                                                                                              #
# Check outputs                                                                               #
# summary(s.clim)                                                                             #  
# data.frame(s.clim)                                                                          #
# ggplot(aes(x = day, y = t90), data = s.clim) + geom_point() + theme_void()                  #
                                                                                              #
## Compare                                                                                    #
#library("ggpubr")                                                                            #
#ggarrange(                                                                                   #
#  ggplot(aes(x = day, y = t90), data = s.clim) + geom_point() + theme_void(),                #
#  ggplot(aes(x = day, y = t90), data = clim) + geom_point() + theme_void()                   #
#)                                                                                            #
###############################################################################################


## Join thresholds with the FULL historical daily temperature data
# Merge the historical 90th percentile thresholds with the full dataset
full_data <- max.temp.ts %>% left_join(clim[,c("day","t90")], by = "day")
# head(full_data)
# summary(full_data)
# unique(full_data$t90)

ggplot(aes(x = t90, y = maxT), data = full_data) +
  geom_point(colour = "grey50", alpha = .1) + 
  geom_abline(xintercept = 1, yintercept = 1,
    color = "red", linetype = "dashed") + 
  xlab("Temperature threshold (90th percentile; °C)") +
  ylab("Daily maximum temperature (°C)") +
  theme_classic()
### --> Points above the threshold would be extremes

# Apply the threshold to identify heat extremes in the full time series (1950-2024)
full_data <- full_data %>% mutate(Heat_Extreme = ifelse(maxT > t90, 1, 0))
# summary(full_data)
# dim(full_data[full_data$Heat_Extreme == 1,])
# nrow(full_data[full_data$Heat_Extreme == 1,])/nrow(full_data) # 15% of the dates are extremes! (that's a LOT) 

### Some visualization
# Plot maximum temperature and thresholds (UGLY plot)
#ggplot(full_data, aes(x = Date)) +
#  geom_line(aes(y = maxT), color = "grey50", size = 1) +
#  geom_point(aes(y = maxT, color = as.factor(Heat_Extreme)), size = 2) +
#  scale_color_manual(values = c("0" = "grey50", "1" = "red")) +
#  labs(title = "Daily Max Temperature and Heat Extremes",
#       x = "Date", y = "Temperature (°C)", color = "Heat Extreme") +
#  theme_minimal()

## To calculate the frequency of heat extremes per year:
freq <- full_data %>% group_by(year) %>% summarise(Count = sum(Heat_Extreme, na.rm = TRUE))

# Plot
p <- ggplot(aes(x = year, y = Count), data = freq) +
  geom_path(colour = "grey") + geom_point() +
  geom_smooth(method = "loess", colour = "red", se = TRUE) + 
  xlab("Year") + ylab("Number of heat extremes") +
  theme_classic() + ggtitle("Number of extreme heat events per year\n(HEG45)")

ggsave(plot = p, filename = "time_series_freq_EHE_single_EP_01.02.25.jpg", dpi = 300, height = 3, width = 6)

### ECEs features and how to compute them:
# - Intensity: The magnitude of the extreme event above the threshold 
# - Duration: The length of consecutive days classified as extremes.
# - Frequency: The number of extreme events within a given period (e.g., annually).

## Step 1: Compute heat extremes and intensity
full_data <- full_data %>% mutate(Intensity = ifelse(Heat_Extreme == 1, maxT - t90, NA))
# summary(full_data[full_data$Heat_Extreme == 1,]) # some venets are > 10°C warmer than their thresolds
# Show mean intensity of heat events per year
intens <- full_data %>% group_by(year) %>% summarise(Intensity = mean(Intensity, na.rm = TRUE))

# Plot
ggplot(aes(x = year, y = Intensity), data = intens) +
  geom_path(colour = "grey") + geom_point() +
  geom_smooth(method = "loess", colour = "red", se = TRUE) + 
  xlab("Year") + ylab("Mean intensity of heat extremes") +
  theme_classic()

## Step 2: Compute duration of consecutive events
duration <- full_data %>%
  mutate(
    Event_ID = rleid(Heat_Extreme) # Give unique IDs for consecutive events
  ) %>%
  group_by(Event_ID) %>%
  mutate(
    Duration = ifelse(Heat_Extreme == 1, n(), NA)
  ) %>%
  ungroup()
# summary(duration) # longest event was 23 days in a row

## Step 3: Summarize event-level statistics
ECEs <- duration %>%
  filter(Heat_Extreme == 1) %>%
  group_by(Event_ID) %>%
  summarise(
    Start_Date = min(Date),
    End_Date = max(Date),
    Duration = n(),
    Mean_Intensity = mean(Intensity, na.rm = TRUE),
    Max_Intensity = max(Intensity, na.rm = TRUE)
  )
# dim(ECEs) # Nb of events
# head(ECEs)
# summary(ECEs)

## Step 4: Compute annual frequency and summaries
annual_frequency <- ECEs %>%
  mutate(Year = year(Start_Date)) %>%
  group_by(Year) %>%
  summarise(
    Number_of_Events = n(),
    Total_Duration = sum(Duration),
    Mean_Duration = mean(Duration, na.rm = TRUE),
    Mean_Intensity = mean(Mean_Intensity, na.rm = TRUE)
  )
# Check
# data.frame(annual_frequency)
# summary(annual_frequency)

# Assemble TS plots in a panel
p1 <- ggplot(aes(x = Year, y = Number_of_Events), data = annual_frequency) +
  geom_vline(xintercept = 1980, linetype = "dotted", color = "grey20") + 
  geom_vline(xintercept = 2008, linetype = "dashed", color = "black") + 
  geom_path(colour = "grey") + geom_point() +
  geom_smooth(method = "loess", colour = "#d53e4f", se = TRUE) +
  xlab("Year") + ylab("Number of Events") +
  theme_classic()

p2 <- ggplot(aes(x = Year, y = Mean_Duration), data = annual_frequency) +
  geom_vline(xintercept = 1980, linetype = "dotted", color = "grey20") +
  geom_vline(xintercept = 2008, linetype = "dashed", color = "black") + 
  geom_path(colour = "grey") + geom_point() +
  geom_smooth(method = "loess", colour = "#d53e4f", se = TRUE) +  
  xlab("Year") + ylab("Mean duration of Events\n(day)") +
  theme_classic()

p3 <- ggplot(aes(x = Year, y = Mean_Intensity), data = annual_frequency) +
  geom_vline(xintercept = 1980, linetype = "dotted", color = "grey20") +
  geom_vline(xintercept = 2008, linetype = "dashed", color = "black") + 
  geom_path(colour = "grey") + geom_point() +
  geom_smooth(method = "loess", colour = "#d53e4f", se = TRUE) + 
  xlab("Year") + ylab("Mean Intensity of Events\n(°C)") +
  theme_classic()

library("ggpubr")
panel <- ggarrange(p1,p2,p3, ncol = 1, nrow = 3)
annotate_figure(panel, top = text_grob("Extreme Heat Events Dynamics at AEG06", color = "black", face = "bold", size = 10))
### Great! Apply to multiple EPs...


### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------