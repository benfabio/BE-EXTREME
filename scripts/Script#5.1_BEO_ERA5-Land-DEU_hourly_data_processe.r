### ------------------------------------------------------------------------------------------------------------

### 11/03/25 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to extract, at every BE EP location, the time series of hourly ERA5 Land data 
### that Adrian Huerta (Klimatologie Group, GIUB) downloaded and deposited on the climcal server: 
### /mnt/climstor2/vol01_ecmwf/era5-land-de/

### R script to: 
### - For every system, region, EP and variable available, extract the hourly time series from the .nc file
###   that Adrian H. downloaded
### - Check quality of data, save as separate .csv or .txt files here: /home/fbenedetti/ERA5-Land-DEU-processed/hourly/ 

### Last update: 23/04/25 (Extracting hourly 'volumetric_soil_water_layer_2' at every EP)

### ------------------------------------------------------------------------------------------------------------

# Libraries
library("dplyr")
library("reshape2")
library("stringr")
library("data.table")
library("lubridate")
library("raster")
library("terra")
library("ncdf4")
library("viridis")
library("scales") 
library("RColorBrewer")
library("parallel")

### ------------------------------------------------------------------------------------------------------------

### 1°) Get plots' metadata (spatial coordinates in WGS84) and correct IDs
setwd("/home/fbenedetti/Exploratories/EP metadata/1000_10_Dataset")
plots <- read.csv("1000_10_data.csv", h = T)
plots <- plots[plots$EP_Plot_ID != "na",]
# Need to add the 0 in the EP id though
BEplotZeros <- function(dat, column, plotnam = "PlotSTD"){
      dat <- as.data.frame(dat)
      funz <- function(x) ifelse((nchar(as.character(x))==4), gsub("(.)$", "0\\1", x), as.character(x)) # eo funz
      dat[,plotnam] <- sapply(dat[,column],funz)
      return(dat)
} # eo FUN
# Apply
plots <- BEplotZeros(plots, "EP_Plot_ID", plotnam = "EP")
ids <- unique(plots$EP)
# ids # good, vector of plot IDs, to use n the master FUN below


### 2°) Write master function to extract ERA5-Land hourly data for each plot and variable

# Vector of variables To extarct as TS
setwd("/mnt/climstor2/vol01_ecmwf/era5-land-de/")
variables <- dir() 
# Remove 'evaporation_from_vegetation_transpiration' for now 
variables <- variables[-c(2)]

# Master FUN: era5_land_extracter()
era5_land_extracter <- function(p) {

    # Useless message 
    message(paste("Extracting ERA5-Land hourly time series of ",var," for plot: ",p,"\n", sep = ""))
    
    # Based on the EP's ID, find if its system (grassland/forest) and region (HND/SCH/SWA)
    # (sanity checks)
    if( grepl("G", p) ) {
        system <- "grassland"
    } else if(grepl("W", p)) {
        system <- "forest"
    } else {
        stop( paste("Plot ",p," is neither a grassland or forest EP; Check plot ID again?", sep = "") )
    } # eo if else loop 

    if( grepl("H", p) ) {
        region <- "HND"
    } else if(grepl("A", p)) {
        region <- "SWA"
     } else if(grepl("S", p)) {
        region <- "SCH"
    } else {
        stop( paste("Plot ",p," does not belong to one of the Exploratories region; Check plot ID again?", sep = "") )
    } # eo if else loop 

    # From 'plots' data.frame, get coordinates of the EP (sanity check)
    coords <- plots[plots$EP == p,c("Longitude","Latitude")] # coords
    # If NA in coords, stop FUN
    if( anyNA(coords) ) {
        stop( paste("Plot ",p," is missing at least one spatial coordinate!", sep = "") )
    } # eo if loop 

    # Go to ERA5-Land-DEU data dir 
    setwd(paste("/mnt/climstor2/vol01_ecmwf/era5-land-de/",var,"/", sep = "")) # dir()
    # Check that dir has at least 50 years/files
    if( length(dir()) < 50 ) {
        stop( paste("Less than 50 years of hourly ERA5 Land data available! Check content of dir again...", sep = "") )
    } # eo if loop 

    # Extract Nb. of years/filenames 
    years <- as.numeric(str_replace_all(dir(),".nc",""))
    
    # In a lapply, open connexion to corresponding .nf file and extract time series
    ts_data <- lapply(years, function(y) {
                
                # Message again
                message(paste("Extracting ERA5-Land ",var," at EP ",p," for ",y, sep = ""))
                
                # Load .nc file with ncdf4 and extarct values based on 'coords'
                nc <- nc_open(paste(y,".nc", sep = ""))
                # Check available variable names
                # print(nc) ; print(nc$var)

                # Extract latitude, longitude, and time variables
                lons <- ncvar_get(nc, "longitude")
                lats <- ncvar_get(nc, "latitude")
                times <- ncvar_get(nc, "valid_time")
                # Convert to Date-Time (POSIXct)
                time_dates <- as.POSIXct(times, origin = "1970-01-01", tz = "UTC")
                # head(time_dates); summary(time_dates)

                # Define target location based on 'coords'
                target_lon <- coords$Longitude
                target_lat <- coords$Latitude

                # Find the closest grid point (indexing)
                lon_idx <- which.min(abs(lons - target_lon))
                lat_idx <- which.min(abs(lats - target_lat))

                # Get name of variable from .nc - should always be the 3rd var in the .nc
                v <- names(nc$var)[3]

                # Extract all time steps at this location
                vals <- ncvar_get(nc, v, start = c(lon_idx, lat_idx, 1), count = c(1,1,-1))
                # The count = c(1, 1, -1) argument in ncvar_get() specifies how many values to
                # extract along each dimension of the NetCDF variable
                
                # Temperatures in ERA5 Land are given in Kelvin --> Convert to °C
                if( var %in% c("2m_temperature","skin_temperature","soil_temperature_level_1",
                    "soil_temperature_level_2","soil_temperature_level_3") ) {
                    vals <- vals - 273.15
                } # eo if loop

                # NOTE: total precipitation is given in meters in ERA5-Land

                # Combine time and values into the same ddf an return
                data <- data.frame(region = region, system = system, plot = p, time = time_dates, var = var, value = vals)

                # Close the NetCDF file
                nc_close(nc)
                # Clean
                rm(vals,lon_idx,lat_idx,target_lat,target_lon,time_dates,lons,lats); gc()
                return(data)

            } # eo FUN

    ) # eo lapply - y
    # Rbind
    ts_ddf <- dplyr::bind_rows(ts_data)
    # Check 
    # dim(ts_ddf); head(ts_ddf); summary(ts_ddf)

    # Save in proper dir as .Rdata
    setwd(paste("/home/fbenedetti/ERA5-Land-DEU-processed/hourly/",system,"/",region,"/", sep = ""))
    save(x = ts_ddf, file = paste("table_hourly_ERA5-Land_data_",var,"_",p,".Rdata", sep = ""))
    rm(ts_ddf); gc()

    # Go back to initial dir
    setwd(paste("/mnt/climstor2/vol01_ecmwf/era5-land-de/", sep = ""))

} # eo FUN

### Apply master FUN in mclapply() to go faster
mc <- 30 # nb of cores to run thz FUN on - be mindful of others!

#for(var in variables) {
var <- "volumetric_soil_water_layer_2"
mclapply(X = ids, FUN = era5_land_extracter, mc.cores = mc)
#} # eo for loop - var in variables

### When finished --> Script#5.2 to derive daily statistics from the hourly data

### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------