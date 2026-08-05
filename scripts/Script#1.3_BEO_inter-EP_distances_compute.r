### ------------------------------------------------------------------------------------------------------------

### 10/06/25 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to compute inter-EPs distances (in km) based on their spatial coordinates (lat./long.).
### We will use these distances for various analyses that may include: accounting for spatial autocorrelation
### in causal linear models, identifying the top 5/10/15/30 closest grassland or forest EPs for 
### measuring microclimatic temperature offsets etc.

### We will use the 'geosphere' R package and the distm() function to compute inter-EPs Haversine distances
### (that integrates the curvature of the Earth - but that does not matter too much on that scale)

### Last update: 05/08/26 (Refreshing to sync with BE-EXTREME)

### ------------------------------------------------------------------------------------------------------------

# Basic libraries 
library("tidyverse")
library("purrr")
library("tidyr")
library("reshape2")
library("lubridate")
library("geosphere")
library("GeNetIt")

### ------------------------------------------------------------------------------------------------------------

### 1°) Load the EPs' metadata (including their spatial coordinates)
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/BEO-UniBern/Data/EP metadata/1000_10_Dataset"); dir()
d <- read.csv("1000_10_data.csv", h = TRUE, sep = ",", dec = ".")
# head(d); str(d)

## Subset EPs (EP_Plot_ID != "na")
d <- d[d$EP_Plot_ID != "na",]
# str(d); head(d); dim(d)
# unique(d$EP_Plot_ID)
# Keep those metadata. But add a '0' in the 'EP_Plot_ID' - @FUN from Marc Beringer
BEplotZeros <- function(dat, column, plotnam = "EP_Plot_ID") {
      dat <- as.data.frame(dat)
      funz <- function(x) ifelse((nchar(as.character(x))==4), gsub("(.)$", "0\\1", x), as.character(x)) # eo funz
      dat[,plotnam] <- sapply(dat[,column],funz)
      return(dat)
} # eo FUN
# Apply
dd <- BEplotZeros(d, "EP_Plot_ID", plotnam = "EP")
str(dd); head(dd); dim(dd)
# Drop 'Plot_ID'
dd <- dd %>% select(-c(Plot_ID))
# Move 'EP' to 2nd position
dd <- dd %>% relocate(EP, .after = EP_Plot_ID)



### 2°) Computing inter-EP distances

# Extract coordinates (longitude first, then latitude)
coords <- as.matrix(dd[,c("Longitude","Latitude")])
# dim(coords); summary(coords) # good

# Compute distance matrix (in meters)
dist_matrix_m <- distm(coords, fun = distHaversine)
# Convert to kilometers
dist_matrix_km <- dist_matrix_m / 1000

# Add rownames/colnames
rownames(dist_matrix_km) <- dd$EP
colnames(dist_matrix_km) <- dd$EP
# head(dist_matrix_km)

# Melt to have rowwise data
dists <- GeNetIt::dmatrix.df(dist_matrix_km)
# head(dists) ; summary(dists)
# unique(dists$from) ; unique(dists$to)
colnames(dists)[c(1,2)] <- c("EP","EP2")

# Merge with 'dd'
df <- merge(dd, dists, by = "EP")
# dim(df) ; head(df) ; summary(df)
# df[8700:8725,]


### BONUS: How to find the 10 closest EP for any EP?

## Specify the ID of your location of interest
#i <- "HEW01"
#N <- 30
## Filter for distances where your location is either 'FROM' or 'TO'
closest <- df %>%
  filter(EP == i | EP2 == i) %>%
  # Remove self-pairs if present
  filter(EP != EP2) %>%
  # Arrange by distance
  arrange(distance) %>%
  # Select the top N closest
  slice_head(n = N)
# Check
# closest


## To find the 10 closest grasslands for one specific forest EP
# top_n <- 5
# x <- df 

find_closest_grass <- function(x, top_n) {
    
    # Identify forest EPs and grassland EPs
    forest_eps <- x %>% filter(Landuse == "Forest") %>% pull(EP) %>% unique()
    grassland_eps <- x %>% filter(Landuse == "Grassland") %>% pull(EP) %>% unique()

    # Filter df to keep only rows where EP is a forest and EP2 is a grassland
    cl <- x %>%
        filter(EP %in% forest_eps, EP2 %in% grassland_eps) %>%
        group_by(EP) %>%
        arrange(distance, .by_group = TRUE) %>%
        slice_head(n = top_n) %>%
        ungroup()

    # Return
    return(cl)

} # eo FUN find_closest_grass

## Apply
closest_grass <- data.frame(find_closest_grass(x = df, top_n = 5))

# Check
# dim(closest) #; head(closest)
# closest[550:580,] # gut

### And to find the reverse: 'top_n' closest forests for each grassland EP

find_closest_forst <- function(x, top_n) {
    
    # Identify forest EPs and grassland EPs
    forest_eps <- x %>% filter(Landuse == "Forest") %>% pull(EP) %>% unique()
    grassland_eps <- x %>% filter(Landuse == "Grassland") %>% pull(EP) %>% unique()

    # Filter df to keep only rows where EP is a forest and EP2 is a grassland
    cl <- x %>%
        filter(EP %in% forest_eps, EP2 %in% grassland_eps) %>%
        group_by(EP2) %>%
        arrange(distance, .by_group = TRUE) %>%
        slice_head(n = top_n) %>%
        ungroup()

    # Return
    return(cl)

} # eo FUN find_closest_forst

## Apply
closest_forst <- data.frame(find_closest_forst(x = df, top_n = 5))
#dim(closest_forst)
#closest_forst[201:231,]

### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
