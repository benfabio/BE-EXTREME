### ------------------------------------------------------------------------------------------------------------

### 24/06/25 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to load dataset #17486 (main tree stand age) and combine it with outputs from R Script#1.3
### (inter-EP distances compute.R) to obtain the base dataset that will later enable you to examine variations
### in temperature (air, skin & soil) offsets associated to microclimatic effects driven by the tree canopy. 

### The purpose here is to create a simple dataset that will inform, for each forest EP (n = 150): 
### - their mean stand age from 2009 to 2024
### - their main tree species (Beech, Spruce, Pine) and whether it is a deciduous or a coniferous forest 
### - the 5 closest grassland EPs and their actual distance to the forest EP (in km)
### - the elevation of each forest and grassland EP (for adiabatic corrections) - dataset #31501

### Based on this, we will be able to examine and model the relationships between mean daily offsets in 
### min/max temperatures and stand age for various types of forest (see R Script#4.4) with MEMs or GAMMs.
### Perhaps look into the effect of forest land use as well (based on the ForMIX dataset #31855). 

### Last update: 24/06/25 (Saving dataset)

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

### 1°) Source R Script#1.3 and find to 5 closest grassland EP for each forest EP
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/BEO-UniBern/Data/Climate data/workspace")
source("Script#1.3_BEO_inter-EP_distances_compute.r")
# ls()
# dim(closest_grass); str(closest_grass)

### To find the 5 closest grasslands for one specific forest EP
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
# dim(closest_grass); str(closest_grass); closest_grass[1:50,] # OK


### 2°) Load dataset #17486 and get the forest EPs' stand age from 2009 to 2024 (sicne you know the mean stand age in 2012) + mean tree species
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/BEO-UniBern/Data/LUI/Stand_age_forest_EPs_2012/17486_2_Dataset")
age <- read.csv("17486_2_data.csv", h = TRUE, sep = ",", dec = ".")
# dim(age); str(age); summary(age); unique(age$EP) # looks good

# For the majority of plots stand age, that is the age of different stand layers or of admixed species was obtained from records
# from the various forest administrations. However, for unmanaged and selection forests in the Hainich-Dün (HAI) stand age,
# assessed as age of the overstorey, was estimated from diameter of the largest 30 trees per ha based on data which Mund (2004)
# sampled in nearby selection forests.

# Notes:
# Fs = European beech - Fagus sylvatica
# Qs = sessile and pedunculate oak - Quercus petrea, Q. robur
# Pa = Norway spruce - Picea abies
# Ps = Scots pine - Pinus sylvestris
# Pm = Douglas fir - Pseudotsuga menziesii
# oHS = other hardwood species (Ash, Maple, Hornbeam) - Fraxinus excelsior, Acer pseudoplatanus, Acer platanoides, Carpinus betulus, and others
# oSS = other softwood broadleaved species (Populus, Prunus, Salix, Sorbus) - Populus, Prunus, Salix, Sorbus

# MTS = main tree species (Hauptbaumart)
# ATS1 = admixed tree species #1 (Mischbaumart 1)
# ATS2 = admixed tree species #2 (Mischbaumart 2)
# age2012_MTS = age of MTS in 2012

## Relative number of MTS?
round((summary(factor(age$MTS))/150)*100, 2)
#  Fs   Fs/oHS   oHS    Pa     Ps     Qs 
# 66.67   1.33   0.67  11.33  14.67   5.33

### From 'age', get the age of each forest EP from 2009 to 2024
# Subset ddf
d <- data.frame(EP = unique(age$EP), Age_2012 = age$age2012_MTS) # d
# Years you want to calculate age for
years <- 2009:2024
# Create long format data frame using expand.grid
long_df <- expand.grid(EP = d$EP, Year = years)
# Merge with 2012 ages
long_df <- merge(long_df, d, by = "EP")
# Calculate age for each year
long_df$Age <- long_df$Age_2012 + (long_df$Year - 2012)
# Drop Age_2012 if you don't need it anymore
long_df$Age_2012 <- NULL
# Sort
long_df <- long_df[order(long_df$EP, long_df$Year), ]
# head(long_df)

## Add main tree species to 'long_df'
long_df$MTS <- NA
for(p in unique(long_df$EP)) { message(p) ; long_df[long_df$EP == p,"MTS"] <- age[age$EP == p,"MTS"] } # eo for loop - EP
# long_df$MTS; summary(factor(long_df$MTS))



### 3°) Load dataset #17706 and get the forest EPs' complete forest type description 
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/BEO-UniBern/Data/LUI/Forest_type_classif_forest_EPs_2008-2014/17706_5_Dataset")
type <- read.csv("17706_5_data.csv", h = TRUE, sep = ",", dec = ".") 
# dim(type); str(type); summary(type)
# unique(type$Forest_type_in_detail)

## Also add to 'long_df'
long_df$Forest_type <- NA
for(p in unique(long_df$EP)) { message(p) ; long_df[long_df$EP == p,"Forest_type"] <- type[type$EP == p,"Forest_type_in_detail"] } # eo for loop - EP
# long_df$Type; summary(factor(long_df$Type))

## And finally, add whether the MTS is deciduous or coniferous
# unique(long_df$MTS)
long_df$MTS_type <- NA
long_df[long_df$MTS == "Fs","MTS_type"] <- "Deciduous"
long_df[long_df$MTS == "Qs","MTS_type"] <- "Deciduous"
long_df[long_df$MTS == "Ps","MTS_type"] <- "Coniferous"
long_df[long_df$MTS == "Fs/oHS","MTS_type"] <- "Deciduous"
long_df[long_df$MTS == "Pa","MTS_type"] <- "Coniferous"
long_df[long_df$MTS == "oHS","MTS_type"] <- "Deciduous"

# Check proportions
(summary(factor(long_df$MTS_type))/2400)*100 # 75% deciduous, 25% coniferous
# summary(long_df)

## Examining distribution of mean stand age as a function of MTS and MTS_type
ggplot(data = long_df, aes(x = factor(MTS), y = Age, fill = factor(MTS))) + geom_violin(colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    scale_fill_brewer(name = "Main tree species", palette = "Paired") + 
    xlab("Main tree species") + ylab("Age (year)") + theme_bw()
### --> Large variability in age for forests made of Fs, Qs, oHS & Ps
### --> Forests made of oHS, Fs/oHS and Pa show lower variability in age. oHS forests are much older (nearly all 150 year old)
###     than Fs/oHS and Pa forests

ggplot(data = long_df, aes(x = factor(MTS_type), y = Age, fill = factor(MTS_type))) + geom_violin(colour = "black") +
    geom_boxplot(fill = "white", colour = "black", width = .1) +
    scale_fill_brewer(name = "Tree type", palette = "Paired") + 
    xlab("Main tree species") + ylab("Age (year)") + theme_bw()
### --> Deciduous forests are on average older than coniferous forests
mean(long_df[long_df$MTS_type == "Coniferous","Age"]) ; sd(long_df[long_df$MTS_type == "Coniferous","Age"])
mean(long_df[long_df$MTS_type == "Deciduous","Age"]) ; sd(long_df[long_df$MTS_type == "Deciduous","Age"])

### 4°) Add the top 5 closest grassland EPs from 'closest_grass'
df <- merge(long_df, closest_grass[,c("EP","Exploratory","EP2","distance")], by = "EP")
# Check
# df[2000:2130,]
colnames(df)[9] <- "Distance_to_EP2_km"

# Move 'Exploratory' to 1st position
df <- df %>% select(Exploratory, everything())
# head(df)


### 5°) In some instances, forests will be located higher up than grasslands. This implies thta adiabatic corrections
### should also be applied (i.e., 0.65°C/100m) when directly comparing the forests' air temperature to the grasslands'
### --> Need to also report each EP altitude! --> dataset #31501
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/BEO-UniBern/Data/EP metadata/31501_6_Dataset")
topo <- read.csv("31501_6_data.csv", h = TRUE, sep = ",", dec = ".")
# dim(topo); str(topo); summary(topo)
# Elevation does range from 15m to > 800m

# unique(topo$EP.ID) # 0 need to be added to get homogeneous EP ID -> BEplotZeros() FUN
topo <- BEplotZeros(topo, "EP.ID", plotnam = "EP")
# unique(topo$EP) # gut gut 

## Add elevation of forest EP
df$EP_elevation <- NA
for(p in unique(df$EP)) { message(p) ; df[df$EP == p,"EP_elevation"] <- topo[topo$EP == p,"Elevation"] } # eo for loop - EP
# summary(df$EP_elevation)

## Add elevation of the top 5 closest grassland EPs
df$EP2_elevation <- NA
for(p in unique(df$EP2)) { message(p) ; df[df$EP2 == p,"EP2_elevation"] <- topo[topo$EP == p,"Elevation"] } # eo for loop - EP
# summary(df$EP2_elevation)
# Check
# head(df); df[1000:1200,]
# summary(df[,c("EP_elevation","EP2_elevation")]) # median elevation is indeed higher for forests 

# Assess covariance of forest age with forest elevation? Are younger forests higher up? (pick one year)
ggplot(data = df[df$Year == 2024,], aes(x = EP_elevation, y = Age, fill = factor(Exploratory))) + 
    geom_point(pch = 21, colour = "black") + scale_fill_brewer(name = "Tree type", palette = "Paired") +
    theme_classic() + facet_wrap(.~factor(Exploratory), scales = "free")
### --> no covariande there

### Save dataet and transfer on cliclam server
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/BEO-UniBern/Data/EP metadata/")
save(x = df, file = "table_stand_age_MTS_distances_elevation_forest_EPs_24_06_25.RData")

### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------