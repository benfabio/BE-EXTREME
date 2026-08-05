### ------------------------------------------------------------------------------------------------------------

### 07/04/26 - ©Fabio Benedetti (Plant Ecology group, IPS, Uni Bern)

### R script to generate the figures (main and suppl.) for the ESSD Data Paper in prep.: 
### "BE-EXTREME: Climate and extreme event reconstructions (1950-2024) for the grasslands and forests of the Biodiversity Exploratories"

### List of Figures to make: 
### - Figure 1: Flow diagram - Made by hand on Keynote - Done
### - Figure 2: Maps of the Biodiversity Exploratories - Done
### - Figure 3: Comparison showing why E-OBS is better than ERA5-Land for air temperature and precipitation - Done
### - Figure 4: Distribution of the evaluation metrics for the 13 target climatic variables in the grasslands - Done
### - Figure 5: Distribution of the evaluation metrics for the 13 target climatic variables in the forests - Done
### - Figure 6: Illustration of the comparison between the fixed-baseline vs. the shifting baseline period - Done
### - Figure 7: Distribution plots of long-term climatic trends (Sen's slopes) for both grasslands and forests - Done
### - Figure 8: Long-term dynamics of extreme climate events for both grasslands and forests - Done
### - Figure 9: Heatmaps illustrating the reconstructed monthly SPEI-6 values for all three regions across all A) grasslands and B) forests - Done

### Potential Suppl. Figures:
### - Fig. S1: Panel of chosen figures illustrating daily offsets distributions per region, tree type and seasons - DONE
### - Fig. S2: Panel GAMM smooth plots - DONE
### - Fig. S3: Long-term dynamics of ECEs for both grasslands and forests - shifting baseline approach - DONE
### - Fig. S4: Heatmaps illustrating the reconstructed monthly SPEI-3 and SPEI-12 values - DONE
### - Fig. S5: Scatterplot evaluating our monthly SPEI reconstructions by comparing them to the SPEI-3/6/12 values of the global SPEI database

### Color codes to use:
# - Yellow grassland: #f9c346
# - Green forest: #3c8e7a

### Last update: 23/07/26 (Fig.S5: Evaluating our monthly SPEI reconstructions by comparing them to the SPEI-3/6/12 values of the global SPEI database)

### ------------------------------------------------------------------------------------------------------------

# Libraries
# install.packages("terra")
# library("tidyverse")
library("dplyr")
library("tidyr")
library("purrr")
library("lubridate")
library("reshape2")
library("pals")
library("maps")
library("ggplot2")
library("ggspatial")
library("ggrepel")
library("scales")
library("ggpubr")
library("RColorBrewer")
library("viridis")

### ------------------------------------------------------------------------------------------------------------

### Figure 2: Maps of the Biodiversity Exploratories

### Basic map of Germany
DEU <- map_data("world", region = "Germany")
map_de <- ggplot(DEU, aes(x = long, y = lat, group = group)) +
    geom_polygon(fill = "grey80", colour = "black") +
    coord_quickmap() + theme_void() + ggtitle("Germany") +
    theme(plot.title = element_text(hjust = 0.5))
# map_de

### Load the shapefile of each region 
# Hainich-Dün
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/BEO-UniBern/Data/Shapefiles exploratories/HND")
hnd_shp <- st_read("ExploHai_Border.shp")
# plot(st_geometry(hnd_shp))
map_HND <- ggplot(data = hnd_shp) + geom_sf(fill = "grey75", color = "black", linewidth = 0.5) +
    theme_void() + labs(title = "Hainich-Dün")

# Schorfheide-Chorin
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/BEO-UniBern/Data/Shapefiles exploratories/SCH")
sch_shp <- st_read("ExploSCH_Border.shp")
# plot(st_geometry(sch_shp))
map_SCH <- ggplot(data = sch_shp) + geom_sf(fill = "grey75", color = "black", linewidth = 0.5) +
    theme_void() + labs(title = "Schorfheide-Chorin")

# Schwäbische-Alb
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/BEO-UniBern/Data/Shapefiles exploratories/SWA")
swa_shp <- st_read("ExploSWA_Border.shp")
# plot(st_geometry(swa_shp))
map_SWA <- ggplot(data = swa_shp) + geom_sf(fill = "grey75", color = "black", linewidth = 0.5) +
    theme_void() + labs(title = "Schwäbische-Alb")

# Arrange in panel
ggarrange(map_SCH,map_HND,map_SWA, ncol = 1, nrow = 3, align = "hv")

### Need to adjust the CRS of the shapefiles so they are WGS84
sch_shp <- st_transform(sch_shp, crs = "+proj=longlat +datum=WGS84 +no_defs")
swa_shp <- st_transform(swa_shp, crs = "+proj=longlat +datum=WGS84 +no_defs")
hnd_shp <- st_transform(hnd_shp, crs = "+proj=longlat +datum=WGS84 +no_defs")

map_DEU_explos <- ggplot() + geom_polygon(data = DEU, aes(x = long, y = lat, group = group), fill = "grey90", colour = "black") +
    geom_sf(data = swa_shp, fill = "#a6d96a", alpha = .5, colour = "#1a9850") + 
    geom_sf(data = hnd_shp, fill = "#a6d96a", alpha = .5, colour = "#1a9850") + 
    geom_sf(data = sch_shp, fill = "#a6d96a", alpha = .5, colour = "#1a9850") + 
    theme_void() + ggtitle("Germany") + theme(plot.title = element_text(hjust = 0.5))

### Now, get the plots' coordinates
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/BEO-UniBern/Maps for project report Spring 2024/1000_9_Dataset")
plots <- read.csv("1000_9_data.csv", h = T, sep = ",", dec = ".")
# str(plots)

### Maps with EPs location
map_EPs_swa <- ggplot(data = swa_shp) + geom_sf(fill = "grey90", color = "black", linewidth = 0.5) +
    geom_point(data = plots[plots$Exploratory == "ALB" & plots$EP_Plot_ID != 'na' & plots$VIP == "yes",],
        aes(x = jitter(Longitude, amount = .003), y = jitter(Latitude, amount = .003), fill = factor(Landuse)), colour = "black", pch = 21) + 
    geom_point(data = plots[plots$Exploratory == "ALB" & plots$EP_Plot_ID != 'na' & plots$VIP == "no",],
        aes(x = jitter(Longitude, amount = .003), y = jitter(Latitude, amount = .003), fill = factor(Landuse)), colour = "black", pch = 21) +
    scale_fill_manual(name = "", values = c("#3c8e7a","#f9c346")) + theme_void() + 
    ggtitle("Schwäbische-Alb") + theme(plot.title = element_text(hjust = 1)) + annotation_scale(location = "br")

map_EPs_sch <- ggplot(data = sch_shp) + geom_sf(fill = "grey90", color = "black", linewidth = 0.5) +
    geom_point(data = plots[plots$Exploratory == "SCH" & plots$EP_Plot_ID != 'na' & plots$VIP == "yes",],
        aes(x = jitter(Longitude, amount = .003), y = jitter(Latitude, amount = .003), fill = factor(Landuse)), colour = "black", pch = 21) + 
    geom_point(data = plots[plots$Exploratory == "SCH" & plots$EP_Plot_ID != 'na' & plots$VIP == "no",],
        aes(x = jitter(Longitude, amount = .003), y = jitter(Latitude, amount = .003), fill = factor(Landuse)), colour = "black", pch = 21) +
    scale_fill_manual(name = "", values = c("#3c8e7a","#f9c346")) + theme_void() + 
    ggtitle("Schorfheide-Chorin") + theme(plot.title = element_text(hjust = 0.5)) + annotation_scale(location = "bl")

map_EPs_hnd <- ggplot(data = hnd_shp) + geom_sf(fill = "grey90", color = "black", linewidth = 0.5) +
    geom_point(data = plots[plots$Exploratory == "HAI" & plots$EP_Plot_ID != 'na' & plots$VIP == "yes",],
        aes(x = jitter(Longitude, amount = .003), y = jitter(Latitude, amount = .003), fill = factor(Landuse)), colour = "black", pch = 21) + 
    geom_point(data = plots[plots$Exploratory == "HAI" & plots$EP_Plot_ID != 'na' & plots$VIP == "no",],
        aes(x = jitter(Longitude, amount = .003), y = jitter(Latitude, amount = .003), fill = factor(Landuse)), colour = "black", pch = 21) +
    scale_fill_manual(name = "", values = c("#3c8e7a","#f9c346")) + theme_void() + 
    ggtitle("Hainich-Dün") + theme(plot.title = element_text(hjust = 0.5)) + annotation_scale(location = "br")

### Save maps in a panel
setwd("/Users/fabiobenedetti/Desktop/")
# Regions maps alone
# panel1 <- ggarrange(plotlist = list(map_EPs_sch,map_EPs_hnd,map_EPs_swa), ncol = 1, nrow = 3, align = "hv", common.legend = TRUE)
# ggsave(plot = panel1, filename = "panel_maps_explos_all_EPs_22.05.25.png", dpi = 300, height = 10, width = 5)
# With map of Germany
panel <- ggarrange(plotlist = list(map_DEU_explos,map_EPs_sch,map_EPs_swa,map_EPs_hnd), ncol = 2, nrow = 2, align = "hv", common.legend = TRUE)
ggsave(plot = panel, filename = "Figure_2.png", dpi = 300, height = 7.5, width = 7.5)
ggsave(plot = panel, filename = "Figure_2.pdf", dpi = 300, height = 7.5, width = 7.5)


### ------------------------------------------------------------------------------------------------------------

### 09/04/26

### Figure 3: Plotting the ERA5-Land and E-OBS biases in daily values (vs. local observations) for every parameter and 
### system (forests and grasslands)

### WARNING: E-OBS only provides air temperature (Ta_200) and precipitation 
### -> ERA5-Land vs. E-OBS biases only feasible for these 2 variables 
### All other variables -> ERA5-Land only
### -> Need to make 2 panels 

# Need to go to climcal server for this 
setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/biases") # ; dir()

### Fig. 3A) Precip. & Ta_200 (ERA5-Land & E-OBS vs. obs)

### Simple helper FUN to read all grassland biases files of one variable
# var = "Ta_200"
# stat = "max"

read_daily_biases <- function(system,var,stat) {
      # Identify the files 
      # Re-write 'stat' argument for precip. data
      if(var == "precipitation") {
          stat <- "total"
      } # eo if loop
      setwd(paste("/home/fbenedetti/ERA5-Land-DEU-processed/daily/biases/", sep = ""))
      files2keep <- dir()[grepl(paste(c(stat,var,system), collapse = "_"), dir())]
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

### Load target variables for both grasslands and forests 
# Load each tabkes and then rbind together
d1 <- na.omit( read_daily_biases(system = "grassland", var = "Ta_200", stat = "max") )
d2 <- na.omit( read_daily_biases(system = "grassland", var = "Ta_200", stat = "min") )
d3 <- na.omit( read_daily_biases(system = "forest", var = "Ta_200", stat = "max") )
d4 <- na.omit( read_daily_biases(system = "forest", var = "Ta_200", stat = "min") )
d5 <- na.omit( read_daily_biases(system = "grassland", var = "precipitation", stat = "total") )
d6 <- na.omit( read_daily_biases(system = "forest", var = "precipitation", stat = "total") )
# Add stat & var for panelling
d1$stat <- "Max" ; d1$var <- "Ta_200" ; d1$system <- "Grassland"
d2$stat <- "Min" ; d2$var <- "Ta_200" ; d2$system <- "Grassland"
d3$stat <- "Max" ; d3$var <- "Ta_200" ; d3$system <- "Forest"
d4$stat <- "Min" ; d4$var <- "Ta_200" ; d4$system <- "Forest"
d5$stat <- "Total" ; d5$var <- "Precipitation" ; d5$system <- "Grassland"
d6$stat <- "Total" ; d6$var <- "Precipitation" ; d6$system <- "Forest"

# Rbind 
biases <- bind_rows(list(d1,d2,d3,d4,d5,d6))
# dim(biases); str(biases)
# Need to melt() to have data source (E-OBS vs. ERA5-Land) in a column
m.biases <- melt(biases, id.vars = c("EP","system","region","Date","month","year","obs","E_OBS","ERA5_Land","stat","var"))
# head(m.biases); str(m.biases)
colnames(m.biases)[c(12,13)] <- c("source","bias")
m.biases$source <- as.character(m.biases$source)
m.biases[m.biases$source == "biases_EOBS","source"] <- "E-OBS"
m.biases[m.biases$source == "biases_ERA5Land","source"] <- "ERA5-Land"
# Re-order so Grasslands are before Forests 
m.biases$system <- factor(m.biases$system, levels = c("Grassland","Forest"))

# Need to make a 2x2 panel of monthly boxplots illustrating the distribution of the biases per month
# Maybe separate it per system first and then assemble with ggarnage
# unique(na.omit(m.biases[m.biases$var == "Ta_200","system"]))

p1 <- ggplot(data = na.omit(m.biases[m.biases$var == "Ta_200",]),
        aes(x = factor(month), y = bias, fill = factor(source))) + geom_boxplot(colour = "black") + 
    scale_fill_manual(name = "Source", values = c("#f9c346","#b94860")) + scale_y_continuous(limits = c(-7.5,7.5)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + xlab("") + ylab("Air temperature bias [°C]") + 
    theme_bw() + facet_wrap(factor(stat)~factor(system))
# Save
ggsave(plot = p1, filename = "Figure_3A.jpg", dpi = 300, width = 10, height = 5.5)

### And now with toal precipitation 
p2 <- ggplot(data = na.omit(m.biases[m.biases$var == "Precipitation",]),
        aes(x = factor(month), y = bias, fill = factor(source))) + geom_boxplot(colour = "black") + 
    scale_fill_manual(name = "Source", values = c("#f9c346","#b94860")) + scale_y_continuous(limits = c(-25,25)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + xlab("Month") + ylab("Total precipitation bias [mm]") + 
    theme_bw() + facet_wrap(.~factor(system))
ggsave(plot = p2, filename = "Figure_3B.jpg", dpi = 300, width = 10, height = 3)

### Save in panel 
panel1 <- ggarrange(p1,p2, align = "hv", nrow = 2, common.legend = TRUE, heights = c(1.8,1))
ggsave(plot = panel1, filename = "Figure_3A+B.jpg", dpi = 300, width = 10, height = 7.5)

### Also make overal boxplots: E-OBS vs. ERA5-Land
# Max Ta_200
b1 <- ggplot(data = na.omit(m.biases[m.biases$var == "Ta_200" & m.biases$stat == "Max",]),
        aes(x = factor(source), y = bias, fill = factor(source))) + geom_boxplot(colour = "black") + 
    scale_fill_manual(name = "Source", values = c("#f9c346","#b94860")) + scale_y_continuous(limits = c(-7.5,7.5)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + ylab("Max air temperature\nbias [°C]") + xlab("") + 
    theme_bw() + guides(fill = FALSE)  

# Min Ta_200
b2 <- ggplot(data = na.omit(m.biases[m.biases$var == "Ta_200" & m.biases$stat == "Min",]),
        aes(x = factor(source), y = bias, fill = factor(source))) + geom_boxplot(colour = "black") + 
    scale_fill_manual(name = "Source", values = c("#f9c346","#b94860")) + scale_y_continuous(limits = c(-7.5,7.5)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + ylab("Min air temperature\nbias [°C]") + xlab("") +
    theme_bw() + guides(fill = FALSE)  

# Total precip.
b3 <- ggplot(data = na.omit(m.biases[m.biases$var == "Precipitation",]),
        aes(x = factor(source), y = bias, fill = factor(source))) + geom_boxplot(colour = "black") + 
    scale_fill_manual(name = "Source", values = c("#f9c346","#b94860")) + scale_y_continuous(limits = c(-25,25)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + ylab("Total precipitation\nbias [mm]") + xlab("") + 
    theme_bw() + guides(fill = FALSE)  

# Save in panel too
panel2 <- ggarrange(b1,b2,b3, align = "hv", nrow = 3, ncol = 1, common.legend = TRUE)
ggsave(plot = panel2, filename = "Figure_3C.jpg", dpi = 300, width = 2.5, height = 7.5)

### Create panel with panel1 + panel2
panel3 <- ggarrange(panel2,panel1, align = "hv", nrow = 1, ncol = 2, common.legend = TRUE, widths = c(1,5))
ggsave(plot = panel3, filename = "Figure_3D.jpg", dpi = 300, width = 12, height = 7.5)


### 13/04/26: Other option: Modelled vs. Observed daily values scatterplot?
### -> Need to melt again to have modelled values as one column instead of 2
test1 <- ggplot(data = na.omit(m.biases[m.biases$var == "Ta_200",]), aes(x = obs, y = E_OBS)) + 
    geom_point(alpha = .5, colour = "grey50") + geom_abline (slope = 1, linetype = "dashed", color = "#b94860") +
    xlab("Observed air temperature") + ylab("Predicted air temperature\n(E-OBS)") + 
    theme_bw() + facet_wrap(factor(system)~factor(stat))

test2 <- ggplot(data = na.omit(m.biases[m.biases$var == "Ta_200",]), aes(x = obs, y = ERA5_Land)) + 
    geom_point(alpha = .5, colour = "grey50") + geom_abline (slope = 1, linetype = "dashed", color = "#b94860") +
    xlab("Observed air temperature") + ylab("Predicted air temperature\n(ERA5-Land)") + 
    theme_bw() + facet_wrap(factor(system)~factor(stat))

# panel_test <- ggarrange(test1,test2, align = "hv", nrow = 2, ncol = 1, common.legend = TRUE)
# ggsave(plot = panel_test, filename = "panel_test.jpg", dpi = 300, width = 6, height = 8)

### Add total daily precip.
test3 <- ggplot(data = na.omit(m.biases[m.biases$var == "Precipitation",]), aes(x = obs, y = E_OBS)) + 
    geom_point(alpha = .5, colour = "grey50") + geom_abline (slope = 1, linetype = "dashed", color = "#b94860") +
    xlab("Observed total precipitation") + ylab("Predicted total precipitation\n(E-OBS)") + 
    theme_bw() + facet_wrap(.~factor(system))

test4 <- ggplot(data = na.omit(m.biases[m.biases$var == "Precipitation",]), aes(x = obs, y = ERA5_Land)) + 
    geom_point(alpha = .5, colour = "grey50") + geom_abline (slope = 1, linetype = "dashed", color = "#b94860") +
    xlab("Observed total precipitation") + ylab("Predicted total precipitation\n(ERA5-Land)") + 
    theme_bw() + facet_wrap(.~factor(system))

panel_test <- ggarrange(test1,test2,test3,test4, align = "hv", nrow = 2, ncol = 2, common.legend = TRUE, heights = c(2,1))
ggsave(plot = panel_test, filename = "Figure_3_scatterplots.jpg", dpi = 300, width = 8.5, height = 7.7)

### For the other variables -> just compute evaluation metrics based on R Script#6.3.2 + Script#6.3.4

### ------------------------------------------------------------------------------------------------------------

### 16/04/26

### Figure 4: Distribution of the five metrics used to evaluate the downscaled and bias-corrected values obtained
### by applying four quantile mapping (QM) strategies to the daily data from E-OBS and ERA5-Land for all 150 grassland EPs

### 16/04/26: Summarizing the mean ± sd values of all 5 evaluation metrics for Fig. 4 and Table SX of the ESSD data paper 

setwd("/home/fbenedetti/ERA5-Land-DEU-processed/daily/quantile_mapping_outputs/evaluation_metrics")

files <- dir()[!grepl("_anoms_",dir())]
files <- files[!grepl("_forests_",files)]
# files # gut
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
# Rbind
table <- bind_rows(res); rm(res); gc()
# dim(table) ; head(table) ; summary(table)


summary_table_raw_values <- data.frame(
    table %>%
    group_by(variable,stat,qm_method) %>%
    summarize(
        MAE = round(mean(mae, na.rm = TRUE),3), MAE_sd = round(sd(mae, na.rm = TRUE),3),
        MBE = round(mean(mbe, na.rm = TRUE),3), MBE_sd = round(sd(mbe, na.rm = TRUE),3),
        Rho = round(mean(corr, na.rm = TRUE),3), Rho_sd = round(sd(corr, na.rm = TRUE),3),
        RMSE = round(mean(rmse, na.rm = TRUE),3), RMSE_sd = round(sd(rmse, na.rm = TRUE),3),
        RIA = round(mean(d1r, na.rm = TRUE),3), RIA_sd = round(sd(d1r, na.rm = TRUE),3),
    )
)
# Check
# head(summary_table_raw_values); dim(summary_table_raw_values)

### 20/04/26: Check scores of individual QM strategies against other
summary(summary_table_raw_values[summary_table_raw_values$qm_method == "global" & summary_table_raw_values$var == "Ta_10", c("MAE","RMSE","RIA")])
summary(summary_table_raw_values[summary_table_raw_values$qm_method == "monthly" & summary_table_raw_values$var == "Ta_10", c("MAE","RMSE","RIA")])
summary(summary_table_raw_values[summary_table_raw_values$qm_method == "anomalies" & summary_table_raw_values$var == "Ta_10", c("MAE","RMSE","RIA")])
summary(summary_table_raw_values[summary_table_raw_values$qm_method == "mw" & summary_table_raw_values$var == "Ta_10", c("MAE","RMSE","RIA")])
# And between 'anomalies' strategy and 'mw' strategy? Who's best based on absolute values? 
summary(summary_table_raw_values[summary_table_raw_values$qm_method == "anomalies", c("MAE","MBE","RMSE","RIA","Rho")])
summary(summary_table_raw_values[summary_table_raw_values$qm_method == "mw", c("MAE","MBE","RMSE","RIA","Rho")])
### Hard to tell whom between 'anomalies' and 'mw' is really best here.
### -> Look for differences in the  anomalies-based metrics


# Merge mean and sd values in one column 
summary_table_raw_values$MAE <- paste(summary_table_raw_values$MAE," ± ", summary_table_raw_values$MAE_sd, sep = "")
summary_table_raw_values$MBE <- paste(summary_table_raw_values$MBE," ± ", summary_table_raw_values$MBE_sd, sep = "")
summary_table_raw_values$Rho <- paste(summary_table_raw_values$Rho," ± ", summary_table_raw_values$Rho_sd, sep = "")
summary_table_raw_values$RMSE <- paste(summary_table_raw_values$RMSE," ± ", summary_table_raw_values$RMSE_sd, sep = "")
summary_table_raw_values$RIA <- paste(summary_table_raw_values$RIA," ± ", summary_table_raw_values$RIA_sd, sep = "")

# Drop '_sd' columns: 5, 7, 9, 11, 13
drop.cols <- colnames(summary_table_raw_values)[c(5,7,9,11,13)]
summary_table_raw_values <- summary_table_raw_values %>% select(-drop.cols)

# Adjust colnames 
colnames(summary_table_raw_values)[c(1,2,3)] <- c("Climatic_parameter","Daily_statistic","QM_strategy")
summary_table_raw_values$System <- "Grasslands"
# Move 'System' to 1st position
summary_table_raw_values <- summary_table_raw_values %>% select(System, everything())
# head(summary_table_raw_values)
# Re-order table based on 'Climatic parameter'
summary_table_raw_values$Climatic_parameter <- factor(summary_table_raw_values$Climatic_parameter,
    levels = c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10","precipitation"))
summary_table_raw_values <- summary_table_raw_values[order(summary_table_raw_values$Climatic_parameter),]
# worked


### 16/04/26: Add forest data
files <- dir()[!grepl("_anoms_",dir())]
files <- files[grepl("_forests_",files)]
# files # gut
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
# Rbind
table <- bind_rows(res); rm(res); gc()
summary_table_raw_values_forests <- data.frame(
    table %>%
    group_by(variable,stat,qm_method) %>%
    summarize(
        MAE = round(mean(mae, na.rm = TRUE),3), MAE_sd = round(sd(mae, na.rm = TRUE),3),
        MBE = round(mean(mbe, na.rm = TRUE),3), MBE_sd = round(sd(mbe, na.rm = TRUE),3),
        Rho = round(mean(corr, na.rm = TRUE),3), Rho_sd = round(sd(corr, na.rm = TRUE),3),
        RMSE = round(mean(rmse, na.rm = TRUE),3), RMSE_sd = round(sd(rmse, na.rm = TRUE),3),
        RIA = round(mean(d1r, na.rm = TRUE),3), RIA_sd = round(sd(d1r, na.rm = TRUE),3),
    )
)
# summary_table_raw_values_forests
summary_table_raw_values_forests$MAE <- paste(summary_table_raw_values_forests$MAE," ± ", summary_table_raw_values_forests$MAE_sd, sep = "")
summary_table_raw_values_forests$MBE <- paste(summary_table_raw_values_forests$MBE," ± ", summary_table_raw_values_forests$MBE_sd, sep = "")
summary_table_raw_values_forests$Rho <- paste(summary_table_raw_values_forests$Rho," ± ", summary_table_raw_values_forests$Rho_sd, sep = "")
summary_table_raw_values_forests$RMSE <- paste(summary_table_raw_values_forests$RMSE," ± ", summary_table_raw_values_forests$RMSE_sd, sep = "")
summary_table_raw_values_forests$RIA <- paste(summary_table_raw_values_forests$RIA," ± ", summary_table_raw_values_forests$RIA_sd, sep = "")
drop.cols <- colnames(summary_table_raw_values_forests)[c(5,7,9,11,13)]
summary_table_raw_values_forests <- summary_table_raw_values_forests %>% select(-drop.cols)
colnames(summary_table_raw_values_forests)[c(1,2,3)] <- c("Climatic_parameter","Daily_statistic","QM_strategy")
summary_table_raw_values_forests$System <- "Forests"
summary_table_raw_values_forests <- summary_table_raw_values_forests %>% select(System, everything())

### Rbind grasslands' table and forests' table 
TableSX.A <- rbind(summary_table_raw_values,summary_table_raw_values_forests)
TableSX.A$Raw_vs_anomalies <- "Raw"


### 16/04/26: Same as above, but for the evaluation metrics based on anomalies to the monthly means
files <- dir()[grepl("_anoms_",dir())]
files <- files[!grepl("_forests_",files)]
# files # gut
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
# Rbind
table <- bind_rows(res); rm(res); gc()
# dim(table) ; head(table) ; summary(table)

summary_table_anoms_values <- data.frame(
    table %>%
    group_by(variable,stat,qm_method) %>%
    summarize(
        MAE = round(mean(mae, na.rm = TRUE),3), MAE_sd = round(sd(mae, na.rm = TRUE),3),
        MBE = round(mean(mbe, na.rm = TRUE),3), MBE_sd = round(sd(mbe, na.rm = TRUE),3),
        Rho = round(mean(corr, na.rm = TRUE),3), Rho_sd = round(sd(corr, na.rm = TRUE),3),
        RMSE = round(mean(rmse, na.rm = TRUE),3), RMSE_sd = round(sd(rmse, na.rm = TRUE),3),
        RIA = round(mean(d1r, na.rm = TRUE),3), RIA_sd = round(sd(d1r, na.rm = TRUE),3),
    )
)
# Check
# head(summary_table_raw_values); dim(summary_table_raw_values)

### 20/04/26: Check scores of individual QM strategies against other -> trying to rank 'mw' and 'anomalies'
#summary(summary_table_anoms_values[summary_table_anoms_values$qm_method == "global", c("MAE","RMSE","RIA")])
#summary(summary_table_anoms_values[summary_table_anoms_values$qm_method == "monthly", c("MAE","RMSE","RIA")])
summary(summary_table_anoms_values[summary_table_anoms_values$qm_method == "anomalies", c("MAE","MBE","RMSE","RIA","Rho")])
summary(summary_table_anoms_values[summary_table_anoms_values$qm_method == "mw", c("MAE","MBE","RMSE","RIA","Rho")])
### Still hard to tell! 

# Merge mean and sd values in one column 
summary_table_anoms_values$MAE <- paste(summary_table_anoms_values$MAE," ± ", summary_table_anoms_values$MAE_sd, sep = "")
summary_table_anoms_values$MBE <- paste(summary_table_anoms_values$MBE," ± ", summary_table_anoms_values$MBE_sd, sep = "")
summary_table_anoms_values$Rho <- paste(summary_table_anoms_values$Rho," ± ", summary_table_anoms_values$Rho_sd, sep = "")
summary_table_anoms_values$RMSE <- paste(summary_table_anoms_values$RMSE," ± ", summary_table_anoms_values$RMSE_sd, sep = "")
summary_table_anoms_values$RIA <- paste(summary_table_anoms_values$RIA," ± ", summary_table_anoms_values$RIA_sd, sep = "")
# Drop '_sd' columns: 5, 7, 9, 11, 13
drop.cols <- colnames(summary_table_anoms_values)[c(5,7,9,11,13)]
summary_table_anoms_values <- summary_table_anoms_values %>% select(-drop.cols)
# Adjust colnames 
colnames(summary_table_anoms_values)[c(1,2,3)] <- c("Climatic_parameter","Daily_statistic","QM_strategy")
summary_table_anoms_values$System <- "Grasslands"
# Move 'System' to 1st position
summary_table_anoms_values <- summary_table_anoms_values %>% select(System, everything())
# head(summary_table_raw_values)
# Re-order table based on 'Climatic parameter'
summary_table_anoms_values$Climatic_parameter <- factor(summary_table_anoms_values$Climatic_parameter,
    levels = c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10","precipitation"))
summary_table_anoms_values <- summary_table_anoms_values[order(summary_table_anoms_values$Climatic_parameter),]
# worked

### 16/04/26: Add forest data
files <- dir()[grepl("_anoms_",dir())]
files <- files[grepl("_forests_",files)]
# files # gut
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
# Rbind
table <- bind_rows(res); rm(res); gc()
summary_table_anoms_forests <- data.frame(
    table %>%
    group_by(variable,stat,qm_method) %>%
    summarize(
        MAE = round(mean(mae, na.rm = TRUE),3), MAE_sd = round(sd(mae, na.rm = TRUE),3),
        MBE = round(mean(mbe, na.rm = TRUE),3), MBE_sd = round(sd(mbe, na.rm = TRUE),3),
        Rho = round(mean(corr, na.rm = TRUE),3), Rho_sd = round(sd(corr, na.rm = TRUE),3),
        RMSE = round(mean(rmse, na.rm = TRUE),3), RMSE_sd = round(sd(rmse, na.rm = TRUE),3),
        RIA = round(mean(d1r, na.rm = TRUE),3), RIA_sd = round(sd(d1r, na.rm = TRUE),3),
    )
)

### 20/04/26: Check between 'anomalies' and 'mw'
summary(summary_table_anoms_forests[summary_table_anoms_forests$qm_method == "anomalies", c("MAE","MBE","RMSE","RIA")])
summary(summary_table_anoms_forests[summary_table_anoms_forests$qm_method == "mw", c("MAE","MBE","RMSE","RIA")])
### Same same...

# summary_table_raw_anoms_forests
summary_table_anoms_forests$MAE <- paste(summary_table_anoms_forests$MAE," ± ", summary_table_anoms_forests$MAE_sd, sep = "")
summary_table_anoms_forests$MBE <- paste(summary_table_anoms_forests$MBE," ± ", summary_table_anoms_forests$MBE_sd, sep = "")
summary_table_anoms_forests$Rho <- paste(summary_table_anoms_forests$Rho," ± ", summary_table_anoms_forests$Rho_sd, sep = "")
summary_table_anoms_forests$RMSE <- paste(summary_table_anoms_forests$RMSE," ± ", summary_table_anoms_forests$RMSE_sd, sep = "")
summary_table_anoms_forests$RIA <- paste(summary_table_anoms_forests$RIA," ± ", summary_table_anoms_forests$RIA_sd, sep = "")
drop.cols <- colnames(summary_table_anoms_forests)[c(5,7,9,11,13)]
summary_table_anoms_forests <- summary_table_anoms_forests %>% select(-drop.cols)
colnames(summary_table_anoms_forests)[c(1,2,3)] <- c("Climatic_parameter","Daily_statistic","QM_strategy")
summary_table_anoms_forests$System <- "Forests"
summary_table_anoms_forests <- summary_table_anoms_forests %>% select(System, everything())

### Rbind grasslands' table and forests' table 
TableSX.B <- rbind(summary_table_anoms_values, summary_table_anoms_forests)
TableSX.B$Raw_vs_anomalies <- "Anomalies"

### Rbind TableSX.A & TableSX.B
TableSX <- rbind(TableSX.A,TableSX.B)
# Check before saving
head(TableSX); tail(TableSX) ; str(TableSX) ; summary(TableSX)
### Save as .csv or .xlx
write.table(TableSX, file = "TableSX_QM_eval_metrics_16.04.26.txt", sep = ";")


### Now, plot distribution of evaluation metrics across QM strategies for grassland EPs only -> making Fig. 4
### For Fig. 4, let's focus on the metrics based on the raw values 
files <- dir()[grepl("_anoms_",dir())]
files <- files[!grepl("_forests_",files)]
# files # gut
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
# Rbind
table <- bind_rows(res); rm(res); gc()
# summary(table)
colnames(table)[c(2:6)] <- c("MBE","MAE","RMSE","Rho","RIA")
# Re-order ddf per 'variable'
table$variable <- factor(table$variable, levels = c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10","precipitation"))
# Re-label 'mw' qm_method
table$qm_method <- as.character(table$qm_method)
table[table$qm_method == "mw","qm_method"] <- "mowing-window"
# Re-order ddf per 'qm_method'
table$qm_method <- factor(table$qm_method, levels = c("global","monthly","mowing-window","anomalies"))

# Fig. 4 should show the distribution of an evaluation metric per QM stragegy, with different facets per stat & vars 
plot <- ggplot(table, aes(x = factor(qm_method), y = RIA, fill = factor(qm_method))) + 
    geom_boxplot(colour = "black") + scale_fill_manual(name = "",
        values = c("#006837","#1a9850","#66bd63","#a6d96a","#d9ef8b")) +
    xlab("") + ylab("Refined index of agreement (RIA) based on anomalies to the monthly mean") + theme_bw() + theme(axis.text.x = element_blank()) + 
    facet_wrap(factor(variable)~factor(stat), scales = "free_y")
# Save
ggsave(plot = plot, filename = "Fig.4_RIA_QM_strat_anoms_free_y_16.04.26.jpg", dpi = 300, width = 7.5, height = 8.5)



### ------------------------------------------------------------------------------------------------------------

### 01/07/26: Figure 5 - Evaluation metrics of the forest microclimate reconstructions
### Keep same style as Fig. 4

setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/prediction_tables/evaluation_tables")

files <- dir()[!grepl("_anoms_",dir())]
# files # gut
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
# Rbind
table <- bind_rows(res); rm(res); gc()
# dim(table) ; head(table) ; summary(table)

### Summarize mean & sd values
summary_table_raw_values <- data.frame(
    table %>%
    group_by(variable,stat) %>%
    summarize(
        MAE = round(mean(mae, na.rm = TRUE),3), MAE_sd = round(sd(mae, na.rm = TRUE),3),
        MBE = round(mean(mbe, na.rm = TRUE),3), MBE_sd = round(sd(mbe, na.rm = TRUE),3),
        Rho = round(mean(corr, na.rm = TRUE),3), Rho_sd = round(sd(corr, na.rm = TRUE),3),
        RMSE = round(mean(rmse, na.rm = TRUE),3), RMSE_sd = round(sd(rmse, na.rm = TRUE),3),
        RIA = round(mean(d1r, na.rm = TRUE),3), RIA_sd = round(sd(d1r, na.rm = TRUE),3),
    )
)
# Check
# head(summary_table_raw_values); dim(summary_table_raw_values)


# Merge mean and sd values in one column 
summary_table_raw_values$MAE <- paste(summary_table_raw_values$MAE," ± ", summary_table_raw_values$MAE_sd, sep = "")
summary_table_raw_values$MBE <- paste(summary_table_raw_values$MBE," ± ", summary_table_raw_values$MBE_sd, sep = "")
summary_table_raw_values$Rho <- paste(summary_table_raw_values$Rho," ± ", summary_table_raw_values$Rho_sd, sep = "")
summary_table_raw_values$RMSE <- paste(summary_table_raw_values$RMSE," ± ", summary_table_raw_values$RMSE_sd, sep = "")
summary_table_raw_values$RIA <- paste(summary_table_raw_values$RIA," ± ", summary_table_raw_values$RIA_sd, sep = "")

# Drop '_sd' columns: 5, 7, 9, 11, 13
drop.cols <- colnames(summary_table_raw_values)[c(4,6,8,10,12)]
summary_table_raw_values <- summary_table_raw_values %>% select(-drop.cols)

# Adjust colnames 
colnames(summary_table_raw_values)[c(1,2)] <- c("Climatic_parameter","Daily_statistic")
summary_table_raw_values$System <- "Forests"
# Move 'System' to 1st position
summary_table_raw_values <- summary_table_raw_values %>% select(System, everything())
# head(summary_table_raw_values)
# Re-order table based on 'Climatic parameter'
summary_table_raw_values$Climatic_parameter <- factor(summary_table_raw_values$Climatic_parameter,
    levels = c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10"))
summary_table_raw_values <- summary_table_raw_values[order(summary_table_raw_values$Climatic_parameter),]
# worked


### Same as above, but for the evaluation metrics based on anomalies to the monthly means
files <- dir()[grepl("_anoms_",dir())]
# files # gut
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
# Rbind
table <- bind_rows(res); rm(res); gc()
# dim(table) ; head(table) ; summary(table)

summary_table_anoms_values <- data.frame(
    table %>%
    group_by(variable,stat) %>%
    summarize(
        MAE = round(mean(mae, na.rm = TRUE),3), MAE_sd = round(sd(mae, na.rm = TRUE),3),
        MBE = round(mean(mbe, na.rm = TRUE),3), MBE_sd = round(sd(mbe, na.rm = TRUE),3),
        Rho = round(mean(corr, na.rm = TRUE),3), Rho_sd = round(sd(corr, na.rm = TRUE),3),
        RMSE = round(mean(rmse, na.rm = TRUE),3), RMSE_sd = round(sd(rmse, na.rm = TRUE),3),
        RIA = round(mean(d1r, na.rm = TRUE),3), RIA_sd = round(sd(d1r, na.rm = TRUE),3),
    )
)
# Check
# head(summary_table_raw_values); dim(summary_table_raw_values)

# Merge mean and sd values in one column 
summary_table_anoms_values$MAE <- paste(summary_table_anoms_values$MAE," ± ", summary_table_anoms_values$MAE_sd, sep = "")
summary_table_anoms_values$MBE <- paste(summary_table_anoms_values$MBE," ± ", summary_table_anoms_values$MBE_sd, sep = "")
summary_table_anoms_values$Rho <- paste(summary_table_anoms_values$Rho," ± ", summary_table_anoms_values$Rho_sd, sep = "")
summary_table_anoms_values$RMSE <- paste(summary_table_anoms_values$RMSE," ± ", summary_table_anoms_values$RMSE_sd, sep = "")
summary_table_anoms_values$RIA <- paste(summary_table_anoms_values$RIA," ± ", summary_table_anoms_values$RIA_sd, sep = "")

# Drop '_sd' columns: 5, 7, 9, 11, 13
drop.cols <- colnames(summary_table_anoms_values)[c(4,6,8,10,12)]
summary_table_anoms_values <- summary_table_anoms_values %>% select(-drop.cols)

# Adjust colnames 
colnames(summary_table_anoms_values)[c(1,2)] <- c("Climatic_parameter","Daily_statistic")
summary_table_anoms_values$System <- "Forests"
# Move 'System' to 1st position
summary_table_anoms_values <- summary_table_anoms_values %>% select(System, everything())
# head(summary_table_raw_values)
# Re-order table based on 'Climatic parameter'
summary_table_anoms_values$Climatic_parameter <- factor(summary_table_anoms_values$Climatic_parameter,
    levels = c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10"))
summary_table_anoms_values <- summary_table_anoms_values[order(summary_table_anoms_values$Climatic_parameter),]
# worked


### Rbind TableSX.A & TableSX.B
TableS4 <- rbind(summary_table_raw_values, summary_table_anoms_values)
# Check before saving
head(TableS4); tail(TableS4) ; str(TableS4) ; summary(TableS4)
### Save as .csv or .xlx
write.table(TableS4, file = "TableS4_forest_recontruct_eval_metrics_01.07.26.txt", sep = ";")


### Now, plot distribution of evaluation metrics across QM strategies for grassland EPs only -> making Fig. 5
### For Fig. 5, let's focus on the metrics based on the raw values 
files <- dir()[grepl("_anoms_",dir())]
# files # gut
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
# Rbind
table <- bind_rows(res)
rm(res); gc()
# summary(table)
colnames(table)[c(3:7)] <- c("MBE","MAE","RMSE","Rho","RIA")
# Re-order ddf per 'variable'
table$variable <- factor(table$variable, levels = c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10"))

# Fig. 5 should show the distribution of an evaluation metric per QM stragegy, with different facets per stat & vars 
plot1 <- ggplot(table, aes(x = factor(variable), y = RIA, fill = factor(variable))) + 
    geom_boxplot(colour = "black") + scale_fill_manual(name = "",
        values = c("#006837","#1a9850","#66bd63","#a6d96a","#d9ef8b","#ffffbf")) +
    xlab("") + ylab("Refined index of agreement (RIA)\nbased on anomalies") + theme_bw() +
    theme(axis.text.x = element_blank()) + facet_wrap(.~factor(stat), scales = "free_y")
# Save
ggsave(plot = plot1, filename = "Fig.5_RIA_anoms_free_y_01.07.26.jpg", dpi = 300, width = 6, height = 3.5)


### Actually, there is space to add the RIA scores based on raw observed values - Make the same figure and 
### combine them in a panel
files <- dir()[!grepl("_anoms_",dir())]
# files # gut
res <- lapply(files, function(f) { d <- get(load(f)); return(d) })
# Rbind
table <- bind_rows(res)
rm(res); gc()
# summary(table)
colnames(table)[c(3:7)] <- c("MBE","MAE","RMSE","Rho","RIA")
# Re-order ddf per 'variable'
table$variable <- factor(table$variable, levels = c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10"))

# Fig. 5 should show the distribution of an evaluation metric per QM stragegy, with different facets per stat & vars 
plot2 <- ggplot(table, aes(x = factor(variable), y = RIA, fill = factor(variable))) + 
    geom_boxplot(colour = "black") + scale_fill_manual(name = "",
        values = c("#006837","#1a9850","#66bd63","#a6d96a","#d9ef8b","#ffffbf")) +
    xlab("") + ylab("Refined index of agreement (RIA)\nbased on raw observed values") + theme_bw() +
    theme(axis.text.x = element_blank()) + facet_wrap(.~factor(stat), scales = "free_y")
# Save
ggsave(plot = plot2, filename = "Fig.5_RIA_raw_free_y_01.07.26.jpg", dpi = 300, width = 6, height = 3.5)

### Combine
fig5 <- ggarrange(plot1, plot2, align = "hv", nrow = 2, ncol = 1, common.legend = TRUE, legend = "right")
ggsave(plot = fig5, filename = "Fig.5_01.07.26.jpg", dpi = 300, width = 5.5, height = 5.75)


### ------------------------------------------------------------------------------------------------------------

### 06/07/26: Making Figure 6: Illustration of the comparison between the fixed-baseline vs. the shifting baseline period
### Using theoretical temperature anomalies and Fig. 2 of Smith et al. (2025) as an example 
### (https://www.sciencedirect.com/science/article/pii/S0079661124002106?via%3Dihub)

library("ggplot2")
library("dplyr")
library("patchwork")

# Critical to reproduce figure below
set.seed(42)

### Create a synthetic time series of annual Ta anomalies
# Warming trend: +1.7 °C over 75 years (linear, representative of Germany)
years <- 1950:2024
n <- length(years)
yr0 <- 1950
yr_range <- 2024 - yr0
# Add a long-term warming trend
trend <- 1.7 * (years - yr0) / yr_range

# AR1 interannual variability (rho = 0.65)
ar1 <- numeric(n)
ar1[1] <- 0
for (i in 2:n) {
  ar1[i] <- 0.65 * ar1[i - 1] + rnorm(1, 0, 0.30)
}

# Annual temperature anomaly (°C, relative to 1950–1980 mean)
anom <- trend + ar1
# Bind a ddf
df_annual <- data.frame(year = years, anom = anom)

# Fixed baseline threshold (1950–1980, 95th percentile)
fixed_thresh <- quantile(df_annual$anom[df_annual$year <= 1980], probs = 0.90) # 0.981 °C

# Shifting baseline (30‑yr trailing window, 95th percentile) ──
shift_thresh <- rep(NA_real_, n)

# Compute shifting thresholds
for(i in seq_len(n)) {
    y <- years[i]
    win <- df_annual$anom[df_annual$year >= (y - 30) & df_annual$year < y]
    if (length(win) >= 30) {
        shift_thresh[i] <- quantile(win, probs = 0.90)
    } # eo if loop 
} # eo for loop 

df_annual$shift_thresh <- shift_thresh


# Identify ECEs based on thresholds
df_annual <- df_annual %>%
    mutate(
        ece_fixed = anom > fixed_thresh,
        ece_shift = !is.na(shift_thresh) & anom > shift_thresh
    )

# Define color palette to use in the plots
col_ts <- "#2c3e50"             # dark navy — time series line/points
col_fixed <- "#c0392b"          # red       — fixed threshold line
col_shift <- "#2980b9"          # blue      — shifting threshold line
col_ece_f <- "#e74c3c"          # red fill  — ECEs under fixed baseline
col_ece_s <- "#3498db"          # blue fill — ECEs under shifting baseline
col_ref <- "#bdc3c7"            # grey      — reference period shading
col_shift_win <- "#aed6f1"      # light blue— shifting window shading

# Add helper columns before plotting
df_annual <- df_annual %>%
    mutate(
        fixed_ymin = ifelse(ece_fixed, fixed_thresh, NA_real_),
        fixed_ymax = ifelse(ece_fixed, anom, NA_real_),
        shift_ymin = ifelse(ece_shift, shift_thresh, NA_real_),
        shift_ymax = ifelse(ece_shift, anom, NA_real_)
  )


### 6. Panel (a): Fixed baseline 
pA <- ggplot(df_annual, aes(x = year)) +
    annotate("rect", xmin = 1950, xmax = 1980, ymin = -Inf, ymax = Inf, fill = col_ref, alpha = 0.35) +
    annotate("text", x = 1965, y = 2.55, label = "Reference period\n(1950–1980)", 
        size = 3, colour = "black", fontface = "italic", lineheight = 0.9) +
    geom_rect(data = df_annual %>% filter(ece_fixed), aes(xmin = year - 0.45, xmax = year + 0.45, ymin = fixed_ymin, ymax = fixed_ymax),
        inherit.aes = FALSE, fill = col_ece_f, alpha = 0.55, colour = "#c0392b") + # ECE shading: ribbon between fixed threshold and anomaly
    geom_line(aes(y = anom), colour = col_ts, linewidth = 0.6, alpha = 0.85) +
    geom_point(aes(y = anom), colour = col_ts, size = 1.2, alpha = 0.7) +
    geom_hline(yintercept = fixed_thresh, colour = col_fixed, linewidth = 0.9, linetype = "dashed") + # Fixed threshold (dashed red)
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4, linetype = "dotted") + # Zero reference line
    scale_x_continuous(breaks = seq(1950, 2020, 10), expand = c(0.01, 0.01)) +
    scale_y_continuous(breaks = seq(-1, 3, 0.5)) +
    coord_cartesian(ylim = c(-1.2, 2.85)) +
    labs(title = "a) Fixed baseline approach (1950–1980)", x = NULL, y = "Temperature anomaly [°C]") +
    theme_classic(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", size = 11, hjust = 0),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 9),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        plot.margin = margin(8, 12, 4, 8)
    )

# Removed: 
# annotate("text", x = 2024, y = fixed_thresh + 0.10, label = paste0("95th percentile (fixed, ", round(fixed_thresh, 2), " °C)"), hjust = 1, size = 2.8, colour = col_fixed, fontface = "italic") +

# Panel (b): Shifting baseline 
shift_thresh_last <- tail(na.omit(df_annual$shift_thresh), 1)

pB <- ggplot(df_annual, aes(x = year)) +
  annotate("rect", xmin = 1950, xmax = 1980, ymin = -Inf, ymax = Inf, fill = col_ref, alpha = 0.35) +
  # Bracket annotations illustrating window concept
  annotate("segment", x = 1950, xend = 1980, y = -1.55, yend = -1.55, colour = col_shift, linewidth = 0.6, 
    arrow = arrow(ends = "both", length = unit(0.12, "cm"), type = "open")) +
  annotate("text", x = 1965, y = -1.64, label = "1950-1980 trailing window", size = 2.6, colour = col_shift, fontface = "italic") +
  annotate("segment", x = 1970, xend = 2000, y = -1.8, yend = -1.8, colour = col_shift, linewidth = 0.6, 
    arrow = arrow(ends = "both", length = unit(0.12, "cm"), type = "open")) +
  annotate("text", x = 1985, y = -1.89, label = "1970-2000 trailing window", size = 2.6, colour = col_shift, fontface = "italic") +
  annotate("segment", x = 1994, xend = 2024, y = -2.05, yend = -2.05, colour = col_shift, linewidth = 0.6,
    arrow = arrow(ends = "both", length = unit(0.12, "cm"), type = "open")) +
  annotate("text", x = 2009, y = -2.14, label = "1994-2024 trailing window", size = 2.6, colour = col_shift, fontface = "italic") +
  geom_rect(data = df_annual %>% filter(ece_shift), aes(xmin = year - 0.45, xmax = year + 0.45,
    ymin = shift_ymin, ymax = shift_ymax), inherit.aes = FALSE,fill = col_ece_s, alpha = 0.55, colour = "#2980b9") +
  geom_line(aes(y = anom), colour = col_ts, linewidth = 0.6, alpha = 0.85) +
  geom_point(aes(y = anom), colour = col_ts, size = 1.2, alpha = 0.7) +
  geom_line(aes(y = shift_thresh), colour = col_shift, linewidth = 0.95, linetype = "dashed", na.rm = TRUE) + # Shifting threshold curve (dashed blue)
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4, linetype = "dotted") +
  scale_x_continuous(breaks = seq(1950, 2020, 10), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(-1, 3.2, 0.5)) +
  coord_cartesian(ylim = c(-2.10, 2.85)) +
  labs(title = "b) Shifting baseline approach (30-year trailing windows)", x = "Year", y = "Temperature anomaly [°C]") +
  theme_classic(base_size = 11) +
  theme(
        plot.title = element_text(face = "bold", size = 11, hjust = 0),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        plot.margin = margin(4,12,8,8)
    )

### 07/07/26: To fix: 
# - Remove background on b) -> FIXED
# - Check why shifting baseline series starts in 1970 instead of 1980... -> FIXED mistake at line 741 
# - Add third 30 year period for shifting approach and complete captions with full period: 
#       1950-1980 window
#       1965-1995 window
#       1994-2024 window

### Combine panels with patchwork 
fig.6 <- pA / pB + plot_annotation(
    theme = theme(
      plot.caption = element_text(
        size = 7.5, colour = "grey45", lineheight = 1.25,
        hjust = 0, margin = margin(t = 6)
      )
    )
  )



### Save 
# Adjust path as needed; here we save in the current working directory
setwd("/Users/fabiobenedetti/Desktop/")
ggsave("Fig.6.jpg", plot = fig.6, width = 6, height = 7.5, dpi = 300)
ggsave("Fig.6.png", plot = fig.6, width = 6, height = 7.5, dpi = 300, bg = "white")

### ------------------------------------------------------------------------------------------------------------

### 13-14/07/26: Figure 7: Distribution plots of long-term climatic trends (Sen's slopes) for both grasslands and forests

### Note: Daily climate time series are strongly autocorrelated, which inflates the Mann-Kendall test statistic and 
### increases false-positive rates. Simply running a standard MK test on raw or naïvely deseasonalized daily data may 
### produce overconfident p-values and should thus be avoided.

### Let's perform a 2-step approach instead: 
# - Step 1) Deseasonalization. 
# Before testing, remove the seasonal cycle from each site-level time series. The standard method is STL decomposition 
# (stats::stl() in R), which extracts a loess-based seasonal component and leaves a residual trend + remainder component.
# Subtracting the seasonal component from the original series yields a seasonality-free time series that retains the long-term signal

# - Step 2) Modified Mann-Kendall + Sen's Slope. 
# Run a Modified Mann-Kendall test that corrects for residual autocorrelation on the deseasonalized daily residuals

# install.packages("modifiedmk")
library("modifiedmk")
library("dplyr")
library("reshape2")
library("tidyr")
library("stringr")
library("zoo")
library("patchwork")
library("parallel")

### 14/07/26: Apply deseasonalize_stl() to each EP and for each parameter

### Helper FUN: Sen's slope 95% CI via Kendall package because modifiedmk::mmkh() does NOT return CIs
# We compute them from the 'S' statistic variance using the standard formula.
sen_ci <- function(x, alpha = 0.05) {
    
    n <- length(x)
    # All pairwise slopes
    idx <- combn(n, 2)
    slopes <- (x[idx[2, ]] - x[idx[1, ]]) / (idx[2, ] - idx[1, ])
    slopes <- sort(slopes)
  
    # Variance of S (Mann-Kendall S statistic), no ties assumed
    var_s <- n * (n - 1) * (2 * n + 5) / 18
    z_a <- qnorm(1 - alpha / 2)
    C_a <- z_a * sqrt(var_s)
  
    M1 <- round((length(slopes) - C_a) / 2)
    M2 <- round((length(slopes) + C_a) / 2) + 1
  
    # Guard against out-of-bounds indices
    M1 <- max(M1, 1)
    M2 <- min(M2, length(slopes))
  
    c(ci_lower = slopes[M1], ci_upper = slopes[M2])

} # eo HELPER - sen_ci

setwd("/home/fbenedetti/BE-EXTREME_files/v1_08.07.26/Reconstructions")
parameters <- c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10","precipitation")

# To test lapply below: 
# param <- "precipitation"

trends <- lapply(parameters, function(param) {

        # Useless message
        message(paste("Loading data for ",param, sep = ""))
        
        # Load the data 
        files <- dir()[grepl(param,dir())]

        if( param == "SM_10" ) {

            d <- read.csv(files[1], header = TRUE)

        } else if( param == "precipitation" ) {

            d_forst <- read.csv(files[1], header = TRUE)
            d2_grass <- read.csv(files[2], header = TRUE)
            
            # Bind and discard the individual files
            # colnames(d_forst) ; colnames(d2_grass)
            names <- colnames(d_forst)[c(2:7,9:13)] # names

            d <- bind_rows(d2_grass[,names], d_forst[,names])

            rm(d_forst, d2_grass); gc()

        } else {

            d_forst <- read.csv(files[1], header = TRUE)
            d2_grass <- read.csv(files[2], header = TRUE)
            
            # Bind and discard the individual files
            names <- colnames(d_forst)[c(2:7,11:15)]

            d <- bind_rows(d2_grass[,names],d_forst[,names])

            rm(d_forst, d2_grass); gc()

        } # eo if else loop - param

        # OK, now, for each EP and Statistic x Parameter -> deseasonalize and compute Sen's slope
        d$Var <- paste(d$Statistic, d$Parameter, sep = "_")
        stats <- unique(d$Var) # stats
        # s <- stats[1]; s
        
        res_stats <- lapply(stats, function(s) {

                message(paste("         Performing deseasonalization and Sen's slope test for ",s, sep = ""))
                subset <- d[d$Var == s,] # dim(subset); head(subset)

                # Perform EP-specific deseasonalization with mclapply
                plots <- unique(subset$EP) # plots
                # p <- plots[4]

                res <- mclapply(plots, function(p) {
                        
                        message(p)

                        sub <- subset[subset$EP == p,]

                        # Use 'sub' object for deseasonalization
                        # Build a ts object: daily with frequency = 365
                        # NOTE: use 365.25 only if leap years cause artefacts; 365 is standard
                        ts_obj <- ts(sub$Final, frequency = 365) # class(ts_obj); str(ts_obj)
  
                        ### STEP 1: STL deseasonalization 
                        # s.window = "periodic" -> seasonal component is constant across years
                        # Alternatively: s.window = 11 (or any odd number ≥ 7) allows the
                        # seasonal shape to evolve slowly over time — more flexible.
                        # robust = TRUE to down-weight outliers during loess fitting
                        
                        if( sum(is.na(ts_obj)) > 0 ) {
                            # Interpolate short gaps; leave long gaps as is if needed
                            ts_obj <- na.approx(ts_obj, na.rm = FALSE)
                            # If any NAs remain (leading/trailing NAs), fill with local mean
                            if( any(is.na(ts_obj)) ) {
                                ts_obj <- na.fill(ts_obj, "extend")
                            } # eo 2nd if loop - NAs
                        } # eo 1st if loop - NAs
                        
                        stl_fit <- stl(ts_obj, s.window  = "periodic", robust = TRUE)
                        # str(stl_fit); head(stl_fit$time.series)
  
                        # Deseasonalized series = original – seasonal component
                        # "We decomposed daily time series into seasonal, trend and remainder components using
                        # STL. For trend detection, we removed the seasonal component and applied a modified 
                        # Mann–Kendall test with Sen’s slope to the deseasonalized daily series."
                        # We prefer to use the "original – seasonal" approach for the MK/Sen analysis 
                        # (instead of using the 'stl_fit$time.series[,"strend"]' directly) because 
                        # it keeps the trend test tied directly to the deseasonalized daily data and avoids extra smoothing assumptions. 
                        deseas <- sub$Final - stl_fit$time.series[,"seasonal"] # summary(deseas)
                        sub$Deseasoned <- as.numeric(deseas)

                        ### STEP 2: Modified Mann-Kendall Test + Sen's Slope
                        ### This applies modifiedmk::mmkh() (Hamed & Rao 1998 variance correction) 
                        ### to each deseasonalized TS and extracts the corrected p-value, Sen's slope, 
                        ### and 95% confidence interval.
                        # See: https://github.com/cran/modifiedmk/blob/master/R/mmkh.R

                        # x <- as.vector(sub$Deseasoned)
                        # ?mmkh
                        ### ISSUE: The modified MK test is doing more than a simple MK: 
                        #   it detrends the series, computes Sen’s slope, constructs a trend-free series,
                        #   then scans autocorrelation across multiple lags and applies the variance correction
                        #   using only significant lags. All of this is implemented in R, not compiled C, 
                        #   and tends to involve multiple passes over the full vector.
                        #   On a ~27k-element series, those passes and the autocorrelation computations
                        #   become noticeably slow...
                        ### -> Run on monthly means

                        # If precip. -> Compute total deseasonalized precip instead of mean! 
                        if( param == "precipitation" ) {
                            
                            df_mon <- sub %>%
                                group_by(Region, EP, Var, Year, Month) %>%
                                summarise(
                                    mean = sum(Deseasoned, na.rm = TRUE),
                                    .groups = "drop"
                                ) %>% arrange(Year, Month)
                            # summary(df_mon)

                        } else {
                            
                            df_mon <- sub %>%
                                group_by(Region, EP, Var, Year, Month) %>%
                                summarise(
                                    mean = mean(Deseasoned, na.rm = TRUE),
                                    .groups = "drop"
                                ) %>% arrange(Year, Month)
                            # summary(df_mon)
                        }
                        
                        x <- df_mon$mean    
  
                        # Skip if too many NAs
                        if( sum(!is.na(x)) < 30 ) {
                            return(
                                tibble(
                                    sen_slope = NA_real_,
                                    ci_lower = NA_real_,
                                    ci_upper = NA_real_,
                                    tau = NA_real_,
                                    p_value = NA_real_,
                                    p_original = NA_real_,
                                    n_eff = NA_real_
                                ) # eo tibble
                            ) # eo return
                        } # eo if loop 

                        # Replace NAs with interpolated values (mmkh does not accept NAs)
                        if( any(is.na(x)) ) {
                            x <- zoo::na.approx(x, na.rm = FALSE)
                            x[is.na(x)] <- median(x, na.rm = TRUE)
                        } # eo if loop 
                        
                        res.test <- mmkh(x, ci = 0.95)
                        ci <- sen_ci(x, alpha = 0.05)
  
                        t <- tibble(
                            EP = p,
                            Parameter = s,
                            sen_slope = res.test["Sen's slope"],
                            ci_lower = ci["ci_lower"],
                            ci_upper = ci["ci_upper"],
                            tau = res.test["Tau"],
                            p_value = res.test["new P-value"],      # autocorrelation-corrected p
                            p_original = res.test["old P.value"],   # original (uncorrected) p for comparison
                            n_eff = res.test["N/N*"],                # ratio of effective to actual sample size
                            slope_decade = res.test["Sen's slope"]*12*10     # Convert slope units per decade
                        ) # eo tibble

                        return(t)

                    }, mc.cores = 25

                ) # eo mclapply - plots

                # Rbind all new objects
                tbl <- bind_rows(res)
                rm(res, subset); gc()
                # Check
                # head(tbl) ; str(tbl) ; summary(tbl)

                return(tbl)

            } # eo FUN - s 
        
        ) # eo lapply - res - stats
        
        # Rbind 
        tbl_stats <- bind_rows(res_stats)
        rm(res_stats); gc()
        # head(tbl_stats) ; str(tbl_stats) ; summary(tbl_stats)
        # unique(tbl_stats$Parameter)

        # Tell whether it's a forest or a grassland plots depending on 'EP' label
        tbl_stats$System <- NA

        for( ep in unique(tbl_stats$EP) ) {
            if( grepl("W",ep) ) {
                tbl_stats[tbl_stats$EP == ep,"System"] <- "Forest"
            } else {
                tbl_stats[tbl_stats$EP == ep,"System"] <- "Grassland"
            }
        } # eo for loop - ep

        # Relocate 
        tbl_stats <- tbl_stats %>% relocate(System) 

        # Return
        return(tbl_stats)

    } # eo 1st FUN - param

) # eo lapply - parameters

### Rbind all
table_trendz_all <- data.frame(bind_rows(trends))
# Check basic stuff 
dim(table_trendz_all)
str(table_trendz_all)
head(table_trendz_all)
summary(table_trendz_all)
unique(table_trendz_all$Parameter) # 13, good
unique(table_trendz_all$EP) # 300, good

# Clean memory
rm(table_trendz_all); gc()


### Add significance labels (four-tier, as in ESSD papers) 
df_trends <- table_trendz_all %>% 
    mutate(
        sig_stars = case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01  ~ "**",
            p_value < 0.05  ~ "*",
            p_value < 0.1   ~ "°",
            TRUE            ~ "ns"
        )
    )
# Check
head(df_trends)
str(df_trends)
summary(factor(df_trends$sig_stars))
#   *   **  ***    °   ns 
# 105   94 2594   68  479 
# Most trends are highly signif.!

# 14/07/26: Save 'df_trends' somewhere so you do not have to re-compute it
# str(df_trends)
setwd("/home/fbenedetti/BE-EXTREME_files/v1_08.07.26")
save(df_trends, file = "table_long-term_trends_all_EPs+params_14_07_26.Rdata")


### Step 3: Make the plots (Fig. 7) to illustrate distribution of trends across systems, regions and parameters

### 14/07/26: Re-load on local machine
setwd("/Users/fabiobenedetti/Desktop/")
df_trends <- get(load("table_long-term_trends_all_EPs+params_14_07_26.Rdata"))
# str(df_trends)

# Define parameter order (adjust to your actual names)
param_order <- c(
    "max_Ta_10", "min_Ta_10",
    "max_Ta_200","min_Ta_200",
    "max_Ts_05", "min_Ts_05",
    "max_Ts_10", "min_Ts_10",
    "max_Ts_20", "min_Ts_20",
    "max_SM_10", "min_SM_10",
    "total_precipitation"
)

# Human-readable labels with units per decade
param_labels <- c(
    max_Ta_10 = "Max Ta_10",
    min_Ta_10 = "Min Ta_10",
    max_Ta_200 = "Max Ta_200",
    min_Ta_200 = "Min Ta_200",
    max_Ts_05 = "Max Ts_05",
    min_Ts_05 = "Min Ts_05",
    max_Ts_10 = "Max Ts_10",
    min_Ts_10 = "Min Ts_10",
    max_Ts_20 = "Max Ts_20",
    min_Ts_20 = "Min Ts_20",
    max_SM_10 = "Max SM_10",
    min_SM_10 = "Min SM_10",
    total_precipitation = "Total precipitation"
)

df_trends2 <- df_trends %>%
    mutate(
        Parameter = factor(Parameter, levels = param_order),
        System = factor(System, levels = c("Grassland","Forest"))
  )
# head(df_trends2)

# We...forgot to save the 'Region' levels! Sorry. Ler's extract it from 'EP' levels then
# unique(df_trends2$EP)
df_trends2 <- df_trends2 %>%
    mutate(
        Region = case_when(
            str_detect(EP, "^HEG") ~ "HND",
            str_detect(EP, "^SEG") ~ "SCH",
            str_detect(EP, "^AEG") ~ "SWA",
            str_detect(EP, "^HEW") ~ "HND",
            str_detect(EP, "^SEW") ~ "SCH",
            str_detect(EP, "^AEW") ~ "SWA"
    ),
    Region = factor(Region, levels = c("SWA","HND","SCH"))
  )
# summary(factor(df_trends2$Region)) # OK


### Panel A – Overview Heatmap of Median Slopes
# Showing the median slope per Parameter × System × Region, with significance indicated via stars.

# Summarise: median slope per Parameter × System × Region
df_heatmap <- df_trends2 %>%
    group_by(Parameter, System, Region) %>%
    summarise(
        med_slope = median(slope_decade, na.rm = TRUE),
        pct_sig = mean(p_value < 0.05, na.rm = TRUE),
        sig_label = case_when(
            mean(p_value < 0.001, na.rm = TRUE) > 0.5 ~ "***",
            mean(p_value < 0.01,  na.rm = TRUE) > 0.5 ~ "**",
            mean(p_value < 0.05,  na.rm = TRUE) > 0.5 ~ "*",
            mean(p_value < 0.10,  na.rm = TRUE) > 0.5 ~ "°",
            TRUE ~ ""
        ), .groups = "drop"
    ) %>%
    mutate(
        sys_region = interaction(System, Region, sep = " · "),
        sys_region = factor(sys_region)
  )
# Check
# head(df_heatmap)
# str(df_heatmap)

# Symmetric color scale around zero
lim <- max(abs(df_heatmap$med_slope), na.rm = TRUE) # lim

p_heatmap <- ggplot(df_heatmap, aes(x = sys_region, y = fct_rev(Parameter), fill = med_slope)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sig_label), size = 3.0, color = "black") +
    scale_fill_gradient2(
        low = "#3288bd",
        mid = "white",
        high = "#d53e4f",
        midpoint = 0,
        limits = c(-lim, lim),
        name = "Median Sen's slope\n(units/decade)"
    ) + scale_y_discrete(labels = param_labels) +
    labs(x = NULL, y = NULL, title = "Summary of the long-term trends across\nsystems and parameters (1950–2024)") +
    theme_bw(base_size = 10) +
    theme(
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        legend.position = "right",
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 10)
    )


### Panel B – Faceted Dot/Forest Plot of Site-Level Slopes
### To show the distribution of slopes and CIs across EPs for each Parameter, separated by System and Region.

# Build a y-axis label combining Region and System
#df_trends_plot <- df_trends2 %>%
#    mutate(
#        row_label = paste0(Region, " · ", System),
#        row_label = factor(row_label)
#    )
# head(df_trends_plot)
# str(df_trends_plot)

#p_forest <- ggplot(df_trends_plot, aes(x = slope_decade, y = row_label, color = System, shape = Region)) +
#    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
#    geom_errorbarh(aes(xmin = ci_lower * 12 * 10, xmax = ci_upper * 12 * 10), 
#        height = 0.25, linewidth = 0.35, alpha = 0.6) +
#    geom_point(size = 1.4, alpha = 0.85) +
#    geom_text(data = df_trends_plot %>% filter(p_value < 0.05), aes(label = sig_stars),
#        hjust = -0.4, size = 2.5, color = "black", show.legend = FALSE
#    ) +
#    facet_wrap(~ Parameter, scales = "free_x", ncol = 4, labeller = labeller(Parameter = param_labels)) +
#    scale_color_manual(name = "System", values = c("Forest" = "#3c8e7a", "Grassland" = "#f9c346")) +
#    scale_shape_manual(name = "Region", values = c("Region1" = 16, "Region2" = 17, "Region3" = 15)) +
#    labs(x = "Sen's slope (units decade⁻¹)", y = NULL, title = "Site-level trend distributions with 95% CI") +
#    theme_bw(base_size = 9) +
#    theme(
#        strip.text = element_text(size = 7.5, face = "bold"),
#        strip.background = element_rect(fill = "grey92", color = NA),
#        axis.text.y = element_text(size = 7.5),
#        axis.text.x = element_text(size = 7),
#        legend.position = "bottom", legend.box = "horizontal",
#        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
#        panel.grid.major.x = element_blank(), panel.grid.minor   = element_blank(),
#        plot.title = element_text(face = "bold", size = 10)
#    )
# Looks good, but WAY TOO busy

### To add signific. test between forests and grasslands
# install.packages("ggsignif")
library("ggsignif")

### Test alternative option: Violin + Boxplot Overlay

p_violin <- ggplot(df_trends2, aes(x = System, y = slope_decade, fill = System)) +
    geom_violin(trim = TRUE, alpha = 0.6, linewidth = 0.3) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    # Pairwise Forest vs. Grassland comparison per facet
    geom_signif(
        comparisons = list(c("Forest", "Grassland")),
        test = "wilcox.test", # Mann-Whitney U
        test.args = list(exact = FALSE),
        map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "°" = 0.1),
        tip_length = 0.01,
        textsize = 3,
        vjust = 1.7 # nudge label below the bracket
    ) +
    facet_wrap(~ Parameter, scales = "free_y", ncol = 4, labeller = labeller(Parameter = param_labels)) +
    scale_fill_manual(values = c("Forest" = "#3c8e7a", "Grassland" = "#f9c346"), guide = "none") +
    labs(x = NULL, y = "Sen's slope (units/decade)", title = "Distribution of long-term trends across sites") +
    theme_bw(base_size = 9) +
    theme(
        strip.text = element_text(size = 7.5, face = "bold"),
        strip.background = element_rect(fill = "grey92", color = NA),
        axis.text.x = element_text(size = 8),
        legend.position = "bottom", panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10)
    )


### Combine into panel
# library("patchwork")
library("ggpubr")

# With ggarrange? 
fig.7 <- ggarrange(p_heatmap, p_violin, nrow = 1, ncol = 2, align = "hv", widths = c(1,1.2), labels = c("A","B"))
fig.7

### Save as .pdf & .png
ggsave("Fig.7.pdf", plot = fig.7, width = 11, height = 6.25, device = cairo_pdf, dpi = 300)
ggsave("Fig.7.jpg", plot = fig.7, width = 11, height = 6.25, dpi = 300)

### Save 'df_trends2' as .csv 
write.csv(df_trends2, file = "table_long-term_trends_all_EPs+params_14.07.26.csv", sep = ";")

### ------------------------------------------------------------------------------------------------------------

### 15/07/26: Making Figure 8: Long-term dynamics of extreme climate events for both grasslands and forests - TO DO

### Fig. 8 aims to illustrate the increasing total number of ECEs per year (all yearly data) that occurred between 
### 1950 and 2024, in the grasslands and the forests (Let's make a panel of 2 bar plots, with A) Grasslands and B) Forests).
### We would like to use bar charts that show the total N of ECEs per year between 1950 and 2024, by stacking the different 
### types of ECEs on the bars and giving them different colors (e.g., shades or reds for the max temperatures/soil moisture
### extremes and shades of blues for min temperature and total precipitation extremes). Here, all the ECEs will be those based on 
### the fixed baseline approach (to focus on those ECEs driven by long-term climatic changes, i.e., driven by 
### anthropogenic climate change). 

### Then, we will make a Suppl. Fig. S3 that will show exactly the same type of information but based on the 
### shifting baseline approach (from 1980 to 2024 then). 

### Got to ECEs tables dir
setwd("/home/fbenedetti/BE-EXTREME_files/v1_08.07.26/ECE_tables")

# Load the tables
t_grass <- read.csv("table_ECEs+features_fixed_baseline_1980-2024_all_parameters+thresholds_grasslands_09.07.26.csv", header = TRUE)
t_forst <- read.csv("table_ECEs+features_fixed_baseline_1980-2024_all_parameters+thresholds_forests_09.07.26.csv", header = TRUE)
# str(t_grass) ; str(t_forst)

### First, let's create an event-type variable from the Statistic × Parameter combination,
### and extract the year from Start_Date (the event year)

# Helper to build an event type label
make_event_type <- function(stat, param) { paste(stat, param, sep = "_") }

# Apply to grasslands data
t_grass_base <- t_grass %>%
    mutate(
        EventType = make_event_type(Statistic, Parameter),
        start_date = as.Date(Start_Date),
        end_date = as.Date(End_Date),
        start_year = year(start_date),
        end_year = year(end_date)
    )
# head(t_grass_base); str(t_grass_base); dim(t_grass_base)

# Let's count the number of events per Year × Region × System × EventType.
counts_start <- t_grass_base %>%
    filter(start_year >= 1950, start_year <= 2024) %>%
    group_by(Region, EP, EventType, Year = start_year) %>%
    summarise( N_events = n(), .groups  = "drop")

# Counts for end year (excluding events that start and end in same year)
counts_end <- t_grass_base %>%
    filter(end_year >= 1950, end_year <= 2024, end_year != start_year) %>%
    group_by(Region, EP, EventType, Year = end_year) %>%
    summarise(N_events = n(), .groups  = "drop")

# Combine start-year and end-year counts
df_counts_grass <- bind_rows(counts_start, counts_end) %>%
    group_by(Region, EP, EventType, Year) %>%
    summarise(N_events = sum(N_events), .groups  = "drop")

# Check
# dim(df_counts_grass) # 144 010 Years x Event Type x EP
# head(df_counts_grass)
# str(df_counts_grass)
# summary(df_counts_grass)
# -> Gives us the stacked-bar heights

rm(counts_ends,counts_start,t_grass_base); gc()

### To make the stacked bar chart, we need to aggregate across EPs, summing N_events over
### all sites within each Region × Year × EventType before plotting
df_counts_grass_plot <- df_counts_grass %>%
    group_by(Region, Year, EventType) %>%
    summarise(N_events = sum(N_events, na.rm = TRUE), .groups = "drop")
# head(df_counts_grass_plot); dim(df_counts_grass_plot)
# summary(df_counts_grass_plot)

# Helper FUN to adjust the ECE type labels
pretty_ece_label_generic <- function(x) {
    
    x <- as.character(x)  # ensure character

    sapply(x, function(s) {
        parts <- strsplit(s, "_")[[1]]
        # Capitalize first part
        parts[1] <- paste0(toupper(substr(parts[1], 1, 1)),
            substr(parts[1], 2, nchar(parts[1])))
        # Rejoin with spaces between first two parts, underscores after
        if( length(parts) > 1 ) {
            paste(parts[1], paste(parts[-1], collapse = "_"))
        } else {
            parts[1]
        } # eo if else loop 
    }, USE.NAMES = FALSE)

} # eo FUN - pretty_ece_label_generic


### Define colour palette for Event Types and order accordingly
# unique(df_counts_grass$EventType) # 13, sanity checked

# Define specific colors per EventType for legend clarity
event_colors <- c(
    "max_Ta_10"            = "#a50026",
    "max_Ta_200"           = "#d73027",
    "max_Ts_05"            = "#fdae61",
    "max_Ts_10"            = "#fee090",
    "max_Ts_20"            = "#ffffbf",
    "min_Ta_10"            = "#313695",
    "min_Ta_200"           = "#4575b4",
    "min_Ts_05"            = "#74add1",
    "min_Ts_10"            = "#abd9e9",
    "min_Ts_20"            = "#e0f3f8",
    "max_SM_10"            = "#35978f",
    "min_SM_10"            = "#dfc27d",
    "total_precipitation"  = "#8e99c1"
)

### Adjust order of 'EventType'
# unique(df_counts_grass_plot$EventType)
event_order <- c("max_Ta_10","max_Ta_200","max_Ts_05","max_Ts_10","max_Ts_20",
    "min_Ta_10","min_Ta_200","min_Ts_05","min_Ts_10","min_Ts_20",
    "max_SM_10","min_SM_10","total_precipitation")

df_counts_grass_plot$EventType <- factor(df_counts_grass_plot$EventType, levels = event_order)

ece_labels <- setNames(pretty_ece_label_generic(event_order), event_order)

p_events <- ggplot(df_counts_grass_plot, aes(x = Year, y = N_events, fill = EventType)) +
    geom_col(position = "stack", color = "grey20", linewidth = 0.1) +
    facet_grid(. ~ factor(Region)) +
    scale_fill_manual(values = event_colors, breaks = names(ece_labels),
        labels = ece_labels, name = "ECE type") +
    labs(x = "Year", y = "Number of ECEs detected by the\nfixed baseline approach") +
    scale_y_continuous(limits = c(0,30500)) + # make sure both grasslands and forest follow the same y axis
    theme_bw(base_size = 9) +
    theme(
        strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "grey92", color = NA),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7), panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "right",
        legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold", size = 10)
    ) + ggtitle("Grassland EPs")

# Print 'p_events'
ggsave("Fig.8A.pdf", plot = p_events, device = "pdf", dpi = 300, height = 4, width = 13)


### Same, but for forests this time! -> Fig.8B and then assemble in panel 
# Apply to forests data
t_forst_base <- t_forst %>%
    mutate(
        EventType = make_event_type(Statistic, Parameter),
        start_date = as.Date(Start_Date),
        end_date = as.Date(End_Date),
        start_year = year(start_date),
        end_year = year(end_date)
    )
# head(t_forst_base); str(t_forst_base); dim(t_forst_base)

# Let's count the number of events per Year × Region × System × EventType.
counts_start <- t_forst_base %>%
    filter(start_year >= 1950, start_year <= 2024) %>%
    group_by(Region, EP, EventType, Year = start_year) %>%
    summarise( N_events = n(), .groups  = "drop")

# Counts for end year (excluding events that start and end in same year)
counts_end <- t_forst_base %>%
    filter(end_year >= 1950, end_year <= 2024, end_year != start_year) %>%
    group_by(Region, EP, EventType, Year = end_year) %>%
    summarise(N_events = n(), .groups  = "drop")

# Combine start-year and end-year counts
df_counts_forst <- bind_rows(counts_start, counts_end) %>%
    group_by(Region, EP, EventType, Year) %>%
    summarise(N_events = sum(N_events), .groups  = "drop")

# Check
# dim(df_counts_forst) # 73 046 Years x Event Type x EP
# head(df_counts_forst)
# str(df_counts_forst)
# summary(df_counts_forst)
# -> Gives us the stacked-bar heights

rm(counts_ends,counts_start,t_forst_base); gc()

df_counts_forst_plot <- df_counts_forst %>%
    group_by(Region, Year, EventType) %>%
    summarise(N_events = sum(N_events, na.rm = TRUE), .groups = "drop")

df_counts_forst_plot$EventType <- factor(df_counts_forst_plot$EventType, levels = event_order)

ece_labels <- setNames(pretty_ece_label_generic(event_order), event_order)

p_events_F <- ggplot(df_counts_forst_plot, aes(x = Year, y = N_events, fill = EventType)) +
    geom_col(position = "stack", color = "grey20", linewidth = 0.1) +
    facet_grid(. ~ factor(Region)) +
    scale_fill_manual(values = event_colors, breaks = names(ece_labels),
        labels = ece_labels, name = "ECE type") +
    labs(x = "Year", y = "Number of ECEs detected by the\nfixed baseline approach") +
    scale_y_continuous(limits = c(0,30500)) + # make sure both grasslands and forest follow the same y axis
    theme_bw(base_size = 9) +
    theme(
        strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "grey92", color = NA),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7), panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "right",
        legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold", size = 10)
    ) + ggtitle("Forest EPs")

# Print 'p_events'
ggsave("Fig.8B.pdf", plot = p_events_F, device = "pdf", dpi = 300, height = 4, width = 13)

### Gather in 2x1 panel
require("ggpubr")
Fig.8 <- ggarrange(p_events, p_events_F, align = "hv", ncol = 1, nrow = 2, labels = c("A","B"), common.legend = TRUE, legend = "bottom")
ggsave("Fig.8.pdf", plot = Fig.8, device = "pdf", dpi = 300, height = 8, width = 10)


### Generate 2 tables (and associated .CSV files) to briefly describe the temporal ECE patterns:
# A) Total N ECEs per types and per year (data used to make the plots above)
# B) Mean annual ECE intensity, duration etc. per ECE type

## A) Total N ECEs per types and per year 
# head(df_counts_grass_plot) ; head(df_counts_forst_plot)
df_counts_grass_plot$Ecosystem <- "Grasslands"
df_counts_forst_plot$Ecosystem <- "Forest"
# Bind and save 
df <- rbind(df_counts_grass_plot, df_counts_forst_plot)
# str(df); dim(df)
# unique(df$Ecosystem); summary(factor(df$Ecosystem))

# Save 
write.csv(x = df, file = "table_N_ECEs_all_types_EPs_15.07.26.csv")
write.table(x = df, file = "table_N_ECEs_all_types_EPs_15.07.26.txt")

rm(df); gc()


## B) Mean annual ECE intensity, duration etc. per ECE type
# str(t_grass) ; str(t_forst)

## Summarize mean ± sd Duration & Mean_Intensity per Year and ECE type and Region
df_ECEs_grass <- data.frame(
    t_grass %>%
    mutate(
        Year = lubridate::year(Start_Date),
        Variable = factor(paste(Statistic,Parameter, sep = "_"))
        ) %>%
    group_by(Region,EP,Variable,Year) %>%
    summarise(
        mean_duration = mean(Duration, na.rm = TRUE),
        sd_duration = sd(Duration, na.rm = TRUE),
        mean_intensity = mean(Mean_Intensity, na.rm = TRUE),
        sd_intensity = sd(Mean_Intensity, na.rm = TRUE),
    ) # eo summarise
) # eo ddf

# Check 
dim(df_ECEs_grass)
head(df_ECEs_grass)
str(df_ECEs_grass)

df_ECEs_grass$Ecosystem <- factor("Grasslands")

## Same for forests
df_ECEs_forst <- data.frame(
    t_forst %>%
    mutate(
        Year = lubridate::year(Start_Date),
        Variable = factor(paste(Statistic,Parameter, sep = "_"))
        ) %>%
    group_by(Region,EP,Variable,Year) %>%
    summarise(
        mean_duration = mean(Duration, na.rm = TRUE),
        sd_duration = sd(Duration, na.rm = TRUE),
        mean_intensity = mean(Mean_Intensity, na.rm = TRUE),
        sd_intensity = sd(Mean_Intensity, na.rm = TRUE),
    ) # eo summarise
) # eo ddf

df_ECEs_forst$Ecosystem <- factor("Forests")

### Rbind
df_ECEs <- rbind(df_ECEs_grass, df_ECEs_forst)
df_ECEs <- df_ECEs %>% relocate(Ecosystem)

# head(df_ECEs); tail(df_ECEs)
# str(df_ECEs)
# summary(df_ECEs)

# Save as .csv
write.csv(x = df_ECEs, file = "table_ECEs_features_all_types_EPs_15.07.26.csv")
rm(df_ECEs); gc()

### ------------------------------------------------------------------------------------------------------------

### 16/07/26: Figure 9 - Heatmaps illustrating the changes in SPEI categories in the grasslands and forests of the 3 Exploratories
### -> Panel of 3 x 2 heatmaps

### Use the local SPEI reconstructions to generate the heatmaps. Re-use R Script#7.4_SPEI_compute.R and adapt based on the 
### SPEI categories from the global SPEI database; https://spei.csic.es/spei_database_2_11/#map_name=spei01#map_position=1487)
### SPEI categories to use: 
#   1.65 - 2.33
#   1.28 - 1.65 
#   0.84 - 1.28
#   0 - 0.84
#   -0.84 - 0
#   -1.28 - -0.84
#   -1.65 - -1.28
#   -2.33 - -1.65

### Go to dthe dir where SPEI data are stored
setwd("/home/fbenedetti/BE-EXTREME_files/v1_08.07.26/SPEI_tables")

# Load the 2 SPEI tables
spei_grass <- read.csv("table_SPEI_monthly_grasslands_09.07.26.csv", header = TRUE)
spei_forst <- read.csv("table_SPEI_monthly_forests_09.07.26.csv", header = TRUE)
# str(spei_grass); str(spei_forst)
# summary(spei_grass$SPEI6) ; summary(spei_forst$SPEI6)
# SPEI-6 ranges from -3.20 to 4.10

## Define SPEI bins for a nice looking and adequate (colorblind friendly) colorbar
breaks <- c(-4.1,-2.33,-1.65,-1.28,-0.84,0,0.84,1.28,1.65,2.33,4.1)
# Give it it adequate labels too
labels <- c("< -2.33","-2.33 to -1.65","-1.65 to -1.28","-1.28 to -0.84","-0.84 to 0","0 to 0.84","0.84 to 1.28","1.28 to 1.65","1.65 to 2.33","> 2.33")

# Add bins to 'ddf' subset
spei_grass <- spei_grass %>% mutate(SPEI3_bin = cut(SPEI3, breaks = breaks, labels = labels, include.lowest = T))
spei_grass <- spei_grass %>% mutate(SPEI6_bin = cut(SPEI6, breaks = breaks, labels = labels, include.lowest = T))
spei_grass <- spei_grass %>% mutate(SPEI12_bin = cut(SPEI12, breaks = breaks, labels = labels, include.lowest = T))
# head(spei_grass)

# Re-order levels of the 'region' factor as follows: SCH > HND > SWA
spei_grass$Region <- factor(spei_grass$Region , levels = c("HND","SCH","SWA"))

# use BrBG palette: 10 bins, manually interpolated from the diverging scale
# RColorBrewer "BrBG" only provides up to 11 colours natively, so we need to pick 10
brbg_cols <- rev(brewer.pal(11, "BrBG"))
brbg_10 <- brbg_cols[c(1:5,7:11)]
rm(brbg_cols); gc()

# Create a Date column from Year + Month (use 1st of each month)
spei_grass <- spei_grass %>% mutate(m = as.Date(paste(Year, Month, "01", sep = "-")))
# head(spei_grass)

### Alternative colorpalettes
# colorspace
library("colorspace")
br3_cols <- diverging_hcl(n = 10, palette = "Blue-Red 3") 

## Plot heatmap
heatmap_grass_spei6 <- ggplot(spei_grass, aes(x = m, y = EP, fill = SPEI6_bin)) + geom_tile(width = 31) + 
    scale_fill_manual(name = "SPEI-6", values = rev(brbg_10), labels = labels, na.value = "grey85", drop = FALSE) +
    # scale_fill_manual(name = "SPEI-6", values = rev(br3_cols), labels = labels, na.value = "grey85", drop = FALSE) + 
    scale_x_date(
        breaks = seq(as.Date("1950-01-01"), as.Date("2024-01-01"), by = "10 years"),
        minor_breaks = seq(as.Date("1950-01-01"), as.Date("2024-01-01"), by = "5 years"),
        date_labels  = "%Y", expand = c(0,0)) +
    facet_wrap(~ Region, nrow = 3, scales = "free_y", strip.position = "left") +
    labs(x = NULL, y = "Experimental plots (EPs)", title = NULL) +
    guides(fill = guide_legend( title.position = "top", nrow  = 2, reverse = FALSE)) +
    theme_bw(base_size = 9) +
    theme(
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        strip.text = element_text(face = "bold", size = 8),
        strip.background = element_blank(), strip.placement  = "outside",
        legend.position = "bottom", legend.key.size  = unit(0.35, "cm"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8, face = "bold"),
        panel.grid = element_blank(), panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4)
    ) + ggtitle("Grasslands")

ggsave("Fig.9Abisbis.pdf", plot = heatmap_grass_spei6, device = "pdf", dpi = 300, height = 4, width = 8)


### Same for forests (Fig.9B -> put both A and B to create panel Fig. 9)
# Add bins to 'ddf' subset
spei_forst <- spei_forst %>% mutate(SPEI3_bin = cut(SPEI3, breaks = breaks, labels = labels, include.lowest = T))
spei_forst <- spei_forst %>% mutate(SPEI6_bin = cut(SPEI6, breaks = breaks, labels = labels, include.lowest = T))
spei_forst <- spei_forst %>% mutate(SPEI12_bin = cut(SPEI12, breaks = breaks, labels = labels, include.lowest = T))
# head(spei_forst)

# Re-order levels of the 'region' factor as follows: SCH > HND > SWA
spei_forst$Region <- factor(spei_forst$Region , levels = c("HND","SCH","SWA"))

# Create a Date column from Year + Month (use 1st of each month)
spei_forst <- spei_forst %>% mutate(m = as.Date(paste(Year, Month, "01", sep = "-")))
# head(spei_forst)

## Plot heatmap
heatmap_forst_spei6 <- ggplot(spei_forst, aes(x = m, y = EP, fill = SPEI6_bin)) + geom_tile(width = 31) + 
    scale_fill_manual(name = "SPEI-6", values = rev(brbg_10), labels = labels, na.value = "grey85", drop = FALSE) +
    # scale_fill_manual(name = "SPEI-6", values = rev(br3_cols), labels = labels, na.value = "grey85", drop = FALSE) + 
    scale_x_date(
        breaks = seq(as.Date("1950-01-01"), as.Date("2024-01-01"), by = "10 years"),
        minor_breaks = seq(as.Date("1950-01-01"), as.Date("2024-01-01"), by = "5 years"),
        date_labels  = "%Y", expand = c(0,0)) +
    facet_wrap(~ Region, nrow = 3, scales = "free_y", strip.position = "left") +
    labs(x = NULL, y = "Experimental plots (EPs)", title = NULL) +
    guides(fill = guide_legend( title.position = "top", nrow  = 2, reverse = FALSE)) +
    theme_bw(base_size = 9) +
    theme(
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        strip.text = element_text(face = "bold", size = 8),
        strip.background = element_blank(), strip.placement  = "outside",
        legend.position = "bottom", legend.key.size  = unit(0.35, "cm"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8, face = "bold"),
        panel.grid = element_blank(), panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4)
    ) + ggtitle("Forests")

ggsave("Fig.9B.pdf", plot = heatmap_forst_spei6, device = "pdf", dpi = 300, height = 4, width = 8)

### Combine in panel 
library("ggpubr")
Fig9 <- ggarrange(heatmap_grass_spei6, heatmap_forst_spei6, align = "hv", ncol = 1, nrow = 2, labels = c("A","B"), common.legend = TRUE, legend = "bottom")
ggsave("Fig.9.pdf", plot = Fig9, device = "pdf", dpi = 300, height = 6, width = 8)
# Nice 

### Compute N severe to extreme wet events and N severe to extreme drought events (using 1.65 and -1.65 as thresholds) 
### for the 1950-1980 period compared to the 1994-2024 period, for grasslands and forests separately.  
### -> Summarize main patterns in Sect.3.4 of the ESSD Data Paper.

# An "event" is defined as a run of consecutive months where SPEI-6 exceeds the threshold in the same direction
# -> Use the rle() function 

# Helper FUN 
count_events <- function(df, threshold, direction = c("wet", "dry")) {
    
    direction <- match.arg(direction)
  
    df <- df %>% arrange(EP, m)
    
    df %>% group_by(EP, Region) %>% 
        mutate(exceeds = if (direction == "wet") SPEI6 > 1.65 else SPEI6 < -1.65) %>%
        summarise(
            n_events = {
                x <- as.logical(na.omit(exceeds))
                if( length(x) == 0 ) {
                    NA_integer_
                } else {
                    r <- rle(x)
                    sum(r$values == TRUE)
                } # eo if else loop 
            },
            
            mean_duration = {
                x <- as.logical(na.omit(exceeds))
                if( length(x) == 0 ) {
                    NA_real_
                } else {
                    r <- rle(x)
                    # lengths of TRUE runs only (i.e. actual events)
                    event_lengths <- r$lengths[r$values == TRUE]
                    
                    if( length(event_lengths) == 0 ) NA_real_ else mean(event_lengths)
                } # eo if else loop
            },
            .groups = "drop"

        ) # eo summarise

} # eo FUN - count_events

### Apply to grasslands and forests

wet_grass_early <- spei_grass %>% filter(Year >= 1950, Year <= 1980) %>% 
    count_events(threshold = 1.65, direction = "wet") %>% 
    mutate(Ecosystem = "Grasslands", Period = "1950-1980", Event_type = "Wet")

wet_grass_late <- spei_grass %>% filter(Year >= 1994, Year <= 2024) %>% 
    count_events(threshold = 1.65, direction = "wet") %>% 
    mutate(Ecosystem = "Grasslands", Period = "1994-2024", Event_type = "Wet")

dry_grass_early <- spei_grass %>% filter(Year >= 1950, Year <= 1980) %>% 
    count_events(threshold = -1.65, direction = "dry") %>% 
    mutate(Ecosystem = "Grasslands", Period = "1950-1980", Event_type = "Dry")

dry_grass_late <- spei_grass %>% filter(Year >= 1994, Year <= 2024) %>% 
    count_events(threshold = -1.65, direction = "dry") %>% 
    mutate(Ecosystem = "Grasslands", Period = "1994-2024", Event_type = "Dry")
# head(wet_grass_early) ; head(wet_grass_late)
# head(dry_grass_early) ; head(dry_grass_late)
summary(dry_grass_early) ; summary(dry_grass_late)

# Repeat for forests (substitute landuse == "forest")
wet_forst_early <- spei_forst %>% filter(Year >= 1950, Year <= 1980) %>% 
    count_events(threshold = 1.65, direction = "wet") %>% 
    mutate(Ecosystem = "Forests", Period = "1950-1980", Event_type = "Wet")

wet_forst_late <- spei_forst %>% filter(Year >= 1994, Year <= 2024) %>% 
    count_events(threshold = 1.65, direction = "wet") %>%
    mutate(Ecosystem = "Forests", Period = "1994-2024", Event_type = "Wet")

dry_forst_early <- spei_forst %>% filter(Year >= 1950, Year <= 1980) %>% 
    count_events(threshold = -1.65, direction = "dry") %>%
    mutate(Ecosystem = "Forests", Period = "1950-1980", Event_type = "Dry")

dry_forst_late <- spei_forst %>% filter(Year >= 1994, Year <= 2024) %>% 
    count_events(threshold = -1.65, direction = "dry") %>% 
    mutate(Ecosystem = "Forests", Period = "1994-2024", Event_type = "Dry")

summary(dry_forst_early) ; summary(dry_forst_late)

### Rbind 
table <- bind_rows(
    wet_grass_early,
    wet_grass_late,
    dry_grass_early,
    dry_grass_late,
    wet_forst_early,
    wet_forst_late,
    dry_forst_early,
    dry_forst_late
)

# Check
dim(table)
head(table); tail(table)
str(table)
summary(table)

# Save as .csv
write.csv(x = table, file = "table_SPEI6_Nevents_duration_all_EPs_16.07.26.csv")
rm(table); gc()


### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------

### 02/07/26: Suppl. Figures

### Fig. S1: Distribution of daily offsets per: Tree Type x Month, Region for the following chosen parameters: 
### min/max Ta_200, min/max Ts_10 and min/max SM_10

### Recycle part of the main FUN of R Script#4.4.3 to plot offsets distributions in a for loop 

# For testing master FUN below while you write it
# var = "Ta_200"
# stat = "max"
# age = 75
# years = 15


plot_offsets_distribution <- function(var, stat, age, years = 15) {

        #' This function takes three arguments and returns a model object of class 'gamm':
        #' @param var The variable to model (character): "Ta_200","Ta_10","Ts_05","Ts_10","Ts_20","SM_10" 
        #' @param stat The daily stat of the associated variable (character): 'max' or 'min'
        #' @param age Maximum stand age to account for in the GAMM approach (integer) 
        #' @param years Number of years of daily data to train the GAMM (integer): 6, 8, 10, 12, 14 or 16
        #' 'age' Should be between 50 and 70 based on R Script#4.4.2
        #' @return A formatted data.frame combining the inputs.
      
        # Useless message
        message(paste("Loading the aggregated data for ",stat," ",var, sep = ""))

        # Read in the data after identifying the corresponding
        setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data") #; dir()

        file <- dir()[grepl(paste("metadata",stat,var,sep = "_"),dir())]
        df <- readRDS(file) # dim(df); str(df) # should be 858'600 daily measurements

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

        # Remove rows with NaN
        combin_df <- combin_df %>% drop_na(StandAge)
        combin_df <- combin_df %>% drop_na(Offset)

        # Create simple vector of years for below
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

        # Define ylabs based on parameters at hand 
        if( var == "Ta_200" ) {
            if( stat == "max" ) {
                lab <- "Offsets in daily max. air temperature\n(max Ta_200) [°C]"
            } else {
                lab <- "Offsets in daily min. air temperature\n(min Ta_200) [°C]"
            } # eo 2nd else if loop
        } else if ( var == "Ts_10" ) {
            if( stat == "max" ) {
                lab <- "Offsets in daily max. soil temperature\n(max Ts_10) [°C]"
            } else {
                lab <- "Offsets in daily min. soil temperature\n(min Ts_10) [°C]"
            } # eo 2nd else if loop
        } else if( var == "SM_10" ) {
            if( stat == "max" ) {
                lab <- "Offsets in daily max. soil moisture\n(max SM_10) [%]"
            } else {
                lab <- "Offsets in daily min. soil moisture\n(min SM_10) [%]"
            } # eo 2nd else if loop
        } # eo 1st else if loop

        # Make plots to illustrate distribution of offsets per: Tree Type x Month, Region
        p1 <- ggplot(data = combin_df, aes(x = factor(Region), y = Offset)) + geom_violin(colour = "black", fill = "grey75") + 
            geom_boxplot(fill = "white", colour = "black", width = .2) + geom_hline(yintercept = 0, linetype = "dashed") + 
            xlab("Region") + ylab(lab) + theme_bw()
        
        p2 <- ggplot(data = combin_df, aes(x = factor(Month), y = Offset, fill = factor(TreeType))) + 
            geom_violin(colour = "black") + geom_boxplot(fill = "white", colour = "black", width = .2) + 
            scale_fill_manual(name = "Tree type", values = c("#238443","#d9f0a3")) + geom_hline(yintercept = 0, linetype = "dashed") + 
            xlab("Months") + ylab(lab) + theme_bw() + theme(legend.position = 'top')

        panel <- ggarrange(p1, p2, align = "hv", ncol = 2, nrow = 1, widths = c(1,2))
        ggsave(plot = panel, filename = paste("panel_distribution_offsets_",stat,"_",var,"_regions_months_02.07.26.jpg", sep = ""), dpi = 300, height = 4.75, width = 9.5)

        rm(panel,p1,p2)
        gc()

} # eo FUN - model_offset_gamm


### Apply plot_offsets_distribution() to chosen cliamtic parameters
for(v in c("Ta_200","Ts_10","SM_10")) {
    for(s in c("max","min")) {
        plot_offsets_distribution(var = v, stat = s, age = 75, years = 15)
    } # eo 2nd for loop
} # eo 1st for loop


### ------------------------------------------------------------------------------------------------------------

### Fig. S2: Panel of GAM's smooth plots with 'ggratia'
### Recycle parts of R Script#4.4.3 to plot the GAMMs' smooth terms with "gratia"
library("gratia")
library("patchwork")

### GAMM objects are stored here: 
setwd("/home/fbenedetti/Exploratories/WP3_Instrumentation_Data/daily_offsets_for_microlimate_modelling/aggregated_data/models") # dir()

# List the GAMM objects of interest (those from early June 2026)
models <- dir()[grep("03.02.26.rds",dir())] # length(models) # should be 12

### Create a function that take the model object based on 'stat' and 'var' and plots the smooths learnt by the GAMM
# stat <- "min"
# var <- "Ta_10"

plot_gamm_smooth <- function(stat, var) {

        #' This function takes three arguments and returns a model object of class 'gamm':
        #' @param var The variable to model (character): "Ta_200","Ta_10","Ts_05","Ts_10","Ts_20","SM_10" 
        #' @param stat The daily stat of the associated variable (character): 'max' or 'min'
        #' @return Smooth plots made with the "gratia" package
      
        # Useless message
        message(paste("Loading the GAMM object for ",stat," ",var, sep = ""))

        model_file <- models[grepl(paste("offsets",stat,var,sep = "_"), models)]
        mod <- readRDS(model_file) 
        rm(model_file)

        # Plots
        pdf(NULL) 
        smooths.plot <- draw(mod$gam) # class(smooths.plot)

        # as.list() gives the underlying list of ggplots
        plots_list <- as.list(smooths.plot) # str(plots_list)
        # To return the names of the smooths (to be adjusted for the paper)
        #   lapply(plots_list, function(p) p$labels$title)
        # Replace specific titles by index
        plots_list[[1]]$labels$title <- "s(Capped stand age)"
        plots_list[[4]]$labels$title <- "s(Anomalies to DOY)"
        # Then rebuild the patchwork
        smooths.plot2 <- wrap_plots(plots_list) # class(smooths.plot2)
        # Combines a list of ggplots back into a patchwork object

        # Add a title
        title <- paste("Modelled smooths for ",stat," ",var, sep = "")
        smooths.plot2 <- smooths.plot2 + plot_annotation(title = title)

        # BONUS: GAMM residuals plots (keep for later maybe)
        resids.plot <- appraise(mod$gam, line_col = "#238443", point_alpha = .01)

        filename1 <- paste("plot_smooths_GAMM",stat,var,"02.07.26.jpg", sep = "_")
        filename2 <- paste("plot_residuals_GAMM",stat,var,"02.07.26.jpg", sep = "_")
        ggsave(plot = smooths.plot2, filename = filename1, dpi = 300, height = 5.75, width = 5.75)
        ggsave(plot = resids.plot, filename = filename2, dpi = 300, height = 7, width = 7.5)

        rm(mod,smooths.plot,resids.plot,smooths.plot2); gc()
        
} # eo FUN - plot_gamm_smooth


for(v in c("Ta_10","Ta_200","Ts_05","Ts_10","Ts_20","SM_10")) {
    for(s in c("max","min")) {
        plot_gamm_smooth(var = v, stat = s)
    } # eo 2nd for loop
} # eo 1st for loop


### ------------------------------------------------------------------------------------------------------------

### 16/07/26: Fig. S3: Long-term dynamics of ECEs for both grasslands and forests - shifting baseline approach
### Same as Fig. 8 but based on the shifting baseline approach - re-use the R code from above

### Got to ECEs tables dir
setwd("/home/fbenedetti/BE-EXTREME_files/v1_08.07.26/ECE_tables")

# Load the tables
t_grass <- read.csv("table_ECEs+features_shifting_baseline_1980-2024_all_parameters+thresholds_grasslands_09.07.26.csv", header = TRUE)
t_forst <- read.csv("table_ECEs+features_shifting_baseline_1980-2024_all_parameters+thresholds_forests_09.07.26.csv", header = TRUE)
# str(t_grass) ; str(t_forst)

### First, let's create an event-type variable from the Statistic × Parameter combination,
### and extract the year from Start_Date (the event year)

# Helper to build an event type label
make_event_type <- function(stat, param) { paste(stat, param, sep = "_") }

# Apply to grasslands data
t_grass_base <- t_grass %>%
    mutate(
        EventType = make_event_type(Statistic, Parameter),
        start_date = as.Date(Start_Date),
        end_date = as.Date(End_Date),
        start_year = year(start_date),
        end_year = year(end_date)
    )
# head(t_grass_base); dim(t_grass_base)
# summary(t_grass_base)

# Let's count the number of events per Year × Region × System × EventType.
counts_start <- t_grass_base %>%
    filter(start_year >= 1980, start_year <= 2024) %>%
    group_by(Region, EP, EventType, Year = start_year) %>%
    summarise( N_events = n(), .groups  = "drop")

# Counts for end year (excluding events that start and end in same year)
counts_end <- t_grass_base %>%
    filter(end_year >= 1980, end_year <= 2024, end_year != start_year) %>%
    group_by(Region, EP, EventType, Year = end_year) %>%
    summarise(N_events = n(), .groups  = "drop")

# Combine start-year and end-year counts
df_counts_grass <- bind_rows(counts_start, counts_end) %>%
    group_by(Region, EP, EventType, Year) %>%
    summarise(N_events = sum(N_events), .groups  = "drop")

# Check
# dim(df_counts_grass) # 86'364 Years x Event Type x EP
# head(df_counts_grass)
# str(df_counts_grass)
# summary(df_counts_grass)
# -> Gives us the stacked-bar heights

rm(counts_ends,counts_start,t_grass_base); gc()

### To make the stacked bar chart, we need to aggregate across EPs, summing N_events over
### all sites within each Region × Year × EventType before plotting
df_counts_grass_plot <- df_counts_grass %>%
    group_by(Region, Year, EventType) %>%
    summarise(N_events = sum(N_events, na.rm = TRUE), .groups = "drop")

# head(df_counts_grass_plot); dim(df_counts_grass_plot)
# summary(df_counts_grass_plot)

# Helper FUN to adjust the ECE type labels
pretty_ece_label_generic <- function(x) {
    
    x <- as.character(x)  # ensure character

    sapply(x, function(s) {
        parts <- strsplit(s, "_")[[1]]
        # Capitalize first part
        parts[1] <- paste0(toupper(substr(parts[1], 1, 1)),
            substr(parts[1], 2, nchar(parts[1])))
        # Rejoin with spaces between first two parts, underscores after
        if( length(parts) > 1 ) {
            paste(parts[1], paste(parts[-1], collapse = "_"))
        } else {
            parts[1]
        } # eo if else loop 
    }, USE.NAMES = FALSE)

} # eo FUN - pretty_ece_label_generic


### Define colour palette for Event Types and order accordingly
# unique(df_counts_grass$EventType) # 13, sanity checked

# Define specific colors per EventType for legend clarity
event_colors <- c(
    "max_Ta_10"            = "#a50026",
    "max_Ta_200"           = "#d73027",
    "max_Ts_05"            = "#fdae61",
    "max_Ts_10"            = "#fee090",
    "max_Ts_20"            = "#ffffbf",
    "min_Ta_10"            = "#313695",
    "min_Ta_200"           = "#4575b4",
    "min_Ts_05"            = "#74add1",
    "min_Ts_10"            = "#abd9e9",
    "min_Ts_20"            = "#e0f3f8",
    "max_SM_10"            = "#35978f",
    "min_SM_10"            = "#dfc27d",
    "total_precipitation"  = "#8e99c1"
)

### Adjust order of 'EventType'
# unique(df_counts_grass_plot$EventType)
event_order <- c("max_Ta_10","max_Ta_200","max_Ts_05","max_Ts_10","max_Ts_20",
    "min_Ta_10","min_Ta_200","min_Ts_05","min_Ts_10","min_Ts_20",
    "max_SM_10","min_SM_10","total_precipitation")

df_counts_grass_plot$EventType <- factor(df_counts_grass_plot$EventType, levels = event_order)

ece_labels <- setNames(pretty_ece_label_generic(event_order), event_order)

p_events_G <- ggplot(df_counts_grass_plot, aes(x = Year, y = N_events, fill = EventType)) +
    geom_col(position = "stack", color = "grey20", linewidth = 0.1) +
    facet_grid(. ~ factor(Region)) +
    scale_fill_manual(values = event_colors, breaks = names(ece_labels),
        labels = ece_labels, name = "ECE type") +
    labs(x = "Year", y = "Number of ECEs detected by the\nshifting baseline approach") +
    scale_y_continuous(limits = c(0,30500)) + # make sure both grasslands and forest follow the same y axis
    theme_bw(base_size = 9) +
    theme(
        strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "grey92", color = NA),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7), panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "right",
        legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold", size = 10)
    ) + ggtitle("Grassland EPs")

# Print 'p_events'
ggsave("Fig.S3A.pdf", plot = p_events_G, device = "pdf", dpi = 300, height = 4, width = 13)


### Same, but for forests
t_forst_base <- t_forst %>%
    mutate(
        EventType = make_event_type(Statistic, Parameter),
        start_date = as.Date(Start_Date),
        end_date = as.Date(End_Date),
        start_year = year(start_date),
        end_year = year(end_date)
    )
# head(t_forst_base); dim(t_forst_base)
# summary(t_forst_base)

# Let's count the number of events per Year × Region × System × EventType.
counts_start <- t_forst_base %>%
    filter(start_year >= 1980, start_year <= 2024) %>%
    group_by(Region, EP, EventType, Year = start_year) %>%
    summarise( N_events = n(), .groups  = "drop")

# Counts for end year (excluding events that start and end in same year)
counts_end <- t_forst_base %>%
    filter(end_year >= 1980, end_year <= 2024, end_year != start_year) %>%
    group_by(Region, EP, EventType, Year = end_year) %>%
    summarise(N_events = n(), .groups  = "drop")

# Combine start-year and end-year counts
df_counts_forst <- bind_rows(counts_start, counts_end) %>%
    group_by(Region, EP, EventType, Year) %>%
    summarise(N_events = sum(N_events), .groups  = "drop")

# Check
# dim(df_counts_forst) # 49'946 Years x Event Type x EP
# head(df_counts_forst)
# str(df_counts_forst)
# summary(df_counts_forst)
# -> Gives us the stacked-bar heights

rm(counts_ends,counts_start,t_forst_base); gc()

df_counts_forst_plot <- df_counts_forst %>%
    group_by(Region, Year, EventType) %>%
    summarise(N_events = sum(N_events, na.rm = TRUE), .groups = "drop")

df_counts_forst_plot$EventType <- factor(df_counts_forst_plot$EventType, levels = event_order)

ece_labels <- setNames(pretty_ece_label_generic(event_order), event_order)

p_events_F <- ggplot(df_counts_forst_plot, aes(x = Year, y = N_events, fill = EventType)) +
    geom_col(position = "stack", color = "grey20", linewidth = 0.1) +
    facet_grid(. ~ factor(Region)) +
    scale_fill_manual(values = event_colors, breaks = names(ece_labels),
        labels = ece_labels, name = "ECE type") +
    labs(x = "Year", y = "Number of ECEs detected by the\nfixed baseline approach") +
    scale_y_continuous(limits = c(0,30500)) + # make sure both grasslands and forest follow the same y axis
    theme_bw(base_size = 9) +
    theme(
        strip.text = element_text(size = 8, face = "bold"),
        strip.background = element_rect(fill = "grey92", color = NA),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7), panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(), legend.position = "right",
        legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        plot.title = element_text(face = "bold", size = 10)
    ) + ggtitle("Forest EPs")

### Gather in 2x1 panel
require("ggpubr")
FigS3 <- ggarrange(p_events_G, p_events_F, align = "hv", ncol = 1, nrow = 2, labels = c("A","B"), common.legend = TRUE, legend = "bottom")
ggsave("Fig.S3.pdf", plot = FigS3, device = "pdf", dpi = 300, height = 8, width = 10)
# Gut gut

### Save N ECEs as .CSV file to briefly sumarise numbers in the Sect. 3.3 of the ESSD data paper 
# head(df_counts_grass_plot) ; head(df_counts_forst_plot)
df_counts_grass_plot$Ecosystem <- "Grasslands"
df_counts_forst_plot$Ecosystem <- "Forest"
# Bind and save 
df <- rbind(df_counts_grass_plot, df_counts_forst_plot)
# str(df); dim(df)
# summary(factor(df$Ecosystem))

# Save 
write.csv(x = df, file = "table_N_ECEs_all_types_EPs_shifting_baseline_16.07.26.csv")
write.table(x = df, file = "table_N_ECEs_all_types_EPs_shifting_baseline_16.07.26.txt")

rm(df); gc()

### ------------------------------------------------------------------------------------------------------------

### 20/07/26: Fig. S4: Heatmaps for SPEI-3 and SPEI-12 (SPEI-6 was shown in main Fig. 9)

### Go to dthe dir where SPEI data are stored
setwd("/home/fbenedetti/BE-EXTREME_files/v1_08.07.26/SPEI_tables")

## Load the 2 SPEI tables
spei_grass <- read.csv("table_SPEI_monthly_grasslands_09.07.26.csv", header = TRUE)
spei_forst <- read.csv("table_SPEI_monthly_forests_09.07.26.csv", header = TRUE)
# str(spei_grass); str(spei_forst)
# summary(spei_grass$SPEI3) ; summary(spei_forst$SPEI3) # scales between -3 and +3
# summary(spei_grass$SPEI12) ; summary(spei_forst$SPEI12) # scales between -3 and +3.3

## Define SPEI bins for a nice looking and adequate (colorblind friendly) colorbar
breaks <- c(-4.1,-2.33,-1.65,-1.28,-0.84,0,0.84,1.28,1.65,2.33,4.1)
breaks_spei3 <- c(-3,-2.33,-1.65,-1.28,-0.84,0,0.84,1.28,1.65,2.33,3)
breaks_spei12 <- c(-3,-2.33,-1.65,-1.28,-0.84,0,0.84,1.28,1.65,2.33,3.3)
labels <- c("< -2.33","-2.33 to -1.65","-1.65 to -1.28","-1.28 to -0.84","-0.84 to 0","0 to 0.84","0.84 to 1.28","1.28 to 1.65","1.65 to 2.33","> 2.33")

# Add bins to 'spei_grass' and 'spei_forst' data.frames
spei_grass <- spei_grass %>% mutate(SPEI3_bin = cut(SPEI3, breaks = breaks_spei3, labels = labels, include.lowest = T))
spei_grass <- spei_grass %>% mutate(SPEI6_bin = cut(SPEI6, breaks = breaks, labels = labels, include.lowest = T))
spei_grass <- spei_grass %>% mutate(SPEI12_bin = cut(SPEI12, breaks = breaks_spei12, labels = labels, include.lowest = T))
# head(spei_grass)

spei_forst <- spei_forst %>% mutate(SPEI3_bin = cut(SPEI3, breaks = breaks_spei3, labels = labels, include.lowest = T))
spei_forst <- spei_forst %>% mutate(SPEI6_bin = cut(SPEI6, breaks = breaks, labels = labels, include.lowest = T))
spei_forst <- spei_forst %>% mutate(SPEI12_bin = cut(SPEI12, breaks = breaks_spei12, labels = labels, include.lowest = T))
# head(spei_forst)

# Re-order levels of the 'region' factor as follows: SCH > HND > SWA
spei_grass$Region <- factor(spei_grass$Region , levels = c("HND","SCH","SWA"))
spei_forst$Region <- factor(spei_forst$Region , levels = c("HND","SCH","SWA"))

# use BrBG palette: 10 bins, manually interpolated from the diverging scale
# RColorBrewer "BrBG" only provides up to 11 colours natively, so we need to pick 10
brbg_cols <- rev(brewer.pal(11, "BrBG"))
brbg_10 <- brbg_cols[c(1:5,7:11)]
rm(brbg_cols); gc()

# Create a Date column from Year + Month (use 1st of each month)
spei_grass <- spei_grass %>% mutate(m = as.Date(paste(Year, Month, "01", sep = "-")))
spei_forst <- spei_forst %>% mutate(m = as.Date(paste(Year, Month, "01", sep = "-")))
# head(spei_grass) ; head(spei_forst)

### Alternative colorpalettes
# colorspace
library("colorspace")
br3_cols <- diverging_hcl(n = 10, palette = "Blue-Red 3") 

## Plot SPEI-3 heatmaps
# Grasslands
heatmap_grass_spei3 <- ggplot(spei_grass, aes(x = m, y = EP, fill = SPEI3_bin)) + geom_tile(width = 31) + 
    scale_fill_manual(name = "SPEI", values = rev(brbg_10), labels = labels, na.value = "grey85", drop = FALSE) +
    scale_x_date(
        breaks = seq(as.Date("1950-01-01"), as.Date("2024-01-01"), by = "10 years"),
        minor_breaks = seq(as.Date("1950-01-01"), as.Date("2024-01-01"), by = "5 years"),
        date_labels  = "%Y", expand = c(0,0)) +
    facet_wrap(~ Region, nrow = 3, scales = "free_y", strip.position = "left") +
    labs(x = NULL, y = "Experimental plots (EPs)", title = NULL) +
    guides(fill = guide_legend( title.position = "top", nrow  = 2, reverse = FALSE)) +
    theme_bw(base_size = 9) +
    theme(
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        strip.text = element_text(face = "bold", size = 8),
        strip.background = element_blank(), strip.placement  = "outside",
        legend.position = "bottom", legend.key.size  = unit(0.35, "cm"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8, face = "bold"),
        panel.grid = element_blank(), panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4)
    ) + ggtitle("Grasslands")

# Forests
heatmap_forst_spei3 <- ggplot(spei_forst, aes(x = m, y = EP, fill = SPEI3_bin)) + geom_tile(width = 31) + 
    scale_fill_manual(name = "SPEI", values = rev(brbg_10), labels = labels, na.value = "grey85", drop = FALSE) +
    scale_x_date(
        breaks = seq(as.Date("1950-01-01"), as.Date("2024-01-01"), by = "10 years"),
        minor_breaks = seq(as.Date("1950-01-01"), as.Date("2024-01-01"), by = "5 years"),
        date_labels  = "%Y", expand = c(0,0)) +
    facet_wrap(~ Region, nrow = 3, scales = "free_y", strip.position = "left") +
    labs(x = NULL, y = "Experimental plots (EPs)", title = NULL) +
    guides(fill = guide_legend( title.position = "top", nrow  = 2, reverse = FALSE)) +
    theme_bw(base_size = 9) +
    theme(
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        strip.text = element_text(face = "bold", size = 8),
        strip.background = element_blank(), strip.placement  = "outside",
        legend.position = "bottom", legend.key.size  = unit(0.35, "cm"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8, face = "bold"),
        panel.grid = element_blank(), panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4)
    ) + ggtitle("Forests")


## Plot SPEI-12 heatmaps
# Grasslands
heatmap_grass_spei12 <- ggplot(spei_grass, aes(x = m, y = EP, fill = SPEI12_bin)) + geom_tile(width = 31) + 
    scale_fill_manual(name = "SPEI", values = rev(brbg_10), labels = labels, na.value = "grey85", drop = FALSE) +
    scale_x_date(
        breaks = seq(as.Date("1950-01-01"), as.Date("2024-01-01"), by = "10 years"),
        minor_breaks = seq(as.Date("1950-01-01"), as.Date("2024-01-01"), by = "5 years"),
        date_labels  = "%Y", expand = c(0,0)) +
    facet_wrap(~ Region, nrow = 3, scales = "free_y", strip.position = "left") +
    labs(x = NULL, y = "Experimental plots (EPs)", title = NULL) +
    guides(fill = guide_legend( title.position = "top", nrow  = 2, reverse = FALSE)) +
    theme_bw(base_size = 9) +
    theme(
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        strip.text = element_text(face = "bold", size = 8),
        strip.background = element_blank(), strip.placement  = "outside",
        legend.position = "bottom", legend.key.size  = unit(0.35, "cm"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8, face = "bold"),
        panel.grid = element_blank(), panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4)
    ) + ggtitle("Grasslands")

# Forests
heatmap_forst_spei12 <- ggplot(spei_forst, aes(x = m, y = EP, fill = SPEI12_bin)) + geom_tile(width = 31) + 
    scale_fill_manual(name = "SPEI", values = rev(brbg_10), labels = labels, na.value = "grey85", drop = FALSE) +
    scale_x_date(
        breaks = seq(as.Date("1950-01-01"), as.Date("2024-01-01"), by = "10 years"),
        minor_breaks = seq(as.Date("1950-01-01"), as.Date("2024-01-01"), by = "5 years"),
        date_labels  = "%Y", expand = c(0,0)) +
    facet_wrap(~ Region, nrow = 3, scales = "free_y", strip.position = "left") +
    labs(x = NULL, y = "Experimental plots (EPs)", title = NULL) +
    guides(fill = guide_legend( title.position = "top", nrow  = 2, reverse = FALSE)) +
    theme_bw(base_size = 9) +
    theme(
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        strip.text = element_text(face = "bold", size = 8),
        strip.background = element_blank(), strip.placement  = "outside",
        legend.position = "bottom", legend.key.size  = unit(0.35, "cm"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8, face = "bold"),
        panel.grid = element_blank(), panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4)
    ) + ggtitle("Forests")


### Combine all 4 heatmaps in Fig. S4 panel 
library("ggpubr")
FigS4A <- ggarrange(heatmap_grass_spei3, heatmap_forst_spei3, heatmap_grass_spei12, heatmap_forst_spei12,
        align = "hv", ncol = 2, nrow = 2, labels = c("A","B","C","D"), common.legend = TRUE, legend = "bottom")
ggsave("FigS4.pdf", plot = FigS4A, device = "pdf", dpi = 300, height = 6.5, width = 10)
# Nice 


### ------------------------------------------------------------------------------------------------------------

### 22-23/07/2026: Fig.S5: Evaluating our monthly SPEI reconstructions by comparing them to the SPEI-3/6/12 values 
### from the global SPEI database (form of technical validation). 

### Rationale: Extract the monthly SPEI-3/6/12 values from the global SPEI database at EP location (use the precise EP coordinates
### or the coordinates of the region) and compare them against our reconstructed values. 
### Gather global and local SPEI values and evaluate similarity with simple corr. coefficient and a bivariate scatterplot

library("terra")
library("sf")

# The EPs coordinates are in the 1000_10 dataset, here
setwd("/home/fbenedetti/Exploratories/EP metadata/1000_10_Dataset")
coords <- read.csv("1000_10_data.csv", header = TRUE, sep = ",")
# str(coords); summary(coords)
# unique(coords$EP_Plot_ID)

# -> Need to adjust EP ids so they match those in our dataset
BEplotZeros <- function(dat, column, plotnam = "EP_Plot_ID") {
      dat <- as.data.frame(dat)
      funz <- function(x) ifelse((nchar(as.character(x))==4), gsub("(.)$", "0\\1", x), as.character(x)) # eo funz
      dat[,plotnam] <- sapply(dat[,column],funz)
      return(dat)
} # eo FUN
# Apply
coords <- BEplotZeros(coords, "EP_Plot_ID", plotnam = "EP")
# str(coords); head(coords)
# unique(coords$EP)

# Drop 'Plot_ID'
coords <- coords %>% select(-c(Plot_ID,EP_Plot_ID))
# Relocate
coords <- coords %>% relocate(EP)
# Subset 
coords <- coords[,c("EP","Longitude","Latitude")]
# Keep only rows with non "na" values in EP 
coords <- coords[coords$EP != "na",]
# dim(coords) # should be 300x3

# If you want the result as points first, make a point object and extract from that.
sites <- vect(coords[,c("Longitude","Latitude")], geom = c("Longitude", "Latitude"), crs = "EPSG:4326")

# The .nc files from the global SPEI database are located here: 
setwd("/home/fbenedetti/SPEI_global_database/raw_netCDF_13.07.26")
# To open the .nc of interest
r_spei3 <- rast("spei03.nc")
r_spei6 <- rast("spei06.nc")
r_spei12 <- rast("spei12.nc")
# r_spei3

# Check CRS
# crs(r_spei3) # crs(r_spei6)
# Looks OK
# Check time dimension
# time(r_spei3)

# Turn the time vector into layer names
names(r_spei3) <- format(time(r_spei3), "%Y-%m-%d")
names(r_spei6) <- format(time(r_spei6), "%Y-%m-%d")
names(r_spei12) <- format(time(r_spei12), "%Y-%m-%d")

# Subset the raster to the first date of interest ("1950-01-01"), then extract from that single layer. S
dates <- which(names(r_spei12) >= "1950-01-16")
r_spei3 <- r_spei3[[dates]]
r_spei6 <- r_spei6[[dates]]
r_spei12 <- r_spei12[[dates]]

# To extract at specififc locations: extract()
# ?terra::extract
table_spei6 <- terra::extract(x = r_spei6, y = sites)
table_spei3 <- terra::extract(x = r_spei3, y = sites)
table_spei12 <- terra::extract(x = r_spei12, y = sites)
# str(table_spei6); head(table_spei6)

# Add EP to these tables, drop 'ID' col from 'table_spei6'
table_spei6 <- cbind(coords,table_spei6[,c(2:length(table_spei6))])
table_spei3 <- cbind(coords,table_spei3[,c(2:length(table_spei3))])
table_spei12 <- cbind(coords,table_spei12[,c(2:length(table_spei12))])

# Melt to have values in a single column and dates and EPs as separate columns
m_table_spei6 <- melt(table_spei6, value.name = "SPEI-6",, id.vars = c(1:3))
m_table_spei3 <- melt(table_spei3, value.name = "SPEI-3", id.vars = c(1:3))
m_table_spei12 <- melt(table_spei12, value.name = "SPEI-12", id.vars = c(1:3))
# head(m_table_spei6); str(m_table_spei6); summary(m_table_spei6)

# Convert "variable" to date and extract Month and Year to match our local SPEI data 
colnames(m_table_spei6)[4] <- "Date"
colnames(m_table_spei3)[4] <- "Date"
colnames(m_table_spei12)[4] <- "Date"
# unique(m_table_spei3$Date)
# head(as.Date(m_table_spei6$Date, format = "%Y-%m-%d")) # seems to work
m_table_spei6$Date <- as.Date(m_table_spei6$Date, format = "%Y-%m-%d")
m_table_spei3$Date <- as.Date(m_table_spei3$Date, format = "%Y-%m-%d")
m_table_spei12$Date <- as.Date(m_table_spei12$Date, format = "%Y-%m-%d")
# Extract 'Month' & 'Year'
m_table_spei6$Month <- lubridate::month(m_table_spei6$Date) ; m_table_spei6$Year <- lubridate::year(m_table_spei6$Date)
m_table_spei3$Month <- lubridate::month(m_table_spei3$Date) ; m_table_spei3$Year <- lubridate::year(m_table_spei3$Date)
m_table_spei12$Month <- lubridate::month(m_table_spei12$Date) ; m_table_spei12$Year <- lubridate::year(m_table_spei12$Date)
# Relocate
m_table_spei6 <- m_table_spei6 %>% relocate("SPEI-6", .after = Year)
m_table_spei3 <- m_table_spei3 %>% relocate("SPEI-3", .after = Year)
m_table_spei12 <- m_table_spei12 %>% relocate("SPEI-12", .after = Year)

# Our local SPEI reconstructions are located here: 
setwd("/home/fbenedetti/BE-EXTREME_files/v1_08.07.26/SPEI_tables")
speiG <- read.csv("table_SPEI_monthly_grasslands_09.07.26.csv", header = TRUE)
speiF <- read.csv("table_SPEI_monthly_forests_09.07.26.csv", header = TRUE)
spei_be <- rbind(speiG,speiF)
rm(speiG,speiF); gc()
# str(spei_be)

# Join local SPEI values and global SPEI values
spei_be <- spei_be %>% left_join(m_table_spei6[,c("EP","Month","Year","SPEI-6")], by = c("EP","Month","Year"))
spei_be <- spei_be %>% left_join(m_table_spei3[,c("EP","Month","Year","SPEI-3")], by = c("EP","Month","Year"))
spei_be <- spei_be %>% left_join(m_table_spei12[,c("EP","Month","Year","SPEI-12")], by = c("EP","Month","Year"))
# Check 
# head(spei_be); str(spei_be)
# summary(spei_be)

### BEWARE: Many EPs fall into the same glibal SPEI grid cell because, for most of Germany,
### the 0.5° × 0.5° resoltuion of the global SPEI database is approximately 56 km × 32–38 km
### As a result, many separate EPs share the same global SPEI values, which will bias our comparison against local SPEI
### -> Perform comparison on the regional monthly mean level: Compute, per region and month, mean SPEI values 
### (both local and global) and run comparison based on that 
### Compute regressions/corr. tests per region

mon_spei <- data.frame(
    spei_be %>% 
    mutate(date = paste(Month,Year, sep = "_")) %>% 
    group_by(Region,date) %>%
    summarise(
        SPEI3_local = mean(SPEI3, na.rm = TRUE),
        SPEI6_local = mean(SPEI6, na.rm = TRUE),
        SPEI12_local = mean(SPEI12, na.rm = TRUE),
        SPEI3_global = mean(get("SPEI-3"), na.rm = TRUE),
        SPEI6_global = mean(get("SPEI-6"), na.rm = TRUE),
        SPEI12_global = mean(get("SPEI-12"), na.rm = TRUE)
    ) # eo summarise
) # eo ddf
# Check
dim(mon_spei); head(mon_spei)
str(mon_spei)
summary(mon_spei)

### Compare global vs. local - Overall
subset <- na.omit(mon_spei) # head(subset)
cor(subset$SPEI3_local, subset$SPEI3_global, method = "spearman")
# 0.873
cor(subset$SPEI6_local, subset$SPEI6_global, method = "spearman")
# 0.875
cor(subset$SPEI12_local, subset$SPEI12_global, method = "spearman")
# 0.884

### Compare regionally 
# HND
subset <- na.omit(mon_spei[mon_spei$Region == "HND",]) # head(subset)
cor(subset$SPEI3_local, subset$SPEI3_global, method = "spearman")
# 0.854
cor(subset$SPEI6_local, subset$SPEI6_global, method = "spearman")
# 0.861
cor(subset$SPEI12_local, subset$SPEI12_global, method = "spearman")
# 0.867

# SCH
subset <- na.omit(mon_spei[mon_spei$Region == "SCH",]) # head(subset)
cor(subset$SPEI3_local, subset$SPEI3_global, method = "spearman")
# 0.878
cor(subset$SPEI6_local, subset$SPEI6_global, method = "spearman")
# 0.881
cor(subset$SPEI12_local, subset$SPEI12_global, method = "spearman")
# 0.890

# SWA
subset <- na.omit(mon_spei[mon_spei$Region == "SWA",]) # head(subset)
cor(subset$SPEI3_local, subset$SPEI3_global, method = "spearman")
# 0.889
cor(subset$SPEI6_local, subset$SPEI6_global, method = "spearman")
# 0.886
cor(subset$SPEI12_local, subset$SPEI12_global, method = "spearman")
# 0.901

### -> All coor coeff > 0.85 and highly signif

### Plot Fig.S5: Facet per Region - Panel with 3 rows: SPEI-3/6/12. Add a dashed 1:1 line in black
### Colour points by time -> check how deviations between local and global change in time. 
### Need to add Date to 'mon_spei' first then 
mon_spei$Date <- as.Date(paste0("01_", mon_spei$date), format = "%d_%m_%Y")
# Actually, just color by Year
mon_spei$Year <- lubridate::year(mon_spei$Date)

# class(mon_spei$Date); head(mon_spei$Date)

FigS5A <- ggplot(data = mon_spei, aes(x = SPEI3_local, y = SPEI3_global, fill = Year)) +
    geom_point(pch = 21, colour = "black") + geom_abline(size = 1, slope = 1, linetype = "dashed", color = "#225ea8") + 
    scale_fill_viridis(name = "Year", option = "rocket") + 
    xlab("Monthly mean SPEI-3 (BE-EXTREME reconstruction)") + ylab("Monthly mean SPEI-3 (Global SPEI)") +
    theme_bw() + facet_wrap(.~factor(Region))

ggsave("FigS5A.pdf", plot = FigS5A, device = "pdf", dpi = 300, height = 3, width = 9)

FigS5B <- ggplot(data = mon_spei, aes(x = SPEI6_local, y = SPEI6_global, fill = Year)) +
    geom_point(pch = 21, colour = "black") + geom_abline(size = 1, slope = 1, linetype = "dashed", color = "#225ea8") + 
    scale_fill_viridis(name = "Year", option = "rocket") + 
    xlab("Monthly mean SPEI-6 (BE-EXTREME reconstruction)") + ylab("Monthly mean SPEI-6 (Global SPEI)") +
    theme_bw() + facet_wrap(.~factor(Region))

ggsave("FigS5B.pdf", plot = FigS5B, device = "pdf", dpi = 300, height = 3, width = 9)


FigS5C <- ggplot(data = mon_spei, aes(x = SPEI12_local, y = SPEI12_global, fill = Year)) +
    geom_point(pch = 21, colour = "black") + geom_abline(size = 1, slope = 1, linetype = "dashed", color = "#225ea8") + 
    scale_fill_viridis(name = "Year", option = "rocket") + 
    xlab("Monthly mean SPEI-12 (BE-EXTREME reconstruction)") + ylab("Monthly mean SPEI-12 (Global SPEI)") +
    theme_bw() + facet_wrap(.~factor(Region))

ggsave("FigS5C.pdf", plot = FigS5C, device = "pdf", dpi = 300, height = 3, width = 9)

require("ggpubr")
FigS5 <- ggarrange(FigS5A,FigS5B,FigS5C, align = "hv", nrow = 3, ncol = 1, labels = c("A","B","C"), common.legend = TRUE, legend = "bottom")
ggsave("FigS5.pdf", plot = FigS5, device = "pdf", dpi = 300, height = 10.5, width = 10.5)


### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------
### ------------------------------------------------------------------------------------------------------------