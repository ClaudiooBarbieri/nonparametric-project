library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

##### dataset ####

covariate <- read.csv("../datasets/LUCAS_final.csv")
shapefile <- st_read("../datasets/shapefile/LUCAS_final.shp")
LC_legend <- read.csv("../datasets/LC_description.csv")
LU_legend <- read.csv("../datasets/LU_description.csv")

# NOTE - if not final adjust before merge (drop duplicated column in shapefile and ensure about pointid names matching)
colnames(covariate)[2] <- "POINTID"
data <- merge(covariate, shapefile, by = "POINTID")
rm(covariate, shapefile)

#drop surely useless columns
data <- subset(data , select = - c(SURVEY_DATE, NUTS_0_2018, NUTS_1_2018, NUTS_2_2018, NUTS_3_2018, TH_LAT, TH_LONG))
attach(data)

# detailed cleansing

table(Depth)
sum(is.na(Depth))
# drop row with depth!=0-20 (and also column related to measure at different depth)
depth <- "0-20 cm"
data <- data[which(Depth=="0-20 cm"), ]
data <- subset(data , select = - c(Depth, OC..20.30.cm., CaCO3..20.30.cm., BD.0.10, BD.10.20, BD.20.30))

# Coarse, Clay, Sand and Silt are 13362 NA; keep or drop??

# convert true type variable
data$Soil_Stones <- as.factor(Soil_Stones)
data$OC_2018 <- as.numeric(OC_2018)
data$CaCO3_2018 <- as.numeric(CaCO3_2018)
data$P_2018 <- as.numeric(P_2018)
data$N_2018 <- as.numeric(N_2018)
data$K_2018 <- as.numeric(K_2018)
data$LC_2018 <- as.factor(LC_2018)
data$LU_2018 <- as.factor(LU_2018)


##### EUROPE #####
# get the UE map to print
world <- ne_countries(scale = "medium", returnclass = "sf")
europe <- world[world$continent == "Europe", ]
cyprus <- world[world$iso_a2 == "CY", ]
europe_with_cyprus <- rbind(europe, cyprus)
lat_range <- c(30, 70)  # latitude from 35°N to 70°N
lon_range <- c(-10, 35) # longitude from 10°W to 35°E
bounding_box <- st_as_sfc(st_bbox(c(xmin = lon_range[1], xmax = lon_range[2],
                                    ymin = lat_range[1], ymax = lat_range[2]), 
                                  crs = st_crs(europe_with_cyprus)))
europe_ue <- st_crop(europe_with_cyprus, bounding_box)
rm(world, europe, cyprus, europe_with_cyprus, bounding_box)

##### PLOTS #####
# plot 
ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") + # Europe map
  geom_sf(data = lucas_2018, aes(colour = '2018'), size = 0.5, alpha = 0.7) + # LUCAS 2018 points
  geom_sf(data = lucas_2015, aes(colour = '2015'), size = 0.5, alpha = 0.7) + # LUCAS 2018 points
  #geom_sf(data = lucas_2018_common, aes(colour = '2015&2018'), size = 0.5, alpha = 0.9) +
  scale_color_manual(name = "Year", values = c("2018" = "red", "2015" = "green", "2015&2018" = 'purple'), aes(size = 2)) +
  guides(color = guide_legend(override.aes = list(size = 4))) + # Enlarge legend dots
  theme_minimal() +
  labs(title = "LUCAS Soil 2018 Points on Europe Map")

# Plot for 2015 point
plot_2015 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = lucas_2015, aes(color = "2015"), size = 1, alpha = 0.7) +
  scale_color_manual(name = "Year", values = c("2018" = "red", "2015" = "green", "2015&2018" = 'purple'), aes(size = 2)) +
  guides(color = guide_legend(override.aes = list(size = 4))) + # Enlarge legend dots
  theme_minimal() +
  labs(title = "LUCAS 2015", color = "Year")

# Plot for 2018 point
plot_2018 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = lucas_2018, aes(color = "2018"), size = 1, alpha = 0.7) +
  scale_color_manual(name = "Year", values = c("2018" = "red", "2015" = "green", "2015&2018" = 'purple'), aes(size = 2)) +
  guides(color = guide_legend(override.aes = list(size = 4))) + # Enlarge legend dots
  theme_minimal() +
  labs(title = "LUCAS 2018", color = "Year")

# plot common points
combined_plot <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = lucas_2018_common, aes(color = '2015&2018'), size = 1, alpha = 0.7) +
  scale_color_manual(name = "Year", values = c("2018" = "red", "2015" = "green", "2015&2018" = 'purple'), aes(size = 2)) +
  guides(color = guide_legend(override.aes = list(size = 4))) + # Enlarge legend dots
  theme_minimal() +
  labs(title = "LUCAS Combined 2015 & 2018", color = "Year")

# Arrange in a three-column layout
(plot_2015 | plot_2018 | combined_plot)