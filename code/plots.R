library(sf)
library(ggplot2)
library(ggridges)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

covariate <- read.csv("../datasets/LUCAS_workon.csv")
shapefile <- st_read("../datasets/shapefile/LUCAS_workon_shape.shp")
LC_legend <- read.csv("../datasets/LC_description.csv")
LU_legend <- read.csv("../datasets/LU_description.csv")
data <- merge(shapefile, covariate, by = "POINTID")
rm(covariate, shapefile)

data$Soil_Stones <- as.factor(data$Soil_Stones)
data$LC_0 <- as.factor(data$LC_0)
data$LC_2018 <- as.factor(data$LC_2018)
data$LU_2018 <- as.factor(data$LU_2018)

table(data$LC_0)
unique(substr(as.character(data$LC_2018), 1, 1))
unique((LC_legend$LC0_Desc_2018))

legend <- data.frame(letter = unique(substr(as.character(LC_legend$LC_2018), 1, 1)), description = unique(LC_legend$LC0_Desc_2018))

data_cropland <- data[which(data$LC_0=='B'),]
data_woodland <- data[which(data$LC_0=='C'),]
data_grassland <- data[which(data$LC_0=='E'),]

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

colors <- c('A'= 'black', 'B'= 'orange', 'C'= 'forestgreen', 'D'= 'deeppink2', 'E'= 'green', 'F'= 'brown', 'G'= 'cyan', 'H'= 'blue')

ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") + # Europe map
  geom_sf(data = data, aes(colour = LC_0), size = 1.5, alpha = 0.7) + 
  scale_color_manual(name = "Land Cover", values = colors, aes(size = 2), labels = legend$description) +
  guides(color = guide_legend(override.aes = list(size = 4))) + # Enlarge legend dots 
  theme_minimal() +
  labs(title = "Points on Europe Map")

#plots for three main cover


plot_cropland <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_cropland, color = colors[2], size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Cropland point")

plot_woodland <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_woodland, color = colors[3], size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Woodland point")

plot_grassland <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_grassland, color = colors[5], size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Grassland point")

# Arrange in a three-column layout
(plot_cropland | plot_woodland | plot_grassland)

datam <- rbind(data_cropland, data_woodland, data_grassland)

boxplot(datam$BDsample_0)
datam$y <- c(rep(1,nrow(data_cropland)),rep(2,nrow(data_woodland)),rep(3,nrow(data_grassland)))

# Create the plot
ggplot(datam, aes(x = BDsample_0, y = y, fill = LC_0)) +
  geom_boxplot(width = 0.15, alpha = 0.5, outlier.shape = NA, color = "black", position = position_nudge(y = -0.2)) + # Boxplot with higher opacity
  geom_jitter(aes(color = LC_0, y = y - 0.2), size = 1, alpha = 0.2, 
              position = position_jitter(height = 0.05), pch = 1) +  # Jittered points
  geom_density_ridges(alpha = 0.6, scale = 0.7, rel_min_height = 0.01) +           # Density ridges
  labs(x = "Bulk density (g cm⁻³)", y =NULL, fill = "Land Cover") +
  theme_minimal() +
  scale_fill_manual(values = colors[c(2,3,5)], labels = legend$description[c(2,3,5)])   +           # Boxplot fill colors
  scale_color_manual(values = colors[c(2,3,5)], guide= 'none') +
  theme(axis.text.y = element_blank())





                       