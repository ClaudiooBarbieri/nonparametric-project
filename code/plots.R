library(sf)
library(ggplot2)
library(ggridges)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(dplyr)

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

# split data according to main covers
data_cropland <- data[which(data$LC_0=='B'),]
data_woodland <- data[which(data$LC_0=='C'),]
data_grassland <- data[which(data$LC_0=='E'),]

# get the UE map to print
world <- ne_countries(scale = "medium", returnclass = "sf")
europe <- world[world$continent == "Europe", ]
cyprus <- world[world$iso_a2 == "CY", ]
europe_with_cyprus <- rbind(europe, cyprus)

# get the NUTS map to print
NUTS <- st_read("../datasets/NUTS/NUTS_RG_20M_2016_4326.shp")
lat_range <- c(30, 70)  # latitude from 35°N to 70°N
lon_range <- c(-10, 35) # longitude from 10°W to 35°E
bounding_box <- st_as_sfc(st_bbox(c(xmin = lon_range[1], xmax = lon_range[2],
                                    ymin = lat_range[1], ymax = lat_range[2]), 
                                  crs = st_crs(NUTS)))
NUTS <- st_crop(NUTS, bounding_box)
europe_ue <- st_crop(europe_with_cyprus, bounding_box)
nuts_level <- 2
NUTS <- NUTS[which(NUTS$LEVL_CODE==nuts_level),]
rm(world, europe, cyprus, europe_with_cyprus, bounding_box, lat_range, lon_range)

# colors for land cover
colors <- c('A'= 'black', 'B'= 'orange', 'C'= 'forestgreen', 'D'= 'deeppink2', 'E'= 'green', 'F'= 'brown', 'G'= 'cyan', 'H'= 'blue')

ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") + # Europe map
  geom_sf(data = data, aes(colour = LC_0), size = 1.5, alpha = 0.7) + 
  scale_color_manual(name = "Land Cover", values = colors, aes(size = 2), labels = legend$description) +
  guides(color = guide_legend(override.aes = list(size = 4))) + # Enlarge legend dots 
  theme_minimal() +
  labs(title = "Points on Europe Map")

#plots for three main cover
colors <- c('cropland'= 'orange', 'woodland'= 'forestgreen', 'grassland' = 'green')

plot_cropland <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_cropland, color = colors[1], size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Cropland point")

plot_woodland <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_woodland, color = colors[2], size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Woodland point")

plot_grassland <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_grassland, color = colors[3], size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Grassland point")

# Arrange in a three-column layout
(plot_cropland | plot_woodland | plot_grassland)

data <- bind_rows(
  data_cropland,
  data_woodland,
  data_grassland,
  .id = "Landcover"
)
data$Landcover <- as.factor(data$Landcover)
levels(data$Landcover) <- c("cropland", "woodland", "grassland")

boxplot(data$BDsample_0)

# Create the plot
ggplot(data, aes(x = BDsample_0, y = 4-as.numeric(Landcover) , fill = Landcover)) +
  geom_boxplot(width = 0.15, alpha = 0.5, outlier.shape = NA, color = "black", position = position_nudge(y = -0.2)) + # Boxplot with higher opacity
  geom_jitter(aes(color = Landcover, y = 4-as.numeric(Landcover) - 0.2), size = 1, alpha = 0.2, 
              position = position_jitter(height = 0.05), pch = 1) +  # Jittered points
  geom_density_ridges(alpha = 0.6, scale = 0.7, rel_min_height = 0.01) +           # Density ridges
  labs(x = "Bulk density (g cm⁻³)", y =NULL, fill = "Land Cover") +
  theme_minimal() +
  scale_fill_manual(values = colors[c(3,2,1)], labels = levels(data$Landcover))   +           # Boxplot fill colors
  scale_color_manual(values = colors[c(3,2,1)], guide= 'none') +
  theme(axis.text.y = element_blank(), legend.text = element_text(size=20), legend.title = element_blank(), legend.key.spacing.y =  unit(10, "cm"))

# by NUTS
my_nuts <- read.csv("../datasets/NUTS/NUTS.csv")
my_nuts <- my_nuts[,c(1,nuts_level+2)]
names(my_nuts)[2] <- "NUTS"
data <- inner_join(data, my_nuts, 'POINTID')
rm(my_nuts)

data_nuts <- NUTS %>%
  st_join(data)

# not complete understood how it works yet
ggplot(data_nuts) +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(aes(fill = pH_H2O), color = "black", size = 0.1) +
  scale_fill_gradientn(
    colors = c("yellow", "green", "cyan", "blue"),
    values = scales::rescale(c(4, 5, 6, 7, 9)),  # Adjust as needed
    name = "pH"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Average pH (H2O)",
       subtitle = "Aggregated at NUTS2 Level")
