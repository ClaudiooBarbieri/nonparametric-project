#load libraries
library(sf)
library(GWmodel)
library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(splines)

# load data
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

data <- data[!is.na(data$EC_2018), ]
coords <- st_coordinates(data)  
n <- dim(data)[1]
dist_mat <- gw.dist(coords, coords, focus=0, p=2, theta=0, longlat=T)


bw_opt <- bw.gwr(BDsample_0 ~ bs(Elev, degree=4) + bs(pH_H2O, degree=4) + bs(EC_2018, degree=4), data, approach="CV", kernel="gaussian",
                 adaptive=FALSE, p=2, theta=0, longlat=T)
fit4 <- gwr.basic(BDsample_0 ~ bs(Elev, degree=4) + bs(pH_H2O, degree=4) + bs(EC_2018, degree=4), data, coords, bw=bw_opt, kernel="gaussian",
                  adaptive=FALSE, p=2, theta=0, longlat=T,dMat=dist_mat)
ss4 <- sum((data$BDsample_0-fit4$lm$fitted.values)^2)


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

res <- fit4$lm$residuals
plot_residuals <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = res), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma") +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density",
    color = "Value Legend"  # Adjust legend title
  )
plot_residuals

fitted <- fit4$lm$fitted.values
plot_fitted <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = fitted), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma") +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density",
    color = "Value Legend"  # Adjust legend title
  )

plot_values <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = BDsample_0), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma") +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density",
    color = "Value Legend"  # Adjust legend title
  )
(plot_fitted | plot_values)
