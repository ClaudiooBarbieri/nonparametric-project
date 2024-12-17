# read data
covariate <- read.csv("../datasets/LUCAS_workon.csv")
LC_legend <- read.csv("../datasets/LC_description.csv")

#consider only BD and LC factor
data <- subset(covariate , select = c(BDsample_0, LC_0))
data$LC_0 <- as.factor(data$LC_0)
head(data)

#boxplots by LC
table(data$LC_0)

plot(data$LC_0, data$BDsample_0, xlab='LC', ylab = "Bulk density", col=rainbow(g),main='Bulk Density by Land Cover')
# split data according to main covers
data_artificial_land <- data[which(data$LC_0=='A'),]
data_cropland <- data[which(data$LC_0=='B'),]
data_woodland <- data[which(data$LC_0=='C'),]
data_shrubland <- data[which(data$LC_0=='D'),]
data_grassland <- data[which(data$LC_0=='E'),]
data_bareland <- data[which(data$LC_0=='F'),]
data_water <- data[which(data$LC_0=='G'),]
data_wetland <- data[which(data$LC_0=='H'),]
data <- bind_rows(
  data_cropland,
  data_woodland,
  data_shrubland,
  data_grassland,
  data_bareland,
  .id = "Landcover"
)
data$Landcover <- as.factor(data$Landcover)
levels(data$Landcover) <- c("cropland",
                            "woodland",
                            "shrubland",
                            "grassland",
                            "bareland")
colors <- c("orange", "darkgreen", "lightgreen", "green", "blue") 
ggplot(data, aes(x = BDsample_0, y = 4-as.numeric(Landcover) , fill = Landcover)) +
  geom_boxplot(width = 0.15, alpha = 0.5, outlier.shape = NA, color = "black", position = position_nudge(y = -0.2)) + # Boxplot with higher opacity
  geom_jitter(aes(color = Landcover, y = 4-as.numeric(Landcover) - 0.2), size = 1, alpha = 0.2, 
              position = position_jitter(height = 0.05), pch = 1) +  # Jittered points
  geom_density_ridges(alpha = 0.6, scale = 0.7, rel_min_height = 0.01) +           # Density ridges
  labs(x = "Bulk density (g cm⁻³)", y =NULL, fill = "Land Cover", color = "Land Cover") +
  theme_minimal() +
  scale_fill_manual(values = colors[c(5,4,3,2,1)])   +           # Boxplot fill colors
  scale_color_manual(values = colors[c(5,4,3,2,1)]) +
  theme(axis.text.y = element_blank(), legend.text = element_text(size = 20),
        legend.title = element_text(size = 20), legend.key.spacing.y =  unit(1, "cm"))

#test gaussianity
# Set up an 8-panel plotting layout (2 rows, 4 columns)
par(mfrow = c(2, 4)) # 2 rows and 4 columns of plots

# Loop through each category and create the Q-Q plots
for (i in c("A", "B", "C", "D", "E", "F", "G", "H")) {
  qqnorm(data$BDsample_0[data$LC_0 == i], 
         main = paste("LC_0 =", i)) # Title for each plot
  qqline(data$BDsample_0[data$LC_0 == i]) # Adds a reference line
}
par(mfrow = c(1, 1))

#shapiro.test
p_values=NULL
for (i in c("A", "B", "C", "D", "E", "F","H")) {
  p_values=c(p_values, shapiro.test(data$BDsample_0[data$LC_0 == i]))
}
p_values

# H0: tau1 = tau2 = tau3 = tau4 = tau5 = tau6 = tau7 = tau8 = 0
# the BDs belong to the same population
g <- nlevels(data$LC_0)
n <- dim(data)[1]
# H1: (H0)^c
# the BDs belong to several different population
# Parametric test:
fit <- aov(BDsample_0 ~ LC_0, data = data)
summary(fit)
# Permutation test:
# Test statistic: F stat
T0 <- summary(fit)[[1]][1,4]
T0

# what happens if we permute the data?
permutazione <- sample(1:n)
weight_perm <- data$BDsample_0[permutazione]
fit_perm <- aov(weight_perm ~ data$Landcover)
summary(fit_perm)

plot(data$LC_0, weight_perm, xlab='treat',col=rainbow(g),main='Permuted Data')


# CMC to estimate the p-value
B <- 1000 # Number of permutations
T_stat <- numeric(B) 

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  weight_perm <- data$BDsample_0[permutation]
  fit_perm <- aov(weight_perm ~ data$LC_0)
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}

layout(1)
hist(T_stat,breaks=30, main = "Histogram of the test_statistic", xlab = "test_statistic", col = "lightblue")
abline(v=T0,col=3,lwd=2)
       
plot(ecdf(T_stat))
abline(v=T0,col=3,lwd=4)

# p-value
p_val <- sum(T_stat>=T0)/B
p_val
# we reject the null hypothesis


library(sf)
library(ggplot2)
library(ggridges)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(dplyr)
library(scales)
#read data
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

#get Europe map
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

#color scale
common_fill_scale <- scale_fill_gradientn(
  colors = c("#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026"),
  values = scales::rescale(c(0,0.5,1,1.5,2)),
  name = "Bulk Density",
  na.value = "darkgray"  # Color for missing data
)

data_cropland <- data[which(data$LC_0=='B'),]
data_woodland <- data[which(data$LC_0=='C'),]
data_grassland <- data[which(data$LC_0=='E'),]

data <- bind_rows(
  data_cropland,
  data_woodland,
  data_grassland,
  .id = "Landcover"
)
data$Landcover <- as.factor(data$Landcover)
levels(data$Landcover) <- c("cropland", "woodland", "grassland")

my_nuts <- read.csv("../datasets/NUTS/NUTS.csv")
my_nuts <- my_nuts[,c(1,nuts_level+2)]
names(my_nuts)[2] <- "NUTS"
data <- inner_join(data, my_nuts, 'POINTID')
rm(my_nuts)

data_nuts <- NUTS %>%
  st_join(data)
data_nuts_c <- data_nuts[data_nuts$Landcover=='cropland',]
data_nuts_w <- data_nuts[data_nuts$Landcover=='woodland',]
data_nuts_g <- data_nuts[data_nuts$Landcover=='grassland',]

#median BD by NUTS2 for each LC (3 main)
plot_by_nuts <- function(v){
  
  aggregated_data <- data_nuts_c %>%
    group_by(NUTS) %>%                  # Group by NUTS region
    summarize(avg_BD = median(BDsample_0))
  
  avg_crop <- ggplot(aggregated_data) +
    geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
    geom_sf(aes(fill = avg_BD), color = "black", size = 0.1) +
    common_fill_scale +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    labs(title = "Median BD Cropland",
         subtitle = "Aggregated at NUTS1 Level")
  
  aggregated_data <- data_nuts_w %>%
    group_by(NUTS) %>%                  # Group by NUTS region
    summarize(avg_BD = median(BDsample_0))
  
  avg_wood <- ggplot(aggregated_data) +
    geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
    geom_sf(aes(fill = avg_BD), color = "black", size = 0.1) +
    common_fill_scale +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    labs(title = "Median BD Woodland",
         subtitle = "Aggregated at NUTS1 Level")
  
  aggregated_data <- data_nuts_g %>%
    group_by(NUTS) %>%                  # Group by NUTS region
    summarize(avg_BD = median(BDsample_0))
  
  avg_grass <- ggplot(aggregated_data) +
    geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
    geom_sf(aes(fill = avg_BD), color = "black", size = 0.1) +
    common_fill_scale +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    ) +
    labs(title = "Median BD Grassland",
         subtitle = "Aggregated at NUTS1 Level")
  
  (avg_crop | avg_wood | avg_grass) + 
    plot_layout(guides = "collect")
}
v <- 'BD'
plot_by_nuts(v)

#median BD by NUTS2
bd_over_europe <-function(v){
  aggregated_data <- data_nuts %>%
    group_by(NUTS) %>%                  # Group by NUTS region
    summarize(avg_BD = median(BDsample_0))
  
  avg_crop <- ggplot(aggregated_data) +
    geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
    geom_sf(aes(fill = avg_BD), color = "black", size = 0.1) +
    common_fill_scale +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    ) +
    labs(title = "Median BD",
         subtitle = "Aggregated at NUTS2 Level")
    avg_crop
} 
v <- 'BD'
bd_over_europe(v)


colors <- c('A'= 'black', 'B'= 'orange', 'C'= 'forestgreen', 'D'= 'deeppink2', 'E'= 'green', 'F'= 'brown', 'G'= 'cyan', 'H'= 'blue')

ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") + # Europe map
  geom_sf(data = data, aes(colour = LC_0), size = 0.1) + 
  scale_color_manual(name = "Land Cover", values = colors, aes(size = 2)) +
  guides(color = guide_legend(override.aes = list(size = 4))) + # Enlarge legend dots 
  theme_minimal() +
  labs(title = "Points on Europe Map")

#plots for three main cover

plot_cropland <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_cropland, aes(color = BDsample_0), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma") +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) +
  labs(
    title = "Bulk density",
    subtitle = "Cropland Point",
    color = "Value Legend"  # Adjust legend title
  )

plot_woodland <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_woodland, aes(color = BDsample_0), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma") +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) +
  labs(
    title = "Bulk density",
    subtitle = "Woodland Point",
    color = "Value Legend"  # Adjust legend title
  )

plot_grassland <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_grassland, aes(color = BDsample_0), size = 1.5) +  # Replace `value_column` with your column name
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
    subtitle = "Grassland Point",
    color = "Value Legend"  # Adjust legend title
  )

# Arrange in a three-column layout
(plot_cropland | plot_woodland | plot_grassland)

#same for all the points
plot_bulk_density <- ggplot() +
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
plot_bulk_density

#look at the covariates
#median Elev by NUTS2
common_fill_scale <- scale_fill_gradientn(
  colors = c("yellow", "orange", "red", "violet", "purple"),
  values = scales::rescale(c(0,0.5,1,1.5,2)),
  name = "Elev",
  na.value = "darkgray"  # Color for missing data
)
Elev_over_europe <-function(v){
  aggregated_data <- data_nuts %>%
    group_by(NUTS) %>%                  # Group by NUTS region
    summarize(avg_Elev = median(Elev))
  
  avg_crop <- ggplot(aggregated_data) +
    geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
    geom_sf(aes(fill = avg_Elev), color = "black", size = 0.1) +
    common_fill_scale +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    ) +
    labs(title = "Median Elevation",
         subtitle = "Aggregated at NUTS2 Level")
  avg_crop
} 
v <- 'elev'
Elev_over_europe(v)

#median pH_H2O by NUTS2
common_fill_scale <- scale_fill_gradientn(
  colors = c("yellow", "orange", "red", "violet", "purple"),
  values = scales::rescale(c(0,0.5,1,1.5,2)),
  name = "pH_H2O",
  na.value = "darkgray"  # Color for missing data
)
pH_H2O_over_europe <-function(v){
  aggregated_data <- data_nuts %>%
    group_by(NUTS) %>%                  # Group by NUTS region
    summarize(avg_pH_H2O = median(pH_H2O))
  
  avg_crop <- ggplot(aggregated_data) +
    geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
    geom_sf(aes(fill = avg_pH_H2O), color = "black", size = 0.1) +
    common_fill_scale +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    ) +
    labs(title = "Median pH_H2O",
         subtitle = "Aggregated at NUTS2 Level")
  avg_crop
} 
v <- 'pH_H2O'
pH_H2O_over_europe(v)

#median OC_2018 by NUTS2
common_fill_scale <- scale_fill_gradientn(
  colors = c("yellow", "orange", "red", "violet", "purple"),
  values = scales::rescale(c(0,0.5,1,1.5,2)),
  name = "OC_2018",
  na.value = "darkgray"  # Color for missing data
)
OC_2018_over_europe <-function(v){
  aggregated_data <- data_nuts %>%
    group_by(NUTS) %>%                  # Group by NUTS region
    summarize(avg_OC_2018 = median(OC_2018))
  
  avg_crop <- ggplot(aggregated_data) +
    geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
    geom_sf(aes(fill = avg_OC_2018), color = "black", size = 0.1) +
    common_fill_scale +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    ) +
    labs(title = "Median OC_2018",
         subtitle = "Aggregated at NUTS2 Level")
  avg_crop
} 
v <- 'OC_2018'
OC_2018_over_europe(v)

#median K_2018 by NUTS2
common_fill_scale <- scale_fill_gradientn(
  colors = c("yellow", "orange", "red", "violet", "purple"),
  values = scales::rescale(c(0,0.5,1,1.5,2)),
  name = "K_2018",
  na.value = "darkgray"  # Color for missing data
)
K_2018_over_europe <-function(v){
  aggregated_data <- data_nuts %>%
    group_by(NUTS) %>%                  # Group by NUTS region
    summarize(avg_K_2018 = median(K_2018))
  
  avg_crop <- ggplot(aggregated_data) +
    geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
    geom_sf(aes(fill = avg_K_2018), color = "black", size = 0.1) +
    common_fill_scale +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    ) +
    labs(title = "Median Potassium (K)",
         subtitle = "Aggregated at NUTS2 Level")
  avg_crop
} 
v <- 'K_2018'
K_2018_over_europe(v)

#1_2_3 quartile of BD by NUTS2
common_fill_scale <- scale_fill_gradientn(
  colors = c("yellow", "orange", "red", "violet", "purple"),
  values = scales::rescale(c(0,0.5,1,1.5,2)),
  name = "Bulk Density",
  na.value = "darkgray"  # Color for missing data
)

plot_quantiles <- function(v){
  
  aggregated_data <- data_nuts %>%
    group_by(NUTS) %>%                  # Group by NUTS region
    summarize(BD_1 = quantile(BDsample_0, 0.25, na.rm = T))
  
  quantile_1 <- ggplot(aggregated_data) +
    geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
    geom_sf(aes(fill = BD_1), color = "black", size = 0.1) +
    common_fill_scale +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    labs(title = "First BD quartile",
         subtitle = "Aggregated at NUTS2 Level")
  
  aggregated_data <- data_nuts %>%
    group_by(NUTS) %>%                  # Group by NUTS region
    summarize(BD_2 = quantile(BDsample_0, 0.5, na.rm = T))
  
  quantile_2 <- ggplot(aggregated_data) +
    geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
    geom_sf(aes(fill = BD_2), color = "black", size = 0.1) +
    common_fill_scale +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    labs(title = "Median BD",
         subtitle = "Aggregated at NUTS2 Level")
  
  aggregated_data <- data_nuts %>%
    group_by(NUTS) %>%                  # Group by NUTS region
    summarize(BD_3 = quantile(BDsample_0, 0.75, na.rm = T))
  
  quantile_3 <- ggplot(aggregated_data) +
    geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
    geom_sf(aes(fill = BD_3), color = "black", size = 0.1) +
    common_fill_scale +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    ) +
    labs(title = "Third BD quartile",
         subtitle = "Aggregated at NUTS2 Level")
  
  (quantile_1 | quantile_2 | quantile_3) + 
    plot_layout(guides = "collect")
}
v <- 'BD'
plot_quantiles(v)
