library(sf)
library(sp)
library(ggplot2)
library(ggridges)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(dplyr)
library(vegan)
library(spdep)
library(units)
library(adespatial)

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

table(data$LC_0)
unique(substr(as.character(data$LC_2018), 1, 1))
unique((LC_legend$LC0_Desc_2018))

legend <- data.frame(letter = unique(substr(as.character(LC_legend$LC_2018), 1, 1)), description = unique(LC_legend$LC0_Desc_2018))
rm(LC_legend, LU_legend)

# distance matrix of location
point_distance <- st_distance(data, data) #Haversine-great circle formula

# difference matrix of bulk density 
#bd_distance <- as.matrix(dist(data$BDsample_0))

# label distances
colnames(bd_distance) <- data$POINTID
colnames(point_distance) <- data$POINTID
rownames(bd_distance) <- data$POINTID
rownames(point_distance) <- data$POINTID

# Mantel tests
#m_test <- mantel(point_distance, bd_distance, permutations = 999) # seems cannot be done for spatial autocorrelation
# Extract the permutation distribution
#perm_dist <- m_test$perm
# Visualize the permutation distribution
#hist(perm_dist, breaks = 30, main = "Permutation Distribution of Mantel Statistic",
#    xlab = "Mantel Statistic", col = "lightblue", border = "black", xlim = c(-0.05,0.11))
#abline(v = m_test$statistic, col = "red", lwd = 2, lty = 2)
#legend("topright", legend = "Observed Statistic", col = "red", lwd = 2, lty = 2)

# Moran's I

# KNN based 
knn <- knearneigh(st_coordinates(data), k = 10) # list of knn
weights <- nb2listw(knn2nb(knn))
nb <- knn2nb(knn) # list of k nn of each point

# distance based
limit <- c(50000,75000,100000)

dnn_50 <- dnearneigh(data, 0, limit[1]/1000) # list of nb within 0-km
dnn_75 <- dnearneigh(data, 0, limit[2]/1000)
dnn_100 <- dnearneigh(data, 0, limit[3]/1000)
dnn_50
dnn_75
dnn_100

dnn <- dnn_100

# Gaussian decay weights
sigma <- 100000  # Set the decay parameter 
dist_matrix <- as.matrix(point_distance)
units(dist_matrix) <- NULL

weights_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))

# compute weights of the neighbors
for (i in 1:length(dnn)) {
  neighbors <- dnn[[i]]
  weights_matrix[i, neighbors] <- exp(-(dist_matrix[i, neighbors]^2) / (2 * sigma^2))
}

listw <- mat2listw(weights_matrix, style = "W", zero.policy = TRUE)

# Perform Moran's I
moran_result <- moran.mc(data$BDsample_0, listw = listw, zero.policy = TRUE, alternative = 'two.sided', nsim = 9999)
moran_result

plot(moran_result, xlab= 'Moran\' I')
abline(v=moran_result$statistic, col = 'red', lwd = 2)
moran_result
morans_scatter <- moran.plot(data$BDsample_0, listw = listw, zero.policy = TRUE, main = "Moran's I Scatter Plot")

#moran_result_np <- moranNP.randtest(data$BDsample_0, listw = listw, zero.policy = TRUE, alter = 'two-sided', nsim = 9999)
#moran_result_np

different_lag <- function(data, dist){
  dnn <- dnearneigh(data, 0, dist)
  
  # Gaussian decay weights
  sigma <- 100000  # Set the decay parameter 
  dist_matrix <- as.matrix(point_distance)
  units(dist_matrix) <- NULL
  
  weights_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))
  
  # compute weights of the neighbors
  for (i in 1:length(dnn)) {
    neighbors <- dnn[[i]]
    weights_matrix[i, neighbors] <- exp(-(dist_matrix[i, neighbors]^2) / (2 * sigma^2))
  }
  
  listw <- mat2listw(weights_matrix, style = "W", zero.policy = TRUE)
  
  # Perform Moran's I
  moran_result <- moran.mc(data$BDsample_0, listw = listw, zero.policy = TRUE, alternative = 'two.sided', nsim = 999)
  moran_result
}
different_lag(data, 50)
different_lag(data, 100)
different_lag(data, 150)
different_lag(data, 200)
different_lag(data, 250)
different_lag(data, 500)
different_lag(data, 1000)
# in the future can try with all covariate too
#moran.randtest(datab[,c(4,5,6,7,8,9,10,11,12,13,14)], listw = listw, zero.policy = TRUE, alternative = 'two.sided', nsim = 999)

# result of the permutation test


local_moran <- localmoran_perm(data$BDsample_0, listw = listw, zero.policy = TRUE, nsim = 9999)
local_moran
data$localmoran <- local_moran[,1]
data$significant <- local_moran[,5] < 0.05
data$p_value <- local_moran[, 5] 


# local moron pval
local_moran_pval <- ggplot(data = data) +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(aes(color = significant), size = 2) +
  scale_color_manual(
    values = c("lightgreen", "darkred"),
    labels = c("Not Significant", "Significant"),
    name = "Significance"
  ) +
  theme_minimal() +
  labs(title = "Significant Clusters (p < 0.05)")

#local moron val
local_moran_val <- ggplot(data) +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(aes(color = localmoran), size = 3) +  # Use Moran's I values
  scale_color_gradient2(
    low = "green",      # Negative values (outliers)
    mid = "white",     # Neutral
    high = "red",      # Positive values (clusters)
    midpoint = 0,      # Center the scale
    name = "Local Moran's I"
  ) +
  theme_minimal() +
  labs(title = "Local Moran's I for Bulk Density",
       subtitle = "Positive: Clusters; Negative: Spatial Outliers")

# Create cluster classifications
data$cluster <- with(data, ifelse(p_value < 0.05,
                                  ifelse(localmoran > 0 & BDsample_0 > mean(BDsample_0), "High-High",
                                         ifelse(localmoran > 0 & BDsample_0 < mean(BDsample_0), "Low-Low",
                                                ifelse(localmoran < 0 & BDsample_0 > mean(BDsample_0), "High-Low", "Low-High"))),
                                  "Non-Significant"))
data$cluster <- factor(data$cluster,
                       levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Non-Significant"))

local_moran_clust <- ggplot(data) +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(aes(color = cluster), size = 2) +
  scale_color_manual(values = c("High-High" = "red", 
                               "Low-Low" = "blue", 
                               "High-Low" = "orange", 
                               "Low-High" = "purple", 
                               "Non-Significant" = "grey"),
                     name = 'Cluster') +
  labs(title = "Local Moran's I Clusters",
       fill = "Cluster Type") +
  theme_minimal()

(local_moran_pval | local_moran_clust)


# Perform Geary's C
geary_result <- geary.mc(data$BDsample_0, listw = nb, zero.policy = TRUE,  nsim=9999)
print(geary_result)
plot(geary_result)

geary_result <- geary.test(data$BDsample_0, listw = nb, zero.policy = TRUE)
geary_stat <- geary_result$statistic

# Add Geary's C results to your spatial object
data$geary <- geary_stat

# Plot Geary's C using ggplot2
ggplot(data) +
  geom_sf(aes(fill = geary)) +  # Color regions by their Geary's C statistic
  scale_fill_viridis_c() +  # Adjust color palette
  theme_minimal() +
  ggtitle("Geary's C Local Spatial Autocorrelation") +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))



local_moran$quadr
