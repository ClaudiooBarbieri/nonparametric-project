library(sf)
library(ggplot2)
library(ggridges)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(dplyr)
library(vegan)
library(spdep)
library(units)


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

point_distance <- st_distance(data, data)
bd_distance <- as.matrix(dist(data$BDsample_0, method = 'manhattan'))

m_test <- mantel(point_distance,bd_distance)

# Extract the permutation distribution
perm_dist <- m_test$perm

# Visualize the permutation distribution
hist(perm_dist, breaks = 30, main = "Permutation Distribution of Mantel Statistic",
     xlab = "Mantel Statistic", col = "lightblue", border = "black", xlim = c(-0.05,0.11))
abline(v = m_test$statistic, col = "red", lwd = 2, lty = 2)
legend("topright", legend = "Observed Statistic", col = "red", lwd = 2, lty = 2)

knn <- knearneigh(st_coordinates(data), k = 5)  # 5 nearest neighbors
weights <- nb2listw(knn2nb(knn))

distance_threshold <- 50000
distance_threshold <- set_units(distance_threshold, m)
weights_matrix <- 1 / point_distance
weights_matrix[point_distance > distance_threshold] <- 0

# Extract numeric values
numeric_values <- as.numeric(weights_matrix)

# Convert to numeric matrix
numeric_matrix <- matrix(numeric_values, nrow = nrow(weights_matrix),
                         ncol = ncol(weights_matrix))
diag(numeric_matrix) <- 0

nb <- mat2listw(numeric_matrix, style = 'M')

# Perform Moran's I
moran_result <- moran.test(data$BDsample_0, listw = nb, zero.policy = TRUE)
print(moran_result)

# Perform Geary's C
geary_result <- geary.test(data$BDsample_0, listw = nb, zero.policy = TRUE)
print(geary_result)




