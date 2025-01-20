library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

##### dataset ####

covariate <- read.csv("../datasets/extra/LUCAS_final.csv")
shapefile <- st_read("../datasets/extra/shp/LUCAS_final.shp")
LC_legend <- read.csv("../datasets/LC_description.csv")
LU_legend <- read.csv("../datasets/LU_description.csv")
                    
# NOTE - if not final adjust before merge (drop duplicated column in shapefile and ensure about pointid names matching)
colnames(covariate)[2] <- "POINTID"
data <- merge(shapefile, covariate, by = "POINTID")
rm(covariate, shapefile)

colSums(!is.na(data))

# save nuts for later eventually 
NUTS <-  subset(data , select = c(POINTID, NUTS_0_2018, NUTS_1_2018, NUTS_2_2018, NUTS_3_2018))
NUTS$geometry <- NULL
write.csv(NUTS, "../datasets/NUTS/NUTS.csv", row.names = FALSE)

#drop surely useless columns
data <- subset(data , select = - c(SURVEY_DATE, NUTS_0_2018, NUTS_1_2018, NUTS_2_2018, NUTS_3_2018, TH_LAT, TH_LONG))

# detailed cleansing

table(data$Depth)
sum(is.na(data$Depth))
# drop row with depth!=0-20 (and also column related to measure at different depth)
depth <- "0-20 cm"
data <- data[which(data$Depth=="0-20 cm"), ]
data <- subset(data , select = - c(Depth, OC..20.30.cm., CaCO3..20.30.cm., BD.0.10, BD.0.20, BD.10.20, BD.20.30))

# convert true type variable
data$Soil_Stones <- as.factor(data$Soil_Stones)
data$OC_2018 <- as.numeric(data$OC_2018)
data$CaCO3_2018 <- as.numeric(data$CaCO3_2018)
data$P_2018 <- as.numeric(data$P_2018)
data$N_2018 <- as.numeric(data$N_2018)
data$K_2018 <- as.numeric(data$K_2018)
data$LC_2018 <- as.factor(data$LC_2018)
data$LU_2018 <- as.factor(data$LU_2018)


boxplot(data$Elev)
#clearly some not possible values -> remove
data <- data[which(data$Elev<=4000),]
data[which.max(data$Elev),]

nrow(data) - sum(is.na(data$Coarse))
nrow(data) - sum(is.na(data$Clay))
nrow(data) - sum(is.na(data$Silt))
nrow(data) - sum(is.na(data$Sand))

# Filter rows where none of the columns are NA
Coarse_Clay_Sand_Silt <- data[!is.na(data$Coarse) & !is.na(data$Clay) & !is.na(data$Silt) & !is.na(data$Sand), ]
sum(is.na(Coarse_Clay_Sand_Silt$BDsample_0))
# 2492 points for which we have Coarse, Clay, Sand and Silt and Bulk density -> drop all column (keep aside)
Coarse_Clay_Sand_Silt <- subset(Coarse_Clay_Sand_Silt , select = c(POINTID, Coarse, Clay, Sand, Silt))
data <- subset(data , select = - c(Coarse, Clay, Sand, Silt))

#colSums(is.na(data))#

# keep just point for which we have bulk density
data <- data[which(!is.na(data$BDsample_0)),]

#set aside nonìt observed values but formula's evaluated
BD_Coarse_estimate <- subset(data , select = c(POINTID, BDfine_0_2, coarse_vol))
data <- subset(data , select = - c(BDfine_0_2, coarse_vol))

data$LC_0 <- factor(substr(as.character(data$LC_2018), 1, 1))
data <- data[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,19,14,15,16,17,18)] #reorder

colSums(is.na(data))

# some NA but lets go on with it, consider later
write.csv(BD_Coarse_estimate, "../datasets/extra/BD_Coarse_estimate.csv", row.names = FALSE)
write.csv(Coarse_Clay_Sand_Silt, "../datasets/extra/Coarse_Clay_Sand_Silt.csv", row.names = FALSE)
rm(BD_Coarse_estimate,Coarse_Clay_Sand_Silt,NUTS)

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

# plot 
ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") + # Europe map
  geom_sf(data = data, color = 'blue', size = 1, alpha = 0.7) + # points
  theme_minimal() +
  labs(title = "LUCAS Soil Points on Europe Map")

shp <- subset(data, select = c(POINTID, geometry))
st_write(shp, "../datasets/shapefile/LUCAS_workon_shape.shp", delete_layer = TRUE)
data$geometry <- NULL
write.csv(data, "../datasets/Lucas_workon.csv", row.names = FALSE)

############################################# 
#this was a trial, not actually used

covariate <- read.csv("../datasets/2018/LUCAS_2018_bd.csv")
shapefile <- st_read("../datasets/2018/shapefile/LUCAS_2018_bd.shp")
shapefile <- subset(shapefile, select = c(POINTID, geometry))
LC_legend <- read.csv("../datasets/2018/LC_description.csv")
LU_legend <- read.csv("../datasets/2018/LU_description.csv")

data <- merge(shapefile, covariate, by = "POINTID")
data <- data[order(data$POINTID),]
rm(covariate, shapefile)

colnames(data)

# save nuts for later eventually 
NUTS <-  subset(data , select = c(POINTID, NUTS_0, NUTS_1, NUTS_2, NUTS_3))

#drop surely useless columns
data <- subset(data , select = - c(SURVEY_DATE, NUTS_0, NUTS_1, NUTS_2, NUTS_3, TH_LAT, TH_LONG))

# detailed cleansing

table(data$Depth)
sum(is.na(data$Depth))
# drop row with depth!=0-20 (and also column related to measure at different depth)
depth <- "0-20 cm"
data <- data[which(Depth=="0-20 cm"), ]
data <- subset(data , select = - c(Depth, OC..20.30.cm., CaCO3..20.30.cm., BD.0.10, BD.0.20, BD.10.20, BD.20.30))

# convert true type variable
# data$Soil_Stones <- as.factor(data$Soil_Stones) -> okay, it should be considered in coarse_mass
data$OC <- as.numeric(data$OC)
data$CaCO3 <- as.numeric(data$CaCO3)
data$P <- as.numeric(data$P)
data$N <- as.numeric(data$N)
data$K <- as.numeric(data$K)
data$LC <- as.factor(data$LC)
data$LU <- as.factor(data$LU)


boxplot(data$Elev)
#clearly some not possible values -> remove
data <- data[which(data$Elev<=4000),]
data[which.max(data$Elev),]

#colSums(is.na(data))#

# keep just point for which we have bulk density
data <- data[which(!is.na(data$BDsample_0)),]

#set aside nonìt observed values but formula's evaluated
BD_Coarse_estimate <- subset(data , select = c(POINTID, BDfine_0_2, coarse_vol))
data <- subset(data , select = - c(BDfine_0_2, coarse_vol))

data$LC_0 <- factor(substr(as.character(data$LC), 1, 1))
data <- data[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,18,14,15,16,17)] #reorder

colSums(is.na(data))

# some NA but lets go on with it, consider later
write.csv(BD_Coarse_estimate, "../datasets/2018/extra/BD_Coarse_estimate.csv", row.names = FALSE)
write.csv(NUTS, "../datasets/2018/extra/NUTS.csv", row.names = FALSE)
rm(BD_Coarse_estimate,Coarse_Clay_Sand_Silt,NUTS)

unique(LC_legend$LC0_Desc)

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

# plot 
ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") + # Europe map
  geom_sf(data = data, color = 'red', size = 1, alpha = 0.7) + # points
  theme_minimal() +
  labs(title = "LUCAS Soil Points on Europe Map")

write.csv(data, "../datasets/2018/Lucas_workon_2018.csv", row.names = FALSE)
