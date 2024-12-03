# load libraries
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(dplyr) 

##### 2015 #####

# read data 2015
lucas_soil_2015 <- read.csv("data/LUCAS2015_topsoildata_20200323/LUCAS_Topsoil_2015_20200323.csv")
head(lucas_soil_2015)

# read shape file 2015
lucas_2015_sf <- st_read("data/LUCAS2015_topsoildata_20200323/LUCAS_Topsoil_2015_20200323-shapefile/LUCAS_Topsoil_2015_20200323.shp")

colnames(lucas_2015_sf)
colnames(lucas_soil_2015)

# spot differences in selected attributes
ids <- c()
for (id in lucas_soil_2015$Point_ID){
  if(lucas_soil_2015[lucas_soil_2015$Point_ID==id,]$LC1_Desc != lucas_2015_sf[lucas_2015_sf$Point_ID==id,]$LC1_Desc){
    ids <- c(ids, id)
  }
}

lucas_soil_2015[lucas_soil_2015$Point_ID==ids[1],]$LC1_Desc
lucas_2015_sf[lucas_2015_sf$Point_ID==ids[1],]$LC1_Desc


# merge them on Point_Id
# lucas_2015 <- merge(lucas_2015_sf, lucas_soil_2015, by = "Point_ID") # no need all feature (almost) already present
lucas_2015 <- lucas_2015_sf

# add missing covariates in the shapefile
lucas_2015$Revisited_point <- lucas_soil_2015$Revisited_point
lucas_2015$Elevation <- lucas_soil_2015$Elevation
lucas_2015$Soil_Stones <- lucas_soil_2015$Soil_Stones

colnames(lucas_2015)[1] <- "POINTID"

# remove useless
rm(id, ids, lucas_soil_2015, lucas_2015_sf)

##### 2018 #####

# read data 2018
lucas_soil_2018 <- read.csv("downloaded/LUCAS-SOIL-2018-data-report-readme-v2/LUCAS-SOIL-2018-v2/LUCAS-SOIL-2018.csv")
head(lucas_soil_2018)

# read shape file 2018
lucas_2018_sf <- st_read("downloaded/LUCAS-SOIL-2018-data-report-readme-v2/LUCAS-SOIL-2018-v2/LUCAS-SOIL-2018 .shp")

# merge them on POINTID
lucas_2018 <- merge(lucas_2018_sf, lucas_soil_2018, by = "POINTID")

# remove useless
rm(lucas_2018_sf, lucas_soil_2018)

#### BULK DENSITY ####

bulk_density_full <- read.csv("downloaded/LUCAS-SOIL-2018-data-report-readme-v2/LUCAS-SOIL-2018-v2/BulkDensity_2018_final-2.csv")
bulk_density_select <- read.csv("downloaded/BD_LUCAS_data_for_paper/BD2018_data_for_paper.csv")
bulk_density_select <- bulk_density_select[-1,-c(7,8)]
bulk_density_select$BDsample_0 <- as.numeric(bulk_density_select$BDsample_0)
bulk_density_select$coarse_mas <- as.numeric(bulk_density_select$coarse_mas)
bulk_density_select$BDfine_0_2 <- as.numeric(bulk_density_select$BDfine_0_2)
bulk_density_select$coarse_vol <- as.numeric(bulk_density_select$coarse_vol)

colnames(bulk_density_full)
colnames(bulk_density_select) # contains bulk density of 0-20cm

# read shape file
bd_shape <- st_read("downloaded/BD_LUCAS_data_for_paper/BD2018_data_for_paper_shapefile/BD2018_fef_landcover.shp") #all covariates present
bulk_density_select <- bd_shape

# read shape file 2018
lucas_2018_sf <- st_read("downloaded/LUCAS-SOIL-2018-data-report-readme-v2/LUCAS-SOIL-2018-v2/LUCAS-SOIL-2018 .shp")
colnames(lucas_2018_sf)[1] <- "POINT_ID"
bulk_density_full <- merge(lucas_2018_sf, bulk_density_full, by = "POINT_ID")

# remove useless
rm(bd_shape, lucas_2018_sf)

##### 2015 + 2018 #####

# check point coherence by ID
common_ids <- intersect(lucas_2018$POINTID, lucas_2015$POINTID) # IDs common to both datasets
common_points <- st_join(lucas_2015, lucas_2018, left=F)

#why does not correspond?
improper_ids <- common_ids[!(common_ids %in% common_points$POINTID.x)]
lucas_2018[lucas_2018$POINTID==improper_ids[1],]$geometry
lucas_2015[lucas_2015$POINTID==improper_ids[1],]$geometry

# slow - DO NOT REPEAT
tot <- 0
for (id in improper_ids){
  tot <- tot + as.numeric(st_distance(
    lucas_2015[lucas_2015$POINTID == improper_ids[1], ],
    lucas_2018[lucas_2018$POINTID == improper_ids[1], ]
  )[1,1])
}
tot/length(improper_ids) # mean of 0.002 meters -> identical

common_points <- st_join(lucas_2015, lucas_2018, left=F, join = st_is_within_distance, dist = 0.01, suffix = c("_2015", "_2018"))

# Filter points that are common small distance problem - > keep 2018 geometry
lucas_2018_common <- lucas_2018[lucas_2018$POINTID %in% common_ids, ]
lucas_2018_common <- lucas_2018_common[order(lucas_2018_common$POINTID),]
lucas_2015_common <- lucas_2015[lucas_2015$POINTID %in% common_ids, ]
lucas_2015_common <- lucas_2015_common[order(lucas_2015_common$POINTID),]

lucas_2015_common$geometry <- lucas_2018_common$geometry

lucas_2015_2018 <- st_join(lucas_2015_common, lucas_2018_common, left=F, suffix = c("_2015", "_2018"))

# remove unused
rm(lucas_2015_common, lucas_2018_common, common_points, common_points_eq, tot, id, common_ids, improper_ids)

##### 2015 + 2018 + BD ####

lucas_2015_bd <- st_join(lucas_2015, bulk_density_full,left = T, join = st_is_within_distance, dist = 0.01)
lucas_2015_bd <- st_join(lucas_2015_bd, bulk_density_select,left = T, join = st_is_within_distance, dist = 0.01)
colnames(lucas_2015_bd)
missing_points_15 <- lucas_2015_bd %>% filter(is.na(POINT_ID.x))

lucas_2018_bd <- st_join(lucas_2018, bulk_density_full,left = T)
lucas_2018_bd <- st_join(lucas_2018_bd, bulk_density_select,left = T)
colnames(lucas_2018_bd)
missing_points_18 <- lucas_2018_bd %>% filter(is.na(POINT_ID.x))

lucas_2015_comm_bd <- st_join(lucas_2015_common, bulk_density_full,left = T)
lucas_2015_comm_bd <- st_join(lucas_2015_comm_bd, bulk_density_select,left = T)
colnames(lucas_2015_comm_bd)
missing_points_15_comm <- lucas_2015_comm_bd %>% filter(is.na(POINT_ID.x))

lucas_2018_comm_bd <- st_join(lucas_2018_common, bulk_density_full,left = T)
lucas_2018_comm_bd <- st_join(lucas_2018_comm_bd, bulk_density_select,left = T)
colnames(lucas_2018_comm_bd)
missing_points_18_comm <- lucas_2018_comm_bd %>% filter(is.na(POINT_ID.x))

# 2015 row
nrow(lucas_2015)
# 2015 common row
nrow(lucas_2015_common)
# 2015 common bd row
nrow(missing_points_15_comm)
# 2015 bd row
nrow(missing_points_15)

# 2018 row
nrow(lucas_2018)
# 2018 common row
nrow(lucas_2018_common)
# 2018 common bd row
nrow(missing_points_18_comm)
# 2018 bd row
nrow(missing_points_18)

# seems all good, remove unused
rm(missing_points_15, missing_points_15_comm, missing_points_18_comm, missing_points_18, lucas_2015_comm_bd, lucas_2018_comm_bd)

# build FULL join
lucas_2015_2018_bd <- st_join(lucas_2015_2018, bulk_density_full, left = T, join = st_is_within_distance, dist = 0.01)
lucas_2015_2018_bd <- st_join(lucas_2015_2018_bd, bulk_density_select,left = T, join = st_is_within_distance, dist = 0.01)

rm(bulk_density_full, bulk_density_select)

#### CLEAN AND EXPORT ####

# 2015
colnames(lucas_2015)
w_lucas_2015 <- lucas_2015[,c(1,24,2,3,4,5,26,25,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)]

lc_legend <- w_lucas_2015 %>%
  st_drop_geometry() %>% 
  select(LC, LC0_Desc, LC1_Desc) %>% 
  distinct() %>%
  arrange(LC)
  
write.csv(lc_legend, "datasets/2015/LC_description.csv", row.names = FALSE)

lu_legend <- w_lucas_2015 %>%
  st_drop_geometry() %>%  
  select(LU, LU1_Desc) %>%  
  distinct() %>%
  arrange(LU)

write.csv(lu_legend, "datasets/2015/LU_description.csv", row.names = FALSE)

w_lucas_2015$LC0_Desc <- NULL
w_lucas_2015$LC1_Desc <- NULL
w_lucas_2015$LU1_Desc <- NULL

st_write(w_lucas_2015, "datasets/2015/shapefile/LUCAS_2015.shp", delete_layer = TRUE)
w_lucas_2015$geometry <- NULL
write.csv(w_lucas_2015, "datasets/2015/LUCAS_2015.csv", row.names = FALSE)
rm(w_lucas_2015)

# 2018
colnames(lucas_2018)
w_lucas_2018 <- lucas_2018[,c(1,21,2,3,4,5,6,11,7,12,8,9,10,13,14,22,15,16,17,18,19,20,23,24,25,26,27,28)]

lc_legend <- w_lucas_2018 %>%
  st_drop_geometry() %>% 
  select(LC, LC0_Desc, LC1_Desc) %>% 
  distinct() %>%
  arrange(LC)

write.csv(lc_legend, "datasets/2018/LC_description.csv", row.names = FALSE)

lu_legend <- w_lucas_2018 %>%
  st_drop_geometry() %>%  
  select(LU, LU1_Desc) %>%  
  distinct() %>%
  arrange(LU)

write.csv(lu_legend, "datasets/2018/LU_description.csv", row.names = FALSE)

w_lucas_2018$LC0_Desc <- NULL
w_lucas_2018$LC1_Desc <- NULL
w_lucas_2018$LU1_Desc <- NULL

st_write(w_lucas_2018, "datasets/2018/shapefile/LUCAS_2018.shp", delete_layer = TRUE)
w_lucas_2018$geometry <- NULL
write.csv(w_lucas_2018, "datasets/2018/LUCAS_2018.csv", row.names = FALSE)
rm(w_lucas_2018)

# 2015 BD
colnames(lucas_2015_bd)
w_lucas_2015_bd <- lucas_2015_bd[,c(1,24,2,3,4,5,25,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,27,28,29,30,33,35,36,37)]

w_lucas_2015_bd$LC0_Desc <- NULL
w_lucas_2015_bd$LC1_Desc <- NULL
w_lucas_2015_bd$LU1_Desc <- NULL

st_write(w_lucas_2015_bd, "datasets/2015/shapefile/LUCAS_2015_bd.shp", delete_layer = TRUE)
w_lucas_2015_bd$geometry <- NULL
write.csv(w_lucas_2015_bd, "datasets/2015/LUCAS_2015_bd.csv", row.names = FALSE)
rm(w_lucas_2015_bd)

# 2018 BD
colnames(lucas_2018_bd)
w_lucas_2018_bd <- lucas_2018_bd[,c(1,21,2,3,4,5,6,11,7,12,8,9,10,13,14,22,15,16,17,18,19,20,23,24,25,26,27,29,30,31,32,35,36,37,38,39)]

w_lucas_2018_bd$LC0_Desc.x <- NULL
w_lucas_2018_bd$LC1_Desc <- NULL
w_lucas_2018_bd$LU1_Desc <- NULL

st_write(w_lucas_2018_bd, "datasets/2018/shapefile/LUCAS_2018_bd.shp", delete_layer = TRUE)
w_lucas_2018_bd$geometry <- NULL
write.csv(w_lucas_2018_bd, "datasets/2018/LUCAS_2018_bd.csv", row.names = FALSE)
rm(w_lucas_2018_bd)

# 2015 2018
colnames(lucas_2015_2018)
w_lucas_2015_2018 <- lucas_2015_2018[,c(1,26,46,23,24,2,3,4,5,6,28,7,29,8,30,9,31,36,10,32,37,11,33,12,34,13,35,14,48,15,49,16,40,17,41,18,42,19,43,25,27,38,39,47,44,45,20,50,21,51,22,52,53)]

length(unique(w_lucas_2015_2018$LU_2015))
length(unique(w_lucas_2015_2018$LU_2018))

lc_legend_15 <- w_lucas_2015_2018 %>%
  st_drop_geometry() %>% 
  select(LC_2015, LC0_Desc_2015, LC1_Desc_2015) %>% 
  distinct() %>%
  arrange(LC_2015)

write.csv(lc_legend_15, "datasets/2015_2018/LC_description_15.csv", row.names = FALSE)

lu_legend_15 <- w_lucas_2015_2018 %>%
  st_drop_geometry() %>%  
  select(LU_2015, LU1_Desc_2015) %>%  
  distinct() %>%
  arrange(LU_2015)

write.csv(lu_legend_15, "datasets/2015_2018/LU_description_15.csv", row.names = FALSE)

lc_legend_18 <- w_lucas_2015_2018 %>%
  st_drop_geometry() %>% 
  select(LC_2018, LC0_Desc_2018, LC1_Desc_2018) %>% 
  distinct() %>%
  arrange(LC_2018)

write.csv(lc_legend_18, "datasets/2015_2018/LC_description_18.csv", row.names = FALSE)

lu_legend_18 <- w_lucas_2015_2018 %>%
  st_drop_geometry() %>%  
  select(LU_2018, LU1_Desc_2018) %>%  
  distinct() %>%
  arrange(LU_2018)

write.csv(lu_legend_18, "datasets/2015_2018/LU_description_18.csv", row.names = FALSE)

w_lucas_2015_2018$LC0_Desc_2015 <- NULL
w_lucas_2015_2018$LC1_Desc_2015 <- NULL
w_lucas_2015_2018$LU1_Desc_2015 <- NULL
w_lucas_2015_2018$LC0_Desc_2018 <- NULL
w_lucas_2015_2018$LC1_Desc_2018 <- NULL
w_lucas_2015_2018$LU1_Desc_2018 <- NULL

w_shp_lucas_2015_2018 <- w_lucas_2015_2018[,c(1,2,47)]
colnames(w_shp_lucas_2015_2018)[1] <- "ID_2015"
colnames(w_shp_lucas_2015_2018)[2] <- "ID_2018"
st_write(w_shp_lucas_2015_2018, "datasets/2015_2018/shapefile/LUCAS_2015_2018.shp", delete_layer = TRUE)
w_lucas_2015_2018$geometry <- NULL
write.csv(w_lucas_2015_2018, "datasets/2015_2018/LUCAS_2015_2018.csv", row.names = FALSE)
rm(w_lucas_2015_2018)

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