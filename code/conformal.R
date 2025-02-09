#load libraries
library(sf)
library(GWmodel)
library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(splines)
library(car)
library(dplyr)
rm(list=ls())

# load data
covariate <- read.csv("../datasets/LUCAS_workon.csv")
shapefile <- st_read("../datasets/shapefile/LUCAS_workon_shape.shp")
LC_legend <- read.csv("../datasets/LC_description.csv")
LU_legend <- read.csv("../datasets/LU_description.csv")
data_full <- merge(shapefile, covariate, by = "POINTID")
rm(covariate, shapefile)

data_full$Soil_Stones <- as.factor(data_full$Soil_Stones)
data_full$LC_0 <- as.factor(data_full$LC_0)
data_full$LC_2018 <- as.factor(data_full$LC_2018)
data_full$LU_2018 <- as.factor(data_full$LU_2018)

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
rm(europe, europe_with_cyprus, cyprus, world, bounding_box, lat_range, lon_range)


# keep only covariate of intrest
data <- data_full %>% select("Elev","OC_2018","N_2018","K_2018","EC_2018","pH_H2O","BDsample_0","LC_0")

data <- data[!is.na(data$OC_2018), ]
data <- data[!is.na(data$N_2018), ]
data <- data[!is.na(data$K_2018), ]
data <- data[!is.na(data$EC_2018), ]
data$Elev <- as.numeric(data$Elev)


# add value for spline regression
data$Elev_reg <- bs(data$Elev, knots = quantile(data$Elev, c(0.5, 0.75, 0.95)), degree=3)
data$EC_2018_reg <- bs(data$EC_2018, knots = quantile(data$EC_2018, c(0.5, 0.75, 0.95)))
data$pH_H2O_reg <- bs(data$pH_H2O, degree=3)
data$OC_2018_reg <- bs(data$OC_2018, degree=3)
data$N_2018_reg <- bs(data$N_2018, degree=3)
data$K_2018_reg <- bs(log(data$K_2018), degree=2)

# function to compute the weights
compute_gaussian_weights <- function(tags, bw) {
  
  #number of tags to consider
  n <- nrow(tags)
  
  #matrix of distances according to tags(coordinates)
  dist_mat <- gw.dist(tags, focus=0, p=2, theta=0, longlat=T)
  
  # Compute Gaussian weights with given bandwidth
  weights <- exp(-0.5 * (dist_mat[-n, n]/bw)^2)
  #weights <- 1/(dist_mat[-n, n])
  
  # Normalize with n+1 weight = 1
  weights <- weights / (sum(weights)+1)
  
  return(c(weights, 1))
}

#prediction names for predict to work (assure same names)
set_pred_names <- function(data_fun){
  
  data_fun$OC_2018_reg1 <- data_fun$OC_2018_reg[,1]
  data_fun$OC_2018_reg2 <- data_fun$OC_2018_reg[,2]
  data_fun$OC_2018_reg3 <- data_fun$OC_2018_reg[,3]
  data_fun$Elev_reg1 <- data_fun$Elev_reg[,1]
  data_fun$Elev_reg2 <- data_fun$Elev_reg[,2]
  data_fun$Elev_reg3 <- data_fun$Elev_reg[,3]
  data_fun$Elev_reg4 <- data_fun$Elev_reg[,4]
  data_fun$Elev_reg5 <- data_fun$Elev_reg[,5]
  data_fun$Elev_reg6 <- data_fun$Elev_reg[,6]
  data_fun$LC_0B <- data_fun$LC_0=="B"
  data_fun$LC_0C <- data_fun$LC_0=="C"
  data_fun$LC_0D <- data_fun$LC_0=="D"
  data_fun$LC_0E <- data_fun$LC_0=="E"
  data_fun$LC_0F <- data_fun$LC_0=="F"
  data_fun$LC_0G <- data_fun$LC_0=="G"
  data_fun$LC_0H <- data_fun$LC_0=="H"
  data_fun$pH_H2O_reg1 <- data_fun$pH_H2O_reg[,1]
  data_fun$pH_H2O_reg2 <- data_fun$pH_H2O_reg[,2]
  data_fun$pH_H2O_reg3 <- data_fun$pH_H2O_reg[,3]
  data_fun$K_2018_reg1 <- data_fun$K_2018_reg[,1]
  data_fun$K_2018_reg2 <- data_fun$K_2018_reg[,2]
  data_fun$N_2018_reg1 <- data_fun$N_2018_reg[,1]
  data_fun$N_2018_reg2 <- data_fun$N_2018_reg[,2]
  data_fun$N_2018_reg3 <- data_fun$N_2018_reg[,3]
  data_fun$EC_2018_reg1 <- data_fun$EC_2018_reg[,1]
  data_fun$EC_2018_reg2 <- data_fun$EC_2018_reg[,2]
  data_fun$EC_2018_reg3 <- data_fun$EC_2018_reg[,3]
  data_fun$EC_2018_reg4 <- data_fun$EC_2018_reg[,4]
  data_fun$EC_2018_reg5 <- data_fun$EC_2018_reg[,5]
  data_fun$EC_2018_reg6 <- data_fun$EC_2018_reg[,6]
  
  return (data_fun)
}

# function that performs non exchangeable split conformal prediction
nexSCP <- function(data, data_predict, alpha = 0.1, calib_percentage = 0.5, bw = 300) {
  
  n <- nrow(data)
  tags <- st_coordinates(data)
  dmat <- gw.dist(tags, focus=0, p=2, theta=0, longlat=T)
  
  calib_indeces <- sample(1:n, size = round(n * calib_percentage))
  
  data_train <- data[-calib_indeces,]
  data_calib <- set_pred_names(data[calib_indeces,])

  dist_mat_train <- gw.dist(st_coordinates(data_train), focus=0, p=2, theta=0, longlat=T)
  dist_mat_train_calib <- gw.dist(st_coordinates(data_train), st_coordinates(data_calib), focus=0, p=2, theta=0, longlat=T)
  
  fit_calib <- gwr.predict(BDsample_0~OC_2018_reg+LC_0+Elev_reg+K_2018_reg+N_2018_reg+pH_H2O_reg+EC_2018_reg, data_train, data_calib, bw=bw, kernel="gaussian",adaptive=FALSE, p=2,
                           theta=0, longlat=F, dMat1=t(dist_mat_train_calib), dMat2=dist_mat_train)
  
  data_predict <- set_pred_names(data_predict)
  
  dist_mat_train_predict <- gw.dist(st_coordinates(data_train), st_coordinates(data_predict), focus=0, p=2, theta=0, longlat=T)
  
  fit_predict <- gwr.predict(BDsample_0~OC_2018_reg+LC_0+Elev_reg+K_2018_reg+N_2018_reg+pH_H2O_reg+EC_2018_reg, data_train, data_predict, bw=bw, kernel="gaussian",adaptive=FALSE, p=2,
                             theta=0, longlat=F, dMat1=t(dist_mat_train_predict), dMat2=dist_mat_train)
  calib_thresh <- c()
  for(i in 1:nrow(data_predict)){
    weights_calib <- compute_gaussian_weights(rbind(tags[calib_indeces,],st_coordinates(data_predict[i,])),bw=1000)
  
    if (sum(weights_calib) >= 1 - alpha) {
      R <- abs(fit_calib$SDF$prediction-data_calib$BDsample_0)
      ord_R <- order(R)
      ind_thresh <- min(which(cumsum(weights_calib[ord_R]) >= 1 - alpha))
      calib_thresh[i] <- sort(R)[ind_thresh]
    } else {
      calib_thresh[i] <- Inf
    }
    
  }
  
  y_PI_df <- data.frame(
    lower = fit_predict$SDF$prediction - calib_thresh,
    prediction = fit_predict$SDF$prediction,
    upper = fit_predict$SDF$prediction + calib_thresh
  )
  
  return(y_PI_df)
}

#import data to predict
data_predict <- st_read("../datasets/2018/shapefile/LUCAS_2018.shp")

#select consistent depth with model
data_predict <- data_predict[which(data_predict$Depth=="0-20 cm"), ]

#keep only not observed point
data_predict <- data_predict %>%
  anti_join(data.frame(POINTID = intersect(data_full$POINTID, data_predict$POINTID)), 
            by = "POINTID")

# keep only covariate of intrest
data_predict <- data_predict %>% select("Elev","OC","N","K","EC","pH_H2O","LC")

data_predict$LC_0 <- factor(substr(as.character(data_predict$LC), 1, 1))
data_predict$OC_2018 <- as.numeric(data_predict$OC)
data_predict$EC_2018 <- as.numeric(data_predict$EC)
data_predict$N_2018 <- as.numeric(data_predict$N)
data_predict$K_2018 <- as.numeric(data_predict$K)

#mtch names
data_predict <- data_predict %>% select("Elev","OC_2018","N_2018","K_2018","EC_2018","pH_H2O","LC_0")


data_predict <- data_predict[!is.na(data_predict$OC_2018), ]
data_predict <- data_predict[!is.na(data_predict$N_2018), ]
data_predict <- data_predict[!is.na(data_predict$K_2018), ]
data_predict <- data_predict[!is.na(data_predict$EC_2018), ]
data_predict$Elev <- as.numeric(data_predict$Elev)

# add value for spline regression
data_predict$Elev_reg <- bs(data_predict$Elev, knots = quantile(data_predict$Elev, c(0.5, 0.75, 0.95)), degree=3)
data_predict$EC_2018_reg <- bs(data_predict$EC_2018, knots = quantile(data_predict$EC_2018, c(0.5, 0.75, 0.95)))
data_predict$pH_H2O_reg <- bs(data_predict$pH_H2O, degree=3)
data_predict$OC_2018_reg <- bs(data_predict$OC_2018, degree=3)
data_predict$N_2018_reg <- bs(data_predict$N_2018, degree=3)
data_predict$K_2018_reg <- bs(log(data_predict$K_2018), degree=2)

result <- nexSCP(data = data,data_predict = data_predict)

n <- nrow(data)
calib_percentage <- 0.5
alpha <- 0.1
bw  <- 300
tags <- st_coordinates(data)
dmat <- gw.dist(tags, focus=0, p=2, theta=0, longlat=T)

calib_indeces <- sample(1:n, size = round(n * calib_percentage))

data_train <- data[-calib_indeces,]
data_calib <- set_pred_names(data[calib_indeces,])

dist_mat_train <- gw.dist(st_coordinates(data_train), focus=0, p=2, theta=0, longlat=T)
dist_mat_train_calib <- gw.dist(st_coordinates(data_train), st_coordinates(data_calib), focus=0, p=2, theta=0, longlat=T)

fit_calib <- gwr.predict(BDsample_0~OC_2018_reg+LC_0+Elev_reg+K_2018_reg+N_2018_reg+pH_H2O_reg+EC_2018_reg, data_train, data_calib, bw=bw, kernel="gaussian",adaptive=FALSE, p=2,
                         theta=0, longlat=F, dMat1=t(dist_mat_train_calib), dMat2=dist_mat_train)

data_predict <- set_pred_names(data_predict)

dist_mat_train_predict <- gw.dist(st_coordinates(data_train), st_coordinates(data_predict), focus=0, p=2, theta=0, longlat=T)

fit_predict <- gwr.predict(BDsample_0~OC_2018_reg+LC_0+Elev_reg+K_2018_reg+N_2018_reg+pH_H2O_reg+EC_2018_reg, data_train, data_predict, bw=bw, kernel="gaussian",adaptive=FALSE, p=2,
                           theta=0, longlat=F, dMat1=t(dist_mat_train_predict), dMat2=dist_mat_train)
calib_thresh <- c()
for(i in 1:nrow(data_predict)){
  weights_calib <- compute_gaussian_weights(rbind(tags[calib_indeces,],st_coordinates(data_predict[i,])),bw=300)
  
  if (sum(weights_calib) >= 1 - alpha) {
    R <- abs(fit_calib$SDF$prediction-data_calib$BDsample_0)
    ord_R <- order(R)
    ind_thresh <- min(which(cumsum(weights_calib[ord_R]) >= 1 - alpha))
    calib_thresh[i] <- sort(R)[ind_thresh]
  } else {
    calib_thresh[i] <- Inf
  }
  
}

y_PI_df <- data.frame(
  lower = fit_predict$SDF$prediction - calib_thresh,
  prediction = fit_predict$SDF$prediction,
  upper = fit_predict$SDF$prediction + calib_thresh
)





