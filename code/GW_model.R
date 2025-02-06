#load libraries
library(sf)
library(GWmodel)
library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(splines)
library(car)

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
                 adaptive=FALSE, p=2, theta=0, dMat=dist_mat, longlat=F)
fit4 <- gwr.basic(BDsample_0 ~ bs(Elev, degree=4) + bs(pH_H2O, degree=4) + bs(EC_2018, degree=4), data, coords, bw=bw_opt, kernel="gaussian",
                  adaptive=FALSE, p=2, theta=0, longlat=T,dMat=dist_mat)
ss4 <- sum((fit4$lm$residuals)^2)
aic4 <- AIC(fit4$lm)


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
    title = "Bulk density true",
    color = "Value Legend"  # Adjust legend title
  )
(plot_fitted | plot_values)

## Covariate selection
data <- data[!is.na(data$OC_2018), ]
data <- data[!is.na(data$N_2018), ]
data <- data[!is.na(data$K_2018), ]
data$Elev <- as.numeric(data$Elev)
pairs(data.frame(data$BDsample_0, data$Elev, data$pH_H2O, data$EC_2018, data$OC_2018, data$N_2018, data$K_2018))

# model for Elev 
data$Elev_reg <- bs(data$Elev, knots = quantile(data$Elev, c(0.5, 0.75, 0.95)), degree=3)
elev_mod <- lm(BDsample_0~bs(Elev, knots = quantile(data$Elev, c(0.5, 0.75, 0.95)), degree=3), data = data)
new_data <-
  with(data, data.frame(
    Elev = seq(range(Elev)[1], range(Elev)[2], by = 0.1)
  ))
preds_elev=predict(elev_mod, new_data, se=T)
se.bands=cbind(preds_elev$fit + 2*preds_elev$se.fit, preds_elev$fit-2*preds_elev$se.fit)
with(data, plot(Elev, BDsample_0, cex =.5, col ="black"))
lines(new_data$Elev, preds_elev$fit ,lwd =2, col ="blue")
matlines(new_data$Elev, se.bands ,lwd =1, col ="red",lty =3)

# model for EC_2018
data$EC_2018_reg <- bs(data$EC_2018, knots = quantile(data$EC_2018, c(0.5, 0.75, 0.95)))
EC_2018_mod <- lm(BDsample_0~bs(EC_2018, knots = quantile(EC_2018, c(0.5, 0.75, 0.95))), data = data)
new_data <-
  with(data, data.frame(
    EC_2018 = seq(range(EC_2018)[1], range(EC_2018)[2], by = 0.1)
  ))
preds_EC_2018=predict(EC_2018_mod, new_data, se=T)
se.bands=cbind(preds_EC_2018$fit + 2*preds_EC_2018$se.fit, preds_EC_2018$fit-2*preds_EC_2018$se.fit)
with(data, plot(EC_2018, BDsample_0, cex =.5, col ="black"))
lines(new_data$EC_2018, preds_EC_2018$fit ,lwd =2, col ="blue")
matlines(new_data$EC_2018, se.bands ,lwd =1, col ="red",lty =3)

# model for pH_H20
data$pH_H2O_reg <- bs(data$pH_H2O, degree=3)
pH_H2O_mod <- lm(BDsample_0~bs(pH_H2O, degree=3), data = data)
new_data <-
  with(data, data.frame(
    pH_H2O = seq(range(pH_H2O)[1], range(pH_H2O)[2], by = 0.1)
  ))
preds_pH_H2O=predict(pH_H2O_mod, new_data, se=T)
se.bands=cbind(preds_pH_H2O$fit + 2*preds_pH_H2O$se.fit, preds_pH_H2O$fit-2*preds_pH_H2O$se.fit)
with(data, plot(pH_H2O, BDsample_0, cex =.5, col ="black"))
lines(new_data$pH_H2O, preds_pH_H2O$fit ,lwd =2, col ="blue")
matlines(new_data$pH_H2O, se.bands ,lwd =1, col ="red",lty =3)

# model for OC_2018
data$OC_2018_reg <- bs(data$OC_2018, degree=3)
OC_2018_mod <- lm(BDsample_0~bs(OC_2018, , degree=3), data = data)
new_data <-
  with(data, data.frame(
    OC_2018 = seq(range(OC_2018)[1], range(OC_2018)[2], by = 0.1)
  ))
preds_OC_2018=predict(OC_2018_mod, new_data, se=T)
se.bands=cbind(preds_OC_2018$fit + 2*preds_OC_2018$se.fit, preds_OC_2018$fit-2*preds_OC_2018$se.fit)
with(data, plot(OC_2018, BDsample_0, cex =.5, col ="black"))
lines(new_data$OC_2018, preds_OC_2018$fit ,lwd =2, col ="blue")
matlines(new_data$OC_2018, se.bands ,lwd =1, col ="red",lty =3)

# model for N_2018
data$N_2018_reg <- bs(data$N_2018, degree=3)
N_2018_mod <- lm(BDsample_0~bs(N_2018, , degree=3), data = data)
new_data <-
  with(data, data.frame(
    N_2018 = seq(range(N_2018)[1], range(N_2018)[2], by = 0.1)
  ))
preds_N_2018=predict(N_2018_mod, new_data, se=T)
se.bands=cbind(preds_N_2018$fit + 2*preds_N_2018$se.fit, preds_N_2018$fit-2*preds_N_2018$se.fit)
with(data, plot(N_2018, BDsample_0, cex =.5, col ="black"))
lines(new_data$N_2018, preds_N_2018$fit ,lwd =2, col ="blue")
matlines(new_data$N_2018, se.bands ,lwd =1, col ="red",lty =3)

# model for K_2018
data$K_2018_reg <- bs(data$K_2018, knots = quantile(data$K_2018, c(0.5, 0.75, 0.95)), degree=4)
K_2018_mod <- lm(BDsample_0~bs(K_2018, knots = quantile(data$K_2018, c(0.5, 0.75, 0.95)), degree=4), data = data)
new_data <-
  with(data, data.frame(
    K_2018 = seq(range(K_2018)[1], range(K_2018)[2], by = 0.1)
  ))
preds_K_2018=predict(K_2018_mod, new_data, se=T)
se.bands=cbind(preds_K_2018$fit + 2*preds_K_2018$se.fit, preds_K_2018$fit-2*preds_K_2018$se.fit)
with(data, plot(K_2018, BDsample_0, cex =.5, col ="black"))
lines(new_data$K_2018, preds_K_2018$fit ,lwd =2, col ="blue")
matlines(new_data$K_2018, se.bands ,lwd =1, col ="red",lty =3)

coords <- st_coordinates(data)  
dist_mat <- gw.dist(coords, coords, focus=0, p=2, theta=0, longlat=T)
bw_opt <- bw.gwr(BDsample_0~OC_2018_reg+pH_H2O_reg+K_2018_reg+EC_2018_reg+Elev_reg, data, approach="CV", kernel="gaussian",
                 adaptive=FALSE, p=2, theta=0, dMat=dist_mat, longlat=F)

model_sel <- gwr.model.selection("BDsample_0",InDeVars=c("Elev_reg", "pH_H2O_reg", "EC_2018_reg", "OC_2018_reg", "N_2018_reg", "K_2018_reg"), data, bw=100,approach="CV",
                    adaptive=F,kernel="gaussian", dMat = dist_mat, p=2, theta=0, longlat=F,
                    parallel.method=F,parallel.arg=NULL)

fit8 <- gwr.basic(BDsample_0~OC_2018_reg+pH_H2O_reg, data, coords, bw=100, kernel="gaussian",
                  adaptive=FALSE, p=2, theta=0, longlat=F,dMat=dist_mat)
fit19 <- gwr.basic(BDsample_0~OC_2018_reg+pH_H2O_reg+K_2018_reg+EC_2018_reg+Elev_reg, data, bw=400, kernel="gaussian",
                   adaptive=FALSE, p=2, theta=0, longlat=F,dMat=dist_mat)

res8 <- fit8$lm$residuals
res19 <- fit19$SDF$residual
common_limits <- range(c(res8, res19), na.rm = TRUE)

plot_residuals8 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = res8), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits=common_limits) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Residuals model8",
    color = "Value Legend"  # Adjust legend title
  )


plot_residuals19 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = res19), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma") +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Residuals model19",
    color = "Value Legend"  # Adjust legend title
  )
(plot_residuals8 | plot_residuals19)


common_limits2 <- range(c(data$BDsample_0, fit8$lm$fitted.values), na.rm = TRUE)
plot_fit8 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = fit8$lm$fitted.values ), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits=common_limits2) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density estimate fit8",
    color = "Value Legend"  # Adjust legend title
  )
(plot_fit8|plot_values)

common_limits3 <- range(c(data$BDsample_0, fit19$SDF[,"yhat"]$yhat), na.rm = TRUE)
plot_fit19 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = fit19$SDF[,"yhat"]$yhat ), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits = common_limits3) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density estimate fit19",
    color = "Value Legend"  # Adjust legend title
  )
(plot_fit19|plot_values)

plot(res8)
plot(res19)
plot(res8, res19)
sum(res8^2)
sum(res19^2)

lm_mod <- lm(BDsample_0~Elev+pH_H2O+EC_2018+OC_2018+N_2018+K_2018, data = data)
sum(lm_mod$residuals^2)

lm_mod2 <- lm(BDsample_0~Elev_reg+pH_H2O_reg+EC_2018_reg+OC_2018_reg+N_2018_reg+K_2018_reg, data = data)
sum(lm_mod2$residuals^2)
