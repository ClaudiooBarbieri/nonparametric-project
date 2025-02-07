#load libraries
library(sf)
library(GWmodel)
library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(splines)
library(car)
library(progress)
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
pH_H2O_mod <- lm(log(BDsample_0)~bs(pH_H2O, degree=3), data = data)
new_data <-
  with(data, data.frame(
    pH_H2O = seq(range(pH_H2O)[1], range(pH_H2O)[2], by = 0.1)
  ))
preds_pH_H2O=predict(pH_H2O_mod, new_data, se=T)
se.bands=cbind(preds_pH_H2O$fit + 2*preds_pH_H2O$se.fit, preds_pH_H2O$fit-2*preds_pH_H2O$se.fit)
with(data, plot(pH_H2O, log(BDsample_0), cex =.5, col ="black"))
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
data$K_2018_reg <- bs(log(data$K_2018), degree=2)
K_2018_mod <- lm(BDsample_0~bs(log(K_2018), degree=2), data = data)
new_data <-
  with(data, data.frame(
    K_2018 = seq(range(log(K_2018))[1], range(log(K_2018))[2], by = 0.1)
  ))
preds_K_2018=predict(K_2018_mod, new_data, se=T)
se.bands=cbind(preds_K_2018$fit + 2*preds_K_2018$se.fit, preds_K_2018$fit-2*preds_K_2018$se.fit)
with(data, plot(log(K_2018), BDsample_0, cex =.5, col ="black"))
lines(new_data$K_2018, preds_K_2018$fit ,lwd =2, col ="blue")
matlines(new_data$K_2018, se.bands ,lwd =1, col ="red",lty =3)

coords <- st_coordinates(data)  
dist_mat <- gw.dist(coords, coords, focus=0, p=2, theta=0, longlat=T)
bw_opt <- bw.gwr(BDsample_0~OC_2018_reg+pH_H2O_reg+K_2018_reg+EC_2018_reg+Elev_reg+N_2018_reg, data, approach="CV", kernel="gaussian",
                 adaptive=FALSE, p=2, theta=0, dMat=dist_mat, longlat=F)

model_sel <- gwr.model.selection("BDsample_0",InDeVars=c("Elev_reg", "pH_H2O_reg", "EC_2018_reg", "OC_2018_reg", "N_2018_reg", "K_2018_reg", "LC_0"), data, bw=bw_opt,approach="CV",
                    adaptive=F,kernel="gaussian", dMat = dist_mat, p=2, theta=0, longlat=F,
                    parallel.method=F,parallel.arg=NULL)

# model with no covariates
fit0 <- gwr.basic(BDsample_0~1, data, bw=bw_opt, kernel="gaussian",
                  adaptive=FALSE, p=2, theta=0, longlat=F,dMat=dist_mat)
#simple but well performing model
fit1 <- gwr.basic(BDsample_0~OC_2018_reg+LC_0+Elev_reg, data, bw=bw_opt, kernel="gaussian",
                  adaptive=FALSE, p=2, theta=0, longlat=F,dMat=dist_mat)
#best model
fit2 <- gwr.basic(BDsample_0~OC_2018_reg+LC_0+Elev_reg+K_2018_reg+N_2018_reg+pH_H2O_reg+EC_2018_reg, data, bw=bw_opt, kernel="gaussian",
                   adaptive=FALSE, p=2, theta=0, longlat=F,dMat=dist_mat)

#basic lm
lm_mod <- lm(BDsample_0~Elev+pH_H2O+EC_2018+OC_2018+N_2018+K_2018+LC_0, data = data)

#GAM
lm_mod2 <- lm(BDsample_0~Elev_reg+pH_H2O_reg+EC_2018_reg+OC_2018_reg+N_2018_reg+K_2018_reg+LC_0, data = data)

semi_mod <- gwr.mixed(BDsample_0~Elev_reg+pH_H2O_reg, data, coords, fixed.vars=c("LC_0"),
          intercept.fixed=FALSE, bw_opt, diagnostic=T, kernel="gaussian",
          adaptive=FALSE, p=2, theta=0, longlat=F,dMat=dist_mat, dMat.rp=dist_mat)

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

res0 <- fit0$SDF$residual
res1 <- fit1$SDF$residual
res2 <- fit2$SDF$residual
common_limits <- range(c(res0, res1, res2), na.rm = TRUE)

plot_residuals0 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = res0), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits=common_limits) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Residuals model0",
    color = "Value Legend"  # Adjust legend title
  )


plot_residuals1 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = res1), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits=common_limits) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Residuals model1",
    color = "Value Legend"  # Adjust legend title
  )

plot_residuals2 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = res2), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits=common_limits) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Residuals model2",
    color = "Value Legend"  # Adjust legend title
  )

(plot_residuals0 | plot_residuals1 | plot_residuals2)


common_limits2 <- range(c(data$BDsample_0, fit0$SDF$yhat), na.rm = TRUE)
plot_fit0 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = fit0$SDF$yhat), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits=common_limits2) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density estimate fit0",
    color = "Value Legend"  # Adjust legend title
  )
(plot_fit0|plot_values)

common_limits3 <- range(c(data$BDsample_0, fit1$SDF$yhat), na.rm = TRUE)
plot_fit1 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = fit1$SDF$yhat ), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits = common_limits3) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density estimate fit1",
    color = "Value Legend"  # Adjust legend title
  )
(plot_fit1|plot_values)

common_limits4 <- range(c(data$BDsample_0, fit2$SDF$yhat), na.rm = TRUE)
plot_fit2 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = fit2$SDF$yhat ), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits = common_limits4) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density estimate fit2",
    color = "Value Legend"  # Adjust legend title
  )
(plot_fit2|plot_values)

common_limits5 <- range(c(data$BDsample_0, lm_mod$fitted.values), na.rm = TRUE)
plot_fit_lm1 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = lm_mod$fitted.values ), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits = common_limits5) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density estimate lm1",
    color = "Value Legend"  # Adjust legend title
  )
(plot_fit_lm1|plot_values)

common_limits6 <- range(c(data$BDsample_0, lm_mod2$fitted.values), na.rm = TRUE)
plot_fit_lm2 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data, aes(color = lm_mod2$fitted.values ), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits = common_limits6) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density estimate lm2",
    color = "Value Legend"  # Adjust legend title
  )
(plot_fit_lm2|plot_values)


#Model comparison
#RSS
sum(lm_mod$residuals^2)
sum(lm_mod2$residuals^2)
fit0$GW.diagnostic$RSS.gw
fit1$GW.diagnostic$RSS.gw
fit2$GW.diagnostic$RSS.gw

#AIC
AIC(lm_mod)
AIC(lm_mod2)
fit0$GW.diagnostic$AIC
fit1$GW.diagnostic$AIC
fit2$GW.diagnostic$AIC

#R2_adj
summary(lm_mod)$adj.r.squared
summary(lm_mod2)$adj.r.squared
fit0$GW.diagnostic$gwR2.adj
fit1$GW.diagnostic$gwR2.adj
fit2$GW.diagnostic$gwR2.adj

## Validation
set.seed(811)  # For reproducibility

# Define percentage for the first subset (e.g., 70%)
split_percentage <- 0.75  
n <- nrow(data)  
sample_indices <- sample(1:n, size = round(n * split_percentage))  

# Create two subsets
data_tr <- data[sample_indices, ]  # 75% of the data
data_val <- data[-sample_indices, ]  # Remaining 25%

coords_tr <- st_coordinates(data_tr)  
dist_mat_tr <- gw.dist(coords_tr, coords_tr, focus=0, p=2, theta=0, longlat=T)
coords_val <- st_coordinates(data_val)  
dist_mat_tr_val <- gw.dist(coords_tr, coords_val, focus=0, p=2, theta=0, longlat=T)

data_val$OC_2018_reg1 <- data_val$OC_2018_reg[,1]
data_val$OC_2018_reg2 <- data_val$OC_2018_reg[,2]
data_val$OC_2018_reg3 <- data_val$OC_2018_reg[,3]
data_val$Elev_reg1 <- data_val$Elev_reg[,1]
data_val$Elev_reg2 <- data_val$Elev_reg[,2]
data_val$Elev_reg3 <- data_val$Elev_reg[,3]
data_val$Elev_reg4 <- data_val$Elev_reg[,4]
data_val$Elev_reg5 <- data_val$Elev_reg[,5]
data_val$Elev_reg6 <- data_val$Elev_reg[,6]
data_val$LC_0B <- data_val$LC_0=="B"
data_val$LC_0C <- data_val$LC_0=="C"
data_val$LC_0D <- data_val$LC_0=="D"
data_val$LC_0E <- data_val$LC_0=="E"
data_val$LC_0F <- data_val$LC_0=="F"
data_val$LC_0G <- data_val$LC_0=="G"
data_val$LC_0H <- data_val$LC_0=="H"
data_val$pH_H2O_reg1 <- data_val$pH_H2O_reg[,1]
data_val$pH_H2O_reg2 <- data_val$pH_H2O_reg[,2]
data_val$pH_H2O_reg3 <- data_val$pH_H2O_reg[,3]
data_val$K_2018_reg1 <- data_val$K_2018_reg[,1]
data_val$K_2018_reg2 <- data_val$K_2018_reg[,2]
data_val$N_2018_reg1 <- data_val$N_2018_reg[,1]
data_val$N_2018_reg2 <- data_val$N_2018_reg[,2]
data_val$N_2018_reg3 <- data_val$N_2018_reg[,3]
data_val$EC_2018_reg1 <- data_val$EC_2018_reg[,1]
data_val$EC_2018_reg2 <- data_val$EC_2018_reg[,2]
data_val$EC_2018_reg3 <- data_val$EC_2018_reg[,3]
data_val$EC_2018_reg4 <- data_val$EC_2018_reg[,4]
data_val$EC_2018_reg5 <- data_val$EC_2018_reg[,5]
data_val$EC_2018_reg6 <- data_val$EC_2018_reg[,6]
# model with no covariates
fit0_pred <- gwr.predict(BDsample_0~1, data_tr, data_val, bw=bw_opt, kernel="gaussian",adaptive=FALSE, p=2,
            theta=0, longlat=F,dMat1=t(dist_mat_tr_val), dMat2=dist_mat_tr)
#simple but well performing model
fit1_pred <- gwr.predict(fit1$GW.arguments$formula, data_tr, data_val, bw=bw_opt, kernel="gaussian",adaptive=FALSE, p=2,
                         theta=0, longlat=F, dMat1=t(dist_mat_tr_val), dMat2=dist_mat_tr)
#best model
fit2_pred <- gwr.predict(fit2$GW.arguments$formula, data_tr, data_val, bw=bw_opt, kernel="gaussian",adaptive=FALSE, p=2,
                         theta=0, longlat=F,dMat1=t(dist_mat_tr_val), dMat2=dist_mat_tr)
#basic lm
lm_mod_tr <- lm(lm_mod$call$formula, data_tr)
lm1_pred <- predict(lm_mod_tr, data_val)
#GAM
lm_mod2_tr <- lm(lm_mod2$call$formula, data_tr)
lm2_pred <- predict(lm_mod2_tr, data_val)

#plots
plot_val <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_val, aes(color = BDsample_0), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma") +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density data_val",
    color = "Value Legend"  # Adjust legend title
  )

common_limits7 <- range(c(data_val$BDsample_0, fit0_pred$SDF$prediction), na.rm = TRUE)
plot_val_fit0 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_val, aes(color = fit0_pred$SDF$prediction), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits = common_limits7) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density prediction fit0",
    color = "Value Legend"  # Adjust legend title
  )
(plot_val_fit0|plot_val)

common_limits8 <- range(c(data_val$BDsample_0, fit1_pred$SDF$prediction), na.rm = TRUE)
plot_val_fit1 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_val, aes(color = fit1_pred$SDF$prediction), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits = common_limits8) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density prediction fit1",
    color = "Value Legend"  # Adjust legend title
  )
(plot_val_fit1|plot_val)

common_limits9 <- range(c(data_val$BDsample_0, fit2_pred$SDF$prediction), na.rm = TRUE)
plot_val_fit2 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_val, aes(color = fit2_pred$SDF$prediction), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits = common_limits9) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density prediction fit2",
    color = "Value Legend"  # Adjust legend title
  )
(plot_val_fit2|plot_val)

common_limits10 <- range(c(data_val$BDsample_0, lm1_pred), na.rm = TRUE)
plot_val_lm1 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_val, aes(color = lm1_pred), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits = common_limits10) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density prediction lm1",
    color = "Value Legend"  # Adjust legend title
  )
(plot_val_lm1|plot_val)

common_limits11 <- range(c(data_val$BDsample_0, lm2_pred), na.rm = TRUE)
plot_val_lm2 <- ggplot() +
  geom_sf(data = europe_ue, fill = "lightgray", color = "black") +
  geom_sf(data = data_val, aes(color = lm2_pred), size = 1.5) +  # Replace `value_column` with your column name
  scale_color_viridis_c(option = "plasma", limits = common_limits11) +  # Use a Viridis color scale, change `option` if needed
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Bulk density prediction lm2",
    color = "Value Legend"  # Adjust legend title
  )
(plot_val_lm2|plot_val)

# compare predictions
#MAE
sum(abs(lm1_pred-data_val$BDsample_0))
sum(abs(lm2_pred-data_val$BDsample_0))
sum(abs(fit0_pred$SDF$prediction-data_val$BDsample_0))
sum(abs(fit1_pred$SDF$prediction-data_val$BDsample_0))
sum(abs(fit2_pred$SDF$prediction-data_val$BDsample_0))

#MSE
sum(abs(lm1_pred-data_val$BDsample_0)^2)
sum(abs(lm2_pred-data_val$BDsample_0)^2)
sum(abs(fit0_pred$SDF$prediction-data_val$BDsample_0)^2)
sum(abs(fit1_pred$SDF$prediction-data_val$BDsample_0)^2)
sum(abs(fit2_pred$SDF$prediction-data_val$BDsample_0)^2)

#Bootstrap
bts0 <- gwr.bootstrap(BDsample_0~1, data_tr, kernel = "gaussian", approach = "AIC",
              R = 1, k.nearneigh = 4, adaptive = TRUE, p = 2,
              theta = 0, longlat = FALSE, dMat=dist_mat_tr, verbose = FALSE,
              parallel.method = FALSE, parallel.arg = NULL)

B <- 10
T.boot <- matrix(nrow = dim(data)[1], ncol = B)
for(b in 1:B){
  response.b <- fit0$SDF$yhat + sample(fit0$SDF$residual, replace = T)
  fm.b <- gwr.basic(response.b~1, data, bw=bw_opt, kernel="gaussian",
                    adaptive=FALSE, p=2, theta=0, longlat=F,dMat=dist_mat)
  T.boot[,b] <- fm.b$SDF$Intercept
}
L.obs <- lm_mod$coefficients
B <- 100
T.boot.L <- matrix(nrow=14, ncol=B)
formula.b <- lm_mod$call$formula

pb <- progress_bar$new(
  format = "  processing [:bar] :percent eta: :eta",
  total = B,
  clear = FALSE)

for (b in 1:B) {
  response.b <- lm_mod$fitted.values + sample(lm_mod$residuals, replace = T)
  fm.b <- lm(response.b~Elev+pH_H2O+EC_2018+OC_2018+N_2018+K_2018+LC_0, data = data)
  sm.b <- summary(fm.b)
  T.boot.L[,b] <- sm.b$coefficients[,1]
  pb$tick()
}
hist(T.boot.L[2,])

alpha <- 0.05
right.quantile.L <- quantile(T.boot.L[2,], 1 - alpha/2)
left.quantile.L <- quantile(T.boot.L[2,], alpha/2)

CI.RP.L <- c(L.obs[2] - (right.quantile.L - L.obs[2]),
             L.obs[2],
             L.obs[2] - (left.quantile.L - L.obs[2]))
names(CI.RP.L) <- c("lwr", "lvl", "upr")
CI.RP.L
