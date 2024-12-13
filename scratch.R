#read data from csv file
data2018 <- read.csv("LUCAS-SOIL-2018-data-report-readme-v2/LUCAS-SOIL-2018-v2/LUCAS-SOIL-2018.csv", head = T)
data2015 <- read.csv("LUCAS2015_topsoildata_20200323/LUCAS_Topsoil_2015_20200323.csv", head = T)
dataBD <- read.csv("BD_LUCAS_data_for_paper/BD2018_data_for_paper.csv", head = T)
head(data2018)
head(data2015)
head(dataBD)
#check data
dim(data2018)
dim(data2015)
dim(dataBD)
names(data2018)
names(data2015)
names(dataBD)
summary(data2018)
summary(data2015)
summary(dataBD)
#convert strings into doubles
#data2018
names(data2018)
sapply(data2018, class)
data2018$OC <- as.numeric(data2018$OC)
data2018$CaCO3 <- as.numeric(data2018$CaCO3)
data2018$P <- as.numeric(data2018$P)
data2018$N <- as.numeric(data2018$N)
data2018$K <- as.numeric(data2018$K)
#data2015
names(data2015)
sapply(data2015, class)
data2015$Coarse <- as.numeric(data2015$Coarse)
data2015$Clay <- as.numeric(data2015$Clay)
data2015$Sand <- as.numeric(data2015$Sand)
data2015$Silt <- as.numeric(data2015$Silt)
data2015$CaCO3 <- as.numeric(data2015$CaCO3)
#dataBD
names(dataBD)
sapply(dataBD, class)
dataBD$BDsample_0 <- as.numeric(dataBD$BDsample_0)
dataBD$coarse_mas <- as.numeric(dataBD$coarse_mas)
dataBD$BDfine_0_2 <- as.numeric(dataBD$BDfine_0_2)
dataBD$coarse_vol <- as.numeric(dataBD$coarse_vol)

#intersection points
p_ID_intersect_15_18 <- intersect(data2015$Point_ID, data2018$POINTID)
p_ID_intersect_15_18_BD <- intersect(p_ID_intersect_15_18, dataBD$POINT_ID)
data2018_inter <- data2018[which(data2018$POINTID %in% p_ID_intersect_15_18_BD),]
data2015_inter <- data2015[which(data2015$Point_ID %in% p_ID_intersect_15_18_BD),]
dataBD_inter <- dataBD[which(dataBD$POINT_ID %in% p_ID_intersect_15_18_BD),]
