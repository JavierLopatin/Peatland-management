
################################################################################
## R-Script: 3_Maps.R                                                         ##
## author: Javier Lopatin                                                     ##
## mail: javierlopatin@gmail.com                                              ##
##                                                                            ##
## Manuscript: Using aboveground vegetation attributes as proxies for mapping ##
## peatland belowground carbon stocks                                         ##
##                                                                            ##
## description: Create the prediction maps                                    ##
##                                                                            ##
################################################################################

library(raster)
library(randomForest)
library(caret)
library(doParallel)

### =============== ######################## =============== ###
### Spatial predicitons of aboveground vegetation attributes ###
### =============== ######################## =============== ###

setwd('/home/javier/Documents/PhD/Peatland/Managements_difference/')

## leave-one-out cross validation
train_control <- trainControl(method="LOOCV")

saveRasters = "rasters.RData"
#load(saveRasters)

# =======================================================
# Prediction with remote sensing data and Random Forests
# =======================================================

# MNF + structural data at the plot location
MNF <- read.table("/media/javier/Elements/Peatland1/MNFData.csv", header=T, sep=",", dec=".")[, 2:42]
names(MNF)

# load RS images
img = stack('/media/javier/Elements/Peatland1/MNFData.tif')
names(img) = colnames(MNF)

# Predict floristic composition (NMDS first axis)
m1data <- data.frame(x = data$NMDS.sp1[-c(33,34,39:44)], MNF)
FC_sp <- train(x ~., data=m1data, trControl=train_control, tuneLength = 20, method="rf")
pred_m1 <- predict(FC_sp)
plot(m1data$x, pred_m1); abline(0,1)
FC_img = predict(img, FC_sp)
plot(FC_img)

# Biomass
# herbaceous
m1data <- data.frame(x = data$Biomasa_herbaceas_kg_m2[-c(33,34,39:44)], MNF)
BM1 <- train(x ~., data=m1data, trControl=train_control, tuneLength = 20, method="rf")
pred_m1 <- predict(BM1)
plot(m1data$x, pred_m1); abline(0,1)
BM1_img = predict(img, BM1)
plot(BM1_img)
# shrubs
m1data <- data.frame(x = data$Biomasa_arbustivas_kg_m2[-c(33,34,39:44)], MNF)
BM2 <- train(x ~., data=m1data, trControl=train_control, tuneLength = 20, method="rf")
pred_m1 <- predict(BM2)
plot(m1data$x, pred_m1); abline(0,1)
BM2_img = predict(img, BM2)
plot(BM2_img)

# Plant diversity
# graminoids
m1data <- data.frame(x = data$gramm_richness[-c(33,34,39:44)], MNF)
DV1 <- train(x ~., data=m1data, trControl=train_control, tuneLength = 20, method="rf")
pred_m1 <- predict(DV1)
plot(m1data$x, pred_m1); abline(0,1)
DV1_img = predict(img, DV1)
plot(DV1_img)
# forbs
m1data <- data.frame(x = data$Herb_richness[-c(33,34,39:44)], MNF)
DV2 <- train(x ~., data=m1data, trControl=train_control, tuneLength = 20, method="rf")
pred_m1 <- predict(DV2)
plot(m1data$x, pred_m1); abline(0,1)
DV2_img = predict(img, DV2)
plot(DV2_img)


# =======================================================
# Make height and NDVI mask
# =======================================================

# load hyperspectral data
hyper <- stack("/media/javier/Elements/Peatland1/hyper_P1_2m.tif")
names(hyper) <- paste0( rep("B", 41), seq(1,41,1) )
hyper[hyper==0]<- NA
plot(hyper[[20]])

# load vegetaiton height data
DCM <- stack("/media/javier/Elements/Peatland1/treesvis/ndsm/DCM_2m.tif")

# resample both images to fit
DCM2 <- resample(DCM, hyper, resample='bilinear')
img <- stack(DCM2, hyper)
img[img$DCM_2m > 2] <- NA; img[img$DCM_2m < 0] <- NA # height mask
plot(img[[2]])

# create NDVI mask
NDVI <- ( hyper[[30]] - hyper[[20]] ) / ( hyper[[30]] + hyper[[20]] )
NDVI[NDVI < 0.3] <- NA
plot(NDVI)
img <- mask(img, NDVI); img <- mask(img, img$DCM_2m)
plot(img[[1]])

# ===============================
# save predictions
# ===============================

out_img = stack(img$DCM_2m, FC_img, BM1_img, BM2_img, DV1_img, DV2_img)
out_img =  mask(out_img, img$DCM_2m)
plot(out_img[[4]])

names(out_img) = PLS$model$gen$mvs_names[1:6]

writeRaster(out_img, 'RF_predictors.tif')

save.image(saveRasters)
