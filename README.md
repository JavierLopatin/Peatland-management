# Peatland-management
#### Here, we explore the above- belowground relationship differences between a peatland with conservation and productive managements in their carbon stocks

##### author: Javier Lopatin
##### mail: javierlopatin@gmail.com

##### Manuscript: Using aboveground vegetation attributes as proxies for mapping
##### peatland belowground carbon stocks



#######################################################################################

##### Prepare data for modeling


```R
# set working directory
setwd("directory")
```

```R
# function to source scripts from GitHub
source_github <- function(u) {
  # load package
  require(RCurl)
  # read script lines from website and evaluate
  script <- getURL(u, ssl.verifypeer = FALSE)
  eval(parse(text = script), envir=.GlobalEnv)
  detach("package:RCurl", unload=TRUE)
}

# Load functions
# Get residuals from PLS-PM object
source_github('https://raw.githubusercontent.com/JavierLopatin/plspmSpatial/master/plspmResiduals.R')
# plspmPredict function for LV and measurement variables prediction
source_github('https://raw.githubusercontent.com/JavierLopatin/plspmSpatial/master/plspmPredict.R')
```

```R
# load general dataset
data <- read.table("data/Peatland1.csv", header=T, sep=",", dec=".")
names(data)

# load Floristic composition and biodiversity
# data produced in https://github.com/JavierLopatin/
ordination <- read.table("data/ordination.csv", header=T, sep=",", dec=".")
PFT <- read.table("data/PFT.csv", header=T, sep=",", dec=".")
biodiv <- read.table("data/diversity.csv", header=T, sep=",", dec=".")
```

################################################################################
###### START TUNING PLS-PM MODEL
################################################################################


```R
library(plspm)
library(plspm.formula)

## set model formula for C prediction
formula1 <-"

  # Outer model
  H =~ Altura_vegetacion_cm
  FC =~ NMDS.sp1
  BM =~ Biomasa_herbaceas_kg_m2 + Biomasa_arbustivas_kg_m2
  DV =~ gramm_richness + Herb_richness
  SD =~ Depth
  C =~ Carbono_Subterraneo_kg_m2 + Carbono_R1_kg_m2 + Carbono_musgo_kg_m2

  # Inner model
  H  ~~ FC
  BM ~~ FC + H
  DV ~~ FC + H
  SD ~~ FC +  BM + DV
  C  ~~ FC + H + BM + DV + SD
  "

# set model formula for soil depth prediction
formula2 <-"

  # Outer model
  H =~ Altura_vegetacion_cm
  FC =~ NMDS.sp1
  BM =~ Biomasa_herbaceas_kg_m2 + Biomasa_arbustivas_kg_m2
  DV =~ gramm_richness + Herb_richness
  SD =~ Depth

  # Inner model
  H  ~~ FC
  BM ~~ FC + H
  DV ~~ FC + H
  SD ~~ FC +  BM + DV
  "

# set the "Modes" of the relations: character vector indicating the type of measurement for each block.
modes = rep("A", 6)
```
```R
# model for belowground C stocks predictions
set.seed(123)
PLS = plspm.formula(formula1, data, modes=modes, maxiter= 500, boot.val = T, br = 500,
             scheme = "factor", scaled = T)
PLS$outer_model
PLS$inner_summary

# model for soul depth predictions
set.seed(123)
soil = plspm.formula(formula2, data, modes=rep('A',5), maxiter= 500, boot.val = T, br = 500,
                     scheme = "factor", scaled = T)
soil$outer_model
soil$inner_summary


# save bootstrapping coefficient path
write.table(PLS$effects, "Managements_difference/effects.csv", sep = ",")
```

################################################################################
###### END TUNING PLS-PM MODEL
################################################################################

################################################################################
##### START SPATIAL AUTOCRRELATION TESTS
################################################################################



```R
# check for spatial autocorrelation
library(ncf) # correlogram
library(ape) # Moran´s I

# get residuals
residuals = plspmResiduals(PLS)
# xy coordinates
xy <- data.frame( x=data$Coordenada_X_WGS84, y=data$Coordenada_Y_WGS84 )
plot(xy)

# obtain distance matrix of XY coordinates
xy.dists <- as.matrix(dist(cbind(xy$x, xy$y)))
xy.dists.inv <- 1/xy.dists # invers
diag(xy.dists.inv) <- 0
xy.dists.inv[1:5, 1:5]# check

# Aboveground biomass
moran.BM <- Moran.I(residuals$inner_residuals[,2], xy.dists.inv); moran.BM
corr.BM <- correlog (xy[,1], xy[,2], z = residuals$inner_residuals[,2],
                     increment = 10, resamp = 500, quiet = T)
plot(corr.BM); grid(); abline(0,0, lty=2)

# Species richness
moran.Rich <- Moran.I(residuals$inner_residuals[,3], xy.dists.inv); moran.Rich
corr.Rich <- correlog (xy[,1], xy[,2], z = residuals$inner_residuals[,3],
                       increment = 10, resamp = 500, quiet = T)
plot(corr.Rich); grid(); abline(0,0, lty=2)

# Soil depht
moran.soil <- Moran.I(residuals$inner_residuals[,4], xy.dists.inv); moran.soil
corr.soil <- correlog (xy[,1], xy[,2], z = residuals$inner_residuals[,4],
                       increment = 10, resamp = 500, quiet = T)
plot(corr.soil); grid(); abline(0,0, lty=2)

# Belowground C stock
moran.C <- Moran.I(residuals$inner_residuals[,5], xy.dists.inv); moran.C
corr.C <- correlog (xy[,1], xy[,2], z = residuals$inner_residuals[,5],
                    increment = 10, resamp = 500, quiet = T)
plot(corr.C); grid(); abline(0,0, lty=2)

save.image("compare.RData")
```

################################################################################
##### END SPATIAL AUTOCRRELATION TESTS
################################################################################

################################################################################
##### START MANAGEMENT DIFFERENCE ASSESSMENT
################################################################################

```R
library(raster)

# raster predictions using RF
img = stack('Managements_difference/RF_predictors.tif')
names(img) = PLS$model$gen$mvs_names[1:6]

# set the bootstrap parameters
N = nrow(data) # N° of observations
B = 500        # N° of bootstrap iterations

# pred C
pred_C = list()
pred_C_cons = list()
pred_C_agri = list()
# plspm.groupsPredict
pred <- list()
pred_cons <- list()
pred_agri <- list()
# store management id
id = list()
# coefficients
coef <- list()
coef_cons <- list()
coef_agri <- list()
# total effects
tot_eff <- matrix(nrow = B, ncol = nrow(PLS$effects)); colnames(tot_eff) = PLS$effects$relationships
tot_eff_cons <- matrix(nrow = B, ncol = nrow(PLS$effects)); colnames(tot_eff_cons) = PLS$effects$relationships
tot_eff_agri <- matrix(nrow = B, ncol = nrow(PLS$effects)); colnames(tot_eff_agri) = PLS$effects$relationships
# accuracies
r2_global <- matrix(ncol=ncol(PLS$scores), nrow=B); colnames(r2_global) = colnames(PLS$scores)
r2_cons <- matrix(ncol=ncol(PLS$scores), nrow=B); colnames(r2_cons) = colnames(PLS$scores)
r2_agri <- matrix(ncol=ncol(PLS$scores), nrow=B); colnames(r2_agri) = colnames(PLS$scores)
#GoF
overall_gof <- 0
cons_gof <- 0
agri_gof <- 0
# maps
SD_img <- list()
SD_img_cons <- list()
SD_img_agri <- list()
C_img <- list()
C_img_cons <- list()
C_img_agri <- list()

for(i in 1:B){

  print(i)

  # if there's a problem during bootstrapping, skip iteration
  skip_to_next <- FALSE

  # create random numbers with replacement to select samples from each group
  set.seed(100+i)
  idx = sample(1:N, N, replace=TRUE)

  # set training data
  train = data[idx,]
  test  = data[-idx,]

  # tune PLSPM
  tryCatch({
    # train C models
    fit = plspm.formula(formula1, train, modes=modes, boot.val = F, scheme = "factor", scaled = T)
    fit_con = plspm.formula(formula1, train[train$Uso == 'Conservacion',], modes=modes, boot.val = F, scheme = "factor", scaled = T)
    fit_agri = plspm.formula(formula1, train[train$Uso == 'Productivo',], modes=modes, boot.val = F, scheme = "factor", scaled = T)

    # coefficients C models
    coef[[i]] <- fit$path_coefs
    coef_cons[[i]] <- fit_con$path_coefs
    coef_agri[[i]] <- fit_agri$path_coefs

    # GoF C models
    overall_gof[i] <- fit$gof
    cons_gof[i] <- fit_con$gof
    agri_gof[i] <- fit_agri$gof

    # R2 C models
    r2_global[i,] <- fit$inner_summary$R2
    r2_cons[i,] <- fit_con$inner_summary$R2
    r2_agri[i,] <- fit_agri$inner_summary$R2

    # Total effects C models
    tot_eff[i,] <- fit$effects$total
    tot_eff_cons[i,] <- fit_con$effects$total
    tot_eff_agri[i,] <- fit_agri$effects$total

    # train soil depth models
    s = plspm.formula(formula2, train, modes=rep('A',5), boot.val = F, scheme = "factor", scaled = T)
    s2 = plspm.formula(formula2, train[train$Uso == 'Conservacion',], modes=rep('A',5), boot.val = F, scheme = "factor", scaled = T)
    s3 = plspm.formula(formula2, train[train$Uso == 'Productivo',], modes=rep('A',5), boot.val = F, scheme = "factor", scaled = T)

    # predict soil depth
    ss = plspmPredict(s, test)
    ss2 = plspmPredict(s2, test)
    ss3 = plspmPredict(s3, test)

    # predict C models using predicted soil depth data
    t1 = test
    t2 = test
    t3 = test
    t1$Depth = ss$mmPredicted[,6]
    t2$Depth = ss2$mmPredicted[,6]
    t3$Depth = ss3$mmPredicted[,6]

    p1 = plspmPredict(fit, t1)
    p2 = plspmPredict(fit_con, t2)
    p3 = plspmPredict(fit_agri, t3)
    pred[[i]] <- p1
    pred_cons[[i]] <- p2
    pred_agri[[i]] <- p3

    # Maps soil depth
    soil_img1 = plspmPredict(s, img)
    soil_img2 = plspmPredict(s2, img)
    soil_img3 = plspmPredict(s3, img)
    SD_img[[i]] <- soil_img1
    SD_img_cons[[i]] <- soil_img2
    SD_img_agri[[i]] <- soil_img3

    # Maps C using the predicted soil depth
    c_img1 = plspmPredict(fit, stack(img, soil_img1$mmPredicted))
    c_img2 = plspmPredict(fit_con, stack(img, soil_img2$mmPredicted))
    c_img3 = plspmPredict(fit_agri, stack(img, soil_img3$mmPredicted))
    C_img[[i]] <- c_img1
    C_img_cons[[i]] <- c_img2
    C_img_agri[[i]] <- c_img3

   }, error = function(e) { skip_to_next <<- TRUE})

}
```

```R
# number of valid iterations
nrow(na.omit(r2_global))
# which are nan
nan_data = which(is.na(r2_global[,1]))

# prepare list of data
overall_data <- list(overall_gof[-nan_data], na.omit(tot_eff[-nan_data,]), r2_global[-nan_data,], coef[-nan_data] )
names(overall_data) <- c('gof','tot_eff', 'r2', 'coeff')

cons_data <- list(cons_gof[-nan_data], na.omit(tot_eff_cons[-nan_data,]), r2_cons[-nan_data,], coef_cons[-nan_data])
names(cons_data) <- c('gof','tot_eff', 'r2', 'coeff')

agri_data <- list(agri_gof[-nan_data], na.omit(tot_eff_agri[-nan_data,]), r2_agri[-nan_data,], coef_agri[-nan_data])
names(agri_data) <- c('gof','tot_eff', 'r2', 'coeff')

# median iterative values
library(robustbase)

median(overall_data$gof)
colMedians(overall_data$tot_eff)
colMedians(overall_data$r2)
apply(simplify2array(overall_data$coeff), 1:2, median)

median(cons_data$gof)
colMedians(cons_data$tot_eff)
colMedians(cons_data$r2)
apply(simplify2array(cons_data$coeff), 1:2, median)

median(agri_data$gof)
colMedians(agri_data$tot_eff)
colMedians(agri_data$r2)
apply(simplify2array(agri_data$coeff), 1:2, median)

```

```R
#------------------------------------------------
## One-sided bootstrapping significance test
#------------------------------------------------

#########################################
# Sig. path coefficients in each model
#########################################

pathSignificance <- function(m1){

  # get median and percentile vaules (alpha=0.05)
  med = apply(simplify2array(m1$coeff), 1:2, median); med[med == 0] = NA
  p05 = apply(simplify2array(m1$coeff), 1:2, quantile, probs=0.05); p05[p05 == 0] = NA
  p95 = apply(simplify2array(m1$coeff), 1:2, quantile, probs=0.95);  p95[p95 == 0] = NA

  # chech significance
  s = sign(sign(p05) == sign(p95))

  # prepare output
  output <- list(med, s==1)
  names(output) = c('median', 'Sig.')

  return(output)
}

# overall
sig_path_over = pathSignificance(overall_data); sig_path_over
# conservation
sig_path_cons = pathSignificance(cons_data); sig_path_cons
# agriculture
sig_path_agri = pathSignificance(agri_data); sig_path_agri

## Function to obtain significant differences between SEM models
Acc_diff_Significance <- function(m1, m2){
  # differences between overall and conservation model
  gof <-  m1$gof - m2$gof
  r2 <-   m1$r2 - m2$r2

  # prepare output
  output <- list(gof, r2[,2], r2[,3],r2[,4],r2[,5],r2[,6])

  # matrix of significances
  a = matrix(nrow = 6, ncol = 3)
  colnames(a) <- c("0.1", "0.05", "0.001")
  rownames(a) <- c('gof', 'r2_H', 'r2_BM', 'r2_SR', 'r2_SD', 'r2_C')
  for (i in 1:nrow(a)){
    # 0.1
    if ( sign(quantile(output[[i]], probs=c(0.1))) == sign(quantile(output[[i]], probs=c(0.9))) ){
      a[i,1] = "True"
    } else{
      a[i,1] = "False"
    }
    # 0.05
    if ( sign(quantile(output[[i]], probs=c(0.05))) == sign(quantile(output[[i]], probs=c(0.95))) ){
      a[i,2] = "True"
    } else{
      a[i,2] = "False"
    }
    # 0.001
    if ( sign(quantile(output[[i]], probs=c(0.001))) == sign(quantile(output[[i]], probs=c(0.995))) ){
      a[i,3] = "True"
    } else{
      a[i,3] = "False"
    }
  }
  return(a)
}

## overall vs conservation
acc_sig_overall_cons <- Acc_diff_Significance(overall_data, cons_data); acc_sig_overall_cons
## overall vs agriculture
acc_sig_overall_agri <- Acc_diff_Significance(overall_data, agri_data); acc_sig_overall_agri
## overall vs conservation
acc_sig_cons_agri <- Acc_diff_Significance(cons_data, cons_data); acc_sig_cons_agri

##########################################################
# Sig. differences in path coefficients between each model
##########################################################

path_diff_Significance <- function(m1,m2){
  # get differences between path coefficients
  diff = list()
  for (i in 1:length(m1$gof)){
    diff[[i]] = m1$coeff[[i]] - m2$coeff[[i]]
  }

  # get percentile vaules permodel
  p01 = apply(simplify2array(diff), 1:2, quantile, probs=0.01); p01[p01 == 0] = NA
  p05 = apply(simplify2array(diff), 1:2, quantile, probs=0.05); p05[p05 == 0] = NA
  p10 = apply(simplify2array(diff), 1:2, quantile, probs=0.10); p10[p10 == 0] = NA
  p90 = apply(simplify2array(diff), 1:2, quantile, probs=0.90);  p90[p90 == 0] = NA
  p95 = apply(simplify2array(diff), 1:2, quantile, probs=0.95);  p95[p95 == 0] = NA
  p99 = apply(simplify2array(diff), 1:2, quantile, probs=0.99);  p99[p99 == 0] = NA

  # chech significance
  s1 = sign(sign(p10) == sign(p90))
  s2 = sign(sign(p05) == sign(p95))
  s3 = sign(sign(p01) == sign(p99))

  output = list(s1,s2,s3)
  names(output) = c('0.1','0.05','0.001')
  return(output)
}

path_diff_over_cons = path_diff_Significance(overall_data, cons_data); path_diff_over_cons
path_diff_over_agri = path_diff_Significance(overall_data, agri_data); path_diff_over_agri
path_diff_cons_agri = path_diff_Significance(cons_data, agri_data); path_diff_cons_agri

##################################################
# Sig. of total effects (direct + indirect) over C
##################################################

getSignificance_tot <- function(m1){
  # total effects only towards C
  FC  = m1$tot_eff[,5]
  H   = m1$tot_eff[,9]
  BM  = m1$tot_eff[,12]
  SR  = m1$tot_eff[,14]
  SD  = m1$tot_eff[,15]

  # prepare output
  output <- list(FC, H, BM, SR, SD)

  # matrix of significances
  a = matrix(nrow = 5, ncol = 3)
  colnames(a) <- c("0.1", "0.05", "0.001")
  rownames(a) <- c('FC', 'H', 'BM', 'SR', 'SD')
  for (i in 1:nrow(a)){
    # 0.1
    if ( sign(quantile(output[[i]], probs=c(0.1))) == sign(quantile(output[[i]], probs=c(0.9))) ){
      a[i,1] = "True"
    } else{
      a[i,1] = "False"
    }
    # 0.05
    if ( sign(quantile(output[[i]], probs=c(0.05))) == sign(quantile(output[[i]], probs=c(0.95))) ){
      a[i,2] = "True"
    } else{
      a[i,2] = "False"
    }
    # 0.001
    if ( sign(quantile(output[[i]], probs=c(0.001))) == sign(quantile(output[[i]], probs=c(0.995))) ){
      a[i,3] = "True"
    } else{
      a[i,3] = "False"
    }
  }

  r = list(a, colMedians( m1$tot_eff)[c(5,9,12,14,15)])
  names(r) = c('Sig','Median')
  return(r)
}

tot_eff2 = getSignificance_tot(overall_data); tot_eff2
tot_eff_cons2 = getSignificance_tot(cons_data); tot_eff_cons2
tot_eff_agri2 = getSignificance_tot(agri_data); tot_eff2

save.image('compare.RData')
```

################################################################################
##### END MANAGEMENT DIFFERENCE ASSESSMENT
################################################################################

################################################################################
##### START PREDICTIVE MAPS
################################################################################
```R
# function to unlist the plspmPredict objects
unlist_raster = function(SD_img, method=1, var=6){
  m = raster(img)
  for (i in seq(1,B)[-nan_data]){
    if (method==1)
      m = addLayer(m, SD_img[[i]]$mmPredicted)
    if (method==2)
      m = addLayer(m, SD_img[[i]]$Scores[[var]])
  }
  return(m)
}

# soil depth
SD_overall_unlist = unlist_raster(SD_img)
SD_cons_unlist = unlist_raster(SD_img_cons)
SD_agri_unlist = unlist_raster(SD_img_agri)

SD_overall_med = calc(SD_overall_unlist, median)
SD_overall_cv = calc(SD_overall_unlist, cv)
plot(SD_overall_med)
plot(SD_overall_cv)

SD_cons_med = calc(SD_cons_unlist, median)
SD_cons_cv = calc(SD_cons_unlist, cv)
plot(SD_cons_med)
plot(SD_cons_cv)

SD_agri_med = calc(SD_agri_unlist, median)
SD_agri_cv = calc(SD_agri_unlist, cv)
plot(SD_agri_med)
plot(SD_agri_cv)

# C
C_overall_unlist = unlist_raster(C_img,2)
C_cons_unlist = unlist_raster(C_img_cons,2)
C_agri_unlist = unlist_raster(C_img_agri,2)

C_overall_med = calc(C_overall_unlist, median)
C_overall_cv = calc(C_overall_unlist, cv)
plot(C_overall_med)
plot(C_overall_cv)

C_cons_med = calc(C_cons_unlist, median)
C_cons_cv = calc(C_cons_unlist, cv)
plot(C_cons_med)
plot(C_cons_cv)

C_agri_med = calc(C_agri_unlist, median)
C_agri_cv = calc(C_agri_unlist, cv)
plot(C_agri_med)
plot(C_agri_cv)

# merge management results
area = shapefile('/media/javier/Elements/Peatland1/area_cut.shp')
area = rasterize(area, raster(img))

# spatial resolution of land cover map
SD_cons_med <- resample(SD_cons_med, area)
SD_cons_cv <- resample(SD_cons_cv, area)
SD_agri_med <- resample(SD_agri_med, area)
SD_agri_cv <- resample(SD_agri_cv, area)

C_cons_med <- resample(C_cons_med, area)
C_cons_cv <- resample(C_cons_cv, area)
C_agri_med <- resample(C_agri_med, area)
C_agri_cv <- resample(C_agri_cv, area)

# convert rasters to dataframes
sd_cons_df_med = values(SD_cons_med)
sd_cons_df_cv = values(SD_cons_cv)
sd_agri_df_med = values(SD_agri_med)
sd_agri_df_cv = values(SD_agri_cv)

C_cons_df_med = values(C_cons_med)
C_cons_df_cv = values(C_cons_cv)
C_agri_df_med = values(C_agri_med)
C_agri_df_cv = values(C_agri_cv)

# convert management areas to vector
fcl_v <- values(area)

# fill up NA values with 0
fcl_v[is.na(fcl_v)] <- 0

# copy vector to have the same dimensions for storing results
res_sd <- fcl_v
res_sd_cv <- fcl_v
res_c <- fcl_v
res_c_cv <- fcl_v

# start the merging
for (i in 1:length(res_v)) {
  print(i)
  # if pixel of landcover map is one of the non relevant classes assign 0
  if (fcl_v[i] == 0){
    res_sd[i] <- 0
    res_sd_cv[i] <- 0
    res_c[i] <- 0
    res_c_cv[i] <- 0
  } else if (fcl_v[i] == 1){
    res_sd[i] <- sd_cons_df_med[i]
    res_sd_cv[i] <- sd_cons_df_cv[i]
    res_c[i] <- C_cons_df_med[i]
    res_c_cv[i] <- C_cons_df_cv[i]

    # if pixel of land-cover is native forest, assign native forest biomass
  } else if (fcl_v[i] == 2) {
    res_sd[i] <- sd_agri_df_med[i]
    res_sd_cv[i] <- sd_agri_df_cv[i]
    res_c[i] <- C_agri_df_med[i]
    res_c_cv[i] <- C_agri_df_cv[i]

  } else {
    print("this cannot be")
  }
}

# write rastes
setwd('/home/javier/Documents/PhD/Peatland/Managements_difference/')

# overall models
writeRaster(SD_overall_med, 'SD_overall_med.tif', overwrite=T)
writeRaster(SD_overall_cv, 'SD_overall_cv.tif', overwrite=T)
writeRaster(C_overall_med, 'C_overall_med.tif', overwrite=T)
writeRaster(C_overall_cv, 'C_overall_cv.tif', overwrite=T)

# write management rasters
writeOut <- function(res_sd, sd_out, out='soildepth_manage.tif'){
  sd_out = raster(area)
  values(sd_out) = res_sd
  sd_out[sd_out==0] <- NA
  writeRaster(sd_out, out)
}

writeOut(res_sd, sd_out, 'soildepth_manage.tif')
writeOut(res_sd_cv, sd_out, 'soildepth_manage_cv.tif')
writeOut(res_c, sd_out, 'C_manage.tif')
writeOut(res_c_cv, sd_out, 'C_manage_cv.tif')
```
################################################################################
##### END PREDICTIVE MAPS
################################################################################
