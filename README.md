# Peatland-management
#### Here, we explore the avove- belowground interdependencies differences between a peatland with conservation and productive managements in their carbon stock

##### author: Javier Lopatin
##### mail: javierlopatin@gmail.com

##### Manuscript: Using aboveground vegetation attributes as proxies for mapping
##### peatland belowground carbon stocks



#######################################################################################

##### Prepare data for modeling


```R
# set working directory
setwd("C:/Users/Lopatin/Dropbox/PhD/Peatland/Managements difference")
```


```R
# load Floristic composition
# data produced in https://github.com/JavierLopatin/Peatland-carbon-stock
ordination <- read.table("data/ordination_sp.csv", header=T, sep=",", dec=".")
```


```R
# new variables to include in the analysis
data$gramm_richness <- data$Reeds_richness + data$Ferns_richness + data$Grass_richness
data$NMDS.sp1 = ordination$NMDS.sp1
```

##### Plot-based models

################################################################################
###### START TUNING PLS-PM MODEL
################################################################################


```R
library(plspm)

# =======================================================
# Set the inner model
# =======================================================

# Set the rows of the inner model matrix.
# 1 = recive an interaction
# 0 = do not recive an interaction
FC = c(0, 0, 0, 0, 0, 0)
H  = c(1, 0, 0, 0, 0, 0)
BM = c(1, 1, 0, 0, 0, 0)
SR = c(1, 1, 0, 0, 0, 0)
SD = c(1, 0, 1, 1, 0, 0)
C  = c(1, 1, 1, 1, 1, 0)

# matrix created by row binding.
inner = rbind(FC, H, BM, SR, SD, C) ; colnames(inner) = rownames(inner)
# plot the inner matrix
innerplot(inner)

# set the "Modes" of the relations: character vector indicating the type of measurement for each block.
modes = rep("A", 6)

# =======================================================
# Set the outer model
# =======================================================

# define which variables are going to be used in each LV.
# selection is done by removing variables with outer correlation (l) below 0.5
# here it is the final selection
outer = list (c("Altura_vegetacion_cm"),                               # heigts
              c("NNMDS.sp1"),                              # FC
              c("Biomasa_herbaceas_kg_m2","Biomasa_arbustivas_kg_m2"), # Biomass
              c("gramm_richness","Herb_richness"),                     # Richness
              c("Depth"),                                              # soil depth
              c("NCarbono_Subterraneo_kg_m2","Carbono_musgo_kg_m2",
                "Carbono_R1_kg_m2"))                                   # C

# Turn some variable into negative values. This is necessary to run PLS-PM when
# variables appear as negative in the outer model. You need to run plspm once
# to notice this. More info at: https://www.gastonsanchez.com/PLS_Path_Modeling_with_R.pdf
data$NNMDS.sp1 <- data$NMDS.sp1 * -1
data$NCarbono_R3_kg_m2 = data$Carbono_R3_kg_m2 * -1
data$NCarbono_Subterraneo_kg_m2 = data$Carbono_Subterraneo_kg_m2 * -1
```


```R
# Run PLSPM for aboveground C stock
set.seed(123)
PLS = plspm(data, inner, outer, modes, maxiter= 500, boot.val = T, br = 500,
            scheme = "factor", scaled = T)
PLS$outer
PLS$inner_summary
PLS$inner_model
PLS$gof
PLS$path_coefs
PLS$boot

# save bootstrapping coefficient path
write.table(PLS$effects, "Managements_difference/effects.csv", sep = ",")

# plot results
innerplot(PLS, arr.pos = 0.35) # inner model

save.image("peatland.RData")
```

################################################################################
###### END TUNING PLS-PM MODEL
################################################################################

################################################################################
###### START MANAGEMENT DIFFERENCE ASSESSMENT
################################################################################


```R
set.seed(123)
compare = plspm.groups(PLS, data$Uso, method = "permutation", reps=100)
compare_data = compare$test

write.csv(compare_data, file="Managements_difference/compare_test.csv")

# significant path coefficients
set.seed(123)
pls1 = plspm(data[data$Uso == 'Conservacion',], inner, outer, modes, maxiter= 500, boot.val = T, br = 500,
              scheme = "factor", scaled = T)

# model scores
Scores <- PLS$scores; colnames(Scores) <- colnames(inner)
# conservation and productive data
obs_cons = grep( "Conservacion", data$Uso )
obs_prod = grep( "Productivo", data$Uso )
Scores1 <- Scores[obs_cons, ]
Scores2 <- Scores[obs_prod, ]

save.image("Managements_difference/compare.RData")

```


```R
### Estimate R²s
# Conservation site
model1 <- compare$group1
colnames(Scores1)
H1     <- model1$H[1] + model1$H[2]*Scores1[,1]
BM1    <- model1$BM[1] + model1$BM[2]*Scores1[,1] + model1$BM[3]*Scores1[,2]
Rich1  <- model1$Rich[1] + model1$Rich[2]*Scores1[,1] + model1$Rich[3]*Scores1[,2]
Depth1 <- model1$Depth[1] + model1$Depth[2]*Scores1[,1] + model1$Depth[3]*Scores1[,3] + model1$Depth[4]*Scores1[,4]
C1     <- model1$C[1] + model1$C[2]*Scores1[,1] + model1$C[3]*Scores1[,2] + model1$C[4]*Scores1[,3] +
          model1$C[5]*Scores1[,4] + model1$C[6]*Scores1[,5]

# obtain fit
r2_H1    = cor(H1, Scores1[,2], method="pearson")^2
r2_BM1   = cor(BM1, Scores1[,3], method="pearson")^2
r2_Rich1 = cor(Rich1, Scores1[,4], method="pearson")^2
r2_Depth1= cor(Depth1, Scores1[,5], method="pearson")^2
r2_C1    = cor(C1, Scores1[,6], method="pearson")^2

# Managed site
model2 <- compare$group2
colnames(Scores2)
H2     <- model2$H[1] + model2$H[2]*Scores2[,1]
BM2    <- model2$BM[1] + model2$BM[2]*Scores2[,1] + model2$BM[3]*Scores2[,2]
Rich2  <- model2$Rich[1] + model2$Rich[2]*Scores2[,1] + model2$Rich[3]*Scores2[,2]
Depth2 <- model2$Depth[1] + model2$Depth[2]*Scores2[,1] + model2$Depth[3]*Scores2[,3] + model2$Depth[4]*Scores2[,4]
C2     <- model2$C[1] + model2$C[2]*Scores2[,1] + model2$C[3]*Scores2[,2] + model2$C[4]*Scores2[,3] +
          model2$C[5]*Scores2[,4] + model2$C[6]*Scores2[,5]

# obtain fit
r2_H2     = cor(H2, Scores2[,2], method="pearson")^2
r2_BM2    = cor(BM2, Scores2[,3], method="pearson")^2
r2_Rich2  = cor(Rich2, Scores2[,4], method="pearson")^2
r2_Depth2 = cor(Depth2, Scores2[,5], method="pearson")^2
r2_C2     = cor(C2, Scores2[,6], method="pearson")^2

### estimate global goodness-of-fit (GoF)
GoF1 <- sqrt( mean(PLS$outer_model$communality) * mean( c(r2_H1, r2_BM1, r2_Rich1, r2_Depth1, r2_C1) ) )
GoF2 <- sqrt( mean(PLS$outer_model$communality) * mean( c(r2_H2, r2_BM2, r2_Rich2, r2_Depth2, r2_C2) ) )

```

################################################################################
###### END MANAGEMENT DIFFERENCE ASSESSMENT
################################################################################
