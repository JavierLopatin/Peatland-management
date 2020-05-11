
### Isomap ordenation

# import data
setwd("/home/javier/Documents/PhD/Peatland/")

library (isopam)
library (vegan)
library (MASS)

# floristic cover data
sp <- read.table("data/Cover_spp.csv", header=T, sep=",", dec=".")
summary(sp2)

# pft cover data
pft <- read.table("data/PFT1.csv", header=T, sep=",", dec=".")

#############################
### Ordination procidure ###
############################

# Selecting the cover of the spp.
ip <- isopam(sp[, 2:50], c.fix=3)
isotab(ip, 2)
ip$flat

# see groups in the space
xy <- data.frame( x=data$Coordenada_X_WGS84, y=data$Coordenada_Y_WGS84 )

svg(file = "Figures/plots_groups.svg", width=5, height=5)
plot(xy[ip$flat==1,], pch=16, ylim=c(min(xy$y),max(xy$y)),xlim=c(min(xy$x),max(xy$x)), cex=1.5,xaxt='n',yaxt='n',ann=FALSE)
points(xy[ip$flat==2,], pch=17, cex=1.5)
points(xy[ip$flat==3,], pch=15, cex=1.5)
dev.off()

##### Apply MDS  #####

# To spp level
# Selecting the number of k (stress value)
nmds.stress <- sapply(1:6, function(x) metaMDS(sp[, 2:50], k=x)$stress)
plot(1:6, nmds.stress)
# Selecting a number of dimensions: compromise between as few dimensions and low stress value
nmds1 <- metaMDS(sp[, 2:(length(sp)-2)], k=4, trymax=100)
nmds1$stress
## plot
ordiplot(nmds1, choices=c(1,2))
ordiplot(nmds1, choices=c(1,3))
ordiplot(nmds1, choices=c(1,4))
ordiplot(nmds1, choices=c(2,3))
ordiplot(nmds1, choices=c(2,4))
ordiplot(nmds1, choices=c(3,4))
# save scores
scor <- scores(nmds1)
NMDS.sp1 <- scor[, 1]
NMDS.sp2 <- scor[, 2]
NMDS.sp3  <- scor[, 3]
NMDS.sp4  <- scor[, 4]

ordination <- data.frame(sp$Plot, NMDS.sp1, NMDS.sp2, NMDS.sp3, NMDS.sp4)
write.table(ordination, file = "data/ordination.csv", sep = ",", col.names = T, row.names = F)

# To PFT level
# Selecting the number of k (stress value)
nmds.stress <- sapply(1:6, function(x) metaMDS(pft[,2:5], k=x)$stress)
plot(1:6, nmds.stress)
# Selecting a number of dimensions: compromise between as few dimensions and low stress value
nmds2 <- metaMDS(pft2, k=3, trymax=100)
nmds2$stress
## plot
ordiplot(nmds2, choices=c(1,2))
ordiplot(nmds2, choices=c(1,3))
ordiplot(nmds2, choices=c(2,3))
# save scores
scor <- scores(nmds2)
NMDS.PFT1 <- scor[, 1]
NMDS.PFT2 <- scor[, 2]
NMDS.PFT3  <- scor[, 3]

PFT <- data.frame(N=pft$Plot, PFT1=NMDS.PFT1, PFT2=NMDS.PFT2, PFT3=NMDS.PFT3)
write.table(PFT, file = "data/PFT.csv", sep = ",", col.names = T, row.names = F)

save.image("ordination.RData")


#-------------------------------------
# Biodiversity
#------------------------------------

# species level
shannon = diversity(sp[, 2:50])
simpson = diversity(sp[, 2:50], index='simpson')

# PFT level
shannon_pft = diversity(pft[, 2:5])
simpson_pft = diversity(pft[, 2:5], index='simpson')

diversity = cbind(sp$Plot, shannon, simpson, shannon_pft, simpson_pft)

write.table(diversity, file = "data/diversity.csv", sep = ",", col.names = T, row.names = F)
