library(vcfR)
library(ape)
library(adegenet)
library(ggplot2)

vcf <- read.vcfR("input/03_366_miss9.recode.vcf", verbose = FALSE)
x <- vcfR2genlight(vcf)
x

#import table
indv_info <- read.table("input/366ALB_China.txt", header = T)
head(indv_info)
x@pop <- as.factor(indv_info$popmap)
x@pop <- as.factor(indv_info$Region)
#-------------------####
#       PCA         ####----->
#-------------------####
pca=glPca(x, parallel = F)
30 #number of axes
#plot pca
myCol <- colorplot(pca$scores, pca$scores, transp = T, cex = 3, col=col)
#with individual names
scatter(pca)
#view coordinates
head(pca$scores)
##for color the plot!
x@strata <- indv_info
#preview info
head(x@strata)
levels(x@strata$popmap)
levels(x@strata$Country)
levels(x@strata$Region)
length(levels(x@strata$popmap))
### plot by region
col <- c("steelblue1","purple1","red2","darkolivegreen4","black")[x@strata$Region]
col
plot(PC2~PC1, data=pca$scores, col=col, pch=(19), cex=1.5)
###try ggplot
...

#-------------------####
#       DAPC        ####----->
#-------------------####
grp=find.clusters(x,max.n.clust = 16, parallel=FALSE)
366#keep all, here 
#choose the least BIC values
6

dapc=dapc(x,grp$grp,parallel=FALSE)
150
5

#a way to estimate the best number of PCs to retain
optim.a.score(dapc)

table(pop(x), grp$grp)

#### new discovery 2021-jan-25, use a prior pop info as group info ####---------------------------------
dapc_a_prior=dapc(x,parallel=FALSE)

scatter(dapc_a_prior,clab=0.6,solid=.7) 
ppi <- 300
dev.copy(png, "output/03_366_miss9_a_prior_0125.png", width=6*ppi, height=6*ppi, res=ppi)
dev.off()

groupe_a_prior=as.data.frame(dapc_a_prior$grp)
write.table(groupe_a_prior, "output/03_366_miss9_a_prior_0125.txt")
#-----------------------------------------------------------------------------

myCol <- c("#7fbf7b","#762a83","#FFDAB9","#ADD8E6","#d9f0d3","#c2a5cf")
# c("#984ea3","#CD2626","#ff7f00","#999999","#377eb8","#4daf4a")
scatter(dapc, col=myCol,clab=0,solid=.7)
#clab is the label

#re-shape the form
#scatter(dapc, scree.da=FALSE, bg="white", pch=20, cell=3, cstar=0, col=myCol, solid=.4,cex=3,clab=1, leg=TRUE, txt.leg=c("1","2","3","4","5","6","7"))
scatter(dapc,  bg="white", 
        pch=20, cell=0, cstar=0, col=myCol, 
        solid=0.9, cex=2,clab=0,leg=FALSE) #scree.da=FALSE,

ppi <- 300
dev.copy(png, "output/03_366_miss9_DAPC_group_0425.png", width=6*ppi, height=6*ppi, res=ppi)
dev.off()

groupe=as.data.frame(dapc$grp)
write.table(groupe, "output/03_366_miss9_whichGroup_0119.txt")
#barplot
compoplot(dapc)

## Loadings
contrib = loadingplot(dapc$var.contr, axis=2, lab.jitter=1)
loadingplot(head(dapc$var.contr, 100), main="Loading plot - first 100 SNPs")
loadingplot(dapc$var.contr, axis=2, lab.jitter=1)
head(dapc$var.contr)

write.csv(contrib$var.values,"output/03_366_miss9_DAPC_contrib_1208.txt")

ppi <- 300
dev.copy(png, "output/03_366_miss9_DAPC_loading_plot_0928.png", width=6*ppi, height=6*ppi, res=ppi)
dev.off()

loading=as.data.frame(dapc$var.contr)
write.table(loading,"output/03_366_miss9_DAPC_loading_table_1208.txt")
which.max(loading$LD1)

#  Plot map and piecharts 
library(maps)
library(mapdata)
library(maptools)
library(mapplots)
library(dplyr)
library(magrittr)
library(tidyr)

coord = read.table("input/coord_16pop.txt",header=T)
coord1<-coord[,c(2,3)]
head(coord)

#-------------------####
#   ploting pie     ####----->
#-------------------####
par(mfrow=c(1,1),oma=c(0,0,0,0))#xiazuoshangyou
plot(coord1, xlim=c(104,130), ylim=c(29.5,46.5),ann = F)
title(xlab= 'Longitude', ylab = 'Latitude', line = 2)
map(add = T,xlim=c(104,130), ylim=c(29.5,46.5))
#map(add = T, col = "white", fill = TRUE,xlim=c(104,130), ylim=c(29,46))
#map('world', fill = TRUE, col = "white",xlim=c(95,135),ylim=c(25,50))

#pie chart on DAPC
#Check that the sample order is the same for the DAPC and phenotypic objects.
table(ifelse(names(dapc$grp) == indv_info$INDV, "y", "n"))

# Make your output object for pie making.
pie_data <- data.frame(sample_id = names(dapc$grp),
                       popmap = indv_info$popmap,
                       cluster = dapc$grp)

# Order your pie making object by population.
#pie_data <- pie_data[order(pie_data$popmap),]

#make PivotTable of pie_data to creat the variance matrix
attach(pie_data)
CC = tapply(cluster,list(popmap,cluster),length) #length indicates the sum number
CC[is.na(CC)]<-0 #replace NA with 0


#map(add = T, col = "grey90", fill = TRUE)

#create loop:z is the variance matrix, x,y are coordinates.
#debug process:Long should be column 1, Lat should be colume 2, but got it the other way around.
for (i in 1:16){
  add.pie(z = as.integer(CC[i,]), x = coord[i,2], y = coord[i,3],labels = '', 
                      col = myCol,radius = 0.6, solid=.5)}
#  %>% text(x = coord[i,2]+1, y = coord[i,3]+1,data.frame(coord$popmap)[i,])}
# add lables except TA
for (i in c(1:4,6:12,14:16)) {
  text(x = coord[i,2]+0.8, y = coord[i,3]+0.9, " ",
                    labels = data.frame(coord$popmap)[i,], cex=0.8,adj=1,pos = NULL)}
#add TA label
text(x = coord[13,2]+1.2, y = coord[13,3]-0.8, " ",
     labels = data.frame(coord$popmap)[13,], cex=0.8,adj=1,pos = NULL)
#add CIX label
text(x = coord[5,2]-0.3, y = coord[5,3]+0.8, " ",
     labels = data.frame(coord$popmap)[5,], cex=0.8,adj=1,pos = NULL)

#label.dist = 1
##### save plot #####
ppi <- 300
dev.copy(png, "output/03_366_miss9_DAPC_0928_2.png", width=6*ppi, height=6*ppi, res=ppi)
dev.off()
#



