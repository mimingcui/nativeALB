#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("LEA")

#below help import input files from the structure format, and to display nice geographic representations of ancestry coefficients with maps.
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

#|---------------------------|#
#|--->  evaluate how many K  |#
#|---------------------------|#####
library(LEA)
genotype=vcf2geno("input/03_366_miss9.recode.vcf","03_366_miss9.geno")
obj.snmfff = snmf(genotype, K = 1:10, entropy = T, repetitions = 10, ploidy = 2, project = "new")
#repetitions = 10 can take long and it doesn't make much difference
tiff("output/03_366_miss9_snmf_cross-entropy_rep10.tiff",res = 300, width = 10,height = 10,units = 'cm')
plot(obj.snmfff,cex = 1, col = "lightblue", pch = 19)
dev.off()
# show the project
show(obj.snmfff)
# summary of the project
summary(obj.snmfff)

#|---------------------------|####
#|--->     first try K = 6   |#
#|---------------------------|####

obj.snmf=snmf("03_366_miss9.geno",K=6,alpha=100,project="new")
qmatrix= Q(obj.snmf, K=6)

#require("RColorBrewer")
#brewer.pal(8, "")
myCol <- c("#762a83","#af8dc3","#e7d4e8","#d9f0d3","#7fbf7b","#1b7837")

#barplot(t(qmatrix),col = myCol,border = NA,space = 0,xlab = "Individuals", ylab = "Admixture coefficients")
# add population info
indv_info <- read.table("input/366ALB_China.txt", sep="\t", header = T)
qmatrix_addpop <- cbind(qmatrix,indv_info)   
head(qmatrix_addpop)
### Import individual order
o <- read.table("input/LEA_barplot_366order.txt",sep="\t",header=F)
colnames(o) <- 'popmap'
head(o)
#qmatrixPO <- merge(o, qmatrix_addpop, by.x = 'popmap', by.y = 'popmap',all.x=T)
qmatrixPO <- qmatrix_addpop[order(match(qmatrix_addpop$popmap, o$popmap)), ] #don't forget o$popmap is as a factor
qmatrixPO1 <- qmatrixPO[,c(1,2,3,4,5,6)] #adjust vector!!

tiff("output/03_366_miss9_snmf_barplot_K=6_0928.tiff", res = 300, width = 20,height = 10,units = 'cm')
barplot(t(qmatrixPO1), col = myCol, border = NA, space = 0,xaxt= "n",
        xlab = "", ylab = "Admixture coefficients")
dev.off()

## ---------- could also try ggplot --------|##
# Convert data to long-form with 'melt' from the reshape2 package. 
library(reshape)
library(ggplot2)
# melt for puting in one column
melt_qmatrix = melt(qmatrixPO, id.vars=c("popmap", "INDV"), 
            measure.vars=c(1,2,3,4,5,6)) 

# I figured it out! NEED to sort first! The issue of 'ggplot2 determines order'.
melt_qmatrix$INDV <- factor(melt_qmatrix$INDV, levels = as.vector(qmatrixPO$INDV))

melt_qmatrix$popmap <- factor(melt_qmatrix$popmap, levels=c("YJ","HRB","CHC","SHY","TOL","CHE","BJ","IMC","SHI","QI","YC","HS","JI","TA","BB"))

p1 <- ggplot(melt_qmatrix, aes(x=INDV, y=value, fill=variable))+
      theme_void() + 
      geom_bar(position="stack", stat="identity",width=2) +
      scale_fill_manual(values=myCol)+
      theme(legend.position='none')+
      facet_grid(. ~ popmap, space="free_x", scales="free_x") 
           #labeller=label_both)
# try interchange INDV and popmap 

p1 <- p1+theme(panel.spacing = unit(-0.1, "lines")) # set width to adjust spaces inbetween the factor bars

k6 <- ggplot(melt_qmatrix, aes(x=INDV, y=value, fill=variable))+
  theme_void() + 
  geom_bar(position="stack", stat="identity",width=2) +
  ylab(paste0("K=6")) +
  scale_fill_manual(values=myCol)+
  theme(legend.position='none',panel.grid= element_blank())+
  theme(axis.title.x=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())

k6

#check the x y range
xrange1<-layer_scales(p1)$x$range$range
xrange2<-layer_scales(p2)$x$range$range

#|---------------------------|####
#|--->   pie chart on a map  |####
#|---------------------------|####
pop=indv_info[6]
K=6
Npop = lengths(unique(pop))#length&lengths--row&column
Npop #be careful!!!should be 16 in this case,but alway show as 1,sometimes is **,but doesn't matter

#create a Structure matrix with 16 rows, K columns
qpop = matrix(NA, ncol = K, nrow = Npop)
#create a Coord matrix with 16 rows, 2 columns
coord.pop = matrix(NA, ncol = 2, nrow = Npop)

coord<-indv_info[,c(2,3)] #pull out long&latitude
coord<-coord[,c(2,1)] #switch order, long should always be on left

#unify data format
is.matrix(qmatrix)
is.matrix(pop)
pop <- as.matrix(pop)
is.data.frame(coord)
qmatrix <- as.data.frame(qmatrix)

for (i in unique(pop)){
  qpop[unique(pop)==i,] = apply(qmatrix[pop == i,], 2, mean)
  coord.pop[unique(pop)==i,] = apply(coord[pop == i,], 2, mean)}

###seperately troubleshoot###
#row.names(qpop)<-c(1:Npop)
#for (i in unique(pop))
#  qpop[unique(pop)==i,] = apply(qmatrix[pop == i,], 2, mean)
#for (i in unique(pop))
#  coord.pop[unique(pop)==i,] = apply(coord[pop == i,,drop=F], 2, mean)
###apply 'drop=F' when error'dim(X) must have a positive length'show up
library(maps)
library(mapplots)
#plot(coord, xlab = "Longitude", ylab = "Latitude", type = "n") 

tiff("output/03_366_miss9_snmf_piechart_K=6_0928.tiff", res = 300, width = 15,height = 15,units = 'cm')
par(mfrow=c(1,1),oma=c(0,0,0,0))#xiazuoshangyou
plot(coord, xlim=c(104,130), ylim=c(29.5,46.5),ann = F,xaxt = "n", yaxt ="n")#ann:whether to show x,y-axis
#title(xlab= 'Longitude', ylab = 'Latitude', line = 2)
map(add = T,xlim=c(104,130), ylim=c(29.5,46.5))
#map(add = T, col = "white", fill = TRUE,xlim=c(90,135), ylim=c(20,55))
for (i in 1:Npop){
  add.pie(z = qpop[i,],x = coord.pop[i,1],y = coord.pop[i,2],labels = "",col = myCol,radius=1)}
###add labels
for (i in c(1:3,6,8:10,12:13,16)) {
  text(x = coord.pop[i,1]-0.2, y = coord.pop[i,2]+1.35, " ",
       labels = data.frame(unique(pop))[i,], cex=1.2,adj=1,pos = NULL)}
for (i in 5) {
  text(x = coord.pop[i,1]+2, y = coord.pop[i,2]-1, " ",
       labels = data.frame(unique(pop))[i,], cex=1.2,pos = NULL,adj=NULL)}
for (i in c(4,11)) {
  text(x = coord.pop[i,1]+2, y = coord.pop[i,2]+1, " ",
       labels = data.frame(unique(pop))[i,], cex=1.2,pos = NULL,adj=NULL)}
dev.off()

#inset
par(mfrow=c(1,1),oma=c(0,0,0,0))#xiazuoshangyou
plot(coord, xlim=c(114,118), ylim=c(36,39),ann = F,bty = "n",xaxt = "n", yaxt ="n")
#%>% map(add = T)#ann:whether to show x,y-axis
map(add = T)
for (i in 1:Npop){
  add.pie(z = qpop[i,],x = coord.pop[i,1],y = coord.pop[i,2],labels = "",col = myCol,radius=0.22)}
for (i in c(14:15)){
  text(x = coord.pop[i,1]-0.3, y = coord.pop[i,2]+0.22, " ",
       labels = data.frame(unique(pop))[i,], cex=3.5,adj=1,pos = NULL)}
for (i in 7){
  text(x = coord.pop[i,1]+0.7, y = coord.pop[i,2]-0.2, " ",
       labels = data.frame(unique(pop))[i,], cex=3.5,adj=1,pos = NULL)}
#|---------------------------|####
#|      plot K=7 K=8         |#-------->
#|---------------------------|####

#------ k=7 --------####
obj.snmf=snmf("03_366_miss9.geno",K=7,alpha=100,project="new")
qmatrix= Q(obj.snmf, K=7)
myCol <- c("#762a83","#af8dc3","#e7d4e8","#d9f0d3","#7fbf7b","#1b7837","#ADD8E6")
qmatrix_addpop <- cbind(qmatrix,indv_info)   
#qmatrixPO <- merge(o, qmatrix_addpop, by.x = 'popmap', by.y = 'popmap',all.x=T)
qmatrixPO <- qmatrix_addpop[order(match(qmatrix_addpop$popmap, o$popmap)), ] #don't forget o$popmap is as a factor
qmatrixPO1 <- qmatrixPO[,c(1,2,3,4,5,6,7)]#adjust vector!!

tiff("output/03_366_miss9_snmf_barplot_K=7_0928.tiff", res = 300, width = 20,height = 10,units = 'cm')
barplot(t(qmatrixPO1), col = myCol, border = NA, space = 0,xaxt= "n",
        xlab = "", ylab = "Admixture coefficients")
dev.off()

### try ggplot,it doesn't follow my order! need to sort
# Convert data to long-form with 'melt' from the reshape2 package. 
# melt for puting in one column
melt_qmatrix = melt(qmatrixPO, id.vars=c("popmap", "INDV"), 
                    measure.vars=c(1,2,3,4,5,6,7)) 

# I figured it out! NEED to sort first! The issue of 'ggplot2 determines order'.
melt_qmatrix$INDV <- factor(melt_qmatrix$INDV, levels = as.vector(qmatrixPO$INDV))

k7 <- ggplot(melt_qmatrix, aes(x=INDV, y=value, fill=variable))+
  theme_void() + 
  geom_bar(position="stack", stat="identity",width=2) +
  ylab(paste0("K=7")) +
  scale_fill_manual(values=myCol)+
  theme(legend.position='none',panel.grid= element_blank())+
  theme(axis.title.x=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
k7

#|---------------------------|####
#|     pie chart on a map    |####----->
#|---------------------------|####
pop=indv_info[6]
K=7
Npop = lengths(unique(pop))#length&lengths--row&column

qpop = matrix(NA, ncol = K, nrow = Npop)
coord.pop = matrix(NA, ncol = 2, nrow = Npop)

coord<-indv_info[,c(2,3)] #pull out long&latitude
coord<-coord[,c(2,1)] #switch order, long should always be on left

#unify data format
is.matrix(qmatrix)
is.matrix(pop)
pop <- as.matrix(pop)
is.data.frame(coord)
qmatrix <- as.data.frame(qmatrix)

for (i in unique(pop)){
  qpop[unique(pop)==i,] = apply(qmatrix[pop == i,], 2, mean)
  coord.pop[unique(pop)==i,] = apply(coord[pop == i,], 2, mean)}

tiff("output/03_366_miss9_snmf_piechart_K=7_0928.tiff", res = 300, width = 15,height = 15,units = 'cm')
par(mfrow=c(1,1),oma=c(0,0,0,0))#xiazuoshangyou
plot(coord, xlim=c(104,130), ylim=c(29.5,46.5),ann = F,xaxt = "n", yaxt ="n")#ann:whether to show x,y-axis
map(add = T,xlim=c(104,130), ylim=c(29.5,46.5))
for (i in 1:Npop){
  add.pie(z = qpop[i,],x = coord.pop[i,1],y = coord.pop[i,2],labels = "",col = myCol,radius=1)}
###add labels
for (i in c(1:3,6,8:10,12:13,16)) {
  text(x = coord.pop[i,1]-0.2, y = coord.pop[i,2]+1.35, " ",
       labels = data.frame(unique(pop))[i,], cex=1.2,adj=1,pos = NULL)}
for (i in 5) {
  text(x = coord.pop[i,1]+2, y = coord.pop[i,2]-1, " ",
       labels = data.frame(unique(pop))[i,], cex=1.2,pos = NULL,adj=NULL)}
for (i in c(4,11)) {
  text(x = coord.pop[i,1]+2, y = coord.pop[i,2]+1, " ",
       labels = data.frame(unique(pop))[i,], cex=1.2,pos = NULL,adj=NULL)}
dev.off()

#inset
par(mfrow=c(1,1),oma=c(0,0,0,0))#xiazuoshangyou
plot(coord, xlim=c(114,118), ylim=c(36,39),ann = F,bty = "n",xaxt = "n", yaxt ="n")
#%>% map(add = T)#ann:whether to show x,y-axis
map(add = T)
for (i in 1:Npop){
  add.pie(z = qpop[i,],x = coord.pop[i,1],y = coord.pop[i,2],labels = "",col = myCol,radius=0.22)}
for (i in c(14:15)){
  text(x = coord.pop[i,1]-0.3, y = coord.pop[i,2]+0.22, " ",
       labels = data.frame(unique(pop))[i,], cex=3.5,adj=1,pos = NULL)}
for (i in 7){
  text(x = coord.pop[i,1]+0.7, y = coord.pop[i,2]-0.2, " ",
       labels = data.frame(unique(pop))[i,], cex=3.5,adj=1,pos = NULL)}

####------ k=8 --------####
obj.snmf=snmf("03_366_miss9.geno",K=8,alpha=100,project="new")
qmatrix= Q(obj.snmf, K=8)
myCol <- c("#762a83","#af8dc3","#e7d4e8","#d9f0d3","#7fbf7b","#1b7837","#ADD8E6","#FFDAB9")
qmatrix_addpop <- cbind(qmatrix,indv_info)   

#qmatrixPO <- merge(o, qmatrix_addpop, by.x = 'popmap', by.y = 'popmap',all.x=T)
qmatrixPO <- qmatrix_addpop[order(match(qmatrix_addpop$popmap, o$popmap)), ] #don't forget o$popmap is as a factor
qmatrixPO1 <- qmatrixPO[,c(1,2,3,4,5,6,7,8)]#adjust vector!!

tiff("output/03_366_miss9_snmf_barplot_K=8_0928.tiff", res = 300, width = 20,height = 10,units = 'cm')
barplot(t(qmatrixPO1), col = myCol, border = NA, space = 0,xaxt= "n",
        xlab = "", ylab = "Admixture coefficients")
dev.off()

#ggplot
melt_qmatrix = melt(qmatrixPO, id.vars=c("popmap", "INDV"), 
                    measure.vars=c(1,2,3,4,5,6,7,8)) 

# I figured it out! NEED to sort first! The issue of 'ggplot2 determines order'.
melt_qmatrix$INDV <- factor(melt_qmatrix$INDV, levels = as.vector(qmatrixPO$INDV))

k8 <- ggplot(melt_qmatrix, aes(x=INDV, y=value, fill=variable))+
  theme_void() + 
  geom_bar(position="stack", stat="identity",width=2) +
  ylab(paste0("K=8")) +
  scale_fill_manual(values=myCol)+
  theme(legend.position='none',panel.grid= element_blank())+
  theme(axis.title.x=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
k8
#|---------------------------|####
#|     pie chart on a map    |####----->
#|---------------------------|####
pop=indv_info[6]
K=8
Npop = lengths(unique(pop))#length&lengths--row&column

qpop = matrix(NA, ncol = K, nrow = Npop)
coord.pop = matrix(NA, ncol = 2, nrow = Npop)

coord<-indv_info[,c(2,3)] #pull out long&latitude
coord<-coord[,c(2,1)] #switch order, long should always be on left

#unify data format
is.matrix(qmatrix)
is.matrix(pop)
pop <- as.matrix(pop)
is.data.frame(coord)
qmatrix <- as.data.frame(qmatrix)

for (i in unique(pop)){
  qpop[unique(pop)==i,] = apply(qmatrix[pop == i,], 2, mean)
  coord.pop[unique(pop)==i,] = apply(coord[pop == i,], 2, mean)}

tiff("output/03_366_miss9_snmf_piechart_K=8_0928.tiff", res = 300, width = 15,height = 15,units = 'cm')
par(mfrow=c(1,1),oma=c(0,0,0,0))#xiazuoshangyou
plot(coord, xlim=c(104,130), ylim=c(29.5,46.5),ann = F,xaxt = "n", yaxt ="n")#ann:whether to show x,y-axis
map(add = T,xlim=c(104,130), ylim=c(29.5,46.5))
for (i in 1:Npop){
  add.pie(z = qpop[i,],x = coord.pop[i,1],y = coord.pop[i,2],labels = "",col = myCol,radius=1)}
###add labels
for (i in c(1:3,6,8:10,12:13,16)) {
  text(x = coord.pop[i,1]-0.2, y = coord.pop[i,2]+1.35, " ",
       labels = data.frame(unique(pop))[i,], cex=1.2,adj=1,pos = NULL)}
for (i in 5) {
  text(x = coord.pop[i,1]+2, y = coord.pop[i,2]-1, " ",
       labels = data.frame(unique(pop))[i,], cex=1.2,pos = NULL,adj=NULL)}
for (i in c(4,11)) {
  text(x = coord.pop[i,1]+2, y = coord.pop[i,2]+1, " ",
       labels = data.frame(unique(pop))[i,], cex=1.2,pos = NULL,adj=NULL)}
dev.off()

#inset
par(mfrow=c(1,1),oma=c(0,0,0,0))#xiazuoshangyou
plot(coord, xlim=c(114,118), ylim=c(36,39),ann = F,bty = "n",xaxt = "n", yaxt ="n")
#%>% map(add = T)#ann:whether to show x,y-axis
map(add = T)
for (i in 1:Npop){
  add.pie(z = qpop[i,],x = coord.pop[i,1],y = coord.pop[i,2],labels = "",col = myCol,radius=0.22)}
for (i in c(14:15)){
  text(x = coord.pop[i,1]-0.3, y = coord.pop[i,2]+0.22, " ",
       labels = data.frame(unique(pop))[i,], cex=3.5,adj=1,pos = NULL)}
for (i in 7){
  text(x = coord.pop[i,1]+0.7, y = coord.pop[i,2]-0.2, " ",
       labels = data.frame(unique(pop))[i,], cex=3.5,adj=1,pos = NULL)}

#### combine multiple plots ####
library(Rmisc)
multiplot(p1,k7,k8, layout = matrix(c(1,2,3), nrow=3, byrow=TRUE))
