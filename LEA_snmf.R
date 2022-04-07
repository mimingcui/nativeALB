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
genotype=vcf2geno("input/alb.vcf","alb.geno")
obj.snmfff = snmf(genotype, K = 1:10, entropy = T, repetitions = 10, ploidy = 2, project = "new")

plot(obj.snmfff,cex = 1, col = "lightblue", pch = 19)
# show the project
show(obj.snmfff)
# summary of the project
summary(obj.snmfff)

#|---------------------------|####
#|     first try K = 6   |#
#|---------------------------|####
obj.snmf=snmf("alb.geno",K=6,alpha=100,project="new")
qmatrix= Q(obj.snmf, K=6)

myCol <- c("#762a83","#af8dc3","#e7d4e8","#d9f0d3","#7fbf7b","#1b7837")

# add population info
indv_info <- read.table("input/pop.txt", sep="\t", header = T)
qmatrix_addpop <- cbind(qmatrix,indv_info)   

### Import individual order
o <- read.table("input/order.txt",sep="\t",header=F)
colnames(o) <- 'popmap'

qmatrixPO <- qmatrix_addpop[order(match(qmatrix_addpop$popmap, o$popmap)), ] #don't forget o$popmap is as a factor
qmatrixPO1 <- qmatrixPO[,c(1,2,3,4,5,6)] #adjust vector!!

## ---------- ggplot --------##
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

library(mapplots)
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

#next to plot K=7 and K=8         
