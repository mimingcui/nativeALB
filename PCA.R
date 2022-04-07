library(vcfR)
library(ape)
library(adegenet)
library(ggplot2)
library(ggrepel)

#--- ---------------- read data
indv_info <- read.table("input/pop.txt", header = T)
vcf <- read.vcfR("input/alb.vcf", verbose = FALSE)

x <- vcfR2genlight(vcf)
x
pca=glPca(x, parallel=T,nf=20) #run PCA, number of axes = 20
str(pca)
eig.perc1 <- round(100*pca$eig/sum(pca$eig))#eigen value
my_pca <- cbind.data.frame(pca$scores,indv_info) #add indv coords to pca


#  -----------plot PCAs               
myShape <- c(16,10,17,1,1,16,17,10,16,4,1,4,16,1,10,17)
myColor <- c("#DE544E","#DE544E","#DE544E","#DE544E","#4daf4a","#4daf4a","#377eb8","#377eb8","#377eb8","#377eb8","#377eb8","#984ea3","#984ea3","#984ea3","#984ea3","#984ea3")
#  -----------re-order populations
my_pca$popmap <- factor(my_pca$popmap, levels=c("CIX","BB","TA","JI","YC","QI","HS","SHI","IMC","BJ","CHE","TOL","SHY","CHC","HRB","YJ"))
ggplot(my_pca, aes(PC1,PC2,col=popmap ,shape=popmap))+
     geom_point(size=2, alpha=0.8)+
     xlab(paste0("PC1 (", eig.perc[1], "% variance)")) +
     ylab(paste0("PC2 (", eig.perc[2], "% variance)")) +
     scale_shape_manual(values=myShape)+
     scale_color_manual(values=myColor)+
     theme_bw()+
     theme(legend.position='none',panel.grid= element_blank())
     #theme(legend.title=element_blank(), legend.key = element_blank())
    # theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# ----------- procrustes +  sampling map 
# procrustes transformation on the pca to minimize the squared distance between the pca and the geographic matrix.
library(MCMCpack)
library(maps)
head(my_pca)
p<-procrustes(as.matrix(my_pca[,1:2]),as.matrix(my_pca[,23:22]),translation=TRUE,dilation=TRUE)
p

# Vegan, another package to do procrustes analysis ####
library(vegan)
pro <- procrustes(as.matrix(my_pca[,1:2]),as.matrix(my_pca[,23:22]), symmetric = FALSE)
pro
plot(pro, kind = 1, type = "text")
protest(as.matrix(my_pca[,1:2]),as.matrix(my_pca[,23:22]),scores = "sites", permutations = 9999) 

colnames(p$X.new) <- c("V1","V2")
my_pca <- cbind(my_pca,p$X.new)
head(my_pca)

#-------> try ggplot
library(ggplot2)
library(ggmap)
library(grid)
library(maptools)
library(rgdal)
library(ggpolypath)

####Create background map####
#Google has a Maps Static API, which is being queried under the hood, as you can see from the error message. You need to enable that API. Go to https://console.developers.google.com/ and search for Maps Static API. Then click the "Enable" button.
ggmap::register_google(key = "AIzaSyB_5TBIdSGNltn050C8sa0yEE1vmj4Tork")
has_google_key()

#Get map of China
map=get_googlemap(center=c(lon=117,lat=36), zoom=5, scale=1,
                  style = 'feature:all|element:labels|visibility:off')
ggmap(map)

ggmap(map)+
geom_point(aes(x=V1,y=V2,color=popmap,shape=popmap),data=my_pca,size=2, alpha=0.8)+
xlab("Longitude")+ylab("Latitude")+
scale_shape_manual(values=c(16,10,17,1,1,16,17,10,16,4,1,4,16,1,10,17))+
scale_color_manual(values=c("#DE544E","#DE544E","#DE544E","#DE544E","#4daf4a","#4daf4a","#377eb8","#377eb8","#377eb8","#377eb8","#377eb8","#984ea3","#984ea3","#984ea3","#984ea3","#984ea3"))+
theme_linedraw() + 
theme(axis.title=element_text(face=1),legend.position="none")+
geom_vline(xintercept = p$tt[1,1], color = 'gray', linetype = 2, size = 0.3) +
geom_hline(yintercept = p$tt[2,1], color = 'gray', linetype = 2, size = 0.3) 
#geom_abline(intercept = 169, slope = p$R[1,2]/p$R[1,1], linetype = 3, size = 0.3) +
#geom_abline(intercept = 100.4, slope = p$R[2,2]/p$R[2,1],linetype = 3, size = 0.3) 

#------- sampling map         
location <- read.table(file = "input/coord_20pop.txt", header = TRUE) 
# re-order populations
location$popmap <- factor(location$popmap, levels=c("CIX","BB","TA","JI","YC","QI","HS","SHI","IMC","BJ","CHE","TOL","SHY","CHC","HRB","YJ","KNA","KOR","ULS","INC"))
myShape1=c(16,10,17,6,1,16,15,10,16,4,1,4,16,1,10,17,16,16,16,16)
myColor1=c("#377eb8","#984ea3","#4daf4a","#DE544E","#00BEC3")
  
ggmap(map)+
  geom_point(aes(x=Longitude,y=Latitude,shape=location$popmap,color=Region,size=Size),data=location)+
  #geom_text_repel(aes(x=Longitude,y=Latitude,label=popmap,color=popmap),data=location,size=2.5,force=3,fontface="bold")+
  geom_text_repel(aes(x=Longitude,y=Latitude,label=popmap),force=4,data=location,size=2.5, fontface="bold")+
  theme_bw() + 
  scale_shape_manual(values=myShape1)+
  scale_color_manual(values=myColor1)+
  theme(axis.title=element_text(face=1))
