library(vcfR)
library(ape)
library(adegenet)
library(ggplot2)
library(ggrepel)
# |-------------------------------|#
# |  01. load data                |#
# |-------------------------------|####

#only China
indv_info <- read.table("input/366ALB_China.txt", header = T)
vcf <- read.vcfR("input/03_366_miss9.recode.vcf", verbose = FALSE)

#with Korea
indv_info1 <- read.table("input/403ALB.txt", header = T)
vcf1 <- read.vcfR("input/380alb.recode.vcf", verbose = FALSE)

indv_info1 <- read.table("input/07_rmkoreansiblings_miss9.txt", header = T)
vcf1 <- read.vcfR("input/07_rmkoreansiblings_miss9.recode.vcf", verbose = FALSE)
# |-------------------------------|#
# |  02. run PCA                  |####
# |-------------------------------|#

#only China
x <- vcfR2genlight(vcf)
x
pca=glPca(x, parallel=T,nf=20) #run PCA, number of axes = 20
str(pca)
eig.perc1 <- round(100*pca$eig/sum(pca$eig))#eigen value
my_pca <- cbind.data.frame(pca$scores,indv_info) #add indv coords to pca

#with Korea
x1 <- vcfR2genlight(vcf1)
x1
pca1=glPca(x1, parallel = T,nf=20)
eig.perc1 <- round(100*pca1$eig/sum(pca1$eig))
my_pca1 <- cbind.data.frame(pca1$scores,indv_info1)

# |-------------------------------|#
# |  03. plot PCAs                |####
# |-------------------------------|#

#with Korea
A <- ggplot(my_pca1, aes(x = PC1,y = PC2,col = Country,shape=plot.shape))+
  geom_point(size=2, alpha=0.7)+
  xlab(paste0("PC1 (", eig.perc1[1], "% variance)")) +
  ylab(paste0("PC2 (", eig.perc1[2], "% variance)")) +
  theme_bw()+
  theme(legend.position=c(0.85,0.995), legend.justification=c(0.005,0.995),legend.key.size=unit(0.1,'cm'),legend.key.width=unit(0.1,'cm'),legend.text = element_text(size = 8))+
  #theme(legend.position='none')
  theme(legend.title=element_blank(), legend.key = element_blank(),panel.grid= element_blank())
  #guides(color = guide_colorbar(order = 1),fill = guide_legend(order = 0))
#legend.background = element_rect(fill = 'white', colour = 'grey'),
A
#A + annotate("text", x=15, y=12,label="5241 SNPs \n 403 samples \n 07_403_miss9_rm_Bayescan_fsthet_Outliers")
ggsave('output/A.tiff', A, width = 5, height = 4.8)


#only China
myShape <- c(16,10,17,1,1,16,17,10,16,4,1,4,16,1,10,17)
myColor <- c("#DE544E","#DE544E","#DE544E","#DE544E","#4daf4a","#4daf4a","#377eb8","#377eb8","#377eb8","#377eb8","#377eb8","#984ea3","#984ea3","#984ea3","#984ea3","#984ea3")
# re-order populations
my_pca$popmap <- factor(my_pca$popmap, levels=c("CIX","BB","TA","JI","YC","QI","HS","SHI","IMC","BJ","CHE","TOL","SHY","CHC","HRB","YJ"))
B <- ggplot(my_pca, aes(PC1,PC2,col=popmap ,shape=popmap))+
     geom_point(size=2, alpha=0.8)+
     xlab(paste0("PC1 (", eig.perc[1], "% variance)")) +
     ylab(paste0("PC2 (", eig.perc[2], "% variance)")) +
     scale_shape_manual(values=myShape)+
     scale_color_manual(values=myColor)+
     theme_bw()+
     theme(legend.position='none',panel.grid= element_blank())
     #theme(legend.title=element_blank(), legend.key = element_blank())
    # theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
B
#B+ annotate("text", x=5, y=10,label="5331 SNPs \n 366 samples \n 03_366_miss9_rm_Bayescan_fsthet_Outliers")
ggsave("output/B.tiff",B,width = 5, height = 4.8)

# |--------------------------------|####
# | 04. procrustes +  sampling map |
# |--------------------------------|####
# procrustes transformation on the pca to minimize the squared distance between the pca and the geographic matrix.
library(MCMCpack)
library(maps)
head(my_pca)
p<-procrustes(as.matrix(my_pca[,1:2]),as.matrix(my_pca[,23:22]),translation=TRUE,dilation=TRUE)
p
#translation must be true, dilation can be false
#map(database="world",xlim=c(95,132),ylim=c(22,52))
#map.axes()
#text(p$X.new,col=c(my_pca$Region),labels=my_pca$popmap,cex=0.15)
#text(my_pca$Longitude,my_pca$Latitude,col=c(my_pca$Region),labels=my_pca$popmap,cex=0.45) 


# Vegan, another package to do procrustes analysis ####
library(vegan)
pro <- procrustes(as.matrix(my_pca[,1:2]),as.matrix(my_pca[,23:22]), symmetric = FALSE)
pro
plot(pro, kind = 1, type = "text")
protest(as.matrix(my_pca[,1:2]),as.matrix(my_pca[,23:22]),scores = "sites", permutations = 9999) 
#-----------------------------------------------------#

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
#some issue related to ggmap
#devtools::session_info()
#devtools::install_github("dkahle/ggmap")

####Create background map####
#Google has a Maps Static API, which is being queried under the hood, as you can see from the error message. You need to enable that API. Go to https://console.developers.google.com/ and search for Maps Static API. Then click the "Enable" button.
ggmap::register_google(key = "AIzaSyB_5TBIdSGNltn050C8sa0yEE1vmj4Tork")
has_google_key()

#Get map of China
map=get_googlemap(center=c(lon=117,lat=36), zoom=5, scale=1,
                  style = 'feature:all|element:labels|visibility:off')
ggmap(map)

C <- 
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
C
ggsave("output/C.tiff",C,width = 5, height = 4.8)

# |--------------------------------|####
# | 0X.       sampling map         |
# |--------------------------------|####
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
D
ggsave("output/D.tiff",D,width = 5, height = 4.8)

####################################################################################
mapworld<-borders("world",colour = "gray50",fill="white") 
location1 <- read.table("input/coord_16pop.txt",header = T)
location1[13,2] <- location1[13,2]-0.15
location1[13,3] <- location1[13,3]-0.15
location1[9,2] <- location1[9,2]+0.1
location1[9,3] <- location1[9,3]+0.15
#china <- readOGR("input/gadm36_CHN_shp/gadm36_CHN_1.shp")
xlim <- c(100,135)
ylim <- c(25,50)

mp <- ggplot()+mapworld+coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE)+
      geom_point(aes(x=Longitude,y=Latitude,size=2,color=pool_seq),data=location1)+
      geom_text_repel(aes(x=Longitude,y=Latitude,label=popmap),data=location1,force=1,size=3)+
      theme_linedraw()+
      xlab("Longitude")+ylab("Latitude")+
      theme_linedraw() + 
      theme(axis.title=element_text(face=2),legend.position="none",legend.title=element_blank(), legend.key = element_blank(),panel.grid= element_blank())+
      annotate("text", x=108, y=26,label="Ten pools for ALB pool-seq")

mp
ggsave("output/pool-seq.tiff",mp,width = 5, height = 4.8)

####################################################################################

####put our five plots together, using multiplot()
layout(matrix(1:4,ncol = 2))
layout.show(3)
library(Rmisc)
multiplot(D,A,C,B,layout = matrix(c(1,1,2,2,3,3,4,4), nrow=2, byrow=TRUE))

pi <- 300
dev.copy(tiff, "output/multi_PCA1.tiff", width=6*ppi, height=6*ppi, res=ppi)
dev.off()


#tiff("output/multi_PCA.tiff")
#multiplot(D,C,A,B,layout = matrix(c(1,2,3,4),ncol = 2))
#dev.off()

#X<-ggarrange(D,A,C,B,labels = c("A", "B", "C","D"),ncol = 2, nrow = 2)
#ggexport(X, filename = "output/multi_PCA1.tiff")

#### -------> plot map only show locations in China ####
location1 <- read.table("input/coord_16pop.txt",header = T)
location1[13,2] <- location1[13,2]-0.15
location1[13,3] <- location1[13,3]-0.1
location1[9,2] <- location1[9,2]+0.15
location1[9,3] <- location1[9,3]+0.1
location1
# re-order populations
location1$popmap <- factor(location1$popmap, levels=c("CIX","BB","TA","JI","YC","QI","HS","SHI","IMC","BJ","CHE","TOL","SHY","CHC","HRB","YJ"))

library(ggrepel)
library(ggmap)
E <- ggmap(map)+
    geom_point(aes(x=Longitude,y=Latitude,color=popmap,shape=popmap),data=location1)+
    geom_label_repel(aes(x=Longitude,y=Latitude,label=popmap,color=popmap),force=5,data=location1,size=2.5, fontface="bold")+
    xlab("Longitude")+ylab("Latitude")+theme_linedraw() + 
    scale_shape_manual(values=c(16,10,17,6,1,16,15,10,16,4,1,4,16,1,10,17))+
    scale_color_manual(values=c("#DE544E","#DE544E","#377eb8","#377eb8","#4daf4a","#4daf4a","#377eb8","#377eb8","#377eb8","#377eb8","#377eb8","#984ea3","#984ea3","#984ea3","#984ea3","#984ea3"))+
    theme(axis.title=element_text(face=1),legend.position="none")
E
ggsave("output/E.tiff",E,width = 5, height = 4.8)


##### Mantel test ######
library(ade4)
geo.dists <- dist(cbind(my_pca$Longitude, my_pca$Latitude))
pca.dists <- dist(my_pca$PC1)

pca.dists <- dist(cbind(my_pca$PC1,my_pca$PC2))


as.matrix(geo.dists)[1:5, 1:5]
as.matrix(pca.dists)[1:5, 1:5]

mantel.rtest(geo.dists, my_pca.dists, nrepet = 9999)

#{Monte-Carlo test
#Call: mantel.rtest(m1 = geo.dists, m2 = my_pca.dists, nrepet = 9999)

#Observation: 0.4185552 

#Based on 9999 replicates
#Simulated p-value: 1e-04 
#Alternative hypothesis: greater 

#Std.Obs  Expectation     Variance 
#2.187664e+01 4.742042e-05 3.659702e-04 }



