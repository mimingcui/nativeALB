  library(dplyr)
  library(magrittr)
  library(tibble)
  library(ggplot2)
  library(adegenet)
  library(vcfR)
  library(vegan)

# load specific useful functions from source
'%ni%' <- Negate('%in%') # reverse of %in%

#-----------------------------------------------------------------------------
#                     01. Load data
#-----------------------------------------------------------------------------
vcf <- vcfR::read.vcfR("input/alb.vcf")

# load population map
strata <- read.table("input/pop.txt", h=T) %>% set_colnames(., c('IND', 'Lat','Long','City','poid','CODE','Region','Country','other')) 

# convert to genind
genind <- vcfR2genind(vcf)
pop(genind) <- as.vector(strata$CODE)

#-----------------------------------------------------------------------------
#                     02. Make connection graph
#-----------------------------------------------------------------------------

#Make df of Coordiantes for  all samples
#use jitter to switch slightly samples positions for connection graph
Coords.samples <- mutate(strata, Lat_jit = jitter(Lat,5), Long_jit = jitter(Long,5)) %>% group_by(CODE) # do(sample_n(., 20)) %>%
Coords.samples <- as.data.frame(ungroup(Coords.samples, CODE))
# ?why group and then followed by ungroup 
head(Coords.samples)
CN <- chooseCN(Coords.samples[,c(10,11)], type=1,d1=NULL,d2=NULL, plot=F, res="listw", ask=F) #Delaunay
CN2 <- chooseCN(Coords.samples[,c(10,11)], type=2,d1=NULL,d2=NULL, plot=F, res="listw", ask=F) #Gabriel graph

par(mar=c(1,1,1,1))
plot(CN, Coords.samples[,c(10,11)],pch=19)
par(font=2, ps=10); text(Coords.samples$CODE, x=Coords.samples$Lat-0.1, y=Coords.samples$Long+0.7, col='blue')

#-----------------------------------------------------------------------------
#                     03. Run sPCA
#-----------------------------------------------------------------------------
CN.spca <- spca(genind, xy=as.matrix(Coords.samples[,10:11]), cn=CN,scannf=FALSE, nfposi=3,nfnega=2)
class(CN.spca)

##---------------------> plots
plot(CN.spca$eig/100,type = "h",lwd = 10,col=rep(c("black", rgb(0, 0, 0, 20, maxColorValue=255)),c(1, 380)),
     bty='l',cex.axis=2,ann = F)
title(xlab="PC", ylab="Eigenvalue", line=2.5,font.lab=2,cex.lab=2)
abline(h = 0, col = "grey")


summary(CN.spca)
plot.spca(CN.spca)
screeplot.spca(x=CN.spca)

# plot s.value with jitter effect for samples
s.value(dfxy = cbind(jitter(Coords.samples[,2],50), jitter(Coords.samples[,2],450)), z = CN.spca$li[,5],
        csize = 0.25, xlim = c(-66,-59), ylim=c(44,50))

# resume individual sPCA info per location
head(CN.spca$li)
sPCA1_pop <- CN.spca1$li %>% mutate(., POP = Coords.samples$CODE) %>% group_by(POP) %>% summarise(., sPC1 = mean(`Axis 1`))
sPCA2_pop <- CN.spca1$li %>% mutate(., POP = Coords.samples$CODE) %>% group_by(POP) %>% summarise(., sPC2 = mean(`Axis 2`))

#load continental shapefile
library(rgdal)
coastline <- readOGR("input/gadm36_CHN_shp/gadm36_CHN_1.shp")
class(coastline)
xlim <- c(90,135)
ylim <- c(20,55)
#area <- crop(x=coastline, y = extent(-66.5,-59,45,49))
plot.new()
par(mar=c(1,1,1,1))
plot(coastline,col="grey95", bg="white", border="grey60")

map('world', fill = TRUE, col = "grey",xlim=c(75,135),ylim=c(15,55))

#plot s.value squares
s.value(dfxy = cbind(unique(Coords.samples[,3]), unique(Coords.samples[,2])), z = sPCA1_pop$sPC1,
        csize = 0.5, add.plot = TRUE, possub = 'top')
par(font=2);text(Coords.samples$CODE, x=(Coords.samples$Long+0.30), y=Coords.samples$Lat, cex=0.7)


#Axis 2
s.value(dfxy = cbind(unique(Coords.samples[,3]), unique(Coords.samples[,2])), z = sPCA2_pop$sPC2,
        csize = 0.5, xlim = c(90,135), ylim=c(20,55), add.plot = TRUE)
#Color plot
colorplot(xy = Coords.samples[,c(3,2)], X=CN.spca1$ls, axes=c(1,2), transp=TRUE, cex=2)#add.plot = TRUE

# |-----------------|####
# |  interpolating  | ===========> 
# |-----------------|####
library(akima)
library(maps)
x <- Coords.samples[,11]
y <- Coords.samples[,10]
temp <- interp(x, y, CN.spca$li[,1])
image(temp, col=azur(100))
points(x,y)

myPal <- colorRampPalette(c("firebrick2", "white", "lightslateblue"))
annot <- function(){
title(main="sPCA - interpolated map of individual scores",xlab = "Longtitude", ylab = "Latitude")
points(x,y,)}

filled.contour(temp$x, temp$y, temp$z, color.pal=myPal, nlev=50,
               xlim = c(105,131),
               ylim = c(29.5,46.5),
               key.title=title("lagged \nscore 1"), plot.title=annot(),
               plot.axes={axis(1);axis(2);points(x,y);map(add=T);
                 for (i in c(1:4,6:16)) {
                   text(x = coord[i,2]-0.1, y = coord[i,3]+0.25, " ",
                        labels = data.frame(coord$popmap)[i,],font=2,cex=0.7,adj=1,pos = NULL)};
                 text(x = coord[5,2]+0.7, y = coord[5,3]+0.3, " ",
                      labels = data.frame(coord$popmap)[5,],font=2,cex=0.7,adj=1,pos = NULL)}
)


#=======> Test significiance sPCA (global and local)
#test for significativity
myGtest <- global.rtest(CN.spca$tab,CN.spca$lw,nperm=999)
myGtest
plot(myGtest)

myLtest <- local.rtest(CN.spca$tab,CN.spca$lw,nperm=999)
myLtest
plot(myLtest)



