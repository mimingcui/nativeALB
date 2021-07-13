#-------------------------------------------------------
#                        01. het
#-------------------------------------------------------
het <- read.table("input/07_403_miss9.het",header = T)
#07_403_miss9_rm_Bayescan_fsthet_Outliers.het
head(het)
summary(het[which(het$Population=='CHC'),]) #summary stats for pops

#box plot, not appropriate when only few indvs
het$O_Het <- 1-(het$O.HOM./het$N_SITES)
het$E_Het <- 1-(het$E.HOM./het$N_SITES)

pop <- read.table("input/403ALB.txt",header = T)
het$Population <- pop$pop_id #in the file pop,there's title 'pop', that's why 'pop' covered 'Population'
het$Population <- factor(het$Population)
het$Country <- pop$Country
plot(O_Het~Population,data = het)


#try ggplot
library(ggplot2)
ggplot(het,aes(x = Population, y = O_Het, col = Population))+
  geom_point()+
  theme_bw()

ggplot(het,aes(x = Population, y = E_Het, col = Population))+
  geom_point()+
  theme_bw()

#re-order the x-axis
positions<- c("KNA","KOR","ULS","INC","BB","CIX","JI","TA","BJ","CHE","HS","IMC","SHI","HRB","CHC","SHY","TOL","YJ","QI","YC")
het$Population<-factor(het$Population,levels = positions)
  
ggplot(het,aes(x = Population, y = O_Het, col = Population))+ #col = Country
     scale_x_discrete(limits=positions)+
     geom_point()+
     theme_bw()+
     theme(axis.text.x=element_text(angle = 30),legend.position = "none")+
     ggtitle("Heterozygosity")

library(forcats)
library(dplyr)
het$Country<-factor(het$Country,levels=c("Korea","China"))
Ho<-
  het %>%
  mutate(Population=fct_reorder(Population,O_Het,mean)) %>% #fct_reorder(Population,O_Het) default by median
  ggplot(aes(x =Population,y = O_Het)) + 
     geom_jitter(col = "grey",width = 0.1)+
     geom_boxplot(notch = F, alpha = 0,outlier.shape = NA)+ 
     #geom_violin(alpha = 0)+
     #scale_x_discrete(limits=positions)+ 
     theme_bw()+
     xlab("")+
     ylab("Heterozygosity")+
     facet_grid(. ~ Country,space="free_x", scales="free_x")+
     theme( panel.grid= element_blank(), axis.text.x=element_text(angle = 30),legend.position = "none")
#same as ggplot(het,aes(x =fct_reorder(Population,O_Het),y = O_Het))+...
ggsave("output/Het_by_mean.png",Ho)
#only draw observed heterozygosity, because dkw expected one seems strange.
#use scale_color_manual(values = c(...)) to set color manually.


{#### 02. plot pi, not done ####
pi_IMC <- read.table("input/07_403_miss9_pi/IMC.windowed.pi",header = T)
head(pi_IMC)
pi_IMC$Population <- c(rep("IMC",1821)
###repeat things

#first_function <- function(x){read.table(paste0(x,".windowed.pi"))}
#extract all_20_files_path_names
x<- list.files("input/07_403_miss9_pi",pattern="*.windowed.pi", full.names=TRUE) #or dir()
all_20_files <- lapply(x,function(x){read.table(x)[-1,]}) #for (i in x){read.table(i)}
#apply (x,header = F)[-1,] to not read the title because when rbind,the title would also count as a row.
#pi_all <- lapply(all_20_files, rbind) #return to list
pi_all <- do.call(rbind,all_20_files) #return to data.frame that we wanted.
#pi_all = do.call(rbind, lapply(files, function(x)read.table(x, stringsAsFactors = FALSE)))
colnames(pi_all)<- c("CHROM","BIN_START","BIN_END","N_VARIANTS","PI")

x #see the order
sapply(all_20_files, nrow)#to see how many rows for each.
pi_all$Population <- c(rep("BB",1840),rep("BJ",1835),rep("CHC",1647),rep("CHE",1876),rep("CIX",1785),rep("HRB",1754),
                       rep("IMC",1821),rep("HS",1889),rep("JI",1815),rep("KNA",384),rep("KOR",609),rep("QI",1891),rep("SHI",1884),
                       rep("SHY",1847),rep("TA",1837),rep("TOL",1868),rep("YC",1897),rep("YJ",1896))
head(pi_all)
#a weird thing happened, seems to relate to the file format. anyway after I write out and read, it works.
write.csv(pi_all,"pi_all_2.csv")
pi_all_2 <- read.csv('pi_all_2.csv')
head(pi_all_2)
B <- ggplot(pi_all_2,aes(x = Population, y = PI, col = Population)) +
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("")
B

B <- ggplot(pi_all_2,aes(x = Population, y = PI, col = Population))+
  geom_jitter(col = "grey",width = 0.1)+ #Add this before the geom_boxplot() to put the points behind the boxes
  geom_boxplot(notch = T, alpha = 0,outlier.shape = NA)+ #alpha makes the background seethrough
  theme_bw()+
  theme(legend.position = "none")+
  xlab("")+
  ylab(expression(pi))
B
ggsave("output/07_403_miss9_pi.png")}
