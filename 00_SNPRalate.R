if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPRelate")
library(gdsfmt)
library(SNPRelate)

#run this each time when re-generate the GDS file
showfile.gds(closeall=TRUE)

#convert vcf to GDS file "test.gds". Vcf file should omit 'Scaffold' in order to run the function snpgdsIBDSelection()
snpgdsVCF2GDS("input/07_403_miss9_rmScaffold_for_snprelate.vcf", "test.gds", method="biallelic.only")
snpgdsSummary("test.gds")

# Open a GDS file
genofile <- snpgdsOpen("test.gds")

# Get sample id
sample.id<- read.gdsn(index.gdsn(genofile, "sample.id"))

# Get population information
poplist <- read.table("input/403ALB.txt", header = T)
head(poplist)
pop_code <- poplist[ ,5]

###extract each populaiton
#CHE.id <- sample.id[pop_code == "CHE"]
#ibd.coeff <- snpgdsIBDSelection(snpgdsIBDMoM(genofile, sample.id=CHE.id, snp.id=genofile$snp.id,
#                                             maf=0.05, missing.rate=0.05, num.thread=2))
#### repeat things using function() ####
x <- unique(pop_code)
x
x <- as.vector(x)

# define a function
#warm up with first_function <- function(x){sample.id[pop_code==x]}
first_function <- function(x){snpgdsIBDSelection(snpgdsIBDMoM(genofile, sample.id=sample.id[pop_code==x], snp.id=genofile$snp.id,
                                                              maf=0.05, missing.rate=0.05, num.thread=2))}
id <- lapply(x, first_function)

#repeating stuff, write.csv(id[],"coeff.csv")
file=list()
for (i in 1:20){
  files = paste0("coeff",i,".csv")
  write.csv(id[i],files)}

library("rlist")
ids <- list.stack(id)
write.csv(ids,"output/all_coeff.csv")

#build network
library(networkD3)
cut <- ids[ids$kinship>=0.15,] #cut by colume 5th(kniship)>=0.22
network <- cut[,1:2]
simpleNetwork(network,fontSize=15)

#build network for poplation KOR
KOR <- id[13]
KOR <- as.data.frame(KOR)
KOR <- KOR[,1:2]
head(KOR)
simpleNetwork(KOR,fontSize=10)



###########33 Estimating IBD Using PLINK method of moments (MoM) ###################
# Estimate IBD coefficients
ibd <- snpgdsIBDMoM(genofile, sample.id=CHE.id, snp.id=genofile$snp.id,
                    maf=0.05, missing.rate=0.05, num.thread=2)
#IBD analysis (PLINK method of moment) on genotypes:
#Excluding 0 SNP on non-autosomes
#Excluding 3,187 SNPs (monomorphic: TRUE, MAF: 0.05, missing rate: 0.05)
#Working space: 15 samples, 537 SNPs
#using 2 (CPU) cores
#PLINK IBD:    the sum of all selected genotypes (0,1,2) = 10299

# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="CHE samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)
#####
######### Estimating IBD Using Maximum Likelihood Estimation (MLE)#################3333
# Estimate IBD coefficients
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
snpset.id <- unlist(unname(snpset))
head(snpset.id)
snp.id=snpset.id

set.seed(100)
snp.id <- sample(snp.id, 1500)  # random 1500 SNPs
ibd <- snpgdsIBDMLE(genofile, sample.id=CHE.id, snp.id=snp.id,
                    maf=0.05, missing.rate=0.05, num.thread=2)
#Identity-By-Descent analysis (MLE) on genotypes:
#Excluding 2,224 SNPs (non-autosomes or non-selection)
#Excluding 1,276 SNPs (monomorphic: TRUE, MAF: 0.05, missing rate: 0.05)
#Working space: 15 samples, 224 SNPs
#using 2 (CPU) cores
#MLE IBD:    the sum of all selected genotypes (0,1,2) = 4274
#MLE IBD:	Mon Jan 13 12:59:53 2020	0%
#MLE IBD:	Mon Jan 13 12:59:54 2020	100%

# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="CHE samples (MLE)")
lines(c(0,1), c(1,0), col="red", lty=2)
#############################################################################


genofile <- snpgdsOpen(snpgdsExampleFileName())
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- read.gdsn(index.gdsn(genofile, path="sample.annot/pop.group"))
YRI.id <- sample.id[pop_code == "YRI"]
########







