###############################################################################
# Part One: Prepare R Environment & Load Phenotypic Data
# (1.0) Clear Global R
rm(list = ls())
# (1.1) Load packages required for analysis
chooseCRANmirror(ind=71)
# Load packages
Pckg.Lst <-c("polyRAD","VariantAnnotation","qqman", "tidyverse")
package.check <- lapply(
  Pckg.Lst,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)}})

# (1.2) Add option so images can be written to file
options(bitmapType='cairo')
###############################################################################
# Part Two: Load VCF, Check it, and filter it.
# (2.0) Clear Global R
mydata <- VCF2RADdata("../Analysis_8/GenotypeFile.F3.vcf.gz",
                      refgenome= "ReferenceGenome/DCARv3d.fna",
                      possiblePloidies=list(2))

# (2.1) Remove High Depth Loci
mydata <- RemoveHighDepthLoci(mydata)
mydata <- MergeRareHaplotypes(mydata, min.ind.with.haplotype = 53)
mydata <- MergeIdenticalHaplotypes(mydata)

# (2.2) Quick Check of Loci Number and Alleles:
nLoci(mydata)
nAlleles(mydata)
nAlleles(mydata)/nLoci(mydata)
# (2.3) run test for mendelian segregation
hh <- HindHe(mydata)
# (2.4) identify total depth (per sample)
TotDepthT <- rowSums(mydata$locDepth)
# (2.5) get the avearge Hind/He per sample
# heterozygosity based QC measured based on individua
hhByInd  <-  rowMeans(hh, na.rm = TRUE)
# (2.6) Make sure there is not a relationship between derpth & He
pdf("Hind.HevDerpth.pdf")
plot(TotDepthT, hhByInd, log= "x", xlab = "Depth", ylab="Hind/He", main = "Samples")
abline(h=0.5, lty =2) #(h=(ploidy - 1)/ploidy)
#below expected would indicate some inbreeding
dev.off()
# Make sure there is not a relationship between depth and genetic diversity
# (2.7) List Those that are high or low heterozygosity
threshold.high <- mean(hhByInd) + 3 * sd(hhByInd)
threshold.low <- mean(hhByInd) - 3 * sd(hhByInd)

Filt.Ind.High <- hhByInd[hhByInd > threshold.high]
Filt.Ind.Low <- hhByInd[hhByInd < threshold.low]
Filter.dat <- c(Filt.Ind.High, Filt.Ind.Low)
write.csv(x=Filter.dat,
          file="Filtered.Individuals.csv", row.names=F)
# (2.8) Filter those unsual heterozygosity

hh <- hh[hhByInd < threshold.high & hhByInd > threshold.low,]
# (2.9) Filter the file!
mydata <- SubsetByTaxon(mydata, rownames(hh))
###############################################################################
# Step Three: Filter by site.
# Step (3.0): Re-calculate He from taxon filtered data
hh2 <- HindHe(mydata)
# Should about match the previous HindHe Calculation
# Step (3.1): Get mean Hind/He
hhByLoc <- colMeans(hh2, na.rm=T) # across sindividual
# Step (3.2) Plot the results
pdf("ByLocus.pdf")
hist(hhByLoc, breaks=30)
dev.off()

print("Estimate Center of Distribution")
Peak <- 0.35 # Change as needed
# (3.3) Calculate mean Hind/He
Est.Inb <- mean(hhByLoc, na.rm=T)
Est.Inb
# (3.4) calculate inbreeding rate
#InbreedingFromHindHe(Est.Inb, 2)
#InbreedingFromHindHe(Peak, 3) # 0.35
# Both at ~ 0.3
#ExpectedHindHe(mydata, inbreeding=0.25)
#ExpectedHindHe(mydata, inbreeding=0.3)
#ExpectedHindHe(mydata, inbreeding=0.35)
#ExpectedHindHe(mydata, inbreeding=0.4)

# 0.3 fits the data fairly well!
#######################################################################
# Step 4: Filter based on histogram data from last step
# (4.0) Plot cut-offs of which loci will be filtered!
png("ByLocus.Thresholds.jpg")
hist(hhByLoc, breaks=30)
abline(v=0.1, lty = 2, col="blue")
abline(v=0.45, lty = 2, col="blue")
dev.off()
# (4.1) Find those loci within the thresholds
keeploci <- names(hhByLoc)[hhByLoc > 0.05 & hhByLoc < 0.50]
# Removes 2752 markers or 11.5% of markers
# (4.2) get locations of "keepers"
keeploci.2 <- which(GetLoci(mydata) %in% keeploci)
# (4.3) Subset data by locus
mydata <- SubsetByLocus(mydata, keeploci.2)
nLoci(mydata)
nAlleles(mydata)

#######################################################################
# (5.0) Estimate Overdispersion for priors
# How good is my allele depth/ratio matching what I expected.
# So if I have a het are my ratios close to 0.50?
# Should be I filtered my vcf file to set to NA w/ ratios >.7/<.3
# Could explain high overdispersion values if I set them to NA.

#overdispersionP <- TestOverdispersion(mydata, to_test = 5:25)
# (5.1) Plot
# my_ovdisp <- overdispersionP$optimal
my_ovdisp <- 16

#######################################################################
# Step 6: Otay time to estimate some genotypes
# (6.0) Set seed for repeatable results
set.seed(1991)
# (6.1) Population structure based genotype estimation
mydataPopStruct <- IteratePopStruct(mydata,
                                    selfing.rate=Est.Inb,
                                    overdispersion = my_ovdisp)
# Can plot Population Structure:
# Give it color based on Eastern/Western/Mixed

# (6.3) Plot Population structure based genotype estimation
png("PopStrhist.jpg")
hist(mydataPopStruct$alleleFreq, breaks = 20)
dev.off()

#######################################################################
# Step 7: Export file formatted for GWA.
dir.create("MergeHaplotypeChange/")
setwd("MergeHaplotypeChange/")
wmgenoPopStruct <- GetWeightedMeanGenotypes(mydataPopStruct)
ProbgenoPopStruct <- GetProbableGenotypes(mydataPopStruct)
write.csv(x=wmgenoPopStruct,file="Weighted_Mean.csv", row.names = F)
write.csv(x=ProbgenoPopStruct$genotypes, file="Probgeno.csv",row.names=T)

# (7.0) Export for GWASpoly as discrete & continuous genotype
Export_GWASpoly(mydataPopStruct, file = "FiltGenoContinuousPoly.csv", naIfZeroReads = TRUE, postmean = TRUE, digits = 3)
Export_GWASpoly(mydataPopStruct, file = "FiltGenoDiscretePoly.csv", naIfZeroReads = TRUE, postmean = FALSE, digits = 3)
# (7.1) Export for GAPIT Continuous
Holder <- ExportGAPIT(mydataPopStruct)
write.csv(x=Holder$GD, file="ContinousGAPITGD.csv", row.names=F)
write.csv(x=Holder$GM, file="ContinousGAPITGM.csv", row.names=F)
# (7.2) Export PCA (i.e., population structure)
mydata <- AddPCA(mydata, nPCsInit=10, maxR2changeratio = 0.03, minPCsOut = 3)
PopStruc <- mydataPopStruct$PCA
write.csv(x=PopStruc, file="PCA.Pop.csv", row.names=F)

#VCF <- RADdata2VCF(mydata, file="polRADvcf.gz")

# The end :)
quit()
