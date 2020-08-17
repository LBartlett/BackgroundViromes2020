##  Analysis Code for Bartlett, L.J. et al 2019.  ##
## Long-term effects of origin and management history on honeybee colony viriomes. ##
## Contact L.J. Bartlett. lewis.bartlett@uga.edu (see online for current affiliation) ##
## LJB ORCID: 0000-0002-4418-8071 ##

## Analysis of ddPCR fluoresence data for understanding long-term effects of
## honeybee colony origin on viriomes

#load libraries [ use install.packages() and then load if not installed ]
library(afex)
library(emmeans)
library(ggbeeswarm)
library(vegan)
library(plyr)
library(tidyr)

# set for reproducibility
set.seed(14000)

# Load 'clean data' file (see data respository associated with publication, currently stored in GitHub Repo)
load(file = "DataTables.RData") 

# Check read-in
head(VirData)

# Strip out NTC wells (we don't want them in downstream analysis)
VirData <- VirData[which(VirData$Treatment != 'Control'), ]

# Get the mean amplitude of the housekeeping gene for each sample
# lazy for-loop it
BAAmp <- as.data.frame(matrix(NA, nrow = length(unique(VirData$Name)), ncol = 3))

colnames(BAAmp) <- c('Name', 'AvAmp', 'Treatment')

BAAmp$Name <- unique(VirData$Name)

for (n in 1:NROW(BAAmp)){
  
  BAAmp$AvAmp[n] <- mean(VirData$Amplitude[which(VirData$Target == 'BetaActinAmel' & VirData$Name == BAAmp$Name[n])])
  
  BAAmp$Treatment[n] <- unique(VirData$Treatment[which(VirData$Name == BAAmp$Name[n])])
  
}

# Test for batch effects of molecular work
# ANOVA of an effect of colony origin on housekeeping gene fluorescence
summary(aov(BAAmp$AvAmp ~ BAAmp$Treatment))

# Evidence of a batch effect. Already taken decisions to use 'relative abundance'
# which will likely control for this quality issue and not introduce a molecular artifact

# Create a data frame of viral amplitudes for each target & origin 
# scaled against mean housekeeping gene to gain an adjusted / relative fluorescence
# again we're lazy for-looping this because apparently I wasn't feeling the apply functions on the day I wrote this
# (I probably wanted an excuse for tea breaks whilst 'code was running')

VirAdj <- VirData[which(VirData$Target != 'BetaActinAmel'),]

AdjAmp <- vector(mode = 'numeric', length = NROW(VirAdj))

for( n in 1:NROW(VirAdj)){
  
  AdjAmp[n] <- ((VirAdj$Amplitude[n]) - (BAAmp$AvAmp[which(BAAmp$Name == VirAdj$Name[n])])) / (BAAmp$AvAmp[which(BAAmp$Name == VirAdj$Name[n])])
  
}

VirAdj$AdjAmp <- AdjAmp

rm(AdjAmp)

## Create data frames of mean scaled viral amplitudes for each target and origin
# this is painfully slow and could be done much, much more quickly with some apply() usage
# but alas, I was young, and lazy, and am now old, and know better, but yet am still lazy

# Create wide version first [will be needed for later analysis]
MedAA <- as.data.frame(matrix(NA, nrow = length(unique(VirAdj$Name)), ncol = length(unique(VirAdj$Target))+2))

colnames(MedAA) <- c('Name', 'Origin', 'AKIV', 'BQCV', 'CBPV', 'DWVA', 'DWVB', 'LSV', 'SBPV', 'SBV')

MedAA$Name <- unique(VirAdj$Name)

for (n in 1:NROW(MedAA)){
  
  MedAA$Origin[n] <- unique(VirAdj$Treatment[which(VirAdj$Name == MedAA$Name[n])]) 
  
  MedAA$AKIV[n] <- mean(VirAdj$AdjAmp[which(VirAdj$Name == MedAA$Name[n] & VirAdj$Target == 'AKIV')])
  
  MedAA$BQCV[n] <- mean(VirAdj$AdjAmp[which(VirAdj$Name == MedAA$Name[n] & VirAdj$Target == 'BQCV')])
  
  MedAA$CBPV[n] <- mean(VirAdj$AdjAmp[which(VirAdj$Name == MedAA$Name[n] & VirAdj$Target == 'CBPV')])
  
  MedAA$DWVA[n] <- mean(VirAdj$AdjAmp[which(VirAdj$Name == MedAA$Name[n] & VirAdj$Target == 'DWVA')])
  
  MedAA$DWVB[n] <- mean(VirAdj$AdjAmp[which(VirAdj$Name == MedAA$Name[n] & VirAdj$Target == 'DWVB')])
  
  MedAA$LSV[n] <- mean(VirAdj$AdjAmp[which(VirAdj$Name == MedAA$Name[n] & VirAdj$Target == 'LSV')])
  
  MedAA$SBPV[n] <- mean(VirAdj$AdjAmp[which(VirAdj$Name == MedAA$Name[n] & VirAdj$Target == 'SBPV')])
  
  MedAA$SBV[n] <- mean(VirAdj$AdjAmp[which(VirAdj$Name == MedAA$Name[n] & VirAdj$Target == 'SBV')])
  
  
}

# create long version

MedAAL <- gather(MedAA, Target, Amplitude, 3:10)


MedAAL$Name <- factor(MedAAL$Name, ordered = FALSE)
MedAAL$Origin <- factor(MedAAL$Origin, ordered = FALSE)
MedAAL$Target <- factor(MedAAL$Target, ordered = FALSE)

## Plot main topic of interest - each viral target's amplitude by colony origin

# Set up plotting

Transpa <- function(color, percent) {
  
  rgb.val <- col2rgb(color)
  
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100)
  
  return(t.col)
  
}

Shades1 <- c('purple', 'green3', 'orange2')

Shades2 <- sapply( X = Shades1, FUN = Transpa, percent = 40)

# Plot mean relative abundance for each colony by Target x Origin

par(mfrow = c(3,3), mar = c(5,3,2,2), bty = "n")

for(VirX in 3:10){
  
  boxplot(MedAA[,VirX] ~ MedAA$Origin, ylim = c(min(MedAA[,VirX])-0.1, max(MedAA[,VirX])+0.1),
          main = colnames(MedAA)[VirX], ylab = NA, yaxt = 'n', xlab = 'Origin', 
          border = 'transparent', cex.axis = 1.3, cex.lab = 1.5, outline = FALSE, lwd = 1.2,
          boxlty = 0, whisklty = 0, staplelty = 0, boxwex = 0.3)
  
  stripchart(MedAA[,VirX] ~ MedAA$Origin,
             vertical = T, add = T, pch = 16, col = Shades2, cex = 1.4, 
             method = 'jitter', lwd = 2)
  
  VMS <- aggregate(MedAA[,VirX] ~ MedAA[,2], FUN = mean)
  stripchart(VMS[,2] ~ VMS[,1],
             vertical = T, add = T, pch = 95, col = Shades1, cex = 4, lwd = 2)
  
  
  
  title(ylab="Relative Abundance", line= 1 , cex.lab=1.5)
  
  
}

plot.new()

legend("topleft", unique(VirAdj$Treatment), fill = Shades1, bty = 'n', 
       border = NA, cex = 2.15)

par(mfrow=c(1,1), mar = c(5,5,2,2))


## Looks like there are some consistent differences between colonies by origin, but that
## esp in the case of feral, this might be viral target dependent.
# Main hypothesis of interest to be further examined

## Same information in different graphical format:

BarMat <- as.matrix(MedAA[order(MedAA$Name),3:10])
rownames(BarMat) <- sort(MedAA$Name)

BarMat <- BarMat + abs(min(BarMat))

SBP <- barplot(t(BarMat), xlab = NA, xaxt = 'n', ylab = 'Relative Viral Abundance', cex.lab = 1.5,
               col = sort(rainbow(8)), border = NA)

axis(side=1,at=SBP[c(7, 21, 35)],labels= c('Feral','High', 'Low'), lwd = 0, cex = 1.5)

axis(side=1,at=SBP[c(14,28)], labels = F, tick = T, lwd = 0 , lwd.ticks = 2, tcl = - 1)

# Please accept my apologies for how utterly disgusting this graph is x

# Make a leged to paste on using software of your choice.
plot.new()

legend("center", rev(colnames(BarMat)), fill = rev(sort(rainbow(8))), bty = 'n', 
       border = NA, cex = 1.2)

# Overall ugly but more intuitive for dimensionality reduction analysis


### Time for main hypothesis testing

## Let's see if and how origin determines viral titre / abundance

# Use the afex package to undertake a factorial mixed anova ('split-plot design')
# uses the lme4 linear mixed modelling engine, but more intuitively packaged for this kind of study

# Set core options for running the afex models

afex_options(return_aov = "afex_aov", method_mixed = "KR", emmeans_model = 'multivariate')
options(contrasts=c('contr.sum','contr.poly'))

# Note this is the 'easy' anova function but works exactly the same as self-specifying error structure
# see ?aov_ez, but overall friendlier. print.formula gives you what it 'creates'
# afex vignette under title '4-Factor mixed design' using the data(obk.long) data set is fully analagous to this
# except this is simply a two-factor

SPANOVA <- aov_ez(id = "Name",dv = "Amplitude", data = MedAAL,
                  between = c("Origin"), within = c("Target"),
                  observed = c("Origin", "Target"),
                  print.formula = TRUE)

# Semi-equivalent model using the mixed function and lme4 
# (less conservative due to a difference in degrees of freedom estimation)

#SPANOVA <- mixed(Amplitude ~ Origin*Target + (1|Name), data = MedAAL)


# Take a look at significance
SPANOVA

# Viral target significant (expected)
# Viral target x origin interaction significant (as speculated from graph/Fig. 1)
# Origin alone not significant - surprising
# Note that origin alone not being significant means our 'batch effect' seen in the earlier anova 
# isn't likely to be a concern, no molecular 'batch' artifacts

# Use further pairwise comparisons to better investigate effects
# emmeans package
PHT <- emmeans(SPANOVA, ~ Origin * Target)

# We're only interested in comparing between the 3 origins within each viral target,
# No interest in comparing between viruses
# Formula here specifies this
# note that without the 'adjust' I don't think an adjustment is made despite package notes.

contrast(PHT, interaction = TRUE, method = "pairwise", by = "Target", adjust = 'fdr')
# some differences. see main manuscript.


## Dimensionality reduction analysis - community ecology, NMDS
# vegan package

# As it expects 'abundance' data and not 'relative abundance data'
# use a + 1 trasform so that it doesn't throw a shit-fit and crash R
# -1 is hypothetically the lowest relative abundance, and given this works on ranks, should be fine

BVComm <- MedAA[,c(-1,-2)] + 1

rownames(BVComm) <- MedAA$Name

# Create Distance Matrix
# Use euclidean as it's very general and this sort of relative abundnace data is a bit unusual...
# ...we don't have 'normal' community ecology data so best be careful and use a very general distance measure

BVDM <- as.matrix((vegdist(BVComm, "euclidean")))

# perform the NMDS with a reduction to 2 dimensions (k = 2)
BeeVirNMDS <- metaMDS(comm = BVDM, k = 2)

# check output
BeeVirNMDS

# stress = 0.060 looks good
#Check if stress plot is monotonically increasing
stressplot(BeeVirNMDS)

# looks good!

# Plot over the k=2 dimensions, colour coded by origin
ShadeV <- MedAA$Origin

for(N in 1: NROW(unique(MedAA$Origin))){
  
  ShadeV[which(ShadeV == unique(MedAA$Origin)[N] )] <- Shades1[N]
  
  
}

BVP<- ordiplot(BeeVirNMDS,type="n",
               xlim = round(c(min(BeeVirNMDS$points[,1]), max(BeeVirNMDS$points[,1]))*1.1, digits = 0), 
               ylim = round(c(min(BeeVirNMDS$points[,2]), max(BeeVirNMDS$points[,2]))*1.1, digits = 2)
)

points(BVP, "sites", pch = 19, col = ShadeV, cex = 1.7)

#Possible clustering (feral occupy smaller / different space?)

adonis(BVDM ~ Origin, data = MedAA, permutations = 10000)

# note that the p value varies as permutation-based. Orbits at around 0.04.

##### Analysis done

