# WGCNA Module Identification of Metabolomics Data

# -------------------- Loading Expression Data ---------------------------------

# set working directory
getwd()
setwd("/Users/dhruv/Desktop/Summer 2022/AH Research/WGCNA")

# load package
library(WGCNA)
options(stringsAsFactors = FALSE)

# note: this is the WGCNA Tutorial method (you can also try preprocessing the data the DEA method - 70% filtering, log2, quatilenorm)
# read in raw expression data
base_exprs <- read.csv("BaseOnly.csv")
rownames(base_exprs) = base_exprs[,1]
base_exprs = base_exprs[,-1]
base_exprs = t(base_exprs)

# -------------------- Filtering Expression Data -------------------------------

# identify good genes/metabolites
gsm = goodSamplesGenes(base_exprs, verbose = 3)
gsm$allOK # 160 removed

if (!gsm$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsm$goodGenes)>0)
    printFlush(paste("Removing metabolites:", paste(names(base_exprs)[!gsm$goodGenes], collapse = ", ")));
  if (sum(!gsm$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(base_exprs)[!gsm$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  base_exprs = base_exprs[gsm$goodSamples, gsm$goodGenes]
}

# clustering
sampleTreeMetabo = hclust(dist(base_exprs), method = "average")
sizeGrWindow(12, 9)
par(cex = 0.6)
par(mar = c())
plot(sampleTreeMetabo, main = "Sample clustering to detect outliers", sub = "", xlab = "",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# no outliers found, so tree cut is not needed

# -------------------- Loading Clinical Data -----------------------------------

# set working directory
getwd()
setwd("/Users/dhruv/Desktop/Summer 2022/AH Research/ORIGINS Multi-Omics : Metabolomics")

# read in clinical data
# group: NGT PreT2D T2D -> 0 1 2
base_clinical = read.csv("clinical_tofix copy.csv")


# make male 1 and females 0
rownames(base_clinical) = base_clinical[,1]
base_clinical = base_clinical[,-1]
gender_row = base_clinical[1,]
for(i in 1:length(gender_row)) {
  if(gender_row[i]=="F") {
    gender_row[i] = 0
  }
  else {
    gender_row[i] = 1
  }
}
base_clinical[1,] = gender_row

base_clinical

# filter only base samples
base_clinical = base_clinical[,-(31:60)]

# transpose
base_clinical = t(base_clinical)

# convert character to numeric
bc <- matrix(as.numeric(base_clinical), ncol = ncol(base_clinical))
rownames(bc) = rownames(base_clinical)
colnames(bc) = colnames(base_clinical)
bc = as.data.frame(bc)


collectGarbage()

# -------------------- Outlier Identification ----------------------------------

# Run 1
# sample network based on Euclidean distance
A = adjacency(t(base_exprs), type = "distance")
# calculate whole network connectivity
k = as.numeric(apply(A,2,sum))-1
# standardize connectivity
Z.k = scale(k)
thresholdZ.k = -5 # or -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
# cluster
sampleTree3Metabo = hclust(as.dist(1-A), method="average")
traitColors3Metabo = data.frame(numbers2colors(bc,signed=FALSE))
dimnames(traitColors3Metabo)[[2]]=paste(names(bc),"C",sep="")
datColors3Metabo = data.frame(outlierC=outlierColor,traitColors3Metabo)
plotDendroAndColors(sampleTree3Metabo, groupLabels = names(datColors3Metabo),
                    colors = datColors3Metabo, main = "Sample dendogram and heat map")
# 0 outliers

# Run 2
# re-cluster data with clinical traits
sampleTree2Metabo = hclust(dist(base_exprs), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColorsMetabo = numbers2colors(bc, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2Metabo, traitColorsMetabo, groupLabels = names(bc),
                    main = "Sample dendogram and heat map")


# -------------------- Identify Soft Threshold ---------------------------------


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(base_exprs, powerVector = powers, verbose = 5)
sft # you want high r^2 and negative slope (using power of 3)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# -------------------- Network Construction and Module Detection ---------------

# One-step network construction and module detection
network = blockwiseModules(base_exprs, corType = "pearson", power = 3,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Base Metabolomics TOM",
                       verbose = 3)

# Lists identified modules and sizes
# 68 modules in descending size (module 0 are unclassified data)
table(network$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColorsMetabo = labels2colors(network$colors)
# Plot dendogram and modules colors underneath
plotDendroAndColors(network$dendrograms[[1]], mergedColorsMetabo[network$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Information saved for later (but needed for later analysis)
moduleLabelsMetbo = network$colors
moduleColorsMetbo = labels2colors(network$colors)
geneTreeMetabo = network$dendrograms[[1]];
MEsMetbo = network$MEs
# optional
save(MEsMetbo, moduleLabelsMetbo, moduleColorsMetbo, geneTreeMetabo,
     file = "SavingModulesPower3.RData")

# -------------------- Relating Modules to trait data --------------------------

# Define numbers of metabolites and samples
nMetabolites = ncol(base_exprs)
nSamples = nrow(base_exprs)

# Recalculate Module Eigengenes (MEs) with color labels
# MEs are summary profiles for each module -> then correlate MEs with clinical
MEs0 = moduleEigengenes(base_exprs, moduleColorsMetbo)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, bc, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(100,100)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(bc),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
















