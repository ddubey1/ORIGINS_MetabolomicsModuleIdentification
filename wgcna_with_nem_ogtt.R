# WGCNA Module Identification of Metabolomics Data

#-------------------------------------------------------------------------------
#-----------------Data Loading, Filtering, and Preprocessing--------------------
#-------------------------------------------------------------------------------

# -------------------- Loading Expression Data ---------------------------------

# set working directory
getwd()
setwd("/Users/dhruv/Desktop/Summer 2022/AH Research/ORIGINS Multi-Omics : Metabolomics")

# load package
library(WGCNA)
options(stringsAsFactors = FALSE)
library(limma)
library(edgeR)

# read in raw expression data
base_nfl_exprs <- read.csv("exprs.csv") # filtered (70% and log2 transformed)
dim(base_nfl_exprs)
names(base_nfl_exprs)

# Better way to filter 70%
zeroes = base_nfl_exprs[,-1] == 0
zeroes_sum = rowSums(zeroes)
idx_keep = which(zeroes_sum < 0.7 * ncol(zeroes))
length(idx_keep)
base_nfl_exprs = base_nfl_exprs[idx_keep,]
dim(base_nfl_exprs)

# log2 transform expression data
base_nfl_exprs[,-1] = log(base_nfl_exprs[,-1] + 1, base = 2)
dim(base_nfl_exprs)

# convert expression data to matrix form
base_nfl_exprs_matrix = as.matrix(base_nfl_exprs, rownames = TRUE)
rownames(base_nfl_exprs_matrix) = base_nfl_exprs_matrix[,1]
base_nfl_exprs_matrix = base_nfl_exprs_matrix[,-1]

# make numeric
base_nfl_exprs_matrix_num <- matrix(as.numeric(base_nfl_exprs_matrix), ncol = ncol(base_nfl_exprs_matrix))
rownames(base_nfl_exprs_matrix_num) = rownames(base_nfl_exprs_matrix)
colnames(base_nfl_exprs_matrix_num) = colnames(base_nfl_exprs_matrix)

# quantile normalization (maybe i can do this in voom to the DGEList???)
bnorm_exprs_matrix = normalizeBetweenArrays(base_nfl_exprs_matrix_num, method="quantile")
bnorm_exprs_matrix = bnorm_exprs_matrix[,-(1:30)] # CHANGED!!!!
bnorm_exprs_matrix = t(bnorm_exprs_matrix)
dim(bnorm_exprs_matrix)

rm(base_nfl_exprs_matrix_num)
rm(base_nfl_exprs_matrix)
rm(base_nfl_exprs)

# expression matrix
bnorm_exprs_matrix 

# remove novel metabolites
colnames(bnorm_exprs_matrix)[907]
bnorm_exprs_matrix = bnorm_exprs_matrix[,1:907]
dim(bnorm_exprs_matrix)

# -------------------- Filtering Expression Data -------------------------------

# identify good genes/metabolites
gsm = goodSamplesGenes(bnorm_exprs_matrix, verbose = 3)
gsm$allOK

if (!gsm$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsm$goodGenes)>0)
    printFlush(paste("Removing metabolites:", paste(names(bnorm_exprs_matrix)[!gsm$goodGenes], collapse = ", ")))
  if (sum(!gsm$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(bnorm_exprs_matrix)[!gsm$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  bnorm_exprs_matrix = bnorm_exprs_matrix[gsm$goodSamples, gsm$goodGenes]
}
dim(bnorm_exprs_matrix)

# clustering
sampleTreeMetabo = hclust(dist(bnorm_exprs_matrix), method = "average")
sizeGrWindow(12, 9)
par(cex = 0.6)
par(mar = c())
plot(sampleTreeMetabo, main = "Sample clustering to detect outliers", sub = "", xlab = "",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# 1 outliers found, not being removed however

# # Removing Outlier
# # Plot a line to show the cut
# abline(h = 70, col = "red");
# # Determine cluster under the line
# clust = cutreeStatic(sampleTreeMetabo, cutHeight = 70, minSize = 10)
# table(clust)
# # clust 1 contains the samples we want to keep.
# keepSamples = (clust==1)
# bnorm_exprs_matrix = bnorm_exprs_matrix[keepSamples, ]
# dim(bnorm_exprs_matrix)

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
for(i in 1:length(gender_row)) { # F is 0, M is 1
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
base_clinical = base_clinical[,-(1:30)]

# transpose
base_clinical = t(base_clinical)

# convert character to numeric
bc <- matrix(as.numeric(base_clinical), ncol = ncol(base_clinical))
rownames(bc) = rownames(base_clinical)
colnames(bc) = colnames(base_clinical)
bc = as.data.frame(bc)

bc[1:10,4] = 0
bc[11:20,4] = 1
bc[21:30,4] = 2

# bc = bc[-c(12),] # removing outliers PreT2D.Base2

collectGarbage()

# clinical is created
bc

# -------------------- Outlier Identification ----------------------------------

# set working directory
getwd()
setwd("/Users/dhruv/Desktop/Summer 2022/AH Research/WGCNA")

# Run 1
# sample network based on Euclidean distance
A = adjacency(t(bnorm_exprs_matrix), type = "distance")
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
# Some outliers
# Not being removed b/c too few samples

# Run 2
# re-cluster data with clinical traits
sampleTree2Metabo = hclust(dist(bnorm_exprs_matrix), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColorsMetabo = numbers2colors(bc, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2Metabo, traitColorsMetabo, groupLabels = names(bc),
                    main = "Sample dendogram and heat map")



#-------------------------------------------------------------------------------
#-----------------Network Construction and Module Identification----------------
#-------------------------------------------------------------------------------

# -------------------- Identify Soft Threshold ---------------------------------


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(bnorm_exprs_matrix, powerVector = powers, verbose = 5)
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

# power estimate
sft$powerEstimate 
# trying 7 for non-novel though

# # Don't use this; use multi-step
# # --------------- One Step Network Construction and Module Detection ----------- 80-90
# 
# # One-step network construction and module detection
# network = blockwiseModules(bnorm_exprs_matrix, corType = "pearson", power = 3,
#                            TOMType = "unsigned", minModuleSize = 20,
#                            reassignThreshold = 0, mergeCutHeight = 0.25,
#                            numericLabels = TRUE, pamRespectsDendro = FALSE,
#                            saveTOMs = TRUE,
#                            saveTOMFileBase = "Base Metabolomics TOM",
#                            verbose = 3)
# 
# # Lists identified modules and sizes
# #  81 modules in descending size (module 0 are unclassified data)
# table(network$colors)
# 
# # open a graphics window
# sizeGrWindow(12, 9)
# # Convert labels to colors for plotting
# mergedColorsMetabo = labels2colors(network$colors)
# # Plot dendogram and modules colors underneath
# plotDendroAndColors(network$dendrograms[[1]], mergedColorsMetabo[network$blockGenes[[1]]],
#                     "Module colors", dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# 
# # Information saved for later (but needed for later analysis)
# moduleLabelsMetbo = network$colors
# moduleColorsMetbo = labels2colors(network$colors)
# geneTreeMetabo = network$dendrograms[[1]];
# MEsMetbo = network$MEs
# # optional
# save(MEsMetbo, moduleLabelsMetbo, moduleColorsMetbo, geneTreeMetabo,
#      file = "SavingModulesPower3.RData")


# BETTER
# --------------- Multi-Step Network Construction and Module Detection --------- 17

# Things that can be tweaked: softPower, minModuleSize, MEDissThresh, deepSplit

# 1. Calculate adjacencies with power
softPower = 7 # sft$powerEstimate 5
adjacencyMat = adjacency(bnorm_exprs_matrix, power = softPower)

# 2. Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacencyMat)
dissTOM = 1-TOM

# 3. Clustering the TOM and identifying modules

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Metabolite clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 10
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
?cutreeDynamic

table(dynamicMods)

# Plotting modules
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Metabolite dendrogram and module colors")


# 4. Merging modules with similar expression profiles (eigengenes)
# Calculate eigengenes
MEList = moduleEigengenes(bnorm_exprs_matrix, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Height cut of 0.25 corresponds to correlation of 0.75
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(bnorm_exprs_matrix, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

# Plot
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

# 24 modules
table(dynamicMods)

#-------------------------------------------------------------------------------
#------------------------Relating Modules to Trait Data-------------------------
#-------------------------------------------------------------------------------

# ---------------------Creating heatmap using trait data------------------------

# Define numbers of metabolites and samples
nMetabolites = ncol(bnorm_exprs_matrix)
nSamples = nrow(bnorm_exprs_matrix)

# Recalculate Module Eigengenes (MEs) with color labels
# MEs are summary profiles for each module -> then correlate MEs with clinical
MEs0 = moduleEigengenes(bnorm_exprs_matrix, moduleColors)$eigengenes
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
png(filename="heatmapClinicalModuleTraitNonNovel OGTT.png",width = 6, height = 10, units = "in", res = 600)
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
               main = paste("Module-trait relationships OGTT"))
dev.off()


# ---------------------Relating Modules to Group--------------------------------

# Define variable weight containing the weight column of datTrait
group = as.data.frame(bc$group)
names(group) = "group"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(bnorm_exprs_matrix, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(bnorm_exprs_matrix, group, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(group), sep="")
names(GSPvalue) = paste("p.GS.", names(group), sep="")


# -------------Metabolites with high GS for Group and MM for Group--------------

# Create descriptions
getwd()
setwd("/Users/dhruv/Desktop/Summer 2022/AH Research/ORIGINS Multi-Omics : Metabolomics")
descrip = read.csv("descrip.csv")

# Extract information from geneModuleMembership and geneTraitSignificance
module = c("magenta","purple","brown","pink")
moduleData = list()
for (i in module) {
  # matching module
  column = match(i, modNames);
  moduleGenes = moduleColors==i;
  
  # plot
  #sizeGrWindow(7, 7);
  png(filename=paste0("moduleTraitPlot_",i,".png"),width = 3, height = 3, units = "in", res = 600)
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", i, "module"),
                     ylab = "Gene significance for group",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = i)
  dev.off()
  
  
  # construct data frame
  moduleData[[i]] = data.frame(metabolite = rownames(geneModuleMembership)[moduleGenes], 
                               module = i, moduleMembership = geneModuleMembership[moduleGenes, column], 
                               groupSignificance = geneTraitSignificance[moduleGenes, 1])
  
  moduleData[[i]]["description"] <- ""
  moddata_rows =  moduleData[[i]][,1]
  
  # adding descriptions to (PreT2D Base vs. NGT Base)
  for(metabolite in 1:length(moddata_rows)) {
    corres_descr = descrip[match(moddata_rows[metabolite], descrip[,1]), 2]
    moduleData[[i]][metabolite, 5] = corres_descr
  }
  
  
  getwd()
  setwd("/Users/dhruv/Desktop/Summer 2022/AH Research/WGCNA/ModuleData")
  write.csv(moduleData[[i]], file = paste0("moduleTraitTable_",i,".csv"))
  
}

#-------------------------------------------------------------------------------







#-------------------------------------------------------------------------------

column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for group",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


rownames(geneModuleMembership)[moduleGenes]
geneTraitSignificance[moduleGenes, 1]
geneModuleMembership[moduleGenes, column]

magentaData = data.frame(metabolite = rownames(geneModuleMembership)[moduleGenes], 
                         module = i, magentaMembership = geneModuleMembership[moduleGenes, column], 
                         groupSignificance = geneTraitSignificance[moduleGenes, 1])




# used to collapse list
#?do.call
# do.call(rbind,listname)

# way to compare matrices
#all(rownames(geneModuleMembership) == rownames(geneTraitSignificance))

# write.csv()

table(dynamicMods)
table(dynamicColors)


#-------------------------------------------------------------------------------




