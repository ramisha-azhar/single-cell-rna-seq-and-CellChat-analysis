# script to perform standard workflow steps to analyze single cell RNA-Seq data and cellchat analysis
# data: 20k Mixture of NSCLC DTCs from 7 donors, 3' v3.1
# data source: https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-controller-1-standard-2-0-0         

setwd("C:/Users/ramis/OneDrive/Desktop/Thesis")
### load libraries
library(Seurat)
library(tidyverse)

### Load the NSCLC dataset
nsclc.sparse.m <- Read10X_h5(filename = '20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')

#Genome matrix has multiple modalities, returning a list of matrices for this genome
# we have other modalities present in it other than gene expression we can see through str() function

str(nsclc.sparse.m)
'''
> str(nsclc.sparse.m)
List of 3
3 means it has three modalities
1)Gene Expression
2)Antibody Capture 
3)Multiplexing Capture
'''
# as we are interested for gene expression data so we use only gene expression data
#gene expression has the counts matrix which we need right now
cts <-  nsclc.sparse.m$`Gene Expression`
dim(cts) #[1]   36601 3862363
cts[1:10,1:10]
class(cts) #"dgCMatrix"

### Initialize the Seurat object with the raw count (non-normalized data).
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)

str(nsclc.seurat.obj)
class(nsclc.seurat.obj) #"Seurat"
nsclc.seurat.obj
# 29552 features across 42081 samples
'''
An object of class Seurat 
29552 features across 42081 samples within 1 assay 
Active assay: RNA (29552 features, 0 variable features)
 1 layer present: counts
'''
### 1. QC -------
View(nsclc.seurat.obj@meta.data)
dim(nsclc.seurat.obj@meta.data) #[1] 42081 3    42081 is the number of cells 

## % MT reads
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)#Here we calculate the percatage of the mitochondrial genes in each cell and add the new column to the metadat
#now the percentage of the mitochondrial gene to each cell has been added

FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')


### 2. Filtering -----------------
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)
nsclc.seurat.obj
#29552 features across 24708 samples after filtering before the samples was 42081
'''
An object of class Seurat 
29552 features across 24708 samples within 1 assay 
Active assay: RNA (29552 features, 0 variable features)
 1 layer present: counts
'''


### 3. Normalize data ----------
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
nsclc.seurat.obj
#2 layers present: counts, data
'''
An object of class Seurat 
29552 features across 24708 samples within 1 assay 
Active assay: RNA (29552 features, 0 variable features)
 2 layers present: counts, data
'''


### 4. Identify highly variable features --------------
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)
nsclc.seurat.ob
'''
> nsclc.seurat.obj
An object of class Seurat 
29552 features across 24708 samples within 1 assay 
Active assay: RNA (29552 features, 2000 variable features)
 2 layers present: counts, data
'''
# Identify the 10 most highly variable genes
#to see our variable features we use this function VariableFeatures() and we select top 10 variable features
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)
top10
'''
 [1] "IGHG1"  "IGKC"   "IGHA1"  "IGHG3"  "IGLC2"  "TPSB2" 
 [7] "TPSAB1" "IGHM"   "JCHAIN" "IGHGP"
'''
# plot variable features with and without labels
#we want to visualize see the features that are variable features
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
#out of 2000 genes which are variable features we can see our top 10 features in our plot with the gene labels 2000 features are in red
LabelPoints(plot = plot1, points = top10, repel = TRUE)


### 5. Scaling -------------
all.genes <- rownames(nsclc.seurat.obj) #29552 all the row names genes are used
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)
#str() function use to see slots in our seurat object 
str(nsclc.seurat.obj)
'''
we have three kind of slots under rna assay
1) counts slot(it has information of raw sparse matrix the count that we read in first step)
2) data slot (when we perform the normalization step the normalize counts are stored in the data slot)
3) scale.data slot (after scaleing all our data stored in this slot )
'''
nsclc.seurat.obj
#3 layers present: counts, data, scale.data
'''
An object of class Seurat 
29552 features across 24708 samples within 1 assay 
Active assay: RNA (29552 features, 2000 variable features)
 3 layers present: counts, data, scale.data
'''


### 6. Perform Linear dimensionality reduction --------------
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))
nsclc.seurat.obj
# 1 dimensional reduction calculated: pca
'''
An object of class Seurat 
29552 features across 24708 samples within 1 assay 
Active assay: RNA (29552 features, 2000 variable features)
 3 layers present: counts, data, scale.data
 1 dimensional reduction calculated: pca
'''
## visualize PCA results
#we will see the top 5 principle component with the top 5 genes of them
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
#we can use heatmap providing the seurat object dims= 1 it is the first PCA with 500 cells
# we can see the heterogenity and we consider only those pca which has heterogeneity
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)


# determine dimensionality of the data
# only choosing those statistical significant pca that capture the majority of the signals in our down stream analysis
# there are all the PCA and it captures the high variance of each PCA 
ElbowPlot(nsclc.seurat.obj)


### 7. Clustering ------------
# we want to cluester the similar cells which has same feature expression pattern we want them to be cluester together
# for that identify the nearest neighbor and we provide the first 15 pca which show the most variation in our data set as we see in elbow plot
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

## understanding resolution
# we want our cells to be assigns to the clusters for that we use function  FindCluster()
# the resolution parameter define the granularity or the resolutions of the clusters
# lower the number few are the clusters higher the number more are the clusters in the resolution parameter
# we can see which resolution works best to distinct out our clusters
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))

##we look at the metadata for each resolution different column has been created in our metadata
# we can see how many clusters we have for each resolution and see what work best for us
View(nsclc.seurat.obj@meta.data)

##we can group the cells by column with this plot we are trying the column 0.1 resolution
# with resolution 0.1 we get 8 clusters and they are distinct well
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE)
 
### 8 . setting identity of clusters------------------------
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)
#we can set the identity of the cells by assigning it to the column or variable of our metadata
#Levels: 0 1 2 3 4 5 6 7

### 9 . non-linear dimensionality reduction --------------

nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)

DimPlot(nsclc.seurat.obj, reduction = "umap")

nsclc.seurat.obj
'''
> nsclc.seurat.obj
An object of class Seurat 
29552 features across 24708 samples within 1 assay 
Active assay: RNA (29552 features, 2000 variable features)
 3 layers present: counts, data, scale.data
 2 dimensional reductions calculated: pca, umap
'''

##########################################################################################
##########################################################################################
##########################################################################################

#Step 1: load libraries
library(CellChat)
library(Seurat) 

#Step 2: Create the CellChat Object

## Extract data and metadata
# our normalize matrix is 29552 24708 store in data slot
data.input <- GetAssayData(nsclc.seurat.obj, assay = "RNA", layer = "data") 
View(data.input)

## Rename our identities so they do not contain 0
#because cellchat do not accept 0 
#it can be anything but not zero
labels <- Idents(nsclc.seurat.obj)
labels <- paste0("cluster_", labels)
meta <- data.frame(group = labels, row.names = names(labels))
View(meta)

##Create the CellChat Object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat
'''
An object of class CellChat created from a single dataset 
 29552 genes.
 24708 cells. 
CellChat analysis of single cell RNA-seq data!
'''
View(cellchat) #we can see all the slots in cellchat object


### Step 3: Set Up the Ligand-Receptor Database

# Set the database to human (since your dataset is NSCLC lung cancer)
# cell chat data bases which is used in our cell cell communication analysis 
# in the cell chat packages

##data base has the four element
#1) human data base has 3233 ligand receptor interactions paires
#2) complex : it is it has 338 ligand receptor complex paires
#3)cofactor : 32 ligand receptor complex pairs cofactor

CellChatDB.human <- CellChatDB.human # Load the human database 
View(CellChatDB.human)

##show the ligand receptor categories
#in human data base 39.6 are signaling pathway in first pi chart which are secreted signaling pathway
showDatabaseCategory(CellChatDB.human)

##Add CellChatDB in our cellchat object
# Set the full human ligand-receptor interaction database
CellChatDB.use <- CellChatDB.human

##Using a subset of CellChatDB for cell-cell communication analysis
#by subsetting we will get 1280 interaction instead of total interaction which is 3233 present in our database
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")

## Add the Secreted Signaling database in the CellChat object
#with View function we can see DB(database) is added to the cellchat object
cellchat@DB <- CellChatDB.use
View(cellchat)


### Step 4:Subset and pre-processing the expression data 
# we have sebset our genes which are in data signalling
# subset the expression data to use less RAM
# 749 genes in our cell chat database are expressed in our cellchat object some of these genes are ligand and some are receptores
# for rest of our analysis we only focus on these 749 genes for our ligand receptor analysis
cellchat <- subsetData(cellchat)

##Pre-processing the expression data
#we identify our expressed genes from 749 signalling gene list

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#The number of highly variable ligand-receptor pairs used for signaling inference is 931 


###Step 5: Compute Communication Probabilities
cellchat <- computeCommunProb(cellchat)


###Step 6: Filter out the cell-cell communication if there are only few number of cells 
cellchat <- filterCommunication(cellchat, min.cells = 10)

###Step 7: Infer the cell-cell communication at a signaling pathway level
# 17 pathway obtained by different cell type
cellchat <- computeCommunProbPathway(cellchat)

###Step 8: Calculate the aggregated cell-cell communication network
# we can plot a network between different cell types
# we aggregated the network cell communication network between different cell type
# count and weight let us see different cell type connecting to each other
cellchat <- aggregateNet(cellchat)
cellchat@net$count
cellchat@net$weight

###Step 9: visualize the aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd=TRUE)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

dev.off()

#  examine the signaling sent from each cell group
#interaction for cluster
mat <- cellchat@net$weight
par(mfrow = c(2, 6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = F, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()

#it will give us the cluster interaction of our interest
mat <- cellchat@net$weight
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:nrow(mat)) {  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))  
  mat2[i, ] <- mat[i, ]  
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])}

#########################################################################################
#########################################################################################
##########################################################################################

#Visualization of Cell-Cell Communication Network

library(CellChat)

#how many significant pathway we identify from cell chat object
# we identify 17 pathway
cellchat@netP[["pathways"]]

'''
[1] "MIF"      "SPP1"     "CypA"     "MK"       "GALECTIN"
[6] "ANNEXIN"  "VISFATIN" "TNF"      "TGFb"     "RESISTIN"
[11] "BAFF"     "APRIL"    "CD70"     "IFN-II"   "CD137"   
[16] "PLAU"     "CSF"
'''
#cellchat has function to extract enrich ligand receptor pairs
#we have 33 ligand receptor pairs  $pairLR
#and we have 56 genes $geneLR
extractEnrichedLR(cellchat, signaling = c(cellchat@netP[["pathways"]]),
                  geneLR.return = TRUE)

# visualize the contribution of each LR pairs to the communication network
# we have 33 ligand receptor pairs
netAnalysis_contribution(cellchat, 
                         signaling = c(cellchat@netP[["pathways"]]), 
                         title = "Contribution of each LR pairs")

#contribution of each LR pairs for first five signaling pathways
netAnalysis_contribution(cellchat, 
                         signaling = c(cellchat@netP[["pathways"]][1:5]), 
                         title = "MIF,SPP1,CypA,MK,GALECTIN")

#to run individual signaling pathway to see how many ligand receptor pair interact
extractEnrichedLR(cellchat, signaling = "TNF", geneLR.return = FALSE)
#two ligand receptor pair contribute to the TNF signaling pathway
'''
 interaction_name
1     TNF_TNFRSF1A
2     TNF_TNFRSF1B
'''

netAnalysis_contribution(cellchat, signaling = "TNF")
#we can see the contribution of each LR pair in TNF signalling pathway

# Circle plot for the individual signalling pathway
#we can plot the aggregated communication network for the TNF signalling pathway
netVisual_aggregate(cellchat, signaling = "TNF", layout = "circle")

#plot individual TNF signalling pathway for each LR pair
netVisual_individual(cellchat, signaling = "TNF", layout = "circle")

#if interested in only one LR pair for one signalling pathway
netVisual_individual(cellchat, signaling = "TNF", 
                     pairLR.use = "TNF_TNFRSF1B",
                     layout = "circle")

# Chord diagram
#parameters for the plot
par(mfrow = c(1, 1), xpd=TRUE)
par(cex = 0.5)

#aggregated network for the TNF signalling pathway

netVisual_aggregate(cellchat, signaling = "TNF", layout = "chord")

netVisual_chord_cell (cellchat, signaling = "TNF")

netVisual_chord_gene (cellchat, signaling = "TNF")


# Heatmap
netVisual_heatmap(cellchat, signaling = "TNF", color.heatmap = "Reds")

# Violin plot 
plotGeneExpression(cellchat, signaling = "TNF")





        