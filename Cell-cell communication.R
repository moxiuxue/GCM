library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(dplyr)
#library(SeuratData)
options(stringsAsFactors = FALSE)

library(stringr)
library(ggpubr)#for many plots in one page
library(data.table)

setwd("./cellchat/")
source('convertR.R', encoding = 'UTF-8')
source('processMeta.R', encoding = 'UTF-8')
source('processData.R', encoding = 'UTF-8')

data1 <- readRDS("GCM_Seurat3_1.rds")
data1@meta.data$typelabels <- data1@active.ident
Idents(data1)<-data1@meta.data$subtype

rGenes <- row.names(data1)
# Basic function to convert Ratto to human gene names
source('convertR.R', encoding = 'UTF-8')
orthologs <- convertR(rGenes)
# Preprocess the meta data
data1 <- processData(data1,orthologs)
data1 <-NormalizeData(data1)

data.input<- data1@assays$RNA@data
meta = data1@meta.data
batch_id = c("GCM") #calculate the cell-cell communication in GCM samples
#batch_id = c("CTL") #calculate the cell-cell communication in CTL samples
meta <- processMeta(meta,batch_id)
cell.use = rownames(meta)
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
meta = data.frame(meta, row.names = colnames(data.input)) # create a dataframe consisting of the cell labels

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "subtype")
#cellchat <- createCellChat(object = data.input, meta = identity, group.by = "labels")
cellchat

#cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "subtype") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)
#dplyr::glimpse(CellChatDB$interaction)
unique(CellChatDB$interaction$annotation)
#CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # use Cell-Cell Contact for cell-cell communication analysis
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use # set the used database in the object

## Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat,thresh.p = 10^(-6))
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

# Part II: Inference of cell-cell communication network
## Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat,raw.use = FALSE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

## Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat)

## Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

## Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

#saveRDS(cellchat, file = "cellchat_CTL.rds")
saveRDS(cellchat, file = "cellchat_GCM.rds")

#Part III: Visualization of cell-cell communication network
cellchat.GCM <- cellchat
#Visualize the aggregated cell-cell communication network. Showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat.GCM@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.GCM@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.GCM@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Examine the signaling sent from each cell group.
mat <- cellchat.GCM@net$weight
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  #iid <- c(14,15,17,19,20,21,22,23,24,25,26,27,30,31,32)
  #for (i in iid) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Heatmap
netVisual_heatmap(cellchat.GCM, font.size = 15, font.size.title = 15)
netVisual_heatmap(cellchat.GCM, font.size = 15, font.size.title = 15, measure = "weight")

pathways.show.all <- cellchat.GCM@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat.GCM@idents)
vertex.receiver = 1:2
for (i in 1:length(pathways.show.all)) {
  #extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
  pairLR <- extractEnrichedLR(cellchat.GCM, signaling = pathways.show.all[i], geneLR.return = FALSE)
  # Visualize communication network associated with both signaling pathway and individual L-R pairs  
  gg <- netAnalysis_contribution(cellchat.GCM, signaling = pathways.show.all[i])
  ggsave(filename=paste0(save_dir,"FB_MT1A_Allpathway_LR/",pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 10, height = 2, units = 'in', dpi = 300)
}



#Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways from some cell groups to other cell groups

set1 <- c("B_Cells", "DC_1", "DC_2", "DC_3", "Endothelial", "Fibroblasts", "Gran_Pro", "ILC", "Lymp_Pro", "Macro_1",
"Macro_2", "Macro_3", "Macro_4", "Macro_5", "NK_1", "NK_2", "T_1", "T_2", "T_3", "T_4")
set2 <- c("Neu_1", "Neu_2", "Neu_3")

#set1<-c("Macro_1","Macro_2", "Macro_3", "Macro_4", "Macro_5")
#Bubble plot
netVisual_bubble(cellchat.GCM, sources.use = set1, targets.use = set2, remove.isolate = FALSE)
netVisual_bubble(cellchat.GCM, sources.use = set2, targets.use = set1, remove.isolate = FALSE)

#Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
netAnalysis_signalingRole_scatter(cellchat.GCM)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.GCM, pattern = "outgoing",width = 10,  height = 20)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.GCM, pattern = "incoming",width = 10,  height = 20)
ht1 + ht2


cellchat.GCM <- computeNetSimilarity(cellchat.GCM, type = "functional")
cellchat.GCM <- netEmbedding(cellchat.GCM, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat.GCM <- netClustering(cellchat.GCM, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat.GCM, type = "functional", label.size = 3.5)


cellchat.GCM <- computeNetSimilarity(cellchat.GCM, type = "structural")
cellchat.GCM <- netEmbedding(cellchat.GCM, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat.GCM <- netClustering(cellchat.GCM, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat.GCM, type = "structural", label.size = 3.5)

##Visualize cell-cell communication mediated by specific signaling pathways from some cell groups to other cell groups
pathways.show <- c("IL1","SPP1","APRIL","BAFF","TWEAK","THBS","THY1","SN")
pathways.show <- c("CXCL","CSF")

# show all the significant interactions (L-R pairs) associated with certain signaling pathways from some cell groups to other cell groups
netVisual_bubble(cellchat.GCM, sources.use = set1, targets.use = set2, signaling = pathways.show, remove.isolate = FALSE)
netVisual_bubble(cellchat.GCM, sources.use = set2, targets.use = set1, signaling = pathways.show, remove.isolate = FALSE)


for (i in 1:length(pathways.show)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual_chord_gene(cellchat.GCM, sources.use = set1, targets.use = set2, signaling = pathways.show[i],legend.pos.x = 10)
  #netVisual_chord_gene(cellchat.GCM, sources.use = set2, targets.use = set1, signaling = pathways.show[i],legend.pos.x = 10)

  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat.GCM, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_L-R_contribution.pdf"), plot=gg, width = 10, height = 3, units = 'in', dpi = 300)
}

#> Comparing communications on a single object
# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat.GCM, signaling = pathways.show)
netVisual_bubble(cellchat.GCM, sources.use = set1, targets.use = set2, pairLR.use = pairLR.use, remove.isolate = TRUE)
#> Comparing communications on a single object
netVisual_chord_gene(cellchat.GCM, sources.use = set1, targets.use = set2, signaling = pathways.show,legend.pos.x = 8)
