
library(Seurat)
library(AUCell)
library(GSEABase)
library(stringr)
library(MAST)
load("F:/Lab_info/wanglab/My_project/crab-eating macaque/mjj/Fw_ crab-eating macaque VISp single nucleus data analysis/sc_maf5.RData")

exprMatrix <- sc@assays$RNA@counts
dim(exprMatrix)
expGenes <- row.names(exprMatrix)



# load PLS genes 
gene1005 <- read.csv("F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/gene1005List.xls",
                     sep="\t")
gene1005 <- gene1005[,1]
load("F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/gene1005_pos_mfas5.RData")
load("F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/gene1005_neg_mfas5.RData")


allgenes <- read.csv("F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/mfas5_ncx_expressed_23613genes.xls",
                     sep="\t")
allgenes_list <- allgenes[,2]

length(intersect(allgenes_list,expGenes)) # 17440

setdiff(expGenes,allgenes)

gene1005_pos_genes <- str_split(gene1005_pos,"_",simplify = T)[,2]
gene1005_neg_genes <- str_split(gene1005_neg,"_",simplify = T)[,2]

length(intersect(gene1005_pos_genes, expGenes)) # 485
length(intersect(gene1005_neg_genes, expGenes)) # 448

gene1005_pos_genes <- intersect(gene1005_pos_genes, expGenes)
gene1005_neg_genes <- intersect(gene1005_neg_genes, expGenes)


extraGeneSets <- GeneSetCollection(
  GeneSet(gene1005_pos_genes,setName="gene1005POS_485"),
  GeneSet(gene1005_neg_genes,setName="gene1005NEG_448"))

setwd("F:/Lab_info/wanglab/My_project/crab-eating macaque/Summay/NC_response/AUCell/")


cells_AUC <- AUCell_run(exprMatrix, extraGeneSets)
save(cells_AUC, file="mfas5_visp_gene1005POSNEG_cells_AUC.RData")

cells_assignment <- AUCell_exploreThresholds(cells_AUC,
                                             plotHist = TRUE,
                                             assignCells = TRUE)
getThresholdSelected(cells_assignment)

cellsTsne <- sc@reductions$umap@cell.embeddings
selectedThresholds<- getThresholdSelected(cells_assignment)
selectedThresholds
AUCell_plotTSNE(tSNE=cellsTsne, exprMat=exprMatrix,
                cex=0.7,borderColor = NA,
                cellsAUC=cells_AUC[1:2,],
                thresholds=selectedThresholds)

par(mfrow=c(2,3))
AUCell_plotTSNE(tSNE=cellsTsne, exprMat=exprMatrix, 
                cex=0.7, 
                alphaOff=0.3,borderColor = NA,
                cellsAUC=cells_AUC[1:2,], 
                thresholds=c(0.05,0.055))


