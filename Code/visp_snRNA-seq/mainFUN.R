if(T){
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(rlist)
  library(ggplotify)
  library(pheatmap)
}

levels = c('Ex_0','Ex_1', 'Ex_2', 'Ex_3', 'Ex_4', 'Ex_5', 'Ex_6', 'Ex_7', 'Ex_8', 
           'Ex_9', 'Ex_10', 'Ex_11','Ex_12', 'Inh_0', 'Inh_1', 'Inh_2', 'AST_0', 'AST_1',
           'Micro', 'oligo_0', 'oligo_1','OPC')


# FUN readdata ------------------------------------------------------------
theme_common = theme_bw() + theme(panel.background = element_blank(), panel.grid = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))

mito = read.csv('/mnt/GLHu/DATA_win/mitoList_macaque.csv')
mito = c(na.omit(mito$gene_name), mito$gene_id)
seuObj = function(pth, samp){
  pth = paste0(pth,samp)
  counts = Read10X(pth)
  fin = CreateSeuratObject(counts,project = samp,min.cells = 3,min.features = 200)
  fin = RenameCells(fin,add.cell.id = samp)
  fin@meta.data$cellID = row.names(fin@meta.data)
  mtExp = intersect(row.names(fin), mito)
  # fin[['pct.mt']] = PercentageFeatureSet(fin , pattern = '^MT')
  fin[['pct.mt']] = PercentageFeatureSet(fin , features = mtExp, assay = 'RNA')
  ## features in mtExp must expressed in the data
  fin[['pct.ribo']] = PercentageFeatureSet(fin , pattern = '^RB[SL]')
  names(fin@meta.data)[1:3] = c('samp', 'nUMI', 'nGene')
  return(fin)
}

featuresPlot = function(sc,plot.featrures = c("nGene", "nUMI", "pct.mt", "pct.ribo"),group = 'samp',pth1,pth2){
  theme.set1 = theme(axis.title.x=element_blank(), 
                     axis.text.x=element_blank(), 
                     axis.ticks.x=element_blank())
  theme.set2 = theme(axis.title.x=element_blank())
  theme_common = theme_bw() + theme(panel.background = element_blank(), panel.grid = element_blank())+
    theme(plot.title = element_text(hjust = 0.5))
  plots = list()
  for(i in seq_along(plot.featrures)){
    plots[[i]] = VlnPlot(sc,pt.size = 0,group.by = group,
                         features = plot.featrures[i]) + theme.set2 + NoLegend()}
  violin <- wrap_plots(plots = plots, nrow=2)    
  ggsave(pth1, plot = violin, width = 10, height = 8) 
  ### plot 
  meta = sc@meta.data
  p1 = ggplot(meta, aes(x = nUMI, y = nGene, color = samp)) + geom_point() + theme_common + NoLegend()
  p2 = ggplot(meta, aes(x = nUMI, y = pct.mt, color = samp)) + geom_point() + theme_common + NoLegend()
  p3 = ggplot(meta, aes(x = nUMI, y = pct.ribo, color = samp)) + geom_point() + theme_common
  p = p1+p2+p3
  ggsave(p, filename = pth2, width =8, height = 10)
}

my_DoubletFinder = function(obj){
  library(DoubletFinder)
  obj = CreateSeuratObject(as.matrix(obj[['RNA']]@counts), meta.data = obj@meta.data[,1:6], 
                           min.cells = 5, min.features = 200)
  pc.num = 1:30 
  obj = obj %>% SCTransform() %>% RunPCA() %>% RunUMAP(dim = pc.num) %>% FindNeighbors() %>% FindClusters()
  sweep.res.list=paramSweep_v3(obj,PCs = pc.num,sct = T,num.core=20) ## sct:Logical representing whether SCTransform was used during original seurat object pre-processing
  ##paramSweep_v3是计算pN和pK参数组合的函数
  sweep.stats=summarizeSweep(sweep.res.list,GT=F)
  opt_pK=find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(opt_pK$pK[which.max(opt_pK$BCmetric)])) ##select the best pK
  homotypic.prop=modelHomotypic(annotations = obj@meta.data$seurat_clusters)
  nExp_poi <- round(0.075*nrow(obj@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  obj <- doubletFinder_v3(obj, PCs = pc.num, pN = 0.2, pK = mpK,
                          nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  return(obj)
}

p_DoubletFinder = function(obj.list, pth, w=7,h=7){
  fin = lapply(obj.list, function(x){
    tar = x@meta.data[,c(1,2,4,5,6,15)]
    names(tar)[6] = 'tag'
    return(tar)
  }) %>% list.rbind()
  p = ggplot(fin, aes(x = samp, y = nUMI,fill = tag)) + theme_common +
    geom_boxplot()
  ggsave(p, filename = pth)
}

my_Heatmap = function(sc,pth){
  tar = sc@reductions$umap@cell.embeddings %>% data.frame()
  tar$cellID = row.names(tar)
  tar = left_join(tar,sc@meta.data, by = 'cellID')
  colours = c('white', 'grey', 'red')
  p_nGene = ggplot(tar,aes(x = UMAP_1,y = UMAP_2, colour = nGene)) + geom_point(size = 0.5, alpha = 0.8)+ scale_color_gradientn(colors = colours)
  p_nUMI = ggplot(tar,aes(x = UMAP_1,y = UMAP_2, colour = nUMI)) + geom_point(size = 0.5, alpha = 0.8)+ scale_color_gradientn(colors = colours)
  p_pct.mt = ggplot(tar,aes(x = UMAP_1,y = UMAP_2, colour = pct.mt)) + geom_point(size = 0.5, alpha = 0.8)+ scale_color_gradientn(colors = colours)
  p_pct.ribo = ggplot(tar,aes(x = UMAP_1,y = UMAP_2, colour = pct.ribo)) + geom_point(size = 0.5, alpha = 0.8)+ scale_color_gradientn(colors = colours)
  p = (p_nGene|p_nUMI)/(p_pct.mt|p_pct.ribo)
  ggsave(p, filename = pth,width = 10, height = 9)
}

if(F){
  monkey_meta = read.csv('/mnt/GLHu/DATA_win/monkey_meta.csv')
  monkey_meta$class_label %>% table
  mouse_meta = read.csv('/mnt/GLHu/DATA_win/mouse_meta.csv')
  mouse_meta$class_label %>% table
  human_meta = read.csv('/mnt/GLHu/DATA_win/human_meta(1).csv')
  human_meta$class_label %>% table

}

my_topMarkers = function(mks,topn = 1:3){
  mks$pct.sig = mks$pct.1 - mks$pct.2
  mks[mks$p_val_adj == 0,]
  fin = lapply(unique(mks$cluster), function(x){
    tar = mks[mks$cluster == x,]
    t1 = tar[order(tar$avg_log2FC,decreasing = T),]
    t1 = t1[topn,]
    t2 = tar[order(tar$pct.sig,decreasing = T ),]
    t2 = t2[topn,]
    fin = unique(rbind(t1,t2))
    return(fin)
  }) %>% list.rbind() %>% data.frame()
  return(fin)
}

my_clusterMarkers = function(mks,pth){
  lapply(unique(mks$cluster) , function(x){
    tar = mks[mks$cluster == x,]
    pth = paste0(pth,x,'.pdf')
    p = VlnPlot(sc, tar$gene, stack = T, combine = F, flip = T,  same.y.lims = T) + NoLegend()
    ggsave(p, filename = pth, width = 7, height = 5)
  })
}











