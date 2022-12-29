ind = dir('/mnt/GLHu/CR/')

sc = lapply(ind,function(x){
  fin = seuObj(pth = '/mnt/GLHu/CR/', samp = x) 
  return(fin)
})
sc = merge(sc[[1]],sc[2])
featuresPlot(sc = sc,pth1 = '/mnt/GLHu/RESULT/features_beforeQC.pdf', 
             pth2 = '/mnt/GLHu/RESULT/pointDistribution_beforeQC.pdf')

# QC ----------------------------------------------------------------------
minGene=500
maxGene=7500
maxnUMI = 30000
pctMT=0.5
sc <- subset(sc, subset = nGene > minGene & 
               nGene < maxGene &
               nUMI< maxnUMI &
               pct.mt < pctMT ) # 8827 cells left
featuresPlot(sc = sc,pth1 = '/mnt/GLHu/afterQC.pdf', 
             pth2 = '/mnt/GLHu/pointDistribution_afterQC.pdf')
sc.list =  SplitObject(sc, 'samp')
tmp = lapply(sc.list, my_DoubletFinder)
sc.list = tmp
p_DoubletFinder(obj.list = sc.list, pth = '/mnt/GLHu/RESULT/DF_boxplot.pdf')

sc.list = lapply(sc.list, function(x){
  names(x@meta.data)[14:15] = c('pANN_1', 'DF_1')
  return(x)
})

sc = merge(sc.list[[1]],sc.list[2])
sc = subset(sc, DF_1 == 'Singlet')

featuresPlot(sc = sc,pth1 = '/mnt/GLHu/RESULT/afterBatchrm.pdf', 
             pth2 = '/mnt/GLHu/RESULT/pointDistribution_afterQC.pdf')
save(sc, file =  '/mnt/GLHu/DATA/sc_maf5.RData')


# cell Cluster ---------------------------------------------------------------

dim = 1:30
sc  = sc %>% SCTransform(vars.to.regress = 'pct.mt') %>% RunPCA()
# p = ElbowPlot(sc, ndims = 50)
# ggsave(p, filename = '/mnt/GLHu/RESULT/p.pdf')
sc = sc %>% RunUMAP(dims = dim) %>% FindNeighbors() %>% FindClusters(resolution = 0.8)
sc@meta.data$seurat_clusters %>% table
sc = sc %>% FindClusters(resolution = 1)

p1 <- DimPlot(sc, reduction = "umap", group.by = "samp")
p2 = DimPlot(sc, reduction = "umap", group.by = 'seurat_clusters',label = T,label.size = 3)
# p2 = DimPlot(sc, reduction = "umap", group.by = 'cellType',label = T,label.size = 3)
p = p1+p2
# ggsave(filename = "/mnt/GLHu/RESULT/UMAPannot.pdf", p, width = 12, height = 5)
ggsave(filename = "/mnt/GLHu/UMAP.pdf", p, width = 12, height = 5)

my_Heatmap(sc,pth = '/mnt/GLHu/HmDistributionFeaatures.pdf')

my_treeSim = function(exp,meta,col){
  foo = lapply(unique(meta[,col]), function(x){
    fin = data.frame(meanExp= rowMeans(exp[,meta[meta[,col] == x,]$cellID]),
                     row.names = row.names(exp))
  })  %>% list.cbind() %>% data.frame()
  names(foo) = paste('cluster', unique(meta[,col]))
  p = pheatmap(foo, show_rownames = F,color = colorRampPalette(colors = c("blue",'blue','white',"white",'white','red',"red"))(100))
  p = as.ggplot(p)
  ggsave(p, filename = '/mnt/GLHu/hclust.pdf', width = 10, height = 6)
}

my_treeSim(exp = sc@assays$SCT@scale.data, meta = sc@meta.data[,c('cellID','seurat_clusters')],col = 'seurat_clusters')




# FindMarkers -------------------------------------------------------------
refMarkers = read.csv('/mnt/GLHu/DATA_win/refMarkers_GLH.csv', row.names = 1)

markersAll = FindAllMarkers(sc, method = 'MAST')
markersAll = markersAll[markersAll$p_val_adj<0.05,] 
write.csv(markersAll, '/mnt/GLHu/data/markersAll.csv')
markersAll_P = markersAll[markersAll$avg_log2FC>0,]
markersAll_P = left_join(markersAll_P,refMarkers, by = 'gene') %>% unique 

write.csv(markersAll_P, '/mnt/GLHu/DATA/markers_p.csv')

my_topMarkers = function(mks,topn = 1:3){
  mks$pct.sig = mks$pct.1 - mks$pct.2
  fin = lapply(unique(mks$cluster), function(x){
    tar = mks[mks$cluster == x,]
    tar = tar[order(tar$avg_log2FC,decreasing = T),]
    t1 = tar[topn,]
    tar = tar[order(mks$pct.sig,decreasing = T ),]
    t2 = tar[topn,]
    tar = unique(rbind(t1,t2))
    return(tar)
  }) %>% list.rbind() %>% data.frame()
  return(fin)
}
tmp = my_topMarkers(mks = unique(markersAll_P[,1:7]), topn = 1:5)
tmp = na.omit(tmp)

p = VlnPlot(sc, tmp$gene, stack = T, combine = F, flip = T,  same.y.lims = T) + NoLegend()
ggsave(p, filename = '/mnt/GLHu/RESULT/markers.pdf', width = 7, height = 10)

my_clusterMarkers(tmp, pth = '/mnt/GLHu/RESULT/')
write.csv(tmp, '/mnt/GLHu/DATA/markers_candidate.csv')


# read markers candidate --------------------------------------------------
markersCandidate = read.csv('/mnt/GLHu/DATA_win/markers_candidate.csv')
markersCandidate$cellType %>% unique 
meta = sc@meta.data[,1:16]
meta$seurat_clusters = as.numeric(as.character(meta$seurat_clusters))
names(markersCandidate)[4] = 'seurat_clusters'
meta = left_join(meta,unique(markersCandidate[,c(4,7)])) 
rownames(meta) = meta$cellID
meta$cellType = factor(meta$cellType, levels = levels)
sc@meta.data = meta
Idents(sc) = sc@meta.data$cellType


p1 <- DimPlot(sc, reduction = "umap", group.by = "cellType", label = T, label.size = 2)
p2 = DimPlot(sc, reduction = "umap", group.by = 'seurat_clusters',label = T,label.size = 3)
# p2 = DimPlot(sc, reduction = "umap", group.by = 'cellType',label = T,label.size = 3)
p = p1+p2
# ggsave(filename = "/mnt/GLHu/RESULT/UMAPannot.pdf", p, width = 12, height = 5)
ggsave(filename = "/mnt/GLHu/UMAP.pdf", p, width = 12, height = 5)

freq =table(Idents(sc)) %>% data.frame()
names(freq) = c('cellType', 'cellNumber')
p = ggplot(freq, aes(x = cellType, y = cellNumber)) + geom_bar(stat = 'identity') + theme_common
ggsave(p, filename = '/mnt/GLHu/RESULT/freq.pdf', width = 10)

sorted = data.frame(cellType = levels, annot_major  = levels)
sorted = left_join(sorted, markersCandidate)


p = VlnPlot(sc, sorted$gene, stack = T, combine = F, flip = T, assay = 'RNA', same.y.lims = F) + NoLegend()
ggsave(p, filename = '/mnt/GLHu/RESULT/markers_all.pdf', width = 7, height = 10)



sum(freq[grep('Ex', freq$cellType),]$cellNumber)/sum(freq$cellNumber)
sum(freq[grep('Inh', freq$cellType),]$cellNumber)/sum(freq$cellNumber)


annot_GLH = data.frame(markers = c('SLC17A7','CADPS2', 'CACNA2D1', 'PHLDB2', 'KCIBIB', 'PDE1A', 'GRB14', 'ND4', 'GAD1', 'PVALB', 'VIP', 'SST', 'LAMP5', 'MAG', 'COL9A1', 'SLC1A2', 'STAB1', 'A2M'))
p = VlnPlot(sc, annot_GLH$markers, stack = T, combine = F, flip = T, same.y.lims = F) + NoLegend()
ggsave(p, filename = '/mnt/GLHu/RESULT/markers_GLH.pdf', width = 7, height = 7)



metaV5 = sc@meta.data
write.csv(metaV5,'/mnt/GLHu/DATA/metaV5.csv')


save(sc, file = '/mnt/GLHu/DATA/sc_maf5.RData')







