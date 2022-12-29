library(reshape2)
library(ggplot2)
library(WGCNA)




# -------------------------------------------------------------------
# 使用左右脑合并的数据构建WGCNA 网络
dir.create("./latest_regionalMean_WGCNA") # 共表达网络的相关结果都存储到该路径下

enableWGCNAThreads(8)
# 确定软阈值
powers = c(c(1:10), seq(from = 12, to=30, by=2))
type = "signed"
# Choose a set of soft-thresholding powers
sft = pickSoftThreshold(t(vsd4.ncx.mat.rmbatch.region),
                        powerVector = powers,
                        networkType = type,
                        corFnc=bicor,
                        verbose = 5)
pdf("./latest_regionalMean_WGCNA/estimate_softpower_vsd4_ncx_rmbatch_23613genes100regions.pdf",width=9.42,height = 4.59)
par(mfrow = c(1,2))
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = "Scale independence");
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()
sft$powerEstimate
# 构建网络
mergingThresh=0.25;softPower=12
vsd4_ncx_rmbatch_100region_net = blockwiseModules(t(vsd4.ncx.mat.rmbatch.region),corType="bicor",
                                                  maxBlockSize=24000,networkType="signed",
                                                  power=softPower,minModuleSize=50,
                                                  mergeCutHeight=mergingThresh,numericLabels=TRUE,
                                                  saveTOMs=TRUE,pamRespectsDendro=FALSE,nThreads=8,
                                                  saveTOMFileBase="./latest_regionalMean_WGCNA/vsd4_ncx_rmbnatch_23613genes100regions_singed_power12_wgcna")

vsd4_ncx_rmbatch_100region_net_mColors=vsd4_ncx_rmbatch_100region_net$colors
vsd4_ncx_rmbatch_100region_net_mColors = labels2colors(vsd4_ncx_rmbatch_100region_net_mColors)
table(vsd4_ncx_rmbatch_100region_net_mColors) # 16 modules
blocknumber=1
datColors=data.frame(vsd4_ncx_rmbatch_100region_net_mColors)[vsd4_ncx_rmbatch_100region_net$blockGenes[[blocknumber]],]
# 画最初始的网络图
plotDendroAndColors(vsd4_ncx_rmbatch_100region_net$dendrograms[[blocknumber]],
                    colors=datColors,
                    dendroLabels=FALSE,
                    hang =0.03,
                    addGuide=TRUE,
                    guideHang=0.05,
                    groupLabels=c("Module Colors"),
                    addTextGuide=TRUE,
                    main="vsd4_ncx_rmbatch_23613genes100regions_net Cluster Dendeogram")

# 合并模块
MEDissThres=0.2 # 对应correlation 0.8
vsd4_ncx_rmbatch_100region_net_merged_net=mergeCloseModules(t(vsd4.ncx.mat.rmbatch.region), 
                                                                 vsd4_ncx_rmbatch_100region_net_mColors,
                                                                 cutHeight = MEDissThres, 
                                                                 verbose = 3)
vsd4_ncx_rmbatch_100region_net_merged_mColors1=vsd4_ncx_rmbatch_100region_net_merged_net$colors;
table(vsd4_ncx_rmbatch_100region_net_merged_mColors1)
vsd4_ncx_rmbatch_100region_net_merged_ME1=vsd4_ncx_rmbatch_100region_net_merged_net$newMEs;

#再使用KME merge modules
vsd4_ncx_rmbatch_100region_net_merged_mColors2=
  moduleMergeUsingKME(t(vsd4.ncx.mat.rmbatch.region),vsd4_ncx_rmbatch_100region_net_merged_mColors1,
  threshPercent=90, mergePercent=40, omitColors=c("grey"))
vsd4_ncx_rmbatch_100region_net_merged_mColors2=vsd4_ncx_rmbatch_100region_net_merged_mColors2$moduleColors;
table(vsd4_ncx_rmbatch_100region_net_merged_mColors2)
plotDendroAndColors(vsd4_ncx_rmbatch_100region_net$dendrograms[[blocknumber]], 
                    cbind(vsd4_ncx_rmbatch_100region_net_mColors, 
                          vsd4_ncx_rmbatch_100region_net_merged_mColors1,
                          vsd4_ncx_rmbatch_100region_net_merged_mColors2),
                    c("Auto", "Merged","KMEmerged"), dendroLabels = FALSE, hang=0.03,
                    main="vsd4_ncx_rmbatch_23613genes100regions_net Cluster Dendeogram ",
                    addGuide = TRUE, guideHang=0.05)


annoRow=data.frame(row.names = ncx_regionOrder,
                   Lobe= ncx_regionOrder_df[ncx_regionOrder,3])
#bk=c(seq(-0.31,0,by=0.01),seq(0.01,0.83,by=0.01))
removeM=which(colnames(vsd4_ncx_rmbatch_100region_net_merged_ME1)=="MEgrey")

pheatmap(t(vsd4_ncx_rmbatch_100region_net_merged_ME1[ncx_regionOrder,-removeM]),
         cluster_cols = FALSE,
         #color = c(colorRampPalette(colors = c("blue","lightyellow"))(40),colorRampPalette(colors = c("lightyellow","red"))(80)),
         #breaks = bk,
         cellwidth = 10, cellheight = 10, 
         color = colorRampPalette(colors = c("blue","lightyellow","red"))(100),
         annotation_col =annoRow)




##---------------------------------------------------------------------------

#--- 尝试2 过滤低表达的基因
vsd4.ncx.mat.rmbatch.region.sd=rowSds(vsd4.ncx.mat.rmbatch.region)
exprSet4.ncx.rowMeans=rowMeans(exprSet4.ncx)
exprSet4.ncx.rowSd=rowSds(exprSet4.ncx)

# 只用蛋白编码加lincRNA，并且去掉线粒体上的基因
keepGenes=(geneAnnot$Gene.type=="protein_coding" | geneAnnot$Gene.type=="lincRNA" | geneAnnot$Gene.type=="miRNA" ) & (geneAnnot$Chromosome.scaffold.name!="MT")
keepGenes=geneAnnot$Gene.stable.ID[keepGenes]
keepGenes=intersect(keepGenes, row.names(vsd4.ncx.mat.rmbatch.region)) # 21023 genes
vsd4.ncx.mat.rmbatch.region.filtered=vsd4.ncx.mat.rmbatch.region[keepGenes,]
vsd4.ncx.mat.rmbatch.region.filtered.rowSd=rowSds(vsd4.ncx.mat.rmbatch.region.filtered)
vsd4.ncx.mat.rmbatch.region.filtered.rowMean=rowMeans(vsd4.ncx.mat.rmbatch.region.filtered)

bottom5percent=quantile(vsd4.ncx.mat.rmbatch.region.filtered.rowMean,probs=seq(0,1,by=0.05))[2]
keepRows=(vsd4.ncx.mat.rmbatch.region.filtered.rowMean>=bottom5percent)
vsd4.ncx.mat.rmbatch.region.filtered=vsd4.ncx.mat.rmbatch.region.filtered[keepRows,] # 19971

# 确定软阈值
powers = c(c(1:10), seq(from = 12, to=30, by=2))
type = "signed"
# Choose a set of soft-thresholding powers
sft = pickSoftThreshold(t(vsd4.ncx.mat.rmbatch.region.filtered),
                        powerVector = powers,
                        networkType = type,
                        corFnc=bicor,
                        verbose = 5)
pdf("./latest_regionalMean_WGCNA/estimate_softpower_vsd4_ncx_rmbatch_19971genes100regions.pdf",width=9.42,height = 4.59)
par(mfrow = c(1,2))
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = "Scale independence");
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()
sft$powerEstimate
# 构建网络
mergingThresh=0.25;softPower=12
vsd4_ncx_rmbatch_100region_filtered_net = blockwiseModules(t(vsd4.ncx.mat.rmbatch.region.filtered),corType="bicor",
                                                  maxBlockSize=24000,networkType="signed",
                                                  power=softPower,minModuleSize=50,
                                                  mergeCutHeight=mergingThresh,numericLabels=TRUE,
                                                  saveTOMs=TRUE,pamRespectsDendro=FALSE,nThreads=8,
                                                  saveTOMFileBase="./latest_regionalMean_WGCNA/vsd4_ncx_rmbnatch_19971genes100regions_singed_power12_wgcna")

vsd4_ncx_rmbatch_100region_filtered_net_mColors=vsd4_ncx_rmbatch_100region_filtered_net$colors
vsd4_ncx_rmbatch_100region_filtered_net_mColors = labels2colors(vsd4_ncx_rmbatch_100region_filtered_net_mColors)
table(vsd4_ncx_rmbatch_100region_filtered_net_mColors) # 16 modules
blocknumber=1
datColors=data.frame(vsd4_ncx_rmbatch_100region_filtered_net_mColors)[vsd4_ncx_rmbatch_100region_filtered_net$blockGenes[[blocknumber]],]
# 画最初始的网络图
plotDendroAndColors(vsd4_ncx_rmbatch_100region_filtered_net$dendrograms[[blocknumber]],
                    colors=datColors,
                    dendroLabels=FALSE,
                    hang =0.03,
                    addGuide=TRUE,
                    guideHang=0.05,
                    groupLabels=c("Module Colors"),
                    addTextGuide=TRUE,
                    main="vsd4_ncx_rmbatch_19971genes100regions_net Cluster Dendeogram")

# 合并模块
MEDissThres=0.2 # 对应correlation 0.8
vsd4_ncx_rmbatch_100region_filtered_net_merged_net=mergeCloseModules(t(vsd4.ncx.mat.rmbatch.region.filtered), 
                                                            vsd4_ncx_rmbatch_100region_filtered_net_mColors,
                                                            cutHeight = MEDissThres, 
                                                            verbose = 3)
vsd4_ncx_rmbatch_100region_filtered_net_merged_mColors1=vsd4_ncx_rmbatch_100region_filtered_net_merged_net$colors;
table(vsd4_ncx_rmbatch_100region_filtered_net_merged_mColors1)
vsd4_ncx_rmbatch_100region_filtered_net_merged_ME1=vsd4_ncx_rmbatch_100region_filtered_net_merged_net$newMEs;

#再使用KME merge modules
vsd4_ncx_rmbatch_100region_filtered_net_merged_mColors2=
  moduleMergeUsingKME(t(vsd4.ncx.mat.rmbatch.region.filtered),vsd4_ncx_rmbatch_100region_filtered_net_merged_mColors1,
                      threshPercent=90, mergePercent=80, omitColors=c("grey"))
vsd4_ncx_rmbatch_100region_filtered_net_merged_mColors2=vsd4_ncx_rmbatch_100region_filtered_net_merged_mColors2$moduleColors;
table(vsd4_ncx_rmbatch_100region_filtered_net_merged_mColors2)
plotDendroAndColors(vsd4_ncx_rmbatch_100region_filtered_net$dendrograms[[blocknumber]], 
                    cbind(vsd4_ncx_rmbatch_100region_filtered_net_mColors, 
                          vsd4_ncx_rmbatch_100region_filtered_net_merged_mColors1,
                          vsd4_ncx_rmbatch_100region_filtered_net_merged_mColors2),
                    c("Auto", "Merged","KMEmerged"), dendroLabels = FALSE, hang=0.03,
                    main="vsd4_ncx_rmbatch_19971genes100regions_net Cluster Dendeogram ",
                    addGuide = TRUE, guideHang=0.05)

annoRow=data.frame(row.names = ncx_regionOrder,
                   Lobe= ncx_regionOrder_df[ncx_regionOrder,3])
#bk=c(seq(-0.31,0,by=0.01),seq(0.01,0.83,by=0.01))
removeM=which(colnames(vsd4_ncx_rmbatch_100region_filtered_net_merged_ME1)=="MEgrey")
pheatmap(t(vsd4_ncx_rmbatch_100region_filtered_net_merged_ME1[ncx_regionOrder,-removeM]),
         cluster_cols = FALSE,
         #color = c(colorRampPalette(colors = c("blue","lightyellow"))(40),colorRampPalette(colors = c("lightyellow","red"))(80)),
         #breaks = bk,
         cellwidth = 10, cellheight = 10, 
         color = colorRampPalette(colors = c("blue","lightyellow","red"))(100),
         annotation_col =annoRow)



##--- 23a, 31, MT都去掉
keepRegions=setdiff(ncx_regionOrder, c("23a","31","MT"))
vsd4.ncx.mat.rmbatch.region.filtered2=vsd4.ncx.mat.rmbatch.region.filtered[,keepRegions]

# 确定软阈值
powers = c(c(1:10), seq(from = 12, to=30, by=2))
type = "signed"
# Choose a set of soft-thresholding powers
sft = pickSoftThreshold(t(vsd4.ncx.mat.rmbatch.region.filtered2),
                        powerVector = powers,
                        networkType = type,
                        corFnc=bicor,
                        verbose = 5)
pdf("./latest_regionalMean_WGCNA/estimate_softpower_vsd4_ncx_rmbatch_19971genes97regions.pdf",width=9.42,height = 4.59)
par(mfrow = c(1,2))
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = "Scale independence");
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()
sft$powerEstimate
# 构建网络
mergingThresh=0.25;softPower=12
vsd4_ncx_rmbatch_97region_filtered_net = blockwiseModules(t(vsd4.ncx.mat.rmbatch.region.filtered2),corType="bicor",
                                                           maxBlockSize=24000,networkType="signed",
                                                           power=softPower,minModuleSize=50,
                                                           mergeCutHeight=mergingThresh,numericLabels=TRUE,
                                                           saveTOMs=TRUE,pamRespectsDendro=FALSE,nThreads=8,
                                                           saveTOMFileBase="./latest_regionalMean_WGCNA/vsd4_ncx_rmbnatch_19971genes97regions_singed_power12_wgcna")

vsd4_ncx_rmbatch_97region_filtered_net_mColors=vsd4_ncx_rmbatch_97region_filtered_net$colors
vsd4_ncx_rmbatch_97region_filtered_net_mColors = labels2colors(vsd4_ncx_rmbatch_97region_filtered_net_mColors)
table(vsd4_ncx_rmbatch_97region_filtered_net_mColors) # 20 modules
blocknumber=1
datColors=data.frame(vsd4_ncx_rmbatch_97region_filtered_net_mColors)[vsd4_ncx_rmbatch_97region_filtered_net$blockGenes[[blocknumber]],]
# 画最初始的网络图
pdf("./latest_regionalMean_WGCNA/vsd4_ncx_rmbnatch_19971genes97regions_singed_power12_Cluster_Dendeogram.pdf", width=9.42,height=6.82)
plotDendroAndColors(vsd4_ncx_rmbatch_97region_filtered_net$dendrograms[[blocknumber]],
                    colors=datColors,
                    dendroLabels=FALSE,
                    hang =0.03,
                    addGuide=TRUE,
                    guideHang=0.05,
                    groupLabels=c("Module Colors"),
                    addTextGuide=TRUE,
                    main="vsd4_ncx_rmbatch_19971genes97regions_net Cluster Dendeogram")
dev.off()
# 合并模块
MEDissThres=0.2 # 对应correlation 0.8
vsd4_ncx_rmbatch_97region_filtered_net_merged_net=mergeCloseModules(t(vsd4.ncx.mat.rmbatch.region.filtered2), 
                                                                     vsd4_ncx_rmbatch_97region_filtered_net_mColors,
                                                                     cutHeight = MEDissThres, 
                                                                     verbose = 3)
vsd4_ncx_rmbatch_97region_filtered_net_merged_mColors1=vsd4_ncx_rmbatch_97region_filtered_net_merged_net$colors;
table(vsd4_ncx_rmbatch_97region_filtered_net_merged_mColors1)
vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1=vsd4_ncx_rmbatch_97region_filtered_net_merged_net$newMEs;
write.table(vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1,file="./latest_regionalMean_WGCNA/vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1.xls",sep="\t", quote=FALSE)

#再使用KME merge modules
vsd4_ncx_rmbatch_97region_filtered_net_merged_mColors2=
  moduleMergeUsingKME(t(vsd4.ncx.mat.rmbatch.region.filtered),vsd4_ncx_rmbatch_97region_filtered_net_merged_mColors1,
                      threshPercent=90, mergePercent=80, omitColors=c("grey"))
vsd4_ncx_rmbatch_97region_filtered_net_merged_mColors2=vsd4_ncx_rmbatch_97region_filtered_net_merged_mColors2$moduleColors;
table(vsd4_ncx_rmbatch_100region_filtered_net_merged_mColors2)
plotDendroAndColors(vsd4_ncx_rmbatch_97region_filtered_net$dendrograms[[blocknumber]], 
                    cbind(vsd4_ncx_rmbatch_97region_filtered_net_mColors, 
                          vsd4_ncx_rmbatch_97region_filtered_net_merged_mColors1,
                          vsd4_ncx_rmbatch_97region_filtered_net_merged_mColors2),
                    c("Auto", "Merged","KMEmerged"), dendroLabels = FALSE, hang=0.03,
                    main="vsd4_ncx_rmbatch_19971genes97regions_net Cluster Dendeogram ",
                    addGuide = TRUE, guideHang=0.05)

# 模块之间的相似性
pheatmap(cor(vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1[,-21]))






annoRow=data.frame(row.names = keepRegions,
                   Lobe= ncx_regionOrder_df[keepRegions,3])
#bk=c(seq(-0.31,0,by=0.01),seq(0.01,0.83,by=0.01))
removeM=which(colnames(vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1)=="MEgrey")
pheatmap(t(vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1[keepRegions,-removeM]),
         cluster_cols = TRUE, 
         #scale="row",
         cellheight = 10,cellwidth = 10,
         color = colorRampPalette(colors = c("blue","lightyellow","red"))(100),
         #filename = "./latest_regionalMean_WGCNA/vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1_hetmap.pdf",
         #width=18,height=6.83,
         annotation_col =annoRow)

# 根据module 特征值之间的相关性得到module的顺序
removeM=which(colnames(vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1)=="MEgrey")
region97_filtered_net_merged_ME1_pheatmap=pheatmap(cor(vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1[keepRegions,-removeM]), display_numbers = T)
region97_filtered_net_merged_moduleOrder=colnames(vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1[keepRegions,-removeM])
region97_filtered_net_merged_moduleOrder=region97_filtered_net_merged_moduleOrder[region97_filtered_net_merged_ME1_pheatmap$tree_col$order]

# 每个lobe单独看, 获得lobe内的顺序
ncx_lobes_order=lobes_order[-7]
lobeRegionOrder=c()
annoRow=data.frame(row.names = keepRegions,
                   Lobe= ncx_regionOrder_df[keepRegions,3])
for(i in 1:length(ncx_lobes_order)){
  tmpLobe=ncx_lobes_order[i]
  tmpLobe_regions=row.names(annoRow)[annoRow[,1]==tmpLobe]
  tmpLobe_tree=pheatmap(t(vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1[tmpLobe_regions,region97_filtered_net_merged_moduleOrder]),
           cluster_cols = TRUE, 
           cluster_rows = FALSE,
           #scale="row",
           cellheight = 10,cellwidth = 10,
           color = colorRampPalette(colors = c("blue","lightyellow","red"))(100),
           #filename = "./latest_regionalMean_WGCNA/vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1_hetmap.pdf",
           #width=18,height=6.83,
           #annotation_col =annoRow
           main=tmpLobe)
  tmpLobe_regions=tmpLobe_regions[tmpLobe_tree$tree_col$order]
  lobeRegionOrder=c(lobeRegionOrder, tmpLobe_regions)
}
lobeRegionOrder

annoRow=data.frame(row.names = lobeRegionOrder,
                   Lobe= ncx_regionOrder_df[lobeRegionOrder,3],
                   Subvison=ncx_regionOrder_df[lobeRegionOrder,5])
pheatmap(t(vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1[lobeRegionOrder,region97_filtered_net_merged_moduleOrder]),
         cluster_cols = FALSE, 
         cluster_rows = FALSE,
         #scale="row",
         cellheight = 10,cellwidth = 10,
         color = colorRampPalette(colors = c("blue","lightyellow","red"))(100),
         filename = "./latest_regionalMean_WGCNA/vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1_hetmap.pdf",
         width=18,height=6.83,
         annotation_col =annoRow)


region97_filtered_net_merged_moduleOrder2=str_split(region97_filtered_net_merged_moduleOrder,"ME",simplify = TRUE)[,2]



##-----将20个模块合并成4个模块，绘制4个模块的MEs热图
manualMerged_mColors2=vsd4_ncx_rmbatch_97region_filtered_net_merged_mColors1

# magenta turquoise grey60 lightgreen greenyellow合并成turquoise
manualMerged_mColors2=gsub("magenta","turquoise",manualMerged_mColors2,fixed=TRUE)
manualMerged_mColors2=gsub("grey60","turquoise",manualMerged_mColors2,fixed=TRUE)
manualMerged_mColors2=gsub("greenyellow","turquoise",manualMerged_mColors2,fixed=TRUE)
manualMerged_mColors2=gsub("lightgreen","turquoise",manualMerged_mColors2,fixed=TRUE)

# black tan  pink salmon brown合并成 brown
manualMerged_mColors2=gsub("tan","brown",manualMerged_mColors2,fixed=TRUE)
manualMerged_mColors2=gsub("pink","brown",manualMerged_mColors2,fixed=TRUE)
manualMerged_mColors2=gsub("salmon","brown",manualMerged_mColors2,fixed=TRUE)
manualMerged_mColors2=gsub("black","brown",manualMerged_mColors2,fixed=TRUE)

# lightcyan cyan green midnightblue purple合并成 cyan
manualMerged_mColors2=gsub("lightcyan","cyan",manualMerged_mColors2,fixed=TRUE)
manualMerged_mColors2=gsub("midnightblue","cyan",manualMerged_mColors2,fixed=TRUE)
manualMerged_mColors2=gsub("green","cyan",manualMerged_mColors2,fixed=TRUE)
manualMerged_mColors2=gsub("purple","cyan",manualMerged_mColors2,fixed=TRUE)

# blue yellow  royalblue lightyellow red 合并成 red
manualMerged_mColors2=gsub("lightyellow","red",manualMerged_mColors2,fixed=TRUE)
manualMerged_mColors2=gsub("yellow","red",manualMerged_mColors2,fixed=TRUE)
manualMerged_mColors2=gsub("royalblue","red",manualMerged_mColors2,fixed=TRUE)
manualMerged_mColors2=gsub("blue","red",manualMerged_mColors2,fixed=TRUE)

table(manualMerged_mColors2)

manualMerged_mColors2_ME=moduleEigengenes(t(vsd4.ncx.mat.rmbatch.region.filtered2), colors = manualMerged_mColors2)
manualMerged_mColors2_ME=manualMerged_mColors2_ME$eigengenes

write.table(manualMerged_mColors2_ME,
            file="vsd4_ncx_rmbatch_97region_filtered_net_manualMerged4_ME.xls",sep="\t", quote=FALSE)

removeM=which(colnames(manualMerged_mColors2_ME)=="MEgrey")
removeM

manualMerged_mColors2_ME_pheatmap=pheatmap(cor(manualMerged_mColors2_ME[keepRegions,-removeM]) )
manualMerged_mColors2_ME_moduleOrder=colnames(manualMerged_mColors2_ME[keepRegions,-removeM])
manualMerged_mColors2_ME_moduleOrder=manualMerged_mColors2_ME_moduleOrder[manualMerged_mColors2_ME_pheatmap$tree_col$order]


ncx_lobes_order=lobes_order[-7]
lobeRegionOrder=c()
annoRow=data.frame(row.names = keepRegions,
                   SaleemNetworks= ncx_regionOrder_df[keepRegions,5],
)
for(i in 1:length(ncx_lobes_order)){
  tmpLobe=ncx_lobes_order[i]
  tmpLobe_regions=row.names(annoRow)[annoRow[,1]==tmpLobe]
  tmpLobe_tree=pheatmap(t(vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1[tmpLobe_regions,region97_filtered_net_merged_moduleOrder]),
                        cluster_cols = TRUE, 
                        cluster_rows = FALSE,
                        #scale="row",
                        cellheight = 10,cellwidth = 10,
                        color = colorRampPalette(colors = c("blue","lightyellow","red"))(100),
                        #filename = "./latest_regionalMean_WGCNA/vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1_hetmap.pdf",
                        #width=18,height=6.83,
                        #annotation_col =annoRow
                        main=tmpLobe)
  tmpLobe_regions=tmpLobe_regions[tmpLobe_tree$tree_col$order]
  lobeRegionOrder=c(lobeRegionOrder, tmpLobe_regions)
}

lobeRegionOrder

annoRow=data.frame(row.names = lobeRegionOrder,
                   Lobe= ncx_regionOrder_df[lobeRegionOrder,3],
                   Subvison=ncx_regionOrder_df[lobeRegionOrder,5])
pheatmap(t(manualMerged_mColors2_ME[lobeRegionOrder,manualMerged_mColors2_ME_moduleOrder]),
         cluster_cols = FALSE, 
         cluster_rows = TRUE,
         #scale="row",
         cellheight = 10,cellwidth = 10,
         color = colorRampPalette(colors = c("blue","lightyellow","red"))(100),
         #filename = "./latest_regionalMean_WGCNA/vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1_hetmap.pdf",
         #width=18,height=6.83,
         annotation_col =annoRow)

manualMerged_mColors2_ME_moduleOrder2=str_split(manualMerged_mColors2_ME_moduleOrder,"ME",simplify = TRUE)[,2]


ncx_lobes_order=lobes_order[-7]
#lobeRegionOrder=c()
annoRow=data.frame(row.names = keepRegions,
                   SaleemNetworks= ncx_regionOrder_df[keepRegions,5],
)
for(i in 1:length(ncx_lobes_order)){
  tmpLobe=ncx_lobes_order[i]
  tmpLobe_regions=row.names(annoRow)[annoRow[,1]==tmpLobe]
  tmpLobe_tree=pheatmap(t(manualMerged_mColors2_ME[tmpLobe_regions,manualMerged_mColors2_ME_moduleOrder]),
                        cluster_cols = TRUE, 
                        cluster_rows = TRUE,
                        #scale="row",
                        cellheight = 10,cellwidth = 10,
                        color = colorRampPalette(colors = c("blue","lightyellow","red"))(100),
                        #filename = "./latest_regionalMean_WGCNA/vsd4_ncx_rmbatch_97region_filtered_net_merged_ME1_hetmap.pdf",
                        #width=18,height=6.83,
                        #annotation_col =annoRow
                        main=tmpLobe)
  tmpLobe_regions=tmpLobe_regions[tmpLobe_tree$tree_col$order]
  #lobeRegionOrder=c(lobeRegionOrder, tmpLobe_regions)
}







