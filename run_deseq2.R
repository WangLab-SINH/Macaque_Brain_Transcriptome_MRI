
#load librarys
library(DESeq2)
library(sva)
library(cqn)
library("BiocParallel")
library(reshape2)
register(MulticoreParam(8))

# set work directory
setwd("./RAnalysis/deseq2/")

# load data
load("crab_819samples_28415genes_counts_phenset.RData")
load("phenSet4_105orderedRegions.RData")

phenSet4_region_df$Lobe[phenSet4_region_df$Region=="HC"]="Subcortical"
row.names(phenSet4_region_df)=phenSet4_region_df[,1]
phenSet4_region_df=phenSet4_region_df[order(phenSet4_region_df$Lobe,phenSet4_region_df$Y),]
regionOrder=phenSet4_region_df[,1]
# define lobe order
lobes_order=c("Frontal", "Parietal", "Temporal", "Cingulate", "Occipital",  "Insula", "Subcortical")
phenSet4_region_df=phenSet4_region_df[order(phenSet4_region_df$Lobe, phenSet4_region_df$Y),]
regionOrder=phenSet4_region_df[,1]
ncx_regionOrder=phenSet4_region_df[phenSet4_region_df$NeocortexRegion=="Yes",1]


# filter genes
keepgenes= unlist( apply(exprSet4,1,function(x){
  state=FALSE
  ratio=sum(x>=1)/length(x)
  if(ratio>=0.1){
    state=TRUE
  }
  return(state)
  
}) )
table(keepgenes) # 23605
exprSet4=exprSet4[keepgenes,]



# get gene annotation
geneAnnot=read.table("../ensembl_mfas5_gene.txt",header=TRUE,sep="\t")
row.names(geneAnnot)=geneAnnot[,1]
for(i in 1:nrow(geneAnnot)){
  if(geneAnnot[i,2]==""){
    geneAnnot[i,2]=geneAnnot[i,1]
  }
}
head(geneAnnot)

geneAnnot$Gtype <- geneAnnot[,3]
smallRNA <- c("Mt_rRNA","Mt_tRNA","miRNA","misc_RNA","rRNA","ribozyme","snRNA","snoRNA","scaRNA","sRNA")
IGorTR <- c("IG_C_gene","IG_V_gene", "TR_C_gene","TR_V_gene")
pseudogene <- c("processed_pseudogene","pseudogene")
mt.genes <- geneAnnot[geneAnnot[,5]=="MT",1]

geneAnnot$Gtype[geneAnnot[,3]%in%smallRNA] <- "smallRNA"
geneAnnot$Gtype[geneAnnot[,3]%in%IGorTR] <- "IGorTR"
geneAnnot$Gtype[geneAnnot[,3]%in%pseudogene] <- "pseudogene"
table(geneAnnot$Gtype)


# get gene GC and length
geneInfo=read.table("../mfas5.0.97.GC.lengths.tsv",header=TRUE,sep="\t")
head(geneInfo) # length and GC content


region_freq=table(phenSet4$Region)
regionList1=names(region_freq)[which(region_freq>=3)] # 102 regions
nRegionList1=length(regionList1)


### regional pairwise DEGs
### -----------------------------------------------------------------
#create dir to store pairwise results
dir.create("deseq2_phenSet4_102regions_pairwise")

# start DESeq2 DEGs
regionList1.pairwise.deseq.res.list=list()
nn=1
##------pair wise comparison
for(i in 1:(nRegionList1-1)){
  r1Name=regionList1[i]
  for(j in (i+1):nRegionList1){
    r2Name=regionList1[j]
    
    ##--- get expression for two group
    r1Index=which(phenSet4$Region==r1Name)
    r2Index=which(phenSet4$Region==r2Name)
    
    r1Sample=row.names(phenSet4)[r1Index]
    r2Sample=row.names(phenSet4)[r2Index]
    
    cat(r1Name, r2Name, length(r1Index), length(r2Index),"\n",sep="\t")
    
    #---get gene count
    r1Count=exprSet4[,r1Index]
    r2Count=exprSet4[,r2Index]
    geneCountTable=cbind(r1Count,r2Count)
    row.names(geneCountTable)=row.names(exprSet4)
    colnames(geneCountTable)=c(r1Sample,r2Sample)
    
    geneCountTable=geneCountTable[rowSums(geneCountTable)>0, ] # Filter non-expressed genes
    
    ##--- mean of count
    # count1=apply(r1Count,1,mean)
    # count2=apply(r2Count,1,mean)
    
    ##---- meta data
    bacth2=c(phenSet4$batch2[r1Index], phenSet4$batch2[r2Index])
    batch1=c(phenSet4$Batch[r1Index], phenSet4$Batch[r2Index])
    
    
    sampDesign = data.frame(row.names = colnames(geneCountTable),
                            batch1=batch1,
                            bacth2=bacth2,
                            condition = c(rep(r1Name,ncol(r1Count)),
                                          rep(r2Name,ncol(r2Count))))
    
    r1Batch2=unique(sampDesign[r1Sample,2])
    r2Batch2=unique(sampDesign[r2Sample,2])
    
    ##-----input from count matrix
    if( (length(r1Batch2)==1) & (length(r2Batch2)==1) ){
      dds = DESeqDataSetFromMatrix(countData = geneCountTable,
                                   colData = sampDesign,
                                   design = ~ condition);
    }
    else if(length(unique(bacth2))>=2){
      dds = DESeqDataSetFromMatrix(countData = geneCountTable,
                                   colData = sampDesign,
                                   design = ~ bacth2+condition);
    }
    else{
      dds = DESeqDataSetFromMatrix(countData = geneCountTable,
                                   colData = sampDesign,
                                   design = ~ condition);
    }
    
    dds$condition = factor(dds$condition,levels=c(r1Name,r2Name));
    
    ##-----cqn normlization
    dds1= dds;
    matchIndex=match(row.names(geneCountTable), row.names(geneInfo))
    cqnObject = cqn(geneCountTable,lengths=as.numeric(geneInfo[matchIndex,1]),
                    x= as.numeric(geneInfo[matchIndex,2]),
                    lengthMethod="smooth",sqn=T);
    cqnOffset = cqnObject$glm.offset;
    cqnNormFactors = exp(cqnOffset);
    cqnNormFactors = cqnNormFactors/mean(cqnNormFactors);
    normalizationFactors(dds1) = cqnNormFactors;
    dds1 = estimateDispersions(dds1);
    
    dds1.normalizedcounts=counts(dds1,normalized=TRUE)
    
    count1=apply(dds1.normalizedcounts[,c(1:ncol(r1Count))],1,mean)
    count2=apply(dds1.normalizedcounts[,c((ncol(r1Count)+1):ncol(dds1.normalizedcounts))],
                 1,mean)
    
    ddsMF = dds1;
    ddsMF = nbinomWaldTest(ddsMF);
    resMF = results(ddsMF);
    
    res = cbind(rownames(resMF),
                count1, count2,
                2**resMF[,2],
                resMF[,2],
                resMF[,5],
                resMF[,6] )
    
    ##---add geen name, gene type and assign change type
    myindex=match(as.character(res[,1]),as.character(geneAnnot[,1]));
    geneName=as.character(geneAnnot[myindex,2])
    geneType=as.character(geneAnnot[myindex,3]);
    res=cbind(res[,1],geneName,geneType,res[,2:ncol(res)]);
    colnames(res)=c("Geneid","geneName","geneType","baseMeanA(count)","baseMeanB(count)",
                    "foldChange","log2FoldChange", "pval","padj");
    if(grepl(pattern = "\\/",r1Name)) {
      r1Name1=gsub(pattern = "\\/","",r1Name)
    }
    else{
      r1Name1=r1Name
    }
    if(grepl(pattern = "\\/",r2Name)) {
      r2Name2=gsub(pattern = "\\/","",r2Name)
    }
    else{
      r2Name2=r2Name
    }
    ##-----output
    outFile = paste("./deseq2_phenSet4_105regions_pairwise_20210308/",r1Name1,"-vs-", r2Name2, ".all.DEX.xls",sep="");
    write.table(res,file=outFile,quote=F,row.names=F, col.names=T,sep="\t");
    
    # select sig results
    filter.res=res[res[,ncol(res)]<0.05,]
    
    regionList1.pairwise.deseq.res.list[[nn]] = filter.res
    names(regionList1.pairwise.deseq.res.list)[nn]=paste(r1Name,r2Name,sep="%vs%")
    nn=nn+1
  }
}
save(regionList1.pairwise.deseq.res.list, file="regionList1.pairwise.deseq.res.list.RData")
