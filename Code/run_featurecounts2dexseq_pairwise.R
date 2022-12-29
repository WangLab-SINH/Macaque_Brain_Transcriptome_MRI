library(DEXSeq)
#library(sva)
#library(stringr)
library(BiocParallel)
library(dplyr)

#setwd("/picb/neurosys/LJ/crab-eating-macaque-transcriptome/RAnalysis/deseq2_20201202update/featurecounts2dexseq/")
setwd("F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/featurecounts2dexseq")

load("featurecounts2dexseq_exoncounts_819samples_20210118_input.RData")
exonCounts.phenSet$Lobe[exonCounts.phenSet$Region=="HC"]="Subcortical"
exonCounts.phenSet$Region[exonCounts.phenSet$Region%in%c("STGr","STGc")]="STG"

#source("/picb/neurosys/LJ/usecode/Subread_to_DEXSeq/load_SubreadOutput.R")
#gtf="/picb/neurosys/LJ/genome/crab-eating-macaque-mfas5/mfas5_ens97_flat.gtf"

source("load_SubreadOutput.R")
gtf="mfas5_ens97_flat.gtf"


regionFreq.table=table(exonCounts.phenSet$Region)
freq2region=names(regionFreq.table)[regionFreq.table>=2] 
regionList <- freq2region
nRegionList <- length(regionList)


dir.create("pairwise")

regionList.pairwise.b1.dexseq.list <- list()
nn<-1
for(i in 1:3){
	r1Name <- regionList[i]
	for(j in (i+1):nRegionList){
		r2Name <- regionList[j]
		cat(nn, "****", r1Name, r2Name,"****\n",sep="\t" )
		outfile=paste0("./pairwise/",r1Name,"vs",r2Name,".exonCountTable.dexseq.results.txt")
	if(file.exists(outfile)){
		dxr=read.csv(outfile, sep="\t", check.names=FALSE, header=TRUE, stringsAsFactors=FALSE)
		regionList.pairwise.b1.dexseq.list[[nn]] <- dxr
	}
	else{
	  r1Index=exonCounts.phenSet[exonCounts.phenSet$Region==r1Name,1]
    r2Index=exonCounts.phenSet[exonCounts.phenSet$Region==r2Name,1]
		
		r1Index=match(r1Index, colnames(exonCounts.dat))
		r2Index=match(r2Index, colnames(exonCounts.dat))			

    cat(r1Name,r2Name, length(r1Index), length(r2Index), "\n", sep="\t")

    exonCountTable <- exonCounts.dat[,c(r1Index,r2Index)]
		# add exon and gene infor
    exonCountTable <- data.frame(exonCounts.dat[,c(1:6)], exonCountTable, stringsAsFactors=FALSE, check.names=FALSE)

    outfile=paste0("./pairwise/",r1Name,"vs",r2Name,".exonCountTable.txt")
    write.table(exonCountTable, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)

    samp <- data.frame(row.names = colnames(exonCounts.dat)[c(r1Index,r2Index)], 
		condition = c( rep(r1Name,length(r1Index)), rep(r2Name, length(r2Index)) ) )
		dxd.fc <- DEXSeqDataSetFromFeatureCounts(outfile,flattenedfile = gtf,sampleData = samp)
		dxr = DEXSeq(dxd.fc,BPPARAM=MulticoreParam(workers=4))
		outfile=paste0("./pairwise/",r1Name,"vs",r2Name,".exonCountTable.dexseq.results.txt")
		write.table(as.data.frame(dxr), file=outfile,sep="\t", quote=FALSE, row.names=FALSE)
		
		regionList.pairwise.b1.dexseq.list[[nn]] <- dxr
		names(regionList.pairwise.b1.dexseq.list)[nn] =paste0(r1Name,"vs",r2Name)
		}
	
	nn <- nn+1
	}
	
}

i=1
nn<-1
  r1Name <- regionList[i]
  for(j in (i+1):nRegionList){
    r2Name <- regionList[j]
    cat(nn, "****", r1Name, r2Name,"****\n",sep="\t" )
    outfile=paste0("./pairwise/",r1Name,"vs",r2Name,".exonCountTable.dexseq.results.txt")
    if(file.exists(outfile)){
      dxr=read.csv(outfile, sep="\t", check.names=FALSE, header=TRUE, stringsAsFactors=FALSE)
      #regionList.pairwise.b1.dexseq.list[[nn]] <- dxr
    }
    else{
      r1Index=exonCounts.phenSet[exonCounts.phenSet$Region==r1Name,1]
      r2Index=exonCounts.phenSet[exonCounts.phenSet$Region==r2Name,1]
      
      r1Index=match(r1Index, colnames(exonCounts.dat))
      r2Index=match(r2Index, colnames(exonCounts.dat))			
      
      cat(r1Name,r2Name, length(r1Index), length(r2Index), "\n", sep="\t")
      
      exonCountTable <- exonCounts.dat[,c(r1Index,r2Index)]
      # add exon and gene infor
      exonCountTable <- data.frame(exonCounts.dat[,c(1:6)], exonCountTable, stringsAsFactors=FALSE, check.names=FALSE)
      
      outfile=paste0("./pairwise/",r1Name,"vs",r2Name,".exonCountTable.txt")
      write.table(exonCountTable, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
      
      samp <- data.frame(row.names = colnames(exonCounts.dat)[c(r1Index,r2Index)], 
                         condition = c( rep(r1Name,length(r1Index)), rep(r2Name, length(r2Index)) ) )
      dxd.fc <- DEXSeqDataSetFromFeatureCounts(outfile,flattenedfile = gtf,sampleData = samp)
      dxr = DEXSeq(dxd.fc,BPPARAM=MulticoreParam(workers=4))
      outfile=paste0("./pairwise/",r1Name,"vs",r2Name,".exonCountTable.dexseq.results.txt")
      write.table(as.data.frame(dxr), file=outfile,sep="\t", quote=FALSE, row.names=FALSE)
      
      #regionList.pairwise.b1.dexseq.list[[nn]] <- dxr
      #names(regionList.pairwise.b1.dexseq.list)[nn] =paste0(r1Name,"vs",r2Name)
    }
    
    nn <- nn+1
  }
  




save(regionList.pairwise.b1.dexseq.list, 
     file="featurecounts2dexseq_exoncounts_819samples_20210118_b1_dexseq_res_list.RData")
save.image("featurecounts2dexseq_exoncounts_819samples_20210118_b1_dexseq_workspace.RData")

