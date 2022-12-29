# identified high expressed genes of each celltye in VISP

# set dirs
setwd("F:/Lab_info/wanglab/My_project/crab-eating macaque/mjj/Fw_ crab-eating macaque VISp single nucleus data analysis")

mfas5_highgenes_df = read.csv("markers_p_maf5.csv")

mfas5_celltypes = read.csv("markers_candidate.csv")


mfas5_highgenes_df$mfas5celltype = mfas5_celltypes$cellType[match(mfas5_highgenes_df$cluster, mfas5_celltypes$cluster)]

table(mfas5_highgenes_df$mfas5celltype)

mfas5_highgenes_df$mfas5BroadCelltype = mfas5_highgenes_df$mfas5celltype
mfas5_highgenes_df$mfas5BroadCelltype[grep("^Ex_",mfas5_highgenes_df$mfas5celltype) ] =  "ExcN"
mfas5_highgenes_df$mfas5BroadCelltype[grep("^Inh_",mfas5_highgenes_df$mfas5celltype) ] =  "IhnN"
mfas5_highgenes_df$mfas5BroadCelltype[grep("^oligo_",mfas5_highgenes_df$mfas5celltype) ] =  "oligo"
mfas5_highgenes_df$mfas5BroadCelltype[grep("^AST_",mfas5_highgenes_df$mfas5celltype) ] =  "AST"

table(mfas5_highgenes_df$mfas5BroadCelltype)

broadCellTypes =unique(mfas5_highgenes_df$mfas5BroadCelltype)
broadCellTypes
broadCellTypes_highExpGenes_list =lapply(broadCellTypes, function(x){
  tmp =  mfas5_highgenes_df[mfas5_highgenes_df$mfas5BroadCelltype==x, ]
  return(unique(tmp$gene[tmp$avg_log2FC>=1]))
})
names(broadCellTypes_highExpGenes_list) = broadCellTypes

mjj_mfas5_broadCellTypes_highExpGenes_list=broadCellTypes_highExpGenes_list

save(mjj_mfas5_broadCellTypes_highExpGenes_list,
     file="mjj_mfas5_broadCellTypes_highExpGenes_list.RData")

mfas5_geneAnno=read.csv("F:/Lab_info/wanglab/My_project/crab-eating macaque/ensembl_mfas5_gene.txt",sep="\t",header=TRUE)
for(i in 1:nrow(mfas5_geneAnno)){
  if(mfas5_geneAnno[i,2]==""){
    mfas5_geneAnno[i,2]=mfas5_geneAnno[i,1]
  }
}

mfas5_highgenes_df$geneID = mfas5_geneAnno$Gene.stable.ID[match(mfas5_highgenes_df$gene, mfas5_geneAnno$Gene.name)]

sub_mfas5_highgenes_df =  unique( mfas5_highgenes_df[,c(9,10,5,1)])
sub_mfas5_highgenes_df = sub_mfas5_highgenes_df[sub_mfas5_highgenes_df$avg_log2FC>=1,]

sub_mfas5_highgenes_df= unique(sub_mfas5_highgenes_df[,c(1,2,3)])
table(sub_mfas5_highgenes_df$mfas5BroadCelltype)

write.table(sub_mfas5_highgenes_df,
            file="mjj_mfas5_broadCellTypes_highExpGenes_list.xls", sep="\t", row.names = F, quote=F)
