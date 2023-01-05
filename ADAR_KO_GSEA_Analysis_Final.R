#ADAR KO RNAseq Analysis
#9/13/21

rm(list=ls())

library(this.path)
library(Biobase)
library(tidyverse)
library(biomaRt)
library(fgsea)
library(pheatmap)
library(DESeq2)

count_cutoff=5
FDR_cutoff=0.05
NES_cutoff=1.8
set.seed(1)

#Load data
setwd(this.path::here())
exp=read.delim('mouse-gene.counts_ja157_ji083021_edit.csv', sep=',')
#Convert names
exp=rename(exp,c("X157.6"="epithelial_KO","X157.5"="endothelial_KO","X157.4"="immune_KO",
                 "X157.3"="epithelial_ctrl","X157.2"="endothelial_ctrl","X157.1"="immune_ctrl"))
#Set gene names as row names
rownames(exp)=exp[,1]
exp[,1]=NULL
#Normalize with DESeq
mat=DESeqDataSetFromMatrix(exp, data.frame(Samples=colnames(exp)), ~1)
mat=estimateSizeFactors(mat)
mat_norm=counts(mat, normalized=TRUE)
#Log2 normalize
exp=log2(mat_norm+1)
#Filter to remove low expressed transcripts (<count_cutoff counts in >50% of samples)
exp=exp[rowSums(exp>=count_cutoff)>ncol(mat_norm)/2,]

#Add gene names back into data frame for conversion
exp=data.frame('mouse_symbol'=rownames(exp),exp)
#Convert from mouse to human symbol
human=useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
mouse=useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org")
human_mouse_lookup=getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = exp$mouse_symbol , mart = mouse,
                          attributesL = c("entrezgene_id"), martL = human, uniqueRows=T)
colnames(human_mouse_lookup)=c('mouse_symbol','human_entrez')
exp=left_join(human_mouse_lookup,exp)
#Remove mouse names
exp$mouse_symbol=NULL
exp$GeneName=NULL
#For multiple human->mouse matches, collapse by average
exp=exp %>% group_by(human_entrez) %>% summarize_all(mean) %>% as.data.frame
#Compute FCs per cell type
FC_table=data.frame(human_entrez=exp$human_entrez, immune=exp$immune_KO-exp$immune_ctrl, endothelial=exp$endothelial_KO-exp$endothelial_ctrl,
                    epithelial=exp$epithelial_KO-exp$epithelial_ctrl)
#Write to file
#write.table(FC_table, 'KOvCtrl_FC_table_human.txt', sep='\t')

#Load KEGG/Reactome/Biocarta modules
KRB=read.delim('Human_ReactomeBiocartaKegg_v4.txt', header=F)
#Convert to list
temp=lapply(1:nrow(KRB), function(x) as.character(na.omit(as.numeric(KRB[x,-(1:2)]))))
names(temp)=KRB$V1
module_list=temp
rm(temp,KRB)
#Filter list to single database, remove database prefix
module_list=module_list[grep('REACTOME_',names(module_list))]
names(module_list)=str_replace_all(names(module_list),'REACTOME_','')
names(module_list)=str_replace_all(names(module_list),'_',' ')

#Remove modules without genes in dataset
module_list=module_list[!sapply(module_list, function(x) is_empty(intersect(x,FC_table$human_entrez)))]

#Run GSEA
comparisons=colnames(FC_table)[-1]
GSEA_res=lapply(comparisons, function(x) fgsea(module_list, setNames(FC_table[,x],FC_table$human_entrez)))
GSEA_NES=sapply(GSEA_res, function(x) {y=x$NES; y})
dimnames(GSEA_NES)=list(GSEA_res[[1]]$pathway, comparisons)
GSEA_p=sapply(GSEA_res, function(x) {y=x$padj; y})
dimnames(GSEA_p)=list(GSEA_res[[1]]$pathway, comparisons)
GSEA_res_sig=lapply(GSEA_res, function(x) {y=x[(abs(x$NES)>NES_cutoff & x$padj<FDR_cutoff),]; z=y[order(-y$NES),]; z})
#Sort all GSEA results
GSEA_res=lapply(GSEA_res, function(x) {y=x[order(-x$NES),]; y})

#Write results to text files
# for (i in 1:length(GSEA_res_sig)) {
#   # #Convert leading edge from list to character vector
#   # GSEA_res_sig[[i]]$leadingEdge=vapply(GSEA_res_sig[[i]]$leadingEdge, paste, collapse = ", ", character(1L))
#   # #Write file
#   # write.table(GSEA_res_sig[[i]], paste0(comparisons[i],'_GSEA_FDR05.txt'), sep='\t')
#   #Convert leading edge from list to character vector
#   GSEA_res[[i]]$leadingEdge=vapply(GSEA_res[[i]]$leadingEdge, paste, collapse = ", ", character(1L))
#   #Write file
#   write.table(GSEA_res[[i]], paste0(comparisons[i],'_KOvCtrl_GSEA_all.txt'), sep='\t')
# }

#Plot significant GSEA results in heatmaps
sig_BTM_all=unique(unlist(lapply(GSEA_res_sig, function(x) x$pathway)))
#Create df of sig BTMs in all comparisons
df=GSEA_NES[match(sig_BTM_all,rownames(GSEA_NES)),]
df_FDR=GSEA_p[match(sig_BTM_all,rownames(GSEA_NES)),]
#Replace nonsig values with 0
#df[df_FDR>=FDR_cutoff | abs(df)<NES_cutoff]=0
df[df_FDR>=FDR_cutoff]=0
#Sort by common/unique
temp=df
temp[temp>0]=1
temp[temp<0]=-1
temp=cbind(temp,rowSums(temp))
temp[temp==0]=-10
df_FDR=df_FDR[do.call(order,-data.frame(temp[,c(4,1,2,3)])),]
#Sig scores only
df=df[do.call(order,-data.frame(temp[,c(4,1,2,3)])),]
#All scores
#df=GSEA_NES[match(sig_BTM_all,rownames(GSEA_NES))[do.call(order,-data.frame(temp[,c(4,1,2,3)]))],]


#Set map limits
ph_lim=ceiling(max(abs(df), na.rm = TRUE)*2)/2
colorLS=colorRampPalette(colors = c("blue", "white", "red"))(n = 100)
breakLS=seq(from = -ph_lim, to = ph_lim, length.out = 100)
#Plotting with pheatmap
result=pheatmap(mat = df,
                color = colorLS,
                breaks = breakLS,
                cluster_cols = FALSE,
                cluster_rows = FALSE,
                show_colnames = TRUE,
                show_rownames = TRUE,
                fontsize = 10,
                cellheight = 15,
                cellwidth = 15,
                filename=paste0('GSEA_Reactome_FDR_',FDR_cutoff,'_NES_',NES_cutoff,'_countcut_',count_cutoff,'.pdf'))
