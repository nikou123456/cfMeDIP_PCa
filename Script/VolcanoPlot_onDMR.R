#==========================================================
# 
#Goal:  Volcano plot 
#
#                 
#==========================================================
rm(list=ls())
options(stringsAsFactors = FALSE)
library(dplyr)
library(ggplot2)
library(DESeq2)
source("./Script/function.R")
out.path <- "./Figure"
if(!dir.exists(out.path)){
  dir.create(out.path)
}


############################
deseq.result <- readRDS("/path/DMR_VPC67_vs_CPC30_with_age_control_raw_deseq.Rdata")
deseq.result <- as.data.frame(results(deseq.result,alpha=0.05))
index1 <- which((deseq.result$log2FoldChange)>=1  & deseq.result$padj<=0.05)
index2 <- which((deseq.result$log2FoldChange)<=(-1)  & deseq.result$padj<=0.05)

plot.volcano<-function(de.data=NULL,p_theme=NULL,padj.cutoff=0.05,log2FC.cutoff=1){
  p<-ggplot(de.data,aes(x=log2FoldChange,y=-log10(padj),col="#999999"))+
    geom_point(color="#999999",alpha=0.5)+
    geom_point(data=subset(de.data,padj<=padj.cutoff & (log2FoldChange)<(-log2FC.cutoff)),
               color="#377EB8",alpha=0.5)+
    geom_point(data=subset(de.data,padj<=padj.cutoff & (log2FoldChange)>log2FC.cutoff),
               color="#E41A1C",alpha=0.5)+p_theme
  p1<- p+geom_hline(yintercept = -log10(0.05), linetype="dotted")+
    geom_vline(xintercept = c(-1, 1), linetype="dotted") +
    labs(x=bquote(~Log[2]~fold~change),y=bquote(~-Log[10]~q~value))+
    p_theme
  return(p1)
}

p<- plot.volcano(de.data=deseq.result,p_theme=p_theme,padj.cutoff=0.05,log2FC.cutoff=1)
tiff("Vocanoplot_DMR_VPC_CPC.tiff")
print(p)
dev.off()