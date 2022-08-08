#########################################################
#
#' @title  Differential gene expression between GR-methylated high and low group
#
#
#########################################################

rm(list=ls())
source("./Script/function.R")
out.path <- "./Figure"
if(!dir.exists(out.path)){
  dir.create(out.path)
}

####################################################################################
#
#Felix Mets GR survival
#
####################################################################################

load("/path/gene_expression.Rdata")
#count.low means gene expression from patients with low GR methylaton
#count.high means gene expression from patients with high GR methylation
mets.de.gene.top10 <- DGE_deseq(count.low = count.low,
                                count.high = count.hight,cutoff = 10,sD1=NULL,
                                cd.label=c("GRLow","GRhigh"),
                                cd.levels= c("GRLow","GRhigh"),
                                out.names="DGE_mets_samples10_GRhigh_vs_GRlow_accroding_ep")


mets.de.gene.top10$gene <- rownames(mets.de.gene.top10)
mets.de.gene.top10.sub <- subset(mets.de.gene.top10,status!="Not Significant")

#-----------------------------------
# Plot figure 
p<- plot.vocano(de.data=mets.de.gene.top10,p_theme=p_theme,padj.cutoff=0.05,log2FC.cutoff=1)
#plot(p)
p1<- p+geom_text_repel(
  data = subset(mets.de.gene.top10, (abs(log2FoldChange)>5) & status!="Not Significant"&(-log10(padj)>3) | (-log10(padj)>8)),
  aes(label=gene),
  size=4,color="black",
  box.padding=unit(0.35,"lines"),
  point.padding=unit(0.3,"lines"))+ guides(color=FALSE)
tiff("Vocanoplot_DGEs_mets_WGBS_sample10_based_ep.tiff")
print(p1)
dev.off()

library(TCGAbiolinks)
Genelist <- as.character(mets.de.gene.top10$gene[which(mets.de.gene.top10$status=="Down regulated")]) 
Genelist <- as.character(mets.de.gene.top10$gene[which(mets.de.gene.top10$status=="Up regulated")])

plot_enrichGO(Genelist=Genelist,TFname=TFname,out.name=out.name,nBar=10)

index <- match(Genelist,wgbs.table$hgnc_symbol)
kegg.resu<- enrichKEGG(gene         = unique(wgbs.table$entrezgene_id[index]),
                       organism     = 'hsa',
                       keyType       = 'ncbi-geneid',  pAdjustMethod = "BH",
                       pvalueCutoff  = 0.5,
                       qvalueCutoff  = 1)

clusterProfiler::dotplot(kegg.resu,showCategory=10)+p_theme
