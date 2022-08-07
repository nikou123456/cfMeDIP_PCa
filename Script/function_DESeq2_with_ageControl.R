rm(list=ls())
options(stringsAsFactors = FALSE)
library(dplyr)
library(BiocParallel)
library(parallel)
library(doParallel)
library(scatterplot3d)
library(DESeq2)
cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)


#===============================================================================
#' @title defining function for Deseq2 with age control
runDESEQ <- function(data=NULL,sD1=NULL,cutoff=10,deseq.result=NULL,out.path=NULL,
                     levels.type=NULL,labels.type=NULL,age.need=NULL){ 
  setwd(out.path)
  if(is.null(sD1)){
    sD1 <- data.frame(
      row.names =colnames(data),
      condition=as.character(sapply(strsplit(colnames(data),"\\_"),"[[",1)) ,
      libType = rep(c("paried-end"), each=ncol(data)),
      age=age.need$age_type
    )
    
    #unique.condition <- unique(sD1$condition)
    sD1$condition <- factor(sD1$condition,
                            levels=levels.type,
                            labels=labels.type)
    sD1$libType <- factor(sD1$libType,levels=unique(sD1$libType))
    
    sD1$age <- factor( sD1$age,levels=unique(sD1$age))
    
    row.sum <- rowSums(data)
    index.zero <- which(row.sum>cutoff)
    data<- data[index.zero,]
    dds1 <- DESeqDataSetFromMatrix(countData = data,
                                   colData = sD1,
                                   design = ~ age + condition)
  }else{
    row.sum <- rowSums(data)
    index.zero <- which(row.sum>cutoff)
    data<- data[index.zero,]
    dds1 <- DESeqDataSetFromMatrix(countData = data,
                                   colData = sD1,
                                   design = ~libType +condition)
  }
  dds_deseq1 <- DESeq(dds1,full=design(dds1))
  out.names <- paste0(deseq.result,"_raw_deseq.Rdata")
  saveRDS(dds_deseq1,file=out.names)
  return(dds_deseq1)
}


deseq2data <- function(deseq.result=NULL,out.name=NULL){
  deseq.result <- as.data.frame(results(deseq.result,alpha=0.05))
  index1 <- which((deseq.result$log2FoldChange)>=1  & deseq.result$padj<=0.05)
  index2 <- which((deseq.result$log2FoldChange)<=(-1)  & deseq.result$padj<=0.05)
  deseq.result$status <- "Not Significant"
  deseq.result$status[index1]<- "hyper methylated"
  deseq.result$status[index2]<-"hypo methylated"
  print(table(deseq.result$status))
  deseq.result$log2FC_state <- "1-2"
  deseq.result$log2FC_state[which(abs(deseq.result$log2FoldChange)>=2  & abs(deseq.result$log2FoldChange)<=3)] <- "2-3"
  deseq.result$log2FC_state[which(abs(deseq.result$log2FoldChange)>=3  & abs(deseq.result$log2FoldChange)<=4)] <- "3-4"
  deseq.result$log2FC_state[which(abs(deseq.result$log2FoldChange)>=4)] <- ">4"
  #table(met.r.j$log2FC_state)
  #table(paste0(met.r.j$status,met.r.j$log2FC_state))
  saveRDS(deseq.result,file = paste0(out.name,"_result_deseq.Rdata"))
  return(deseq.result)
}


#####################################################################################################
#' @title  DMR between CPC and CPC with age control
#' @description
#'
data <- readRDS("cfmedip_profile.Rdata")
age <- readRDS("age_information")
index <- match(gsub("^Mets_|^Local_","",colnames(data)),age$sampleID)
age.need <- age[index,]
age.need$age_type <-as.character(cut(age.need$age,5))

#*******************************************************************************************###
out.path <- "./DMR_DESeq2"
setwd(out.path)
deseq.cpc.vpc <- runDESEQ(data=data,cutoff=10,sD1=NULL,
                          deseq.result = "vpc_vs_cpc"
                          ,out.path = out.path,
                          levels.type=c("Local","Mets"),
                          labels.type=c("CPC","VPC"),
                          age.need=age.need)
dmr.cpc.vpc <-deseq2data(deseq.result=deseq.cpc.vpc,
                         out.name="vpc_vs_cpc")
parallel::stopCluster(cl)