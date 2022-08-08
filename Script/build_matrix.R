library(dplyr)
library(BiocParallel)
library(parallel)
library(doParallel)
###modified from /cluster/projects/hansengroup/wye/5mc/OICR/Sort_Bamfiles_redo/Build_count_rpkm.R
args <- commandArgs(trailingOnly = TRUE);
#data.path <- "/cluster/projects/hansengroup/sujunc/methylation/5mC/output/wig"
data.path <- args[1];
myname <- args[2];
setwd(data.path)
validation.files <- list.files(data.path,"*.Counts.wig.txt")
count.data <- data.frame(index=as.character(c(1:10318935)))
rpkm.data <- data.frame(index=as.character(c(1:10318935)))
for (file in validation.files){
  data <- read.table(file,sep="\t",header=FALSE)
  data$V1 <- as.integer(data$V1)
  data.sum <- summarise(data,totl.sum=sum(V1))
  data$V2 <- round((data$V1*10^9)/(300*data.sum$totl.sum),digits = 2)
  count.temp <- data.frame(value=data$V1)
  colnames(count.temp)<- sapply(strsplit(file,"\\_"),"[[",2 )
  rpkm.temp <- data.frame(value=data$V2)
  colnames(rpkm.temp)<- sapply(strsplit(file,"\\_"),"[[",2 )
  count.data <- cbind(count.data,count.temp)
  rpkm.data <- cbind(rpkm.data,rpkm.temp)
}
rownames(count.data) <- as.character(count.data$index)
rownames(rpkm.data) <- as.character(rpkm.data$index)
count.data$index <- NULL
rpkm.data$index <- NULL
saveRDS(count.data,file=paste0(Sys.Date() ,'_count_', myname, '.rds'))
saveRDS(rpkm.data,file=paste0(Sys.Date() ,'_rpkm_', myname, '.rds'))

