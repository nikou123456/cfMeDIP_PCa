rm(list=ls())
options(stringsAsFactors = FALSE)
library(CaSpER)
library(knitr)
library(BiocParallel)
library(parallel)
library(doParallel)
load("/path/CasPER/Prepared_CNVdata_update.RData")

data("hg19_cytoband")

annotation.bin$Chr <- annotation.bin$seqnames
annotation.bin$Chr <- gsub("chr","",annotation.bin$Chr)
rownames(annotation.bin) <- annotation.bin$Gene
loh <- readRDS("/path/CaSpER/loh.Rdata")
index <- match(colnames(data.jv)[grep("^V",colnames(data.jv))],names(loh))
loh.jv <- loh[index]
loh.name.mapping.jv <- data.frame(loh.name=names(loh.jv),
                                  sample.name=names(loh.jv))

control.sample.ids <- colnames(data.jv)[grep("^G|^J",colnames(data.jv))]

unique.chr <- unique(annotation.bin$seqnames)
for (i in seq(1,24,2)) {
  cl <- parallel::makeCluster(20)
  doParallel::registerDoParallel(cl)
  index.chr <- which(annotation.bin$seqnames %in% unique.chr[i:(i+1)])
  print(table(annotation.bin$seqnames[index.chr])) 
  object.jv <- CreateCasperObject(raw.data=data.jv[index.chr,],
                                  loh.name.mapping=loh.name.mapping.jv, 
                                  sequencing.type="bulk", 
                                  cnv.scale=3, loh.scale=3, 
                                  matrix.type="normalized", 
                                  expr.cutoff=4.5,
                                  log.transformed=TRUE,
                                  annotation=annotation.bin[index.chr,], 
                                  method="iterative", 
                                  loh=loh.jv,
                                  filter="median",
                                  genomeVersion="hg19",
                                  control.sample.ids=control.sample.ids, 
                                  cytoband=cytoband)
  
  final.objects.jv <- runCaSpER(object.jv, removeCentromere=T,
                                cytoband=cytoband, method="iterative")
  setwd("/path/CNV/CaSpER/")
  saveRDS(final.objects.jv,file=paste0("final.objects.jv.",unique.chr[i],"-",unique.chr[(i+1)] ,".Rdata"))
  parallel::stopCluster(cl)
}

