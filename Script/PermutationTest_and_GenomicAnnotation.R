#==========================================================
# 
#Goal:  Permutation analysis and annotation on DMR
#
#                   
#                 
#==========================================================
rm(list=ls())
options(stringsAsFactors = FALSE)
library(dplyr)
library(BiocParallel)
library(parallel)
library(doParallel)
source("./Script/function.R")
out.path <- "./Figure"
if(!dir.exists(out.path)){
  dir.create(out.path)
}



##############################################################################################
#
#' @title Figure supplement Permutation test
#' @description Annotation plot
#
################################################################################

#===============================================================================
# loadoing genomic feature and background bins
load("/path/annotation_Granges_genome.Rdata")
annotations <- readRDS("/path/annotation.Rdata")
back.region <- GRanges.genome
back.region.reduce <- reduce(back.region)

#===============================================================================
#' @tile perfor permustation test
library(regioneR)


#===============================================================================

dmr.vpc.cpc <- readRDS("/path/DMR_VPC67_vs_CPC30_with_age_control_result_deseq.Rdata")
dmr.vpc.cpc.$status <- factor(dmr.vpc.cpc.$status,
                                   levels = c("hypo methylated","Not Significant","hyper methylated"),
                                   labels=c("hypo methylated","Not Significant","hyper methylated"))
dmr.vpc.cpc.$status <- as.character(dmr.vpc.cpc.$status)
pe.vpc.hypo <- subset(dmr.vpc.cpc.,status=="hypo methylated")
pe.vpc.hyper <- subset(dmr.vpc.cpc.,status=="hyper methylated")
hypo.region <- GRanges.genome[as.numeric(rownames(pe.vpc.hypo))]
hyper.region <- GRanges.genome[as.numeric(rownames(pe.vpc.hyper))]
hypo.region.reduce <- reduce(GRanges.genome[as.numeric(rownames(pe.vpc.hypo))])
hyper.region.reduce <- reduce(GRanges.genome[as.numeric(rownames(pe.vpc.hyper))])

library(annotatr)
annotations.list <- annotations
#-------HyperDMR after mering adjacent regions
hyper.reduce <-permutation(special=hyper.region.reduce,
                           back.region=back.region.reduce,
                           annotations.list=annotations.list,
                           annotation.region=annotation.region) 
plot_permutation(hyper.reduce)
ggsave(file="F2_D.png")

#-------HypoDMR after mering adjacent regions
hypo.reduce <-permutation(special=hypo.region.reduce,
                          back.region=back.region.reduce,
                          annotations.list=annotations.list,
                          annotation.region=annotation.region) 
plot_permutation(hypo.reduce)
ggsave(file="F2_E.png")


#-------HyperDMR on TR----------------------------------
TF.list<- readRDS("/path/MF-Chip-seq-26-region.Rdata")
tf.hyper.reduce <-permutation(special=hyper.region.reduce,
                              back.region=back.region.reduce,
                              annotations.list=TF.list,
                              annotation.region=TF.region)
plot_permutation.tf(tf.hyper.reduce)
ggsave(file="SF2_G.png")

#-------HypoDMR on TR----------------------------------
tf.hypo.reduce <-permutation(special=hypo.region.reduce,
                             back.region=back.region.reduce,
                             annotations.list=TF.list,
                             annotation.region=TF.region) 
plot_permutation.tf(tf.hypo.reduce)
ggsave(file="SF2_H.png")

