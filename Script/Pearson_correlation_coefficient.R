
#==========================================================
# 
#Goal: Cor matched patients and no-matched patients
#      
#      
#                   Wenbin Ye
#                 -- 20210707
#==========================================================
rm(list=ls())
options(stringsAsFactors = FALSE)
library(extrafont)
library(ggplot2)
library(export)
library(ggpubr)
source("./Script/function.R")
out.path <- "./Figure"
if(!dir.exists(out.path)){
  dir.create(out.path)
}


melted_cormat <- readRDS("./Data/Correlation_bet_tissue_liquid_VPC.Rdata")
my_comparisons <- list(c("No","Yes"))

ggplot(melted_cormat,aes(x=same, y=value))+
  geom_boxplot( outlier.shape = NA)+geom_jitter(position=position_jitter(0.2),color="gray")+
  labs(x=NULL,y="Pearson correlation coefficient")+ 
  p_theme+ guides(color=guide_legend(title=NULL))+theme(legend.position  ="none")+
  stat_n_text()+stat_mean_sd_text()+
  stat_compare_means(comparisons = my_comparisons,#label = "p.signif",
                     method = "wilcox.test",hide.ns=TRUE,size=4)
ggsave(file="S1_A.png",width=4,height=4)


melted_cormat.sub  <- readRDS("./Data/Correlation_bet_tissue_liquid_WCDT.Rdata")
my_comparisons <- list(c("No","Yes"))
ggplot(melted_cormat.sub,aes(x=same, y=value))+
  geom_boxplot( outlier.shape = NA)+geom_jitter(position=position_jitter(0.2),color="gray")+
  labs(x=NULL,y="Pearson correlation coefficient")+ 
  p_theme+ guides(color=guide_legend(title=NULL))+theme(legend.position  ="none")+
  stat_compare_means(comparisons = my_comparisons,#label = "p.signif",
                     method = "wilcox.test",hide.ns=TRUE)+
  stat_n_text()+stat_mean_sd_text()

ggsave(file="F1_B.png",width=4,height=4)






