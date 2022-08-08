#==========================================================
# 
#Goal:  GR infromation 
#
#                   Wenbin Ye
#                 -- 20210629
#==========================================================
rm(list=ls())
options(stringsAsFactors = FALSE)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(export)
library(extrafont)
colors <- brewer.pal(8,"Dark2")
source("./Script/function.R")
out.path <- "./Figure"
if(!dir.exists(out.path)){
  dir.create(out.path)
}


col.data <- readRDS("./Data/sample_information.Rdata")
vpc.gr <- readRDS("./Data/GR_VPC.Rdata")
#===============================================================================
#Only check VPC based on RPKM
ggplot(vpc.gr, 
       aes(x=ctDNA,y=log2(GR_RPKM+1)))+geom_point()+theme_classic()+
  p_theme+
  geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="% ctDNA",y=bquote(~Log[2]~RPKM))+guides(color=guide_legend(title=NULL))+
  theme( strip.text.x = element_text(size=12),
         legend.background = element_blank(),
         legend.position = c(0.9,0.9))+
  geom_text_repel(
    data = subset(col.data.rpkm, log2(GR_RPKM+1)>1),
    aes(label=simpleName),
    size=4,
    box.padding=unit(0.35,"lines"),
    point.padding=unit(0.3,"lines"))
ggsave(file="F3_B.png")


#================================================================================
# Survival analysis only for 67 CPV

vpc.clinical <- readRDS("./Data/VPC_clinical.Rdata")
vpc.clinical$GR_RPKM <- vpc.gr$GR_RPKM[match(vpc.clinical$Sample,vpc.gr$Sample)]
survival.data.sub <- subset(vpc.clinical,duplicated ==0)
survival.data.sub$OS.time <- as.numeric(survival.data.sub$`Published_Days to death or last followup`)
survival.data.sub$OS <- survival.data.sub$Published_Censored...14
survival.data.sub$OS[which(survival.data.sub$OS=="Yes")] <- 1
survival.data.sub$OS[which(survival.data.sub$OS=="No")] <- 2
survival.data.sub$OS <- as.numeric(survival.data.sub$OS)
survival.data.sub$RFS.time <- as.numeric(survival.data.sub$`Published_Days to progression or last followup`)
survival.data.sub$RFS <- survival.data.sub$Published_Censored...12
table(survival.data.sub$RFS )
survival.data.sub$RFS[which(survival.data.sub$RFS=="Yes")] <- 1
survival.data.sub$RFS[which(survival.data.sub$RFS=="No")] <- 2
survival.data.sub$RFS <- as.numeric(survival.data.sub$RFS)

survival.data.sub$GR_RPKM_log <- log2(survival.data.sub$GR_RPKM+1)
survival.data.sub$cluster <- "Low"
survival.data.sub$cluster[which(survival.data.sub$GR_RPKM_log>0)] <- "High"
#survival.data.sub$cluster[which(survival.data.sub$GR_RPKM_log<quantile(survival.data.sub$GR_RPKM_log,0.5))] <- "Low"

survival.data.sub$cluster <- factor(survival.data.sub$cluster,levels=c("High","Low"))
table(survival.data.sub$cluster)

library(survival)
library(survminer)
OS.fit <- survfit(Surv(OS.time, OS) ~ cluster, 
                  data = survival.data.sub)
ggsurvplot(OS.fit,size = 1,# change line size
           palette = c("#2E9FDF", "#E7B800"), # custom color palette
           conf.int = FALSE, # Add confidence interval
           legend.title = "GR",
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col ="strata", # Risk table color by groups
           ylab="Overall survival",
           xlab="Time (days)",
           legend.labs=c("High","Low"),
           ggtheme = p_theme# Change ggplot2 theme
)
ggsave(file="S3_B.png",width=4,height=5)


RFS.fit <- survfit(Surv(RFS.time, RFS) ~ cluster, 
                   data = survival.data.sub)
ggsurvplot(RFS.fit,size = 1,# change line size
           palette = c("#2E9FDF", "#E7B800"), # custom color palette
           conf.int = FALSE, # Add confidence interval
           legend.title = "GR",
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col ="strata", # Risk table color by groups
           ylab="Progression-free survival",
           xlab="Time (days)",
           legend.labs=c("High","Low"),
           ggtheme = p_theme# Change ggplot2 theme
)
ggsave(file="S3_C.png",width=4,height=5)



###########################################################################
#
#'@title Figure 3C-E Tissue and GR information
#
###########################################################################
####################################################################################
#
#Felix Mets GR survival_gene expression
#
####################################################################################
WCDT.clinical <- readRDS("./Data/GR_WCDT.Rdata")
ggplot(WCDT.clinical , 
       aes(x=beta,y=log2(tpm+1)))+geom_point()+theme_classic()+
  p_theme+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="Beta",y="Expression")+guides(color=guide_legend(title=NULL))+
  theme( strip.text.x = element_text(size=12),
         legend.background = element_blank(),
         legend.position = c(0.9,0.9))+ylim(0,8.1)+xlim(0,0.8)
ggsave(file="S3_G.png",width = 4, height = 4)


WCDT.clinical.ep <- WCDT.clinical[which(WCDT.clinical$tpm>=0),]
WCDT.clinical <- WCDT.clinical[which(!is.na(WCDT.clinical$raw.name)),]
WCDT.clinical$timepoint <- as.character(sapply(strsplit(WCDT.clinical$raw.name,"-"),"[[",3))
WCDT.clinical$cluster.beta<-"Mid"
WCDT.clinical$cluster.beta[which(WCDT.clinical$beta >as.numeric(quantile(WCDT.clinical$beta,0.5)))] <- "High"
WCDT.clinical$cluster.beta[which(WCDT.clinical$beta <=as.numeric(quantile(WCDT.clinical$beta,0.5)))] <- "Low"
WCDT.clinical$cluster.beta <- factor(WCDT.clinical$cluster.beta,level=c("High","Low"))

WCDT.clinical.ep$cluster.tpm<-"Mid"
WCDT.clinical.ep$cluster.tpm[which(WCDT.clinical.ep$tpm >as.numeric(quantile(WCDT.clinical.ep$tpm,0.5)))] <- "High"
WCDT.clinical.ep$cluster.tpm[which(WCDT.clinical.ep$tpm <=as.numeric(quantile(WCDT.clinical.ep$tpm,0.5)))] <- "Low"
WCDT.clinical.ep$cluster.tpm <- factor(WCDT.clinical.ep$cluster.tpm,level=c("High","Low"))

os.fit <- survfit(Surv(OS_time, OS) ~cluster.beta,
                  data = WCDT.clinical)
ggsurvplot(os.fit,size = 1,# change line size
           palette = c("#2E9FDF","#E7B800"), # custom color palette
           conf.int = FALSE, # Add confidence interval
           legend.title = "Beta",
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col ="strata",
           ylab="Overall survival",
           xlab="Time (days)",
           legend.labs=c("High","Low"),
           ggtheme = p_theme# Change ggplot2 theme
)
ggsave(file="F3_E.png",width = 4, height = 5)

os.fit <- survfit(Surv(OS_time, OS) ~cluster.tpm,
                  data = WCDT.clinical.ep)
ggsurvplot(os.fit,size = 1,# change line size
           palette = c("#2E9FDF","#E7B800"), # custom color palette
           conf.int = FALSE, # Add confidence interval
           legend.title = "Expression",
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col ="strata",
           ylab="Overall survival",
           xlab="Time (days)",
           legend.labs=c("High","Low"),
           ggtheme = p_theme# Change ggplot2 theme
)
ggsave(file="S3.A OS_WCDT_WGBS_GR_Gene_Expressionn.pdf",width = 4, height = 5)
