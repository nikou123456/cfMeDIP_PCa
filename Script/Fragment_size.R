#===============================================================================
#Fragment size based on cfMeDIP-seq of prostate cancer
#                           
#                            
#
#==================================================================================
#Fragment size distribution
rm(list=ls())
options(stringsAsFactors = FALSE)
library(dplyr)
library(scatterplot3d)
library(tidyr)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(extrafont)
source("./Script/function.R")
out.path <- "./Figure"
if(!dir.exists(out.path)){
  dir.create(out.path)
}


#===============================================================================
#' @title Loading fragment data for VPC cohort
col.data <- readRDS("./Data/col.datale_information.Rdata")
inster.v <- readRDS("./Data/fragment_vpc.Rdata")
index <- which(col.data$group.cp=="VPC")
vpc.clinical <- col.data[index,]
q0.25 <- as.numeric(quantile(vpc.clinical$ctDNA,0.25)) 
q0.5<- as.numeric(quantile(vpc.clinical$ctDNA,0.5)) 
q0.75<- as.numeric(quantile(vpc.clinical$ctDNA,0.75)) 
q1<- as.numeric(quantile(vpc.clinical$ctDNA,1)) 
vpc.clinical$ctDNA_type <- "QT0-25"
vpc.clinical$ctDNA_type[which(vpc.clinical$ctDNA>q0.25 & vpc.clinical$ctDNA<=q0.5)]<-"QT25-50"
vpc.clinical$ctDNA_type[which(vpc.clinical$ctDNA>q0.5 & vpc.clinical$ctDNA<=q0.75)]<-"QT50-75"
vpc.clinical$ctDNA_type[which(vpc.clinical$ctDNA>q0.75 & vpc.clinical$ctDNA<=q1)]<-"QT75-100"
index <- match(inster.v$col.dataleID,vpc.clinical$col.dataleID)
inster.v$ctDNA_type <- vpc.clinical$ctDNA_type[index]
identical(inster.v$col.dataleID,vpc.clinical$col.dataleID[index])


#===============================================================================
#' @title Plot Fig SF1
gd <- inster.v %>%
  group_by(ctDNA_type,insert_size) %>%
  summarise(proportation = mean(proportation))
col.datale.names <- unique(gd$ctDNA_type)


ggplot(data = dplyr::filter(.data =inster.v, 
                            ctDNA_type %in% c( "QT0-25",  "QT25-50",  "QT50-75", "QT75-100" )),
       aes(x=insert_size, y=proportation,color=ctDNA_type ))+
  labs(title=NULL,x="Insert size (bp)", y = "% of total reads",color="Type")+
  guides(color=guide_legend(title=NULL))+xlim(100,225)+p_theme+
  theme(legend.position=c(0.8,0.7),legend.background = element_blank())+
  stat_smooth(method="gam",formula=y~s(x,k=20),size=1,se=TRUE)+
  scale_color_manual(values = rev(RColorBrewer::brewer.pal(5,'Purples')))
ggsave(file="SF1_J.png",width=4,height=4)


#===============================================================================
#' @title Plotting Figure 1G
fg.data <- readRDS("./Data/fragment_ratio.Rdata")
my_comparisons <- list( 
  c("CPC","Barrier"),
  c("CPC","VPC"),
  c("CPC","WCDT"))
library(RColorBrewer)
fg.data$group.cp <- factor(fg.data$group.cp ,levels=c("CPC","Barrier","VPC","WCDT"))
colors <- brewer.pal(8,"Dark2")
ggplot(fg.data,
       aes(x=group.cp, y=proportation,color=group.cp))+
  geom_boxplot( outlier.shape = NA)+geom_jitter(position=position_jitter(0.2))+
  scale_color_manual(values=colors[c(1:3,7)])+ylim(0,0.75)+
  labs(x=NULL,y="Longer fragment ratio")+
  p_theme+ guides(color=guide_legend(title=NULL))+theme(legend.position  ="none")+
  stat_compare_means(comparisons = my_comparisons,#label = "p.signif",
                     method = "wilcox.test",hide.ns=TRUE)#+ # Add pairwise comparisons 
ggsave(file="F1_G.png",width=4,height=4)


#===============================================================================
#' @title Plotting Figure 1H
fg.data.v <- subset(fg.data,group.cp %in% c("VPC"))
index <- match(fg.data.v$filename,vpc.clinical$simpleName)
fg.data.v$ctDNA <- vpc.clinical$ctDNA[index]
fg.data.v$ctDNA_type <- vpc.clinical$ctDNA_type[index]
fg.data.v$ctDNA_type <- factor(fg.data.v$ctDNA_type,
                                  levels=c("QT0-25","QT25-50","QT50-75","QT75-100"))
my_comparisons <- list( c("QT0-25","QT75-100"),
                        c("QT0-25","QT50-75"),
                        c("QT0-25","QT25-50"),
                        c("QT25-50","QT75-100"))

ggplot(fg.data.v , 
       aes(x=ctDNA,y=proportation,color=group.cp))+geom_point()+#theme_classic()+
  p_theme+scale_color_manual(values=colors[c(3)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="% ctDNA",y="Longer fragment ratio")+
  theme( strip.text.x = element_text(size=12),
         legend.background = element_blank(),
         legend.position = c(0.9,0.9))+guides(color=FALSE)+ylim(0,0.6)
ggsave(file="F1_H.pdf",width=4,height=4)

#==================================================================================
#' @title  survival analysis for Fig I and Fig J
#==================================================================================
vpc.col <- readRDS("./Data/VPC_clinical.Rdata")
survival.data.sub <- subset(vpc.col,duplicated ==0)
survival.data.sub$OS.time <- as.numeric(survival.data.sub$`Published_Days to death or last followup`)
survival.data.sub$OS <- survival.data.sub$Published_Censored...14
survival.data.sub$OS[which(survival.data.sub$OS=="Yes")] <- 1
survival.data.sub$OS[which(survival.data.sub$OS=="No")] <- 2
survival.data.sub$OS <- as.numeric(survival.data.sub$OS)
survival.data.sub$RFS.time <- as.numeric(survival.data.sub$`Published_Days to progression or last followup`)
survival.data.sub$RFS <- survival.data.sub$Published_Censored...12
survival.data.sub$RFS[which(survival.data.sub$RFS=="Yes")] <- 1
survival.data.sub$RFS[which(survival.data.sub$RFS=="No")] <- 2
survival.data.sub$RFS <- as.numeric(survival.data.sub$RFS)

index <- match(survival.data.sub$simpleName,fg.data.v$filename)
survival.data.sub$fragment <- fg.data.v$proportation[index]
survival.data.sub$fra.cluster <- NA
survival.data.sub$fra.cluster[which(survival.data.sub$fragment>=quantile(survival.data.sub$fragment,0.5))] <- "Longer"
survival.data.sub$fra.cluster[which(survival.data.sub$fragment<quantile(survival.data.sub$fragment,0.5))] <- "Shorter"
survival.data.sub$fra.cluster <- factor(survival.data.sub$fra.cluster,levels=c("Shorter","Longer"))

library(survival)
library(survminer)
OS.fit <- survfit(Surv(OS.time, OS) ~ fra.cluster, 
                  data = survival.data.sub)
ggsurvplot(OS.fit,size = 1,# change line size
           palette = c("#2E9FDF", "#E7B800"), # custom color palette
           conf.int = FALSE, # Add confidence interval
           legend.title = "Insert",
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col ="strata", # Risk table color by groups
           ylab="Overall survival",
           xlab="Time (days)",
           legend.labs=c("Shorter","Longer"),
           ggtheme = p_theme# Change ggplot2 theme
)

ggsave(file="F1_I.png",width=4,height=5)


RFS.fit <- survfit(Surv(RFS.time, RFS) ~ fra.cluster, 
                   data = survival.data.sub)
ggsurvplot(RFS.fit,size = 1,# change line size
           palette = c("#2E9FDF", "#E7B800"), # custom color palette
           conf.int = FALSE, # Add confidence interval
           legend.title = "Insert",
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col ="strata", # Risk table color by groups
           ylab="Progression-free survival",
           xlab="Time (days)",
           legend.labs=c("Shorter","Longer"),
           ggtheme = p_theme# Change ggplot2 theme
)
ggsave(file="F1_J.png",width=4,height=5)


#==================================================================================
#' @title  Fragment ratio Standard deviation - F1 K
#==================================================================================
sd.plot <- readRDS("./Data/fragment_sd.Rdata")
sd.plot <- sd.plot[!sd.plot$group%in%c('Benign', 'HC', 'OCIR'), ];
sd.plot$group <- factor(sd.plot$group, levels = c('CPC', 'Barrier', 'VPC', 'WCDT'));
sd.plot <- unique(sd.plot);
ggplot(sd.plot,
       aes(x=group, y=mean,color=group))+
  geom_boxplot( outlier.shape = NA)+geom_jitter(position=position_jitter(0.2))+
  scale_color_manual(values=colors[c(1:3,7)])+ylim(0,0.06)+
  labs(x=NULL,y="Ratio standard deviation")+
  p_theme+ guides(color=guide_legend(title=NULL))+theme(legend.position  ="none")+
  stat_compare_means(comparisons = my_comparisons,#label = "p.signif",
                     method = "wilcox.test",hide.ns=TRUE)#+ # Add pairwise comparisons p-value
ggsave(file="F1_K.png",width=4,height=4)


