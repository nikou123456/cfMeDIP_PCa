#==========================================================
#' @title Plot Pincipal component for CPC and Barrier VPC and WCDT
#'        based on normalzation data on CPC/barrier/vpc/wcdt using deseq2::combet
#' 
#==========================================================
rm(list=ls())
options(stringsAsFactors = FALSE)
library(ggplot2)
library(export)
library(scatterplot3d)
library(extrafont)
source("./Script/function.R")
out.path <- "./Figure"
if(!dir.exists(out.path)){
  dir.create(out.path)
}


col.data <- readRDS("./Data/col.datale_information.Rdata")

#===============================================================================
#****************Principal component analysis************************
#log.norm.counts <- readRDS("log.norm.counts.Rdata")
#pac.raw <-  prcomp(t(log.norm.counts))
#saveRDS(list(x=pac.raw$x,pac.raw$sdev),file="./Data/principal_components_analysis.Rdata")
pac.raw <- readRDS("./Data/principal_components_analysis.Rdata")
names(pac.raw) <- c("x","sdev")
percentVar <- pac.raw$sdev^2 / sum( pac.raw$sdev^2 )
explain.pac <- data.frame(Var=percentVar*100,
                          cumm=as.numeric(cumsum(percentVar)*100) ,
                          nameID=colnames(pac.raw$x))
                        
library(RColorBrewer)
colors <- brewer.pal(8,"Dark2")
explain.pac$nameID <- factor(explain.pac$nameID,levels=explain.pac$nameID)
ggplot(explain.pac[1:10,],aes(x=nameID,group=1)) + 
  geom_bar(aes(y=Var), stat = "identity",fill="#3B9AB2") + #
  scale_x_discrete(labels=  explain.pac$nameID[1:20])+
  geom_line(aes(y=cumm)) + geom_point(aes(y=cumm)) +
  scale_x_discrete(name = "Principal component", labels = 1:20) +
  scale_y_continuous(name = "Percentage of explained variance", 
                     sec.axis = sec_axis(~./1, name = "Cumulative explained variance")) +p_theme+
  guides(fill=FALSE)

ggsave(file="F1_Percentage_of_explained_variance.png",width=4,height=4)

#===============================================================================
#****************Plot 3D and 2D plot************************
pac.3d <- as.data.frame(pac.raw$x)
pac.3d$sample <- rownames(pac.3d)
pac.3d$ctDNA <- col.data$ctDNA
pac.3d$group <- col.data$group
pac.3d$sampleID <- col.data$sampleID
pac.3d$source <- col.data$source
pac.3d$total.reads <- col.data$total_reads_medips
pac.3d$label <- paste0(pac.3d$group,"",pac.3d$ctDNA)


library(RColorBrewer)
type <- c("CPC", "Barrier", "VPC" , "WCDT" )
colors <- brewer.pal(8,"Dark2")
colors <- colors[c(1:3,7)]
use.colors <- colors[match( pac.3d$group,type )]

pac.3d$group <- factor(pac.3d$group,levels=c("CPC","Barrier","VPC","WCDT"))

colors <- brewer.pal(8,"Dark2")
ggplot(pac.3d, aes(PC1, PC2, color=group)) +
  geom_point(size=2.2) + scale_color_manual(values=colors[c(1:3,7)])+
  labs(title=NULL,x=paste0("PC1: ",round(percentVar[1] * 100),"% variance"),
       y=paste0("PC2: ",round(percentVar[2] * 100),"% variance"))+
  p_theme+    guides(color=FALSE) +
  geom_text_repel(
    data = subset(pac.3d, ctDNA<7),
    aes(label=ctDNA),
    size=4,
    box.padding=unit(0.35,"lines"),
    point.padding=unit(0.3,"lines"))+
  geom_text_repel(
    data = subset(pac.3d, PC2>200),
    aes(label=sampleID),
    size=4,
    box.padding=unit(0.35,"lines"),
    point.padding=unit(0.3,"lines"))

ggsave(file="S1_H.png",width=4,height=4)


ggplot(pac.3d, aes(PC1, PC3, color=group)) +
  geom_point(size=2.2) + scale_color_manual(values=colors[c(1:3,7)])+
  labs(title=NULL,x=paste0("PC1: ",round(percentVar[1] * 100),"% variance"),
       y=paste0("PC3: ",round(percentVar[3] * 100),"% variance"))+
  p_theme+    guides(color=FALSE)+
  geom_text_repel(
    data = subset(pac.3d, ctDNA<4.4),
    aes(label=ctDNA),
    size=4,
    box.padding=unit(0.35,"lines"),
    point.padding=unit(0.3,"lines"))+
  geom_text_repel(
    data = subset(pac.3d, PC3>300),
    aes(label=sampleID),
    size=4,
    box.padding=unit(0.35,"lines"),
    point.padding=unit(0.3,"lines"))
ggsave(file="S1_HH.png",width=4,height=4)

#
#
ggplot(pac.3d, aes(PC2, PC3, color=group)) +
  geom_point(size=2.2) + scale_color_manual(values=colors[c(1:3,7)])+
  labs(title=NULL,x=paste0("PC2: ",round(percentVar[2] * 100),"% variance"),
       y=paste0("PC3: ",round(percentVar[3] * 100),"% variance"))+
  p_theme+    guides(color=FALSE) +
  geom_text_repel(
    data = subset(pac.3d, PC3>200|PC2>(300)),
    aes(label=sampleID),
    size=4,
    box.padding=unit(0.35,"lines"),
    point.padding=unit(0.3,"lines"))
ggsave(file="S1_I.png",width=4,height=4)


#-------------------Legend--------------------------------------
p<-ggplot(pac.3d, aes(PC1, PC2, color=group)) +
  geom_point(size=2.2) + scale_color_manual(values=colors[c(1:3,7)])+
  labs(title=NULL,x=paste0("PC1: ",round(percentVar[1] * 100),"% variance"),
       y=paste0("PC2: ",round(percentVar[2] * 100),"% variance"))+
  p_theme+theme(legend.position = "top")+guides(color=guide_legend(title=NULL))
leg <- get_legend(p)
# Convert to a ggplot and print
as_ggplot(leg)
ggsave(file="F1_PCA_legend.pdf",width=4,height=2)





#============================================================================
# plot ctDNA 
ggplot(subset(pac.3d,group %in% c("CPC","VPC")), aes(PC1, PC2, color=ctDNA)) +
  geom_point(size=2.2) + scale_colour_gradient(high = "#DA141E", low = "#FFEFA4", na.value = "gray")+
  # scale_color_gradient2(values=colorRampPalette( brewer.pal( 9 , "YlOrRd" )))+
  labs(title=NULL,x=paste0("PC1: ",round(percentVar[1] * 100),"% variance"),
       y=paste0("PC2: ",round(percentVar[2] * 100),"% variance"))+
  p_theme+guides(color=FALSE)# + #+theme(legend.position = "top")#+


ggplot(subset(pac.3d,group %in% c("CPC","VPC")), aes(PC1, PC3, color=ctDNA)) +
  geom_point(size=2.2) + scale_colour_gradient(high = "#DA141E", low = "#FFEFA4", na.value = "gray")+
  # scale_color_gradient2(values=colorRampPalette( brewer.pal( 9 , "YlOrRd" )))+
  labs(title=NULL,x=paste0("PC1: ",round(percentVar[1] * 100),"% variance"),
       y=paste0("PC3: ",round(percentVar[3] * 100),"% variance"))+
  p_theme+guides(color=FALSE)# + #+theme(legend.position = "top")#+

ggplot(subset(pac.3d,group %in% c("CPC","VPC")), aes(PC2, PC3, color=ctDNA)) +
  geom_point(size=2.2) + scale_colour_gradient(high = "#DA141E", low = "#FFEFA4", na.value = "gray")+
  # scale_color_gradient2(values=colorRampPalette( brewer.pal( 9 , "YlOrRd" )))+
  labs(title=NULL,x=paste0("PC2: ",round(percentVar[2] * 100),"% variance"),
       y=paste0("PC3: ",round(percentVar[3] * 100),"% variance"))+
  p_theme+guides(color=FALSE)# + #+theme(legend.position = "top")#+



p<-ggplot(subset(pac.3d,group %in% c("CPC","VPC")), aes(PC1, PC2, color=ctDNA)) +
  geom_point(size=2.2) + scale_colour_gradient(high = "#DA141E", low = "#FFEFA4", na.value = "gray")+
  # scale_color_gradient2(values=colorRampPalette( brewer.pal( 9 , "YlOrRd" )))+
  labs(title=NULL,x=paste0("PC1: ",round(percentVar[1] * 100),"% variance"),
       y=paste0("PC2: ",round(percentVar[2] * 100),"% variance"))+
  p_theme
leg <- get_legend(p)
# Convert to a ggplot and print
as_ggplot(leg)



intervals   = seq(0,0.9,1/10)
x.5      = findInterval(pac.3d$ctDNA/100, intervals)
gradient <- colorRampPalette( brewer.pal( 9 , "YlOrRd" ))
colours <- gradient(length(intervals))
x.6     = colours[x.5]
index <- which(is.na(x.6))
final.colors.raw <- x.6 
final.colors.raw[index] <- "gray"

index<- which(pac.3d$group %in% c("CPC","VPC"))
temp <- subset(pac.3d,group %in% c("CPC","VPC"))
scatterplot3d(x=temp$PC1,y=temp$PC2,z=temp$PC3,
              color=final.colors.raw[index] ,
              pch=16,
              xlab="PC1",
              ylab="PC2",
              zlab="PC3",grid=TRUE)  
ggsave(file="F1_PCA_3D_ctDNA.png",width=5,height=5)



library(scatterplot3d)

scatterplot3d(x=pac.3d$PC1,y=pac.3d$PC2,z=pac.3d$PC3,
              color=use.colors ,
              pch=16,
              xlab="PC1",
              ylab="PC2",
              zlab="PC3",grid=TRUE)  
ggsave(file="F1_F.png",width=5,height=5)
