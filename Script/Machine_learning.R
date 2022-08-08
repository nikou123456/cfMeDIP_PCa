#==========================================================
# 
#Goal:  Building model by discovery CPC and VPC data
#       Selction feacture: Top150 HyerDMR and HypoDMR
#       Remaing sample as testing dataset
#=============================
rm(list=ls())
options(stringsAsFactors = FALSE)
library(DESeq2)
library(RWeka) 
library(pROC)
library(ROCR)
library(randomForest)
library(dplyr)
library(extrafont)
set.seed(123)
source("./Script/function.R")
out.path <- "./Figure"
if(!dir.exists(out.path)){
  dir.create(out.path)
}


#===============================================================================
#' @title M Mechine learing
#' @param pe.vpc DESeq2 result of  21VPC/OICR VPC versus 21CPCGENE/Localized
#' @param pe_hyper  RPKM table based on 26630 DMRs

#==================Loading  normalied count table===============================
norm.count.new <- readRDS("norm.count")
col.data <- readRDS("col_data.Rdata")
log.norm.counts <- log2(norm.count.new +1)
rm(norm.count.new)
if(!identical(col.data$simpleName,colnames(log.norm.counts))){
  stop("sample order error")
}

col.data$cancer_type <- "Localized" 
col.data$cancer_type[grep("Barrier|VPC|WCDT",col.data$group)] <-"Mets"
col.data$cancer_type <-as.character(col.data$cancer_type)



log.norm.counts.raw <- log.norm.counts
log.norm.counts.raw <- as.data.frame(log.norm.counts.raw)

colnames(log.norm.counts.raw) <- gsub("AIX_","AIX-",colnames(log.norm.counts.raw))


#################Loading RandomForest###########################################
#
#               Loading DMR result
#
################################################################################
input.path <- "/path/DMR_traning"
col.data.raw <- readRDS("col_data.Rdata")


setwd(input.path)
for (i in 1:50) {
  #i<-2
  print("#############Runing Time:")
  print(i)
  DMR.file <- paste0("Combined_Unique_Lolicazed_vs_VPC_with_age_control_time_",i,"_result_deseq.Rdata")
  DMR.result  <- readRDS(DMR.file)
  DMR.result$status <- gsub("\\s+regulated","",DMR.result$status)
  DMR.result <- DMR.result[-grep("Not Significant",DMR.result$status),]
  print(table(DMR.result$status))
  common.dmr <- intersect(rownames(DMR.result),rownames(log.norm.counts.raw))
  DMR.result <-DMR.result[common.dmr,]
  print(table(DMR.result$status))
  DMR.result$rank_score <- as.numeric(-log10(DMR.result$padj)*DMR.result$log2FoldChange)
  DMR.result  <-  DMR.result[order( DMR.result$rank_score,decreasing = TRUE),]
  DMR.result.top <- DMR.result)
  log.norm.counts <- log.norm.counts.raw[rownames(DMR.result.top),]
  rownames(log.norm.counts) <-paste0(DMR.result.top$status,"_",rownames(DMR.result.top))
  
  train_path <- paste0("TrainID_Combined_Unique_Lolicazed_vs_VPC_with_age_control_time_",i,".Rdata")
  train_id <- readRDS(train_path)
  train_id <- gsub("Local_|Mets_","",train_id)
  train_id <- gsub("AIX_","AIX-",train_id)
  index <- which(colnames(log.norm.counts) %in% train_id)
  data.rf <- log.norm.counts[,index]
  data.vd.raw <- log.norm.counts[,-index]
  
  rf.vpc<- colnames(data.rf)[grep("^V|^\\d+",colnames(data.rf))]
  vd.vpc <- colnames(data.vd.raw)[grep("^V|^\\d+$",colnames(data.vd.raw ))]
  
  rf.patient.vpc <- col.data.raw$patient[match(rf.vpc,col.data.raw$simpleName)]
  #col.data.raw$simpleName[match(colnames(data.rf),col.data.raw$simpleName)]
  vd.patient.vpc <- col.data.raw$patient[match(vd.vpc,col.data.raw$simpleName)]
  #identical(col.data.raw$simpleName[match(vd.vpc,col.data.raw$simpleName)],vd.vpc)
  
  overlap.pateint.id <- which(vd.patient.vpc %in% rf.patient.vpc)
  #intersect(vd.patient.vpc[overlap.pateint.id],rf.patient.vpc)
  overlap.vpc.id <- vd.vpc[overlap.pateint.id]  #need to delete
  data.vd <- data.vd.raw[,-which(colnames(data.vd.raw) %in% overlap.vpc.id)]
  
  

  index1 <- match(as.character(colnames(data.rf))   ,col.data$simpleName)
  index2 <- match(as.character(colnames(data.vd))   ,col.data$simpleName)
  train <- data.frame(t(data.rf),class= col.data$cancer_type[index1] )
  test <- data.frame(t(data.vd ),class= col.data$cancer_type[index2] )
  train$class <- as.factor(train$class)
  test$class <- as.factor(test$class) 
  setwd("/path/Output")
  #infor.information <- InfoGainAttributeEval(class ~.,data=train)
  
  infor.information <- readRDS(paste0("information_gain_list_time_",i,".Rdata"))
  index.hyper <- names(infor.information)[grep("^hyper",names(infor.information))][1:150]
  index.hypo <- names(infor.information)[grep("^hypo",names(infor.information))][1:150]
  trainSet <- train[,c(index.hyper,index.hypo,"class")]
  testSet <- test[,c(index.hyper,index.hypo,"class")]
  set.seed(123)

  
  setwd("/path/Output")
  # model <- randomForest(class ~ ., data = trainSet,
  #                       proximity = TRUE,
  #                       type = "classification",
  #                       importance = TRUE)
  model <- readRDS(paste0("RF_model_time_",i,".Rdata"))
  print("Have done construction of RF model")
  pred <- predict(model, testSet[, -ncol(testSet)])
  proba <- as.data.frame(predict(model, testSet[, -ncol(testSet)], type = "prob"))
  
  ROC <- plotROC(proba, testSet$class)
  result <- as.matrix(table(pred = pred, true = testSet$class))
  accuracy <- sum(diag(result))/sum(result)
  tp <- diag(result)
  fp <- rowSums(result) - diag(result)
  fn <- colSums(result) - diag(result)
  tn <- sum(result) - tp - fn - fp
  sp <- tn/(tn + fp)
  precision <- diag(result)/rowSums(result)
  recall <- (diag(result)/colSums(result))
  f1 <- 2 * precision * recall/(precision + recall)
  temp <- rbind(data.frame(precision, sp, recall, f1))
  data.temp <- temp
  data.temp$calss <- rownames(data.temp)
  data.temp <- reshape2::melt(data.temp)
  data.temp$variable<- factor(data.temp$variable,
                              levels = c("precision" , "sp", "recall", "f1"),
                              labels =c("Precision" , "Specificity", "Sensitivity" ,"F1 score") )
  
  data.temp$calss <- factor(data.temp$calss,levels = c("Localized","Mets"),
                            labels=c("Localized","Metastatic"))
  setwd("/path/Output")
  saveRDS(list(result=result,performance=data.temp,prediction=pred,
               accuracy=accuracy,ROC=ROC,probability=proba),
          file = paste0("Prediction_time_Unique_Patient_with_Timepoint",i,".Rdata"))
  setwd(input.path)
}



#================================================================================
#
#' @title Ploting Probability 
setwd("/path/Output")
files <- list.files("./","Prediction_time_Unique_Patient_with_Timepoint\\S+")  
proba <- data.frame()
for(i in 1:length(files)){
  temp <- readRDS(files[i])
  prob.temp <- temp$probability
  prob.temp$time <- i
  prob.temp$sample <- rownames(prob.temp)
  proba <- rbind(proba,prob.temp)
}
table(proba$time)
table(proba$sample)


subset.justion <- subset(proba,type %in%  c("CPC","PARD","HC") )
table(subset.justion$type)
subset.justion$type <- factor(subset.justion$type ,levels=c("CPC","PARD","HC"),
                              labels=c("CPC","OHS","Healthy"))
unique(subset.justion$sample)

my_comparisons <- list( 
  c("CPC","Healthy"),
  c("OHS","Healthy"))
table(subset.justion$type)

library(ggpubr)
ggplot(subset.justion,aes(x=type,y=Localized,fill=type))+geom_boxplot(outlier.shape = NA)+p_theme+
  #geom_violin(trim=FALSE)+ geom_boxplot(width=0.1,fill="white")+
  geom_jitter(position=position_jitter(0.2),color="#999999")+
  labs(x=NULL,y="Probability of localized PCa")+guides(fill="none")+
  scale_fill_brewer(palette="Set3")+ylim(0,1.19)+
  stat_compare_means(comparisons = my_comparisons,#label = "p.signif",
                     method = "wilcox.test",hide.ns=TRUE)#+stat_n_text(size=4)
ggsave(file="Prediciton_Probability_HC.pdf",width=4,height=4)
