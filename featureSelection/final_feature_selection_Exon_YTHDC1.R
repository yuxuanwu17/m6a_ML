# library(doParallel)
# registerDoParallel(makePSOCKcluster(7))

library(S4Vectors)
library(m6ALogisticModel)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fitCons.UCSC.hg19)
library(phastCons100way.UCSC.hg19)
library(caret)
library(ROCR)
library(pROC)

set.seed(42)

source("~/m6A reader/featureSelection/F_score_fn.R")
source("/home/kunqi/m7G/method/class1.R")
source("/home/kunqi/m7G/method/class2.R")
source("/home/kunqi/m7G/method/class3.R")
source("/home/kunqi/m7G/method/class4.R")
source("/home/kunqi/m7G/method/class5.R")
source("/home/kunqi/m7G/method/class6.R")


Additional_features_hg19 = list(
  HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
  miR_targeted_genes = miR_targeted_genes_grl,
  TargetScan = TargetScan_hg19_gr,
  Verified_miRtargets = verified_targets_gr,
  METTL3_TREW = METTL3_TREW,
  METTL14_TREW = METTL14_TREW,
  WTAP_TREW = WTAP_TREW,
  METTL16_CLIP = METTL16_CLIP,
  ALKBH5_PARCLIP = ALKBH5_PARCLIP,
  FTO_CLIP = FTO_CLIP,
  FTO_eCLIP = FTO_eCLIP
)

# get data
ExonYTHDC1_p <- readRDS("/home/kunqi/m6A reader/YTHDC1/Exon/postive.rds")
read_neg <- readRDS(paste("/home/kunqi/m6A reader/YTHDC1/Exon/negative",1,".rds",sep = ""))

matureSE <- SummarizedExperiment()
mcols(ExonYTHDC1_p) <- NULL
rowRanges(matureSE) <- c(ExonYTHDC1_p,read_neg-2) # -2!
mature_FE <- predictors_annot(se = matureSE,
                              txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                              bsgnm = Hsapiens,
                              fc = fitCons.UCSC.hg19,
                              motif = c(),
                              pc = phastCons100way.UCSC.hg19,
                              struct_hybridize = Struc_hg19,
                              feature_lst = Additional_features_hg19,
                              hk_genes_list = HK_hg19_eids,
                              genes_ambiguity_method = "average")
matrixResults <- mcols(mature_FE)
predf <- data.frame(matrixResults)

# normalization
predf <- (predf-min(predf))/(max(predf)-min(predf))

# F score
F_score <- F_score_fn(predf)

# add labels for training sample
label_train <- c(rep("P",(length(ExonYTHDC1_p))),rep("N",(length(read_neg))))
label_train <- factor(label_train , labels = c("N","P"))
predf$class<-label_train

# normalization for training only
imp_normal <- (F_score-min(F_score))/(max(F_score)-min(F_score))

# Varmatrix 
VarMatrix <- matrix(NA,nrow = length(rownames(imp_normal)),ncol = 2)
VarMatrix[,1]  <- imp_normal[,1][order(imp_normal[,1],decreasing = T)]
VarMatrix[,2]<- colnames(predf)[order(imp_normal[,1],decreasing = T)] 
colnames(VarMatrix) <- c("F_score","features")
VarMatrix <- as.data.frame(VarMatrix)
VarMatrix$F_score <- as.numeric(as.character(VarMatrix$F_score))
VarMatrix$features <- factor(VarMatrix$features,levels = VarMatrix$features)

# plot
Varplot <- ggplot(VarMatrix,aes(x=features,y=F_score))+
  geom_histogram(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90))
Varplot
ggsave("/home/yuxuan.wu/m6A reader/auroc/AUROC_nomal_f_score_YTHDC1_Exon.pdf",Varplot,width = 8.63, height = 4.7)


# svm model
fitControl <- trainControl(method = "cv",
                           number = 5,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary)

# create a matrix to find the optimal combination 
ROC_max <- matrix(NA, nrow = length(F_score[,1]), ncol = 1)
index <- createDataPartition(y = predf$class, p = .75, list = FALSE)
training <- predf[ index,]
testing  <- predf[-index,]

ROC_max <- matrix(NA, nrow = length(F_score[,1]), ncol = 1)
for(j in 1:length(F_score[,1])){
  top_selected <- which(colnames(predf)%in%VarMatrix[,2][1:j]==TRUE)
  reModel <- train(class ~ ., data = predf[,c(top_selected,ncol(predf))], 
                   method = 'svmRadial',  
                   preProc = c("center", "scale"),
                   trControl = fitControl,
                   metric = "ROC")
  ROC_max[j,]<-max(reModel$results$ROC)
}
top_selected_final <- which(colnames(predf)%in%VarMatrix[,2][1:which(ROC_max==max(ROC_max))]==TRUE)
setwd("/home/yuxuan.wu/m6A reader/top_selected_final")
saveRDS(top_selected_final, "feature_YTHDC1_Exon.rds") # you have to save it !
saveRDS(ROC_max,"ROC_max_YTHDC1_Exon.rds")

# plot auc
VarM <- matrix(NA,nrow = (ncol(predf)-1),ncol = 2)
VarM[,1] <- ROC_max[,1]
VarM[,2] <- colnames(predf)[order(F_score[,1],decreasing = T)] 
colnames(VarM) <- c("AUROC","top_N")
VarM <- as.data.frame(VarM)
VarM$AUROC <- as.numeric(as.character(VarM$AUROC))
VarM$top_N <- c(1:58)
# dev.off()
VarMresult <- ggplot(VarM,aes(x=top_N,y=AUROC))+geom_point(size=1)
VarMresult <- VarMresult+ylab("AUC")+xlab("Performance of site prediction\n(TopN Features Used)")
VarMresult

# save the plot
setwd("/home/yuxuan.wu/m6A reader/auroc")
ggsave("auroc_Exon_YTHDC1_mix.pdf",VarMresult,width = 8.63, height = 4.7)

