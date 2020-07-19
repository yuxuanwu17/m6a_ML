# library(doParallel)
# registerDoParallel(makePSOCKcluster(13))

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
FullYTHDF3_p <- readRDS("/home/kunqi/m6A reader/YTHDF3/Full/postive.rds")
p_index <- grep(1,FullYTHDF3_p$GSE86214)
training_FullYTHDF3_p <- FullYTHDF3_p[-p_index]

read_neg <- readRDS(paste("/home/kunqi/m6A reader/YTHDF3/Full/negative",1,".rds",sep = ""))
training_FullYTHDF3_n <- read_neg[-p_index]

matureSE <- SummarizedExperiment()
mcols(training_FullYTHDF3_p) <- NULL
rowRanges(matureSE) <- c(training_FullYTHDF3_p,training_FullYTHDF3_n-2) # -2!
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
label_train <- c(rep("P",(length(training_FullYTHDF3_p))),rep("N",(length(training_FullYTHDF3_n))))
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
# saveRDS(VarMatrix, "varMatrix86214.rds")

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
saveRDS(top_selected_final, "feature_YTHDF3_Full.rds") # you have to save it !

# plot auc
VarM <- matrix(NA,nrow = (ncol(predf)-1),ncol = 2)
VarM[,1] <- ROC_max[,1]
VarM[,2] <- colnames(predf)[order(F_score[,1],decreasing = T)] 
colnames(VarM) <- c("AUROC","top_N")
VarM <- as.data.frame(VarM)
VarM$AUROC <- as.numeric(as.character(VarM$AUROC))
VarM$top_N <- c(1:56)
# dev.off()
VarMresult <- ggplot(VarM,aes(x=top_N,y=AUROC))+geom_point(size=1)
VarMresult

# save the plot
setwd("/home/yuxuan.wu/m6A reader/auroc")
ggsave("auroc_Full_YTHDF3.pdf",VarMresult)


#--------------------------here we begin our desired loop

# positive_Full
FullYTHDF3_Pos <- readRDS("/home/kunqi/m6A reader/YTHDF3/Full/postive.rds")

# testing_sample
index<- grep(1,FullYTHDF3_Pos$GSE86214)
testT_testing <- as.character(DNAStringSet(Views(Hsapiens,FullYTHDF3_Pos[index]+20)))
testT_testing_genome <- FullYTHDF3_Pos[index]

# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,FullYTHDF3_Pos[-index]+20)))
testT_genome <- FullYTHDF3_Pos[-index]

# genome method
genome_method <- c("genome_chemiPro","genome")
format=""

# create the matrix to store the auc and roc value 
aucMatrix<- matrix (NA,ncol = 1,nrow = 10)
storeNegMatrix <- matrix(NA, nrow=10, ncol = 1)
FullYTHDF3_StoreMatrix <- matrix(NA, nrow=9, ncol=2) 
rownames(FullYTHDF3_StoreMatrix) <- c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")
colnames(FullYTHDF3_StoreMatrix) <- c("crossValid","indepTest")

# begin the loop
for (method in genome_method){
  for(i in 1:10){
    storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/YTHDF3/Full/negative",i,".rds",sep = "")  
    testN_read <- readRDS(storeNegMatrix[i,])
    
    # get the desired length of negative training sample
    testN <- sample(testN_read,size =length(index_train) )
    
    #prepare the sample used in the genome encoding 
    testN_genome <- testN-2
    
    #prepare the sample used in other encoding method as a character
    testN <- as.character(DNAStringSet(Views(Hsapiens,testN+18)))
    
    
    testN_testing_genome<- sample(testN_read,size=length(index_testing))-2
    testN_testing <- testN[index_testing]
    
    
    testAll <-c(testN,testT) 
    
    # analysis the test sample
    
    testAll_testing <- c(testN_testing,testT_testing)
    
    format=""
    if (method=="genome"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/genome.R")
    }
    
    if (method=="genome_chemiPro"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/genome_chemiPro.R")
    }
    
    training <- as.data.frame(matrixResults)
    testing <- as.data.frame(testingResults)
    
    
    # add labels for training sample
    label_train <- c(rep("Neg",(length(testN_genome))), rep("Pos",(length(testT_genome))))
    label_train <- factor(label_train , labels = c("Neg","Pos") )
    
    label_testing <- c(rep("Neg",(length(testN_testing_genome))), rep("Pos",(length(testT_testing_genome))))
    label_testing <- factor(label_testing , labels = c("Neg","Pos") )
    
    training$class<-label_train
    
    testing$class <- label_testing
    
    fitControl <- trainControl(method = "cv",
                               number = 5,
                               classProbs = TRUE,
                               summaryFunction = twoClassSummary)
    
    reModel2 <- train(class ~ ., data = training[,c(top_selected_final,ncol(training))], 
                      method = 'svmRadial',  
                      preProc = c("center", "scale"),
                      trControl = fitControl,
                      metric = "ROC")
    
    
    try(testing$pred_class2 <- predict(reModel2, testing[,c(top_selected_final,(ncol(training)))]))
    try(pred_result2 <- prediction(as.numeric(testing$pred_class2),testing$class))
    try(pred_result_auc2 <- performance(pred_result2,"auc"))
    try(aucMatrix[i,] <-attr(pred_result_auc2,"y.values")[[1]][1])
  }
  avgAUC <- apply(aucMatrix,MARGIN = 2, FUN = mean, na.rm=TRUE)
  indepTest <- avgAUC
  indepTest
  
  #caculate the max ROC of the prediction model
  maxROC <- max(ROC_max)
  crossValid <- maxROC
  crossValid
  
  results <- c(crossValid,indepTest)
  FullYTHDF3_StoreMatrix[method,] <- results
}

setwd("/home/yuxuan.wu/m6A reader/storeMatrix")
saveRDS(FullYTHDF3_StoreMatrix,"FullYTHDF3_StoreMatrix_genome_86214.rds")
