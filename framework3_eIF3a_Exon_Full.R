library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(ROCR)
library(pROC)
set.seed(520)

# positive_Full
FulleIF3a_Pos <- readRDS("/home/kunqi/m6A reader/eIF3a/Full/postive.rds")
readRDS("/home/kunqi/m6A reader/eIF3a/Full/negative1.rds")
index_train <- grep(1,FulleIF3a_Pos$GSE73405)
FulleIF3a_Pos_indep <- FulleIF3a_Pos[index_train]

# testing_sample
index_testing <- grep(1,FulleIF3a_Pos$GSE65004)
testT_testing <- as.character(DNAStringSet(Views(Hsapiens,FulleIF3a_Pos[index_testing]+20)))
testT_testing_genome <- FulleIF3a_Pos[index_testing]+20

# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,FulleIF3a_Pos_indep+20)))
testT_genome <- FulleIF3a_Pos_indep

library(caret)

# adding the sources
source("/home/kunqi/m7G/method/class1.R")
source("/home/kunqi/m7G/method/class2.R")
source("/home/kunqi/m7G/method/class3.R")
source("/home/kunqi/m7G/method/class4.R")
source("/home/kunqi/m7G/method/class5.R")
source("/home/kunqi/m7G/method/class6.R")

storeNegMatrix <- matrix(NA, nrow=10, ncol = 1)
aucMatrix <- matrix(NA, nrow = 10, ncol = 1)
rocMatrix <- matrix(NA, nrow = 10, ncol = 1)

# create a matrix to store the value
# rename the matrix's column and row 
FulleIF3a_StoreMatrix <- matrix(NA, nrow=6, ncol=2) 
rownames(FulleIF3a_StoreMatrix) <- c("CONPOSITION","Chemiprop+SeqFeature+one-hot","EIIP","AutoCo+PseKNC","AutoCo","PSNP")
colnames(FulleIF3a_StoreMatrix) <- c("crossValid","indepTest")

num=0
# for (method in c("CONPOSITION","Chemiprop+SeqFeature+one-hot","EIIP","AutoCo+PseKNC","AutoCo","PSNP")){
# num <- num+1
method <- "onehot"
  for (i in 1:10){
    # get the negative sample 
    storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/eIF3a/Exon/negative",i,".rds",sep = "")  
    testN_read <- readRDS(storeNegMatrix[i,])
    testN <- sample(testN_read,size =length(index_train) )
    testN_genome <- testN-2
    testN <- as.character(DNAStringSet(Views(Hsapiens,testN+18)))
    testAll <-c(testN,testT) 
    
    # analysis the test sample
    testN_testing <- testN[index_testing]
    testAll_testing <- c(testN_testing,testT_testing)
    
    if (method=="CONPOSITION"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/CONPOSITION.R")
    }
    
    if (method =="Chemiprop+SeqFeature+one-hot"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/Chemiprop+SeqFeature+one-hot.R")
    }
    
    if (method == "EIIP"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/EIIP.R")
    }
    
    if (method =="AutoCo"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/AutoCo.R")
    }
    
    if (method == "AutoCo+PseKNC"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/AutoCo+PseKNC.R")
    }
    
    if (method=="PSNP"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/PSNP.R")
    }
    
    if (method=="onehot"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/onehot.R")
    }
    
    #convert it to the dataframe format
    dataframeR<-as.data.frame(matrixResults)
    testing <- as.data.frame(testingResults)
    
    # label the data 
    label_train <- c(rep("Neg",(length(testN))), rep("Pos",(length(testT))))
    label_train <- factor(label_train , labels = c("Neg","Pos") )
    
    label_testing <- c(rep("Neg",(length(testN_testing))), rep("Pos",(length(testT_testing))))
    label_testing <- factor(label_testing , labels = c("Neg","Pos") )
    
    # add the label to the dataframe
    dataframeR$class<-label_train
    testing$class <- label_testing
    # import the independent testing sample
    
    training <- dataframeR
    
    fitControl <- trainControl(method = "cv",
                               number = 5,
                               classProbs = TRUE,
                               summaryFunction = twoClassSummary)
    #svm prediction
    Cmotifvsnon_pred <- train(class ~ ., data = training, 
                              method = 'svmRadial', 
                              preProc = c("center", "scale"),
                              trControl = fitControl,
                              metric = "ROC")
    
    testing$pred_class <- predict(Cmotifvsnon_pred, testing)
    BIOmotifvsnon_testppred <- prediction(as.numeric(testing$pred_class),testing$class)
    BIOmotifvsnon_testpppred_auc <- performance(BIOmotifvsnon_testppred,"auc")
    aucMatrix[i,]<-attr(BIOmotifvsnon_testpppred_auc,"y.values")[[1]][1]
    rocMatrix[i,] <- max(Cmotifvsnon_pred$results$ROC)
  }             #end the first for loop (10 times)
  
  #caculate the average AUC after 10 times 
  avgAUC <- apply(aucMatrix,MARGIN = 2, FUN = mean)
  avgAUC
  
  #caculate the max ROC of the prediction model
  avgROC <- apply(rocMatrix,MARGIN = 2, FUN = mean)
  avgROC
  
  results <- c(avgROC,avgAUC)
  
  FulleIF3a_StoreMatrix[num,] <- results
# }

# setwd("/home/yuxuan.wu/m6A reader/storeMatrix")
# saveRDS(FulleIF3a_StoreMatrix,"storeMatrix_FulleIF3a.rds")

# saveRDS(FulleIF3a_StoreMatrix,"storeMatrix.rds")
# FulleIF3a_StoreMatrix<- readRDS("storeMatrix.rds")
# readRDS("storeMatrix_ExoneIF3a.rds")
