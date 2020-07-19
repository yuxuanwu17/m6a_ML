library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(ROCR)
library(pROC)
library(kernlab)
library(caret)

# library(doParallel)
# cl <- makePSOCKcluster(10)
# registerDoParallel(cl)

set.seed(2123)

# positive_Full
FullYTHDF3_Pos <- readRDS("/home/kunqi/m6A reader/YTHDF3/Full/postive.rds")

# testing_sample
index<- grep(1,FullYTHDF3_Pos$GSE83438)
testT_testing <- as.character(DNAStringSet(Views(Hsapiens,FullYTHDF3_Pos[index]+20)))
testT_testing_genome <- FullYTHDF3_Pos[index]

# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,FullYTHDF3_Pos[-index]+20)))
testT_genome <- FullYTHDF3_Pos[-index]

# adding the sources
source("/home/kunqi/m7G/method/class1.R")
source("/home/kunqi/m7G/method/class2.R")
source("/home/kunqi/m7G/method/class3.R")
source("/home/kunqi/m7G/method/class4.R")
source("/home/kunqi/m7G/method/class5.R")
source("/home/kunqi/m7G/method/class6.R")

storeNegMatrix <- matrix(NA, nrow=10, ncol = 1)
aucMatrix <- matrix(NA, nrow = 10, ncol = 1)
aucMatrix2 <- matrix(NA, nrow = 10, ncol = 1)
rocMatrix <- matrix(NA, nrow = 10, ncol = 1)

# create a matrix to store the value
# rename the matrix's column and row 
FullYTHDF3_StoreMatrix <- matrix(NA, nrow=9, ncol=2) 
rownames(FullYTHDF3_StoreMatrix) <- c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")
colnames(FullYTHDF3_StoreMatrix) <- c("crossValid","indepTest")

#encoding method as a list
encodingMethod <- c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")
genomeMethod <- c("CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","onehot")
tst1 <- c("CONPOSITION","Chemiprop_SeqFeature")
tst2 <- c("EIIP","AutoCo_PseKNC")
tst3 <- c("AutoCo","PSNP","onehot")
tst_all <- c(tst1,tst2,tst3)

format=""
for (method in tst_all){
  for (i in 1:10){
    # get the negative sample 
    storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/YTHDF3/Full/negative",i,".rds",sep = "")  
    testN_read <- readRDS(storeNegMatrix[i,])
    
    #negative training sample
    testN <- testN_read[-index]
    testN_genome <- testN-2
    
    #prepare the sample used in other encoding method as a character
    testN <- as.character(DNAStringSet(Views(Hsapiens,testN+18)))
    
    testN_testing_genome<- testN_read[index]-2
    testN_testing <- as.character(DNAStringSet(Views(Hsapiens,testN_read[index]+18)))
    
    testAll <-c(testN,testT) 
    testAll_testing <- c(testN_testing,testT_testing)
    
    if (method=="genome_chemiPro"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/genome_chemiPro.R")
    }
    
    if (method=="CONPOSITION"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/CONPOSITION.R")
    }
    
    if (method =="Chemiprop_SeqFeature"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/Chemiprop_SeqFeature.R")
    }
    
    if (method == "EIIP"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/EIIP.R")
    }
    
    if (method =="AutoCo"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/AutoCo.R")
    }
    
    if (method == "AutoCo_PseKNC"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/AutoCo_PseKNC.R")
    }
    
    if (method=="PSNP"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/PSNP.R")
    }
    
    if (method=="genome"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/genome.R")
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
    testing$pred_class <- NA
    BIOmotifvsnon_testppred <- NA
    BIOmotifvsnon_testpppred_auc <- NA
    
    try(testing$pred_class <- predict(Cmotifvsnon_pred, testing,type="prob"))
    try(BIOmotifvsnon_testppred <- prediction(testing$pred_class$Pos,testing$class))
    try(BIOmotifvsnon_testpppred_auc <- performance(BIOmotifvsnon_testppred,"auc"))
    try(aucMatrix[i,]<-attr(BIOmotifvsnon_testpppred_auc,"y.values")[[1]][1])
    try(rocMatrix[i,] <- max(Cmotifvsnon_pred$results$ROC))
  } #end the 10 loops
  
  #store the matrix as names 
  setwd("~/testing")
  aucName <- paste("aucMatrix",method,"-tstallFull83438.rds",sep = "")  
  saveRDS(aucMatrix,aucName)
  
  #caculate the average AUC after 10 times 
  avgAUC <- apply(aucMatrix,MARGIN = 2, FUN = mean,na.rm=T)
  avgAUC
  
  #caculate the max ROC of the prediction model
  avgROC <- apply(rocMatrix,MARGIN = 2, FUN = mean,na.rm=T)
  avgROC
  
  # generate the results
  results <- c(avgROC,avgAUC)
  
  FullYTHDF3_StoreMatrix[method,] <- results
}

setwd("/home/yuxuan.wu/m6A reader/storeMatrix")
saveRDS(FullYTHDF3_StoreMatrix,"storeMatrix_FullYTHDF3_tstall83438.rds")

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

setwd("/home/yuxuan.wu/m6A reader/top_selected_final")
top_selected_final <- readRDS("feature_YTHDF3_Full.rds") # you have to save it !
ROC_max <- readRDS("ROC_max_YTHDF3_Full.rds") 

# positive_Full
FullYTHDF3_Pos <- readRDS("/home/kunqi/m6A reader/YTHDF3/Full/postive.rds")
readRDS("/home/kunqi/m6A reader/YTHDF3/Full/negative1.rds")
index_train <- grep(1,FullYTHDF3_Pos$GSE83438)
FullYTHDF3_Pos_indep <- FullYTHDF3_Pos[index_train]

# testing_sample
index_testing <- grep(1,FullYTHDF3_Pos$GSE83438)
testT_testing <- as.character(DNAStringSet(Views(Hsapiens,FullYTHDF3_Pos[index_testing]+20)))
testT_testing_genome <- FullYTHDF3_Pos[index_testing]

# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,FullYTHDF3_Pos_indep+20)))
testT_genome <- FullYTHDF3_Pos_indep

# genome method
genome_method <- c("genome","genome_chemiPro")
format=""

# create the matrix to store the auc and roc value 
aucMatrix<- matrix (NA,ncol = 1,nrow = 10)
storeNegMatrix <- matrix(NA, nrow=10, ncol = 1)
FullYTHDF3_StoreMatrix <- matrix(NA, nrow=9, ncol=2) 
rownames(FullYTHDF3_StoreMatrix) <- c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")
colnames(FullYTHDF3_StoreMatrix) <- c("crossValid","indepTest")
FullYTHDF3_StoreMatrix <-  readRDS("/home/yuxuan.wu/m6A reader/storeMatrix/storeMatrix_FullYTHDF3_tstall83438.rds")

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
    
    if (method=="genome"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/genome.R")
      training <- as.data.frame(matrixResults)
      testing <- as.data.frame(testingResults)
      genomeMatrix_train <- as.data.frame(matrixResults)
      genomeMatrix_test <- as.data.frame(testingResults)
      
      
      # label the data 
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
      
      
      try(testing$pred_class2 <- predict(reModel2, testing[,c(top_selected_final,(ncol(training)))],type="prob"))
      try(pred_result2 <- prediction(testing$pred_class2$Pos,testing$class))
      try(pred_result_auc2 <- performance(pred_result2,"auc"))
      try(aucMatrix[i,] <-attr(pred_result_auc2,"y.values")[[1]][1])
    }
    
    if (method=="genome_chemiPro"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/Chemiprop.R")
      training <- as.data.frame(matrixResults)
      testing <- as.data.frame(testingResults)
      
      chemi_index <- 57:ncol(matrixResults)
      training <- cbind(genomeMatrix_train,training)
      testing <- cbind(genomeMatrix_test,testing)
      
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
      
      reModel2 <- train(class ~ ., data = training[,c(top_selected_final,chemi_index,ncol(training))], 
                        method = 'svmRadial',  
                        preProc = c("center", "scale"),
                        trControl = fitControl,
                        metric = "ROC")
      
      
      try(testing$pred_class2 <- predict(reModel2, testing[,c(top_selected_final,chemi_index,(ncol(testing)))],type="prob"))
      try(pred_result2 <- prediction(testing$pred_class2$Pos,testing$class))
      try(pred_result_auc2 <- performance(pred_result2,"auc"))
      try(aucMatrix[i,] <-attr(pred_result_auc2,"y.values")[[1]][1])
    }
    
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
saveRDS(FullYTHDF3_StoreMatrix,"Final_FullYTHDF3_StoreMatrix83438.rds")
