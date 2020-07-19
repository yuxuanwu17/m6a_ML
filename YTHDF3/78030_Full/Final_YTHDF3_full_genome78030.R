library(doParallel)
registerDoParallel(makePSOCKcluster(12))

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

# testing_sample
index<- grep(1,FullYTHDF3_Pos$GSE78030)
testT_testing <- as.character(DNAStringSet(Views(Hsapiens,FullYTHDF3_Pos[index]+20)))
testT_testing_genome <- FullYTHDF3_Pos[index]

# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,FullYTHDF3_Pos[-index]+20)))
testT_genome <- FullYTHDF3_Pos[-index]

# genome method
genome_method <- c("genome","genome_chemiPro")
format=""

# create the matrix to store the auc and roc value 
aucMatrix<- matrix (NA,ncol = 1,nrow = 10)
storeNegMatrix <- matrix(NA, nrow=10, ncol = 1)
FullYTHDF3_StoreMatrix <- readRDS("/home/yuxuan.wu/m6A reader/storeMatrix/storeMatrix_FullYTHDF3_tstall78030.rds")

# begin the loop
for (method in genome_method){
  for(i in 1:10){
    storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/YTHDF3/Full/negative",i,".rds",sep = "")  
    testN_read <- readRDS(storeNegMatrix[i,])
    
    # get the desired length of negative training sample
    testN <- sample(testN_read,size =length(testT) )
    
    #prepare the sample used in the genome encoding 
    testN_genome <- testN-2
    
    #prepare the sample used in other encoding method as a character
    testN <- as.character(DNAStringSet(Views(Hsapiens,testN+18)))
    
    
    testN_testing_genome<- sample(testN_read,size=length(testT_testing))-2
    testN_testing <- testN[index]
    
    
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
      
      
      try(testing$pred_class2 <- predict(reModel2, testing[,c(top_selected_final,(ncol(training)))]))
      try(pred_result2 <- prediction(as.numeric(testing$pred_class2),testing$class))
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
      
      
      try(testing$pred_class2 <- predict(reModel2, testing[,c(top_selected_final,chemi_index,(ncol(testing)))]))
      try(pred_result2 <- prediction(as.numeric(testing$pred_class2),testing$class))
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
saveRDS(FullYTHDF3_StoreMatrix,"Final_FullYTHDF3_StoreMatrix78030.rds")
