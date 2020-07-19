# library(doParallel)
# registerDoParallel(makePSOCKcluster(20))

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
top_selected_final <- readRDS("feature_YTHDF3_Exon.rds") # you have to save it !
ROC_max <- readRDS("ROC_max_YTHDF3_Exon.rds") 

ExonYTHDF3_Pos <- readRDS("/home/kunqi/m6A reader/YTHDF3/Exon/postive.rds")
# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,ExonYTHDF3_Pos+20)))
testT_genome <- ExonYTHDF3_Pos


# genome method
genome_method <- c("genome","genome_chemiPro")
format="mixup"

# create the matrix to store the auc and roc value 
aucMatrix<- matrix (NA,ncol = 1,nrow = 10)
storeNegMatrix <- matrix(NA, nrow=10, ncol = 1)
ExonYTHDF3_StoreMatrix <- readRDS("/home/yuxuan.wu/m6A reader/storeMatrix/ExonYTHDF3_StoreMatrix_mix.rds")

# begin the loop
for (method in genome_method){
  for(i in 1:10){
    # get the negative sample 
    storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/YTHDF3/Exon/negative",i,".rds",sep = "")  
    testN <- readRDS(storeNegMatrix[i,])
    
    #prepare the sample used in the genome encoding 
    testN_genome <- testN-2
    
    #prepare the sample used in other encoding method as a character
    testN <- as.character(DNAStringSet(Views(Hsapiens,testN+18)))
    # combine the two positive and negative sample together
    testAll <-c(testN,testT) 
    
    if (method=="genome"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/genome.R")
      genomeMatrix <- as.data.frame(matrixResults)
      # #convert it to the dataframe format
      dataframeR<-as.data.frame(matrixResults)
      
      # label the data 
      label_train <- c(rep("Neg",(length(testN))), rep("Pos",(length(testT))))
      label_train <- factor(label_train , labels = c("Neg","Pos") )
      
      #add the label to the dataframe
      dataframeR$class<-label_train
      
      #split the data
      intrain <- createDataPartition(y = dataframeR$class, p= 0.8, list = FALSE)
      training <- dataframeR[intrain,]
      testing <- dataframeR[-intrain,]
      
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
      chemi_index <- 58:ncol(matrixResults)
      
      dataframeR <- cbind(genomeMatrix,matrixResults)
      
      # label the data 
      label_train <- c(rep("Neg",(length(testN))), rep("Pos",(length(testT))))
      label_train <- factor(label_train , labels = c("Neg","Pos") )
      
      #add the label to the dataframe
      dataframeR$class<-label_train
      
      #split the data
      intrain <- createDataPartition(y = dataframeR$class, p= 0.8, list = FALSE)
      training <- dataframeR[intrain,]
      testing <- dataframeR[-intrain,]
      
      fitControl <- trainControl(method = "cv",
                                 number = 5,
                                 classProbs = TRUE,
                                 summaryFunction = twoClassSummary)
      
      reModel2 <- train(class ~ ., data = training[,c(top_selected_final,chemi_index,ncol(training))], 
                        method = 'svmRadial',  
                        preProc = c("center", "scale"),
                        trControl = fitControl,
                        metric = "ROC")
      
      try(testing$pred_class2 <- predict(reModel2, testing[,c(top_selected_final,chemi_index,(ncol(testing)))], type="prob"))
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
  ExonYTHDF3_StoreMatrix[method,] <- results
}

setwd("/home/yuxuan.wu/m6A reader/storeMatrix")
saveRDS(ExonYTHDF3_StoreMatrix,"Final_ExonYTHDF3_StoreMatrix_mixup.rds")
