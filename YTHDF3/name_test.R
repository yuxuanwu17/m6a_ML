library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(ROCR)
library(pROC)
library(kernlab)
# library(doParallel)
# cl <- makePSOCKcluster(5)
# registerDoParallel(cl)

set.seed(5123)

# positive_Full
FulleIF3a_Pos <- readRDS("/home/kunqi/m6A reader/eIF3a/Full/postive.rds")
readRDS("/home/kunqi/m6A reader/eIF3a/Full/negative1.rds")
index_train <- grep(1,FulleIF3a_Pos$GSE73405)
FulleIF3a_Pos_indep <- FulleIF3a_Pos[index_train]

# testing_sample
index_testing <- grep(1,FulleIF3a_Pos$GSE65004)
testT_testing <- as.character(DNAStringSet(Views(Hsapiens,FulleIF3a_Pos[index_testing]+20)))
testT_testing_genome <- FulleIF3a_Pos[index_testing]

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
aucMatrix2 <- matrix(NA, nrow = 10, ncol = 1)
rocMatrix <- matrix(NA, nrow = 10, ncol = 1)

# create a matrix to store the value
# rename the matrix's column and row 
FulleIF3a_StoreMatrix <- matrix(NA, nrow=9, ncol=2) 
rownames(FulleIF3a_StoreMatrix) <- c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")
colnames(FulleIF3a_StoreMatrix) <- c("crossValid","indepTest")
encodingMethod <- c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")
ecdMethod_except_genome <- c("CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","onehot")

format=""
for (method in ecdMethod_except_genome){
  for (i in 1:2){
    # get the negative sample 
    storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/eIF3a/Full/negative",i,".rds",sep = "")  
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
      
      Cmotifvsnon_pred <- train(class ~ ., data = training, 
                                method = 'svmRadial', 
                                preProc = c("center", "scale"),
                                trControl = fitControl,
                                metric = "ROC")
      
      try(testing$pred_class <- predict(Cmotifvsnon_pred, testing))
      try(BIOmotifvsnon_testppred <- prediction(as.numeric(testing$pred_class),testing$class))
      try(BIOmotifvsnon_testpppred_auc <- performance(BIOmotifvsnon_testppred,"auc"))
      try(aucMatrix[i,]<-attr(BIOmotifvsnon_testpppred_auc,"y.values")[[1]][1])
      try(rocMatrix[i,] <- max(Cmotifvsnon_pred$results$ROC))
  } #end the 10 loops
  
  #store the matrix as names 
  setwd("~/testing")
  aucName <- paste("aucMatrix",method,".rds",sep = "")  
  saveRDS(aucMatrix,aucName)
  
  
  #caculate the average AUC after 10 times 
  avgAUC <- apply(aucMatrix,MARGIN = 2, FUN = mean,na.rm=T)
  avgAUC
  
  #caculate the max ROC of the prediction model
  avgROC <- apply(rocMatrix,MARGIN = 2, FUN = mean,na.rm=T)
  avgROC
  
  # generate the results
  results <- c(avgROC,avgAUC)
  # FulleIF3a_StoreMatrix <- rbind(results,FulleIF3a_StoreMatrix)
  FulleIF3a_StoreMatrix[method,] <- results
}

setwd("/home/yuxuan.wu/m6A reader/storeMatrix")
# saveRDS(FulleIF3a_StoreMatrix,"storeMatrix_FulleIF3a_09_23.rds")

# saveRDS(FulleIF3a_StoreMatrix,"storeMatrix.rds")
# FulleIF3a_StoreMatrix<- readRDS("storeMatrix.rds")
# FulleIF3a_StoreMatrix <- readRDS("storeMatrix_FulleIF3a_09_23.rds")
