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
top_selected_final <- readRDS("feature_eIF3a_Exon.rds") # you have to save it !
ROC_max <- readRDS("ROC_max_eIF3a_Exon.rds") 

# positive_Exon
ExoneIF3a_Pos <- readRDS("/home/kunqi/m6A reader/eIF3a/Exon/postive.rds")
readRDS("/home/kunqi/m6A reader/eIF3a/Exon/negative1.rds")
index_train <- grep(1,ExoneIF3a_Pos$GSE73405)
ExoneIF3a_Pos_indep <- ExoneIF3a_Pos[index_train]

# testing_sample
index_testing <- grep(1,ExoneIF3a_Pos$GSE65004)
testT_testing <- as.character(DNAStringSet(Views(Hsapiens,ExoneIF3a_Pos[index_testing]+20)))
testT_testing_genome <- ExoneIF3a_Pos[index_testing]

# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,ExoneIF3a_Pos_indep+20)))
testT_genome <- ExoneIF3a_Pos_indep

# genome method
genome_method <- c("genome","genome_chemiPro")
format=""

# create the matrix to store the auc and roc value 
aucMatrix<- matrix (NA,ncol = 1,nrow = 10)
storeNegMatrix <- matrix(NA, nrow=10, ncol = 1)
ExoneIF3a_StoreMatrix <- matrix(NA, nrow=9, ncol=2) 
rownames(ExoneIF3a_StoreMatrix) <- c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")
colnames(ExoneIF3a_StoreMatrix) <- c("crossValid","indepTest")
ExoneIF3a_StoreMatrix <-  readRDS("/home/yuxuan.wu/m6A reader/storeMatrix/eIF3a/storeMatrix_ExoneIF3a.rds")

# begin the loop
for (method in genome_method){
    i=1
    storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/eIF3a/Exon/negative",i,".rds",sep = "")  
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
      
      
      try(testing$pred_class2 <- predict(reModel2, testing[,c(top_selected_final,(ncol(training)))], type = "prob"))
      try(pred_result2 <- prediction(testing$pred_class2,testing$class))
      try(pred_result_auc2 <- performance(pred_result2,"auc"))
      try(aucMatrix[i,] <-attr(pred_result_auc2,"y.values")[[1]][1])
    }
    
    if (method=="genome_chemiPro"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/Chemiprop.R")
      training <- as.data.frame(matrixResults)
      testing <- as.data.frame(testingResults)
      
      chemi_index <- 59:ncol(matrixResults)
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
      
      
      try(testing$pred_class2 <- predict(reModel2, testing[,c(top_selected_final,chemi_index,(ncol(testing)))],type = "prob"))
      try(pred_result2 <- prediction(testing$pred_class2,testing$class))
      try(pred_result_auc2 <- performance(pred_result2,"auc"))
      try(aucMatrix[i,] <-attr(pred_result_auc2,"y.values")[[1]][1])
      saveRDS(reModel2,"/home/yuxuan.wu/m6A reader/model/eIF3a_Exon.rds")
    }
    
  
  avgAUC <- apply(aucMatrix,MARGIN = 2, FUN = mean, na.rm=TRUE)
  indepTest <- avgAUC
  indepTest
  
  #caculate the max ROC of the prediction model
  maxROC <- max(ROC_max)
  crossValid <- maxROC
  crossValid
  
  results <- c(crossValid,indepTest)
  ExoneIF3a_StoreMatrix[method,] <- results
}

# setwd("/home/yuxuan.wu/m6A reader/storeMatrix")
# saveRDS(ExoneIF3a_StoreMatrix,"Final_ExoneIF3a_StoreMatrix.rds")
