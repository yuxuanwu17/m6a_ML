library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
set.seed(520)

# positive_Exon
ExonIF3a_Pos <- readRDS("/home/kunqi/m6A reader/eIF3a/Exon/postive.rds")
readRDS("/home/kunqi/m6A reader/eIF3a/Exon/negative1.rds")
index_train <- grep(1,ExonIF3a_Pos$GSE73405)
ExonIF3a_Pos_indep <- ExonIF3a_Pos[index_train]

# testing_sample
index_testing <- grep(1,ExonIF3a_Pos$GSE65004)
testT_testing <- as.character(DNAStringSet(Views(Hsapiens,ExonIF3a_Pos[index_testing]+20)))

# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,ExonIF3a_Pos_indep+20)))

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
ExonIF3a_StoreMatrix <- matrix(NA, nrow=6, ncol=2) 
rownames(ExonIF3a_StoreMatrix) <- c("CONPOSITION","Chemiprop+SeqFeature+one-hot","EIIP","AutoCo+PseKNC","AutoCo","PSNP")
colnames(ExonIF3a_StoreMatrix) <- c("crossValid","indepTest")

num=0
for (method in c("CONPOSITION","Chemiprop+SeqFeature+one-hot","EIIP","AutoCo+PseKNC","AutoCo","PSNP")){
num <- num+1
for (i in 1:10){
  # get the negative sample 
  storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/eIF3a/Exon/negative",i,".rds",sep = "")  
  testN <- readRDS(storeNegMatrix[i,])
  testN <- as.character(DNAStringSet(Views(Hsapiens,testN+18)))
  testAll <-c(testN,testT) 
  
  # analysis the test sample
  testN_testing <- testN[index_testing]
  testAll_testing <- c(testN_testing,testT_testing)
  
  if (method=="CONPOSITION"){
  matrixResults <- CONPOSITION(testAll,NI=3,NTYPE = "DNA",Freq = 2)
  testingResults <- CONPOSITION(testAll_testing,NI=3,NTYPE = "DNA", Freq = 2)
  }
  
  if (method =="Chemiprop+SeqFeature+one-hot"){
  matrixResults_CP <- ChemicalProperty(testAll)
  matrixResults_SQ <- ChemicalProperty(testAll)
  matrixResults_onehot <- ChemicalProperty(testAll)
  matrixResults <- cbind(matrixResults_CP,matrixResults_SQ,matrixResults_onehot)
  testingResults_CP <- ChemicalProperty(testAll_testing)
  testingResults_SQ<- ChemicalProperty(testAll_testing)
  testingResults_onehot<- ChemicalProperty(testAll_testing)
  testingResults <- cbind(testingResults_CP,testingResults_SQ,testingResults_onehot)
  }

  if (method == "EIIP"){
    matrixResults <- electronIonInteraction(testAll)
    testingResults <- electronIonInteraction(testAll_testing)
  }

  if (method =="AutoCo"){
    for(i in 1:ncol(RNA_SCORE)){
      if(i==1){
        matrixResults <- autoCovariance(testAll,convert=1,descriptor = colnames(RNA_SCORE)[i])
      }
      matrixResults <- cbind(matrixResults,autoCovariance(testAll,convert=1,descriptor = colnames(RNA_SCORE)[i]))
    }
    
    for(i in 1:ncol(RNA_SCORE)){
      if(i==1){
        testingResults <- autoCovariance(testAll_testing,convert=1,descriptor = colnames(RNA_SCORE)[i])
      }
      testingResults <- cbind(testingResults,autoCovariance(testAll_testing,convert=1,descriptor = colnames(RNA_SCORE)[i]))
    }
  }
  
  if (method == "AutoCo+PseKNC"){
    for(i in 1:3){
      if(i==1){
      matrixResults_ACC <- autoCovariance(testAll,convert=1,descriptor = colnames(RNA_SCORE)[10-i])
      }
      matrixResults_ACC <- cbind(matrixResults_ACC,autoCovariance(testAll,convert=1,descriptor = colnames(RNA_SCORE)[i]))
    }
    
    for(i in 1:3){
      if(i==1){
      testingResults_ACC <- autoCovariance(testAll_testing,convert=1,descriptor = colnames(RNA_SCORE)[10-i])
      }
      testingResults_ACC <- cbind(testingResults_ACC,autoCovariance(testAll_testing,convert=1,descriptor = colnames(RNA_SCORE)[i]))
    }
    
    matrixResults_pse <- PseKNC(testAll, NI=2,NTYPE = "DNA")
    testingResults_pse <- PseKNC(testAll_testing, NI = 2, NTYPE = "DNA")
    
    matrixResults <- cbind(matrixResults_ACC,matrixResults_pse)
    testingResults <- cbind(testingResults_ACC,testingResults_pse)
  }
  
  if (method=="PSNP"){
    matrixResults<- PSNP(testN,testT,NI=3,NTYPE = "DNA")
    testingResults <- PSNP(testN_testing,testT_testing,NI=3,NTYPE = "DNA")
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
  
  library(ROCR)
  library(pROC)
  
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

ExonIF3a_StoreMatrix[num,] <- results
}

# saveRDS(ExonIF3a_StoreMatrix,"storeMatrix_ExonIF3a.rds")

# saveRDS(ExonIF3a_StoreMatrix,"storeMatrix.rds")
# ExonIF3a_StoreMatrix<- readRDS("storeMatrix.rds")
