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

# positive_Exon
ExonYTHDF3_Pos <- readRDS("/home/kunqi/m6A reader/YTHDF3/Exon/postive.rds")

# testing_sample
index<- grep(1,ExonYTHDF3_Pos$GSE83438)
testT_testing <- as.character(DNAStringSet(Views(Hsapiens,ExonYTHDF3_Pos[index]+20)))
testT_testing_genome <- ExonYTHDF3_Pos[index]

# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,ExonYTHDF3_Pos[-index]+20)))
testT_genome <- ExonYTHDF3_Pos[-index]

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
ExonYTHDF3_StoreMatrix <- matrix(NA, nrow=9, ncol=2) 
rownames(ExonYTHDF3_StoreMatrix) <- c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")
colnames(ExonYTHDF3_StoreMatrix) <- c("crossValid","indepTest")

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
    storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/YTHDF3/Exon/negative",i,".rds",sep = "")  
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
    
    try(testing$pred_class <- predict(Cmotifvsnon_pred, testing))
    try(BIOmotifvsnon_testppred <- prediction(as.numeric(testing$pred_class),testing$class))
    try(BIOmotifvsnon_testpppred_auc <- performance(BIOmotifvsnon_testppred,"auc"))
    try(aucMatrix[i,]<-attr(BIOmotifvsnon_testpppred_auc,"y.values")[[1]][1])
    try(rocMatrix[i,] <- max(Cmotifvsnon_pred$results$ROC))
  } #end the 10 loops
  
  #store the matrix as names 
  setwd("~/testing")
  aucName <- paste("aucMatrix",method,"-tstallexon83438.rds",sep = "")  
  saveRDS(aucMatrix,aucName)
  
  #caculate the average AUC after 10 times 
  avgAUC <- apply(aucMatrix,MARGIN = 2, FUN = mean,na.rm=T)
  avgAUC
  
  #caculate the max ROC of the prediction model
  avgROC <- apply(rocMatrix,MARGIN = 2, FUN = mean,na.rm=T)
  avgROC
  
  # generate the results
  results <- c(avgROC,avgAUC)
  
  ExonYTHDF3_StoreMatrix[method,] <- results
}

setwd("/home/yuxuan.wu/m6A reader/storeMatrix")
saveRDS(ExonYTHDF3_StoreMatrix,"storeMatrix_ExonYTHDF3_tstall83438.rds")
