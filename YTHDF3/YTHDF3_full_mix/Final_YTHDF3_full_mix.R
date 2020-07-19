library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(ROCR)
library(pROC)

# library(doParallel)
# cl <- makePSOCKcluster(5)
# registerDoParallel(cl)

set.seed(520)

# positive_Full
FullYTHDF3_Pos <- readRDS("/home/kunqi/m6A reader/YTHDF3/Full/postive.rds")
# show the first negative sample 
# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,FullYTHDF3_Pos+20)))
testT_genome <- FullYTHDF3_Pos

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
aucMatrix2 <- matrix(NA, nrow = 10, ncol = 1)
# create a matrix to store the value
# rename the matrix's column and row 
FullYTHDF3_StoreMatrix <- matrix(NA, nrow=9, ncol=2) 
rownames(FullYTHDF3_StoreMatrix) <- c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")
colnames(FullYTHDF3_StoreMatrix) <- c("crossValid","indepTest")


format="mixup"
encodingMethod <- c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")
genomeMethod <- c("CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","onehot")
tst1 <- c("CONPOSITION","Chemiprop_SeqFeature")
tst2 <- c("EIIP","AutoCo_PseKNC")
tst3 <- c("AutoCo","PSNP","onehot")
tst_all <- c(tst1,tst2,tst3)

for(method in tst_all){
  for (i in 1:10){
    # get the negative sample 
    storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/YTHDF3/Full/negative",i,".rds",sep = "")  
    testN <- readRDS(storeNegMatrix[i,])
    
    #prepare the sample used in the genome encoding 
    testN_genome <- testN-2
    
    #prepare the sample used in other encoding method as a character
    testN <- as.character(DNAStringSet(Views(Hsapiens,testN+18)))
    # combine the two positive and negative sample together
    testAll <-c(testN,testT) 
    
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
    
    # #convert it to the dataframe format
    dataframeR<-as.data.frame(matrixResults)
    
    # label the data 
    label_train <- c(rep("Neg",(length(testN))), rep("Pos",(length(testT))))
    label_train <- factor(label_train , labels = c("Neg","Pos") )
    
    #add the label to the dataframe
    dataframeR$class<-label_train
    
    #split the data
    intrain <- createDataPartition(y = dataframeR$class, p= 0.75, list = FALSE)
    training <- dataframeR[intrain,]
    testing <- dataframeR[-intrain,]
    
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
      
      try(testing$pred_class <- predict(Cmotifvsnon_pred, testing, type="prob"))
      try(BIOmotifvsnon_testppred <- prediction(testing$pred_class$Pos,testing$class))
      try(BIOmotifvsnon_testpppred_auc <- performance(BIOmotifvsnon_testppred,"auc"))
      try(aucMatrix[i,]<-attr(BIOmotifvsnon_testpppred_auc,"y.values")[[1]][1])
      try(rocMatrix[i,] <- max(Cmotifvsnon_pred$results$ROC))
  }
      saveRDS(Cmotifvsnon_pred,"/home/yuxuan.wu/m6A reader/model/YTHDF3_Full_mix.rds")
      #caculate the average AUC after 10 times 
      avgAUC <- apply(aucMatrix,MARGIN = 2, FUN = mean,na.rm=T)
      avgAUC
      
      #caculate the max ROC of the prediction model
      avgROC <- apply(rocMatrix,MARGIN = 2, FUN = mean,na.rm=T)
      avgROC
  
  results <- c(avgROC,avgAUC)
  FullYTHDF3_StoreMatrix[method,] <- results
}

setwd("/home/yuxuan.wu/m6A reader/storeMatrix")
saveRDS(FullYTHDF3_StoreMatrix,"FullYTHDF3_StoreMatrix_mix.rds")

