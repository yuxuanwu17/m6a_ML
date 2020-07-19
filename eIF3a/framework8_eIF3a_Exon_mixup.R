library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(ROCR)
library(pROC)
# library(doParallel)
# cl <- makePSOCKcluster(8)
# registerDoParallel(cl)

set.seed(520)

# positive_Exon
ExoneIF3a_Pos <- readRDS("/home/kunqi/m6A reader/eIF3a/Exon/postive.rds")
# show the first negative sample 
readRDS("/home/kunqi/m6A reader/eIF3a/Exon/negative1.rds")


# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,ExoneIF3a_Pos+20)))
testT_genome <- ExoneIF3a_Pos

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
ExoneIF3a_StoreMatrix <- matrix(NA, nrow=9, ncol=2) 
rownames(ExoneIF3a_StoreMatrix) <- c("genome_chemiPro","CONPOSITION","Chemiprop+SeqFeature","EIIP","AutoCo+PseKNC","AutoCo","PSNP","genome","onehot")
colnames(ExoneIF3a_StoreMatrix) <- c("crossValid","indepTest")

num=0
format="mixup"
for (method in c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")){
  num <- num+1
  for (i in 1:10){
    # get the negative sample 
    storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/eIF3a/Exon/negative",i,".rds",sep = "")  
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
    intrain <- createDataPartition(y = dataframeR$class, p= 0.8, list = FALSE)
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
  # ExoneIF3a_StoreMatrix <- rbind(results,ExoneIF3a_StoreMatrix)
  ExoneIF3a_StoreMatrix[num,] <- results
}

setwd("/home/yuxuan.wu/m6A reader/storeMatrix")
saveRDS(ExoneIF3a_StoreMatrix,"storeMatrix_ExoneIF3a_0915_mix.rds")

# saveRDS(ExoneIF3a_StoreMatrix,"storeMatrix.rds")
# ExoneIF3a_StoreMatrix<- readRDS("storeMatrix.rds")
# ExoneIF3a_StoreMatrix <- readRDS("storeMatrix_ExoneIF3a_0915_mix.rds")
