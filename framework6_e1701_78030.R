library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(ROCR)
library(pROC) 
# library(doParallel)
# registerDoParallel(makePSOCKcluster(10))
set.seed(520)

# positive_Exon
ExonYTHDF3_Pos <- readRDS("/home/kunqi/m6A reader/YTHDF3/Exon/postive.rds")

# readRDS("/home/kunqi/m6A reader/YTHDF3/Full/postive.rds")
readRDS("/home/kunqi/m6A reader/YTHDF3/Exon/negative1.rds")

# select  as independent-test 
index_train <- grep(0,ExonYTHDF3_Pos$GSE78030)   ##### this is the training set as positive
ExonYTHDF3_Pos_indep <- ExonYTHDF3_Pos[index_train]

# testing_sample 
index_testing <- grep(1,ExonYTHDF3_Pos$GSE78030) #### this is the testing set as positive
testT_testing <- as.character(DNAStringSet(Views(Hsapiens,ExonYTHDF3_Pos[index_testing]+20)))

testT_testing_genome <- ExonYTHDF3_Pos[index_testing]+20
# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,ExonYTHDF3_Pos_indep+20)))
testT_genome <- ExonYTHDF3_Pos_indep+20

library(caret)
library(e1071)

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
ExonYTHDF3_StoreMatrix <- matrix(NA, nrow=9, ncol=1) 

rownames(ExonYTHDF3_StoreMatrix) <- c("genome_chemiPro","CONPOSITION","Chemiprop+SeqFeature","EIIP","AutoCo+PseKNC","AutoCo","PSNP","genome","onehot")
colnames(ExonYTHDF3_StoreMatrix) <- "indepTest"

format=""
num=0
for (method in c("CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","genome_chemiPro","onehot")){
num <- num+1
  for (i in 1:10){
  i=1
  # get the negative sample 
  storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/YTHDF3/Exon/negative",i,".rds",sep = "")  
  testN_read <- readRDS(storeNegMatrix[i,])
  testN <- sample(testN_read,size=length(index_train))
  testN_genome <- testN-2
  testN <- as.character(DNAStringSet(Views(Hsapiens,testN+18)))
  testAll <-c(testN,testT) 
  
  # analysis the test sample
  testN_testing <- sample(testN, size = length(index_testing))
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
  
  #svm model
  Cmotifvsnon_pred <- svm(class ~ ., data = training, 
                            type="C",
                            cross=5,
                            probability = TRUE,
                            kernel="radial")
  
  testing$pred_class <- predict(Cmotifvsnon_pred, testing)
  
  # ROC curve & AUC value
  svm_roc <- roc(testing$class,as.numeric(testing$pred_class))
  
  BIOmotifvsnon_testppred <- prediction(as.numeric(testing$pred_class),testing$class)
  BIOmotifvsnon_testpppred_auc <- performance(BIOmotifvsnon_testppred,"auc")
  aucMatrix[i,]<-attr(BIOmotifvsnon_testpppred_auc,"y.values")[[1]][1]
}             #end the first for loop (10 times)

#caculate the average AUC after 10 times 
avgAUC <- apply(aucMatrix,MARGIN = 2, FUN = mean)
indepTest <- avgAUC
indepTest

#caculate the max ROC of the prediction model

ExonYTHDF3_StoreMatrix[num,] <- indepTest
}

setwd("/home/yuxuan.wu/m6A reader/storeMatrix")
saveRDS(ExonYTHDF3_StoreMatrix,"storeMatrix_ExonYTHDF3.rds")
# ExonYTHDF3_StoreMatrix<- readRDS("storeMatrix_ExonYTHDF3.rds")
