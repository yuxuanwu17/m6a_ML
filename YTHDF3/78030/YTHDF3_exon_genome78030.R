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

# positive_Exon
ExonYTHDF3_Pos <- readRDS("/home/kunqi/m6A reader/YTHDF3/Exon/postive.rds")

# testing_sample
index<- grep(1,ExonYTHDF3_Pos$GSE78030)
testT_testing <- as.character(DNAStringSet(Views(Hsapiens,ExonYTHDF3_Pos[index]+20)))
testT_testing_genome <- ExonYTHDF3_Pos[index]

# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,ExonYTHDF3_Pos[-index]+20)))
testT_genome <- ExonYTHDF3_Pos[-index]

# genome method
genome_method <- c("genome","genome_chemiPro")
format=""

# create the matrix to store the auc and roc value 
aucMatrix<- matrix (NA,ncol = 1,nrow = 10)
storeNegMatrix <- matrix(NA, nrow=10, ncol = 1)
ExonYTHDF3_StoreMatrix <- matrix(NA, nrow=9, ncol=2) 
rownames(ExonYTHDF3_StoreMatrix) <- c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")
colnames(ExonYTHDF3_StoreMatrix) <- c("crossValid","indepTest")
ExonYTHDF3_StoreMatrix <- readRDS("/home/yuxuan.wu/m6A reader/storeMatrix/storeMatrix_Exon_YTHDF3_mix.rds")

# begin the loop
for (method in genome_method){
  for(i in 1:10){
    storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/YTHDF3/Exon/negative",i,".rds",sep = "")  
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
    
    format=""
    if (method=="genome"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/genome.R")
    }
    
    if (method=="genome_chemiPro"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/genome_chemiPro.R")
    }
    
    training <- as.data.frame(matrixResults)
    testing <- as.data.frame(testingResults)
    
    
    # add labels for training sample
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
saveRDS(ExonYTHDF3_StoreMatrix,"ExonYTHDF3_StoreMatrix_genome_78030.rds")
