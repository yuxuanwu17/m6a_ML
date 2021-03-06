# library(doParallel)
# registerDoParallel(makePSOCKcluster(5))

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

library(MLmetrics)

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

ExoneIF3a_Pos <- readRDS("/home/kunqi/m6A reader/eIF3a/Exon/postive.rds")
# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,ExoneIF3a_Pos+20)))
testT_genome <- ExoneIF3a_Pos

# genome method
genome_method <- c("genome","genome_chemiPro")
format="mixup"

# create the matrix to store the auc and roc value 
aucMatrix<- matrix (NA,ncol = 1,nrow = 10)
prauc <- matrix(NA,ncol = 1,nrow=10)
accMatrix <- matrix(NA,ncol = 1,nrow=10)
mattMatrix <- matrix(NA,ncol = 1,nrow=10)

storeNegMatrix <- matrix(NA, nrow=10, ncol = 1)
ExoneIF3a_StoreMatrix <- readRDS("/home/yuxuan.wu/m6A reader/storeMatrix/ExoneIF3a_StoreMatrix_mix.rds")


# begin the loop
method <- c("genome_chemiPro")
    i = 1 
    # get the negative sample 
    storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/eIF3a/Exon/negative",i,".rds",sep = "")  
    testN <- readRDS(storeNegMatrix[i,])
    
    #prepare the sample used in the genome encoding 
    testN_genome <- testN
    
    #prepare the sample used in other encoding method as a character
    testN <- as.character(DNAStringSet(Views(Hsapiens,testN+18)))
    # combine the two positive and negative sample together
    testAll <-c(testN,testT) 
    
    rest_genome <-readRDS("/home/yuxuan.wu/rest.rds")
    rest_seq <- as.character(DNAStringSet(Views(Hsapiens,rest_genome+18)))
    
    if (method=="genome_chemiPro"){
      source("/home/yuxuan.wu/m6A reader/encoding_method/genome.R")
      genomeMatrix <- as.data.frame(matrixResults)
      resMatrix <- as.data.frame(resResults)
      source("/home/yuxuan.wu/m6A reader/encoding_method/Chemiprop.R")
      chemi_index <- 58:ncol(matrixResults)
      resResult <- as.data.frame(resResults)
      
      
      dataframeR <- cbind(genomeMatrix,matrixResults)
      resCbd <- cbind(resMatrix,resResult)
      
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
      
      finalResults <- predict(reModel2, testing[,c(top_selected_final,chemi_index,(ncol(testing)))], type="prob")
      
}


saveRDS(finalResults,"/home/yuxuan.wu/m6A reader/vennPlot/ExoneIF3a_mixup.rds")


