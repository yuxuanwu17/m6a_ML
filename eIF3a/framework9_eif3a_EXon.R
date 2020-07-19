library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(ROCR)
library(pROC)
library(kernlab)
# library(doParallel)
# cl <- makePSOCKcluster(20)
# registerDoParallel(cl)

set.seed(5123)

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
ExoneIF3a_StoreMatrix <- matrix(NA, nrow=9, ncol=2) 
rownames(ExoneIF3a_StoreMatrix) <- c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")
colnames(ExoneIF3a_StoreMatrix) <- c("crossValid","indepTest")

encodingMethod <- c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")
genomeMethod <- c("CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","onehot")
tst1 <- c("CONPOSITION","Chemiprop_SeqFeature")
tst2 <- c("EIIP","AutoCo_PseKNC")
tst3 <- c("AutoCo","PSNP","onehot")
tst_all <- c(tst1,tst2,tst3)
format=""
# for (method in c("genome_chemiPro","CONPOSITION","Chemiprop_SeqFeature","EIIP","AutoCo_PseKNC","AutoCo","PSNP","genome","onehot")){
for(method in tst_all) {
  for (i in 1:10){
    # get the negative sample 
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
    
    if(method=="genome" || method=="genome_chemiPro"){
      fitControl <- trainControl(method = "cv",
                                 number = 5,
                                 classProbs = TRUE,
                                 summaryFunction = twoClassSummary)
      
      myModel <- train(class ~ ., data = training, 
                       method = 'svmRadial',  
                       preProc = c("center", "scale"),
                       trControl = fitControl,
                       metric = "ROC")
      
      # prediction
      testing$pred_class <- predict(myModel, testing)
      pred_result <- prediction(as.numeric(testing$pred_class),testing$class)
      pred_result_auc <- performance(pred_result,"auc")
      aucMatrix[i,]<-attr(pred_result_auc,"y.values")[[1]][1]
      
      # importance
      var_genome_imp <- varImp(myModel)$importance
      
      # normalization
      var_genome_imp_normal <- (var_genome_imp-min(var_genome_imp))/(max(var_genome_imp)-min(var_genome_imp))
      
      # store in a matrix
      VarMatrix <- matrix(NA,nrow = length(rownames(var_genome_imp_normal)),ncol = 2)
      VarMatrix[,1]  <- var_genome_imp_normal[,1][order(var_genome_imp_normal[,1],decreasing = T)]
      VarMatrix[,2]<- colnames(training)[order(var_genome_imp_normal[,1],decreasing = T)] 
      colnames(VarMatrix) <- c("importance","features")
      VarMatrix <- as.data.frame(VarMatrix)
      VarMatrix$importance <- as.numeric(as.character(VarMatrix$importance))
      VarMatrix$features <- factor(VarMatrix$features,levels = VarMatrix$features)
      
      # find the optimal combination of features
      ROC_max <- matrix(NA, nrow = length(var_genome_imp[,1]), ncol = 1)
          for(j in 1:length(var_genome_imp[,1])){
        # for
          top_selected <- which(colnames(training)%in%VarMatrix[,2][1:j]==TRUE)
          reModel <- train(class ~ ., data = training[,c(top_selected,ncol(training))], 
                         method = 'svmRadial',  
                         preProc = c("center", "scale"),
                         trControl = fitControl,
                         metric = "ROC")
          ROC_max[j,]<-max(reModel$results$ROC)
      }
      
          top_selected_final <- which(colnames(training)%in%VarMatrix[,2][1:which(ROC_max==max(ROC_max))]==TRUE)
          reModel2 <- train(class ~ ., data = training[,c(top_selected_final,ncol(training))], 
                            method = 'svmRadial',  
                            preProc = c("center", "scale"),
                            trControl = fitControl,
                            metric = "ROC")
          
          # prediction
          testing$pred_class2 <- predict(reModel2, testing)
          pred_result2 <- prediction(as.numeric(testing$pred_class2),testing$class)
          pred_result_auc2 <- performance(pred_result2,"auc")
          aucMatrix2[i,]<-attr(pred_result_auc2,"y.values")[[1]][1]
          aucMatrix <- aucMatrix2
          avgAUC <- apply(aucMatrix,MARGIN = 2, FUN = mean)
          avgAUC
          
          avgROC <- max(ROC_max)
          avgROC
    }
    else{
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

    #caculate the average AUC after 10 times 
    avgAUC <- apply(aucMatrix,MARGIN = 2, FUN = mean)
    avgAUC
    
    #caculate the max ROC of the prediction model
    avgROC <- apply(rocMatrix,MARGIN = 2, FUN = mean)
    avgROC
    }
  }             #end the first for loop (10 times)
  
  results <- c(avgROC,avgAUC)

  ExoneIF3a_StoreMatrix[method,] <- results
}

setwd("/home/yuxuan.wu/m6A reader/storeMatrix")
saveRDS(ExoneIF3a_StoreMatrix,"storeMatrix_ExoneIF3a_test.rds")

# saveRDS(ExoneIF3a_StoreMatrix,"storeMatrix.rds")
# ExoneIF3a_StoreMatrix<- readRDS("storeMatrix.rds")
# ExoneIF3a_StoreMatrix <- readRDS("storeMatrix_ExoneIF3a_09_22.rds")
