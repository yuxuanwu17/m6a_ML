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
library(m6ALogisticModel)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fitCons.UCSC.hg19)
library(phastCons100way.UCSC.hg19)
set.seed(42)

# matrix
storeNegMatrix <- matrix(NA, nrow=10, ncol = 1)
aucMatrix <- matrix(NA, nrow = 10, ncol = 1)
aucMatrix2 <- matrix(NA, nrow = 10, ncol = 1)


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

format=""
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


# generate genome features
source("/home/yuxuan.wu/m6A reader/encoding_method/genome.R")
  
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

# svm model
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
# saveRDS(VarMatrix, "/home/zhendi/2019-fall/result/varmatrix.rds")

# plot
# Varplot <- ggplot(VarMatrix,aes(x=features,y=importance))+
#   geom_histogram(stat = "identity")+
#   theme(axis.text.x = element_text(angle = 90))
# Varplot
# ggsave("/home/zhendi/2019-fall/result/Feature_select_varplot.pdf",Varplot)

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
  # saveRDS(ROC_max,"/home/zhendi/2019-fall/result/feature_selection_genome_auc_max.rds")
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
 }

maxROC <- max(ROC_max)
avgAUC <- apply(aucMatrix2,MARGIN = 2, FUN = mean, na.rm=TRUE)
results <- c(maxROC,avgAUC)

setwd("/home/yuxuan.wu/m6A reader/storeMatrix")
# ExoneIF3a_StoreMatrix <- readRDS("storeMatrix_ExoneIF3a_0912.rds")
# ExoneIF3a_StoreMatrix["genome",] <- results
# saveRDS(ExoneIF3a_StoreMatrix,"storeMatrix_ExoneIF3a_09_22.rds")
