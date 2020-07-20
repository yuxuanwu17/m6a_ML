library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)
library(ROCR)
library(pROC)
library(kernlab)
library(caret)
library(plyr)

# library(doParallel)
# cl <- makePSOCKcluster(10)
# registerDoParallel(cl)

set.seed(2123)

# positive_Full
FulleIF3a_Pos <- readRDS("/home/kunqi/m6A reader/eIF3a/Full/postive.rds")

# testing_sample
index<- grep(1,FulleIF3a_Pos$GSE65004)
testT_testing <- as.character(DNAStringSet(Views(Hsapiens,FulleIF3a_Pos[index]+20)))


# positive sample (all)
testT<-as.character(DNAStringSet(Views(Hsapiens,FulleIF3a_Pos[-index]+20)))



storeNegMatrix <- matrix(NA, nrow=10, ncol = 1)

storeMatrix <- matrix(NA, ncol = 20, nrow = 1500)

list1 <- list()

i <- 1

    # get the negative sample 
    storeNegMatrix[i,] <- paste("/home/kunqi/m6A reader/eIF3a/Full/negative",i,".rds",sep = "")  
    testN_read <- readRDS(storeNegMatrix[i,])
    
    #negative training sample
    testN <- testN_read[-index]
    
    #prepare the sample used in other encoding method as a character
    testN <- as.character(DNAStringSet(Views(Hsapiens,testN+20)))
    testN_testing <- as.character(DNAStringSet(Views(Hsapiens,testN_read[index]+20)))
    
    train_All <-c(testT,testN) 
    testing_All<- c(testT_testing,testN_testing)
    
    num <- c(1:length(train_All))
    list1[[(2*i-1)]] <- data.frame(num,train_All)
    
    num <- c(1:length(testing_All))
    list1[[(2*i)]] <- data.frame(num,testing_All)
    
    x <- merge(x = list1[[(2*i-1)]],y=list1[[(2*i)]],by = "num",all.x =T)
    z <- as.data.frame(x)

write.csv(z, file = "/home/yuxuan.wu/deep_learning/sequence_data/eif3a_full_40.csv")

  # for (i in 1:10) {
  #   if (i==1) {
  #     x <- merge(x = list1[[(2*i-1)]],y=list1[[(2*i)]],by = "num",all.x =T)
  #     z <- as.data.frame(x)
  #   }else{
  #     x <- merge(x = list1[[(2*i-1)]],y=list1[[(2*i)]],by = "num",all.x =T)
  #     y <- as.data.frame(x)
  #     z <- cbind(z,y)
  #   }
  # }



