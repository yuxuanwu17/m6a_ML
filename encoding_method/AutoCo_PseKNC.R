if(format!="mixup"){
for(i in 1:3){
  if(i==1){
    matrixResults_ACC <- autoCovariance(testAll,convert=1,descriptor = colnames(RNA_SCORE)[11-i])
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
}else{
  for(i in 1:3){
    if(i==1){
      matrixResults_ACC <- autoCovariance(testAll,convert=1,descriptor = colnames(RNA_SCORE)[11-i])
    }
    matrixResults_ACC <- cbind(matrixResults_ACC,autoCovariance(testAll,convert=1,descriptor = colnames(RNA_SCORE)[i]))
  }
  matrixResults_pse <- PseKNC(testAll, NI=2,NTYPE = "DNA")
  matrixResults <- cbind(matrixResults_ACC,matrixResults_pse)
}