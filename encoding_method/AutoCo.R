if (format =="mixup"){
for(i in 1:ncol(RNA_SCORE)){
  if(i==1){
    matrixResults <- autoCovariance(testAll,convert=1,descriptor = colnames(RNA_SCORE)[i])
  }
  matrixResults <- cbind(matrixResults,autoCovariance(testAll,convert=1,descriptor = colnames(RNA_SCORE)[i]))
}
}else{
########################## two groups training & testing
for(i in 1:ncol(RNA_SCORE)){
  if(i==1){
    matrixResults <- autoCovariance(testAll,convert=1,descriptor = colnames(RNA_SCORE)[i])
  }
  matrixResults <- cbind(matrixResults,autoCovariance(testAll,convert=1,descriptor = colnames(RNA_SCORE)[i]))
}

for(i in 1:ncol(RNA_SCORE)){
  if(i==1){
    testingResults <- autoCovariance(testAll_testing,convert=1,descriptor = colnames(RNA_SCORE)[i])
  }
  testingResults <- cbind(testingResults,autoCovariance(testAll_testing,convert=1,descriptor = colnames(RNA_SCORE)[i]))
}
}