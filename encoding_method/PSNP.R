if(format!="mixup"){
matrixResults<- PSNP(testN,testT,NI=3,NTYPE = "DNA")
testingResults <- PSNP(testN_testing,testT_testing,NI=3,NTYPE = "DNA")
}else{
  matrixResults<- PSNP(testN,testT,NI=3,NTYPE = "DNA")
}