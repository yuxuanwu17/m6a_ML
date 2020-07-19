if(format!="mixup"){
matrixResults <- CONPOSITION(testAll,NI=3,NTYPE = "DNA",Freq = 2)
testingResults <- CONPOSITION(testAll_testing,NI=3,NTYPE = "DNA", Freq = 2)
}else{
  matrixResults <- CONPOSITION(testAll,NI=3,NTYPE = "DNA",Freq = 2)
}