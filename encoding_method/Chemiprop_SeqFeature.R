if(format!="mixup"){
  matrixResults_CP<- ChemicalProperty(testAll)
  matrixResults_SF <- sequenceFeatures(testAll,NTYPE = "DNA")
  testingResults_CP <- ChemicalProperty(testAll_testing)
  testingResults_SF <- sequenceFeatures(testAll_testing,NTYPE = "DNA")
  matrixResults <- cbind(matrixResults_CP,matrixResults_SF)
  testingResults <- cbind(testingResults_CP,testingResults_SF)
  
}else{
  matrixResults_CP<- ChemicalProperty(testAll)
  matrixResults_SF <- sequenceFeatures(testAll,NTYPE = "DNA")
  matrixResults <- cbind(matrixResults_CP,matrixResults_SF)
}