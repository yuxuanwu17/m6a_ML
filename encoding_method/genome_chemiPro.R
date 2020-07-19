library(m6ALogisticModel)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fitCons.UCSC.hg19)
library(phastCons100way.UCSC.hg19)

Additional_features_hg19 = list(
  HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
  miR_targeted_genes = miR_targeted_genes_grl,
  TargetScan = TargetScan_hg19_gr,
  Verified_miRtargets = verified_targets_gr,
  METTL3_TREW = METTL3_TREW,
  METTL14_TREW = METTL14_TREW,
  WTAP_TREW = WTAP_TREW,
  METTL16_CLIP = METTL16_CLIP,
  ALKBH5_PARCLIP = ALKBH5_PARCLIP,
  FTO_CLIP = FTO_CLIP,
  FTO_eCLIP = FTO_eCLIP
)

if(format!="mixup"){
matureSE <- SummarizedExperiment()
mcols(testT_genome) <- NULL
rowRanges(matureSE) <- c(testT_genome,testN_genome)
mature_FE <- predictors_annot(se = matureSE,
                              txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                              bsgnm = Hsapiens,
                              fc = fitCons.UCSC.hg19,
                              motif = c(),
                              pc = phastCons100way.UCSC.hg19,
                              struct_hybridize = Struc_hg19,
                              feature_lst = Additional_features_hg19,
                              hk_genes_list = HK_hg19_eids,
                              genes_ambiguity_method = "average")
matrixResults_genome <- mcols(mature_FE)

matureSE <- SummarizedExperiment()
mcols(testT_testing_genome) <- NULL
rowRanges(matureSE) <- c(testT_testing_genome,testN_testing_genome)
mature_FE <- predictors_annot(se = matureSE,
                              txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                              bsgnm = Hsapiens,
                              fc = fitCons.UCSC.hg19,
                              motif = c(),
                              pc = phastCons100way.UCSC.hg19,
                              struct_hybridize = Struc_hg19,
                              feature_lst = Additional_features_hg19,
                              hk_genes_list = HK_hg19_eids,
                              genes_ambiguity_method = "average")
testingResults_genome <- mcols(mature_FE)
matrixResults_chemi <- ChemicalProperty(testAll)
testingResults_chemi <- ChemicalProperty(testAll_testing)

matrixResults <- cbind(matrixResults_genome,matrixResults_chemi)
testingResults <- cbind(testingResults_genome,testingResults_chemi)
}else{
  matureSE <- SummarizedExperiment()
  mcols(testT_genome) <- NULL
  rowRanges(matureSE) <- c(testT_genome,testN_genome)
  mature_FE <- predictors_annot(se = matureSE,
                                txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                bsgnm = Hsapiens,
                                fc = fitCons.UCSC.hg19,
                                motif = c(),
                                pc = phastCons100way.UCSC.hg19,
                                struct_hybridize = Struc_hg19,
                                feature_lst = Additional_features_hg19,
                                hk_genes_list = HK_hg19_eids,
                                genes_ambiguity_method = "average")
  matrixResults_genome <- mcols(mature_FE)
  
  matrixResults_chemi <- ChemicalProperty(testAll)
  matrixResults <- cbind(matrixResults_genome,matrixResults_chemi)
}