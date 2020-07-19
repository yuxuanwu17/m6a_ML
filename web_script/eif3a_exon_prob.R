library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg19)
library(fitCons.UCSC.hg19)
library(phastCons100way.UCSC.hg19)
library(caret)
library(ROCR)
library(pROC)
library(m6ALogisticModel)
source("/home/kunqi/m7G/method/class1.R")
source("/home/kunqi/m7G/method/class2.R")
source("/home/kunqi/m7G/method/class3.R")
source("/home/kunqi/m7G/method/class4.R")
source("/home/kunqi/m7G/method/class5.R")
source("/home/kunqi/m7G/method/class6.R")


# 

GenoFgenreation <- function(data){
  Notmethylated_sample <- readRDS("/home/kunqi/m6A reader/eIF3a/Exon/negative1.rds")
  mature_sample <- readRDS("/home/kunqi/m6A reader/eIF3a/Exon/postive.rds")
  analysis_data <- data

  mcols(mature_sample) <- NULL
  analysis_data <- c(analysis_data-2,mature_sample,Notmethylated_sample-2)
  matureSE <- SummarizedExperiment()
  rowRanges(matureSE) <- analysis_data 

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
  matureVSnon  <- mcols(mature_FE)
  top_selected_final <- readRDS("/home/yuxuan.wu/m6A reader/top_selected_final/feature_eIF3a_Exon.rds")
  matureVSnon <- matureVSnon[,top_selected_final]
  return(matureVSnon[1:length(data),])
}

SeqFgeneration <- function(data){
  analysisData <- data
  analysisData<-as.character(DNAStringSet(Views(Hsapiens,analysisData+20)))
  CP <- ChemicalProperty(analysisData)
  return(CP)
}

SitePredictionCollection <- function(data){
  GF <- GenoFgenreation(data)
  SeqFR <- SeqFgeneration(data)
  BothFR <- cbind(GF,SeqFR)
  reModel2 <- readRDS("/home/yuxuan.wu/m6A reader/model/eIF3a_Exon.rds")
  results <- predict(reModel2, BothFR ,type="prob")
  results <- results
  return(results)
}



