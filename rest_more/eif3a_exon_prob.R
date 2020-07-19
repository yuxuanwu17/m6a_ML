source("/home/yuxuan.wu/m6A reader/web_script/eif3a_exon_prob.R")
testPos <- readRDS("/home/yuxuan.wu/rest.rds")
testPos <- testPos+2
posProb <- SitePredictionCollection(testPos)
saveRDS(posProb,"/home/yuxuan.wu/m6A reader/rest_more/eif3a_exon.rds")