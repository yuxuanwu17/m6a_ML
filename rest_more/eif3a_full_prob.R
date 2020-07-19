source("/home/yuxuan.wu/m6A reader/web_script/eif3a_full_prob.R")
testPos <- readRDS("/home/yuxuan.wu/m6A reader/rest_more.rds")
testPos <- testPos+2
posProb <- SitePredictionCollection(testPos)
saveRDS(posProb,"/home/yuxuan.wu/m6A reader/rest_more/eif3a_full.rds")
