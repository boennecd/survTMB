# set.seed(1)
# eortc <- subset(coxme::eortc, center %in% head(unique(center), 10))
# eortc <- eortc[sample(NROW(eortc)), ]
# saveRDS(eortc, file.path("tests", "testthat", "eortc.RDS"))
eortc <- if(file.exists("eortc.RDS"))
  readRDS("eortc.RDS") else
    readRDS(file.path("tests", "testthat", "eortc.RDS"))
