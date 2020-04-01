# set.seed(1)
# eortc <- subset(coxme::eortc, center %in% head(unique(center), 10))
# eortc <- eortc[sample(NROW(eortc)), ]
# saveRDS(eortc, file.path("tests", "testthat", "eortc.RDS"))
eortc <- if(file.exists("eortc.RDS"))
  readRDS("eortc.RDS") else
    readRDS(file.path("tests", "testthat", "eortc.RDS"))

formals(expect_known_value) $update <- FALSE
formals(expect_known_output)$update <- FALSE
options(width = 80, digits = 4, useFancyQuotes = FALSE,
        setWidthOnResize = FALSE)
