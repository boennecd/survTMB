library(testthat)
library(survTMB)
library(survival)

get_test_file_name <- function(f_name){
  f_name <- file.path("test-data", paste0(f_name))
  f_name_use <- file.path("..", "..", "survTMB", f_name)

  out <- if(file.exists(f_name_use))
    f_name_use
  else {
    f_name_use <- file.path("..", "..", "inst", f_name)
    if(file.exists(f_name_use))
      f_name_use
    else
      file.path("inst", f_name)
  }

  stopifnot(file.exists(out))
  out
}

# set.seed(1)
# eortc <- subset(coxme::eortc, center %in% head(unique(center), 10))
# eortc <- eortc[sample(NROW(eortc)), ]
# saveRDS(eortc, file.path("tests", "testthat", "eortc.RDS"))
eortc <- if(file.exists("eortc.RDS"))
  readRDS("eortc.RDS") else
    readRDS(file.path("tests", "testthat", "eortc.RDS"))

formals(expect_known_value) $update <- FALSE
formals(expect_known_output)$update <- FALSE
options(width = 80, digits = 3, useFancyQuotes = FALSE,
        setWidthOnResize = FALSE)

test_res_dir <- if(!dir.exists("test-res"))
  file.path("tests", "testthat", "test-res") else
    "test-res"

# if these change then update the man page!
.MGSM_ADFun_members <- c("laplace", "gva", "snva", "y", "event",
                         "X", "XD", "Z", "grp", "terms", "link", "cl",
                         "opt_func", "dense_hess", "sparse_hess")
.MGSM_fit_members <- c("params", "va_params", "link", "ADFun_cl", "fit_cl",
                       "method", "optim", "is_va", "fix_names", "rng_names")
