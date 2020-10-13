#' @importFrom TMB config
.onLoad <- function(libname, pkgname){
  fix_atomic_seqfault()
  # see https://groups.google.com/g/tmb-users/c/4EZSIOJQF6Q
  config(tape.parallel = FALSE, DLL = "survTMB")
}
