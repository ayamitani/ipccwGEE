#' @export

#----------------------------------
# Print summary
#----------------------------------

print.summary.ipccwgee <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients)
  cat("\nNumber of iterations:\n")
  print(x$niter)
}
