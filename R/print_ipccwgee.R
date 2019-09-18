#' @export

#----------------------------------
# Print output
#----------------------------------

print.ipccwgee <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}
