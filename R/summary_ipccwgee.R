#' @export

#------------------------------------
# Summary
#------------------------------------

summary.ipccwgee <- function(object, ...) {

  TAB <- cbind(Estimate = coef(object),
               Robust_S.E. = object$robust.se,
               Wald_Chisq = object$wald.chisq,
               p.value = object$p.value)
  TAB <- round(TAB, 5)
  colnames(TAB) <- c("Estimate", "Std.Error", "Wald.Chisq", "Pr > Chisq")
  rownames(TAB) <- object$coef.names
  res <- list(coefficients = TAB,
              niter = object$niter,
              call = object$call)
  class(res) <- "summary.ipccwgee"

  res
}
