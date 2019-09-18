#' Fits the drop-model or missingness-model using logistic regression for estimating and computing inverse probability censoring weights.
#'
#' \code{dropmodel} fits a logitistic regression model via \code{\link[stats]{glm}} to predict drop-out or missingness
#' to estimate and compute inverse probability censoring weights (IPCWs).
#' The output contains
#' 1. Coeffficient estimates and inference,
#' 2. A data set that contains all
#' original variables and computed IPCWs (variable name \code{ipcw}), and
#' 3. The cluster-level score functions that will  be used to compute the robust variance in the final model.
#' Users can proceed to \code{\link{ipccwGEE}} to fit a marginal model via inverse probability censoring cluster-weighted generalized estimating equations.
#'
#' The \code{data} must be provided in case level or equivalently in `long' format.
#'
#' @param formula a formula expression for the drop-out model. The outcome should be the missingness indicator (binary).
#' @param data an optional data frame containing the variables provided in
#' \code{formula}, \code{cluster.var}, \code{unit.var} and \code{time.var}.
#' @param cluster.var a vector that identifies the clusters.
#' @param unit.var a vector that identifies the unit within a cluster.
#' @param time.var a vector that identifies the repeated observation of a unit.
#'
#' @return Returns a list with components:
#' \item{call}{the matched call.}
#' \item{coefficients}{the estimated regression parameters, standard errors, z-values, and  p-values of the logistic regression model, fitted via \code{glm}.}
#' \item{outdata}{the original data set used in the function plus computed IPCW (variable name "ipcw") as a new variable.}
#' \item{Slist}{the cluster-level score functions that will be used to compute the robust variance in the final model.}
#' @author Aya Mitani
#' @examples
#' data(dental)
#' outdrop <- dropmodel(toothstat ~ prevmaxcal5mm + basenumteeth, data = dental,
#' cluster.var = subject, unit.var = tooth, time.var = visit)
#' outdat <- outdrop$outdata
#' Slist <- outdrop$Slist
#' @export



dropmodel <- function(formula, data, cluster.var, unit.var, time.var){

  thenames <- names(data)

  names(data)[names(data) == deparse(substitute(cluster.var))] <- "cluster.var"
  names(data)[names(data) == deparse(substitute(unit.var))] <- "unit.var"
  names(data)[names(data) == deparse(substitute(time.var))] <- "time.var"

  yname <- all.vars(formula)[1]
  xname <- all.vars(formula)[-1]

  outmod <- glm(formula, family = binomial("logit"), data = data)
  summod <- summary(outmod)
  lambdavec <- predict(outmod, data, "response")
  lambdavec[is.na(lambdavec)] <- 1

  ### compute weights
  newid <- paste(data$cluster.var, data$unit.var)
  wtveclist <- vector("list", length(unique(newid)))
  for(ij in unique(newid)){

    lambdavecij <- lambdavec[newid == ij]
    cumlambdavecij <- lambdavecij
    for (jj in 2:length(lambdavecij)){
      cumlambdavecij[jj] <- cumlambdavecij[jj-1] * cumlambdavecij[jj]
    }
    wtvecij <- 1/cumlambdavecij
    wtveclist[[which(unique(newid) == ij)]] <- wtvecij
  }


  ### score function for each cluster
  dat <- data[c("cluster.var", "unit.var", "time.var", yname, xname)]
  dat2 <- dat[-which(dat$time.var==1),]
  yvar <- as.matrix(dat2[yname])
  xvar <- as.matrix(dat2[xname])
  estalphavec <- as.numeric(outmod$coefficients)
  Y <- yvar
  X <- cbind(rep(1, length(Y)), xvar)
  cluster <- as.numeric(dat2$cluster.var)


  Slist <- list()
  ii <- 0

  for (i in unique(cluster)){

    ii <- ii + 1

    y <- as.matrix( Y[cluster==i] )
    x <- as.matrix( X[cluster==i,] )

    u <- as.matrix( exp( x %*% estalphavec ) / ( 1 + exp( x %*% estalphavec ) ) )

    s <- t(x) %*% as.matrix( y - u )

    row.names(s) <- c()

    Slist[[ii]] <- s

  }

  names(data) <- thenames
  data$ipcw <- unlist(wtveclist)

  result <- list()
  result$call <- match.call()
  result$coefficients <- summod$coefficients
  result$outdata <- data
  result$Slist <- Slist
  result

}



