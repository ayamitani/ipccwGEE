#' Inverse probability censoring cluster weighted GEE for clustered longitudinal data with informative cluster size and informative drop-out.
#'
#' Solves the cluster-weighted generalized estimating equations
#' with inverse probability censoring weights for correlated binary
#' responses in clustered longitudinal data with informative cluster size and informative drop-out
#' assuming a logit model for the marginal probabilities.
#'
#' The \code{data} must be provided in case level or equivalently in `long' format.
#'
#' @param formula a formula expression as for the regression model.
#' @param data an optional data frame containing the variables provided in
#' \code{formula}, \code{cluster.var}, \code{unit.var} and \code{ipcw}.
#' @param cluster.var a vector that identifies the clusters.
#' @param unit.var a vector that identifies the unit within a cluster.
#' @param ipcw a vector that contains the inverse probability censoring weights for each observation, computing from \code{dropmodel}.
#' @param Slist a list that contains the cluster-level score functions from the drop-out model fitted using \code{dropmodel}.
#' @param show.iter logical. If \code{'TRUE'}, displays the iteration number during the estimation process.
#'
#' @return Returns an object of the class \code{"ipccwgee"}. This has components:
#' \item{call}{the matched call.}
#' \item{coefficients}{the estimated regression parameter vector of the marginal model.}
#' \item{coef.names}{the variable name of the coefficients.}
#' \item{robust.variance}{the estimated "robust" covariance matrix.}
#' \item{robust.se}{the estimated "robust" standard errors.}
#' \item{wald.chisq}{the Wald Chi-square test statistic for coefficient estimates.}
#' \item{p.value}{the p-value based on a Wald Chi-square test statistic that no covariates are statistically significant.}
#' \item{niter}{the number of iterations the model took to converge.}

#' @author Aya Mitani
#' @examples
#' data(dental)
#' ## First fit the drop model or missingness model.
#' outdrop <- dropmodel(toothstat ~ prevmaxcal5mm + basenumteeth, data = dental,
#' cluster.var = subject, unit.var = tooth, time.var = visit)
#' ## Save the outputs.
#' outdat <- outdrop$outdata
#' Slist <- outdrop$Slist
#' ## Restrict data set to include visits where the unit (tooth) was observed.
#' subdental <- subset(outdat, outdat$toothstat == 1)
#' ## Fit marginal model with cluster-weights and IPCWs.
#' ipccwGEEout <- ipccwGEE(maxcal5mm ~ baseage + edu, data = subdental, cluster = subject, unit = tooth,
#' ipcw = ipcw, Slist = Slist, show.iter = TRUE)
#' summary(ipccwGEEout)
#' @export


ipccwGEE <- function(formula, data, cluster.var, unit.var, ipcw, Slist, show.iter = FALSE){

  call <- match.call()
  mcall <- match.call(expand.dots = FALSE)
  mf <- match(c("formula", "data", "cluster.var", "unit.var", "ipcw"), names(mcall), 0L)
  m <- mcall[c(1L, mf)]
  if (is.null(m$cluster.var))
    m$cluster.var <- as.name("id")
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")

  cluster <- model.extract(m, "cluster.var")
  unit <- model.extract(m, "unit.var")
  ipcw <- model.extract(m, "ipcw")

  mterms <- attr(m, "terms")
  xvars <- as.character(attr(mterms, "variables"))[-c(1:2)]
  #names(m) <- c("y", xvars, "cluster", "unit", "ipcw")

  p <- length(xvars) + 1
  q <- dim(Slist[[1]])[1]
  X <- model.matrix(formula, m)  ### NEED MODEL MATRIX
  Y <- model.response(m, "numeric")

  ### Get initial estimates of beta using GLM
  ### initialize
  beta <- rep(0,p)
  bdiff <- 1
  iter <- 0

  while (sum(abs(bdiff))>.000001){

    ### initial estimation of betas using identity correlation structure

    DWZ <- matrix(0, nrow = p)
    DWD <- matrix(0, ncol = p, nrow = p)
    DWZZWD <- matrix(0, ncol = p, nrow = p)

    for (i in unique(cluster)){

      y <- as.matrix( Y[cluster==i] )
      x <- as.matrix( X[cluster==i,] )

      u <- exp( x %*% beta ) / ( 1 + exp( x %*% beta ) )
      dudb <- exp( x %*% beta ) / ( 1 + exp( x %*% beta ) ) ^ 2

      D <- x[,1] * dudb
      for (pp in 2:ncol(x)){
        D <- cbind(D, x[,pp] * dudb)
      }

      ### use identity corr structure
      R <- matrix( nrow = length(y), ncol = length(y), 0 )
      diag(R) <-  1

      vv <- u * ( 1 - u )
      V <- matrix(ncol = length(vv), nrow = length(vv), 0)
      diag(V) <- vv
      W <- V ^ (1/2) %*% R %*% V ^ (1/2)

      invW <- solve(W)

      DWZ <- DWZ + t(D) %*% invW %*% ( y - u )
      DWD <- DWD + t(D) %*% invW %*% D

      DWZZWD <- DWZZWD + t(D) %*% invW %*% ( y - u ) %*% t( y - u )  %*% invW %*% D

    }

    invDWD <- solve(DWD)
    bdiff <- invDWD %*% DWZ
    beta <- beta + bdiff

    covbeta <- invDWD %*% DWZZWD %*% invDWD

    expXb <- exp( X %*% beta )
    U <- expXb / ( 1 + expXb )
    Z <- ( Y - U ) / sqrt( U * ( 1 - U ) )

    iter <- iter + 1

  }

  beta0 <- beta

  bdiff <- 1
  iter <- 0

  while (sum(abs(bdiff))>.000001){

    DWZ <- matrix(0, nrow = p)
    DWD <- matrix(0, ncol = p, nrow = p)
    DWZZWD <- matrix(0, ncol = p, nrow = p)
    SS <- matrix(0, ncol = q, nrow = q)
    US <- matrix(0, nrow = p, ncol = q)
    EE <- matrix(0, ncol = p, nrow = p)
    Ulist <- list()

    ii <- 0

    for (i in unique(cluster)){

      uniti <- unit[cluster==i]
      csizei <- length( table(uniti) )
      wi <- 1/csizei
      yi <- matrix( Y[cluster==i])
      xi <- matrix( X[cluster==i,] , byrow = FALSE, ncol = p)
      ipcwi <- matrix( ipcw[cluster==i])

      DWZj <- matrix(0, nrow = p)
      DWDj <- matrix(0, ncol = p, nrow = p)

      for(j in unique(uniti)){

        yij <- matrix( yi[uniti==j] )
        xij <- matrix( xi[uniti == j,], ncol = p)
        ipcwij <- matrix( ipcwi[uniti==j] )

        u <- exp( xij %*% beta ) / ( 1 + exp( xij %*% beta ) )
        dudb <- exp( xij %*% beta ) / ( 1 + exp( xij %*% beta ) ) ^ 2

        if(ncol(xij) > 1){
          D <- xij[,1] * dudb
          for (pp in 2:ncol(xij)){
            D <- cbind(D, xij[,pp] * dudb)
          }
        }else{
          D <- as.matrix( xij[,1] * dudb )
        }

        vv <- u * ( 1 - u )
        V <- matrix(ncol = length(vv), nrow = length(vv), 0)
        diag(V) <- vv
        R <- matrix(ncol = length(vv), nrow = length(vv), 0)
        diag(R) <- 1
        W <- V ^ (1/2) %*% R %*% V ^ (1/2)

        invW <- solve(W)

        Delta <- matrix(ncol = length(vv), nrow = length(vv), 0)
        diag(Delta) <- ipcwij


        DWZj <- DWZj + t(D) %*% invW %*% Delta %*% ( yij - u )
        DWDj <- DWDj + t(D) %*% invW %*% Delta %*% D

      }

      wDWZj <- wi * DWZj
      ii <- ii + 1
      Sj <- Slist[[ii]]
      SjSj <- Sj %*% t( Sj )
      UjSj <- wDWZj %*% t( Sj )


      DWZ <- DWZ + wDWZj
      DWD <- DWD + wi * DWDj
      DWZZWD <- DWZZWD + ( wDWZj ) %*% t( wDWZj )

      US <- US + UjSj
      SS <- SS + SjSj
      Ulist[[ii]] <- wDWZj

    }

    invDWD <- solve(DWD)
    invSS <- solve(SS)
    bdiff <- invDWD %*% DWZ
    beta <- beta + bdiff

    U <- exp(X %*% beta) / (1 + exp(X %*% beta))
    var <- U * (1 - U)
    Z <- (Y - U) / sqrt(var)

    iter <- iter + 1

    if(show.iter == "TRUE"){
      print(iter)
    }

  }

  for (i in 1:length(unique(cluster))){

    Ui <- Ulist[[i]]
    Si <- Slist[[i]]
    Ei <- Ui - US %*% invSS %*% Si
    EE <- EE + Ei %*% t( Ei )

  }

  ### Compute variance
  covbeta <- invDWD %*% EE %*% invDWD

  #covbeta <- invDWD %*% DWZZWD %*% invDWD
  robust.se <- sqrt(diag(covbeta))

  ### Wald test
  wald.chisq <- (beta/sqrt(diag(covbeta)))^2
  p.value <- 1 - pchisq(wald.chisq, df = 1)

  coef.names <- c("int", xvars)

  result <- list()
  result$call <- call
  result$coefficients <- beta
  result$coef.names <- coef.names
  result$robust.variance <- covbeta
  result$robust.se <- robust.se
  result$wald.chisq <- wald.chisq
  result$p.value <- p.value
  result$niter <- iter
  class(result) <- "ipccwgee"
  result

}



