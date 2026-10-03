#' Multi Split conformal prediction intervals with Multivariate Response
#'
#' Compute prediction intervals using Multi Split conformal inference for a
#' multivariate response.
#'
#' @param x Feature matrix of dimension n x p.
#' @param y Response matrix of dimension n x q.
#' @param x0 New points to evaluate, matrix of dimension n0 x p.
#' @param train.fun Function to perform model training, producing an estimator of E(Y|X).
#'   Input arguments: x (features), y (responses).
#' @param predict.fun Function to predict responses at new feature values.
#'   Input arguments: out (output from train.fun), newx (new features).
#' @param alpha Miscoverage level for prediction intervals. Default 0.1.
#' @param split Indices defining the training split. Default NULL (random split).
#' @param seed Integer seed for random split. Ignored if split is provided. Default FALSE.
#' @param randomized Logical, whether to use the randomized approach. Default FALSE.
#' @param verbose Logical, print progress? Default FALSE.
#' @param training_size Proportion of data used for training. Default 0.5.
#' @param s_type Type of modulation function: "identity", "st-dev", or "alpha-max". Default "st-dev".
#' @param B Number of repetitions. Default 100.
#' @param lambda Smoothing parameter. Default 0.
#' @param tau Smoothing parameter for intersection method:
#'   \describe{
#'     \item{tau = 1 - 1/B}{Bonferroni intersection method.}
#'     \item{tau = 0}{Unadjusted intersection.}
#'   }
#'   Default 1 - (B + 1)/(2 * B).
#' @param seed_beta Seed for the randomized version. Default FALSE.
#' @param score Nonconformity measure to use for the split conformal function.
#'
#' @return A list with components x0, lo, and up. lo and up are matrices of dimension n0 x q.
#'
#' @details This function extends the univariate Multi Split conformal approach to the multivariate case.
#'   Parallelization is performed via the \code{\link[future.apply]{future_sapply}} function.
#'
#' @references Solari, Djordjilovic (2021), "Multi Split Conformal Prediction" (baseline for univariate case)
#'
#' @example inst/examples/ex.msplit.R
#' @export conformal.multidim.msplit




conformal.multidim.msplit = function(x,y, x0, train.fun, predict.fun, alpha=0.1,
                       split=NULL, seed=FALSE, randomized=FALSE,seed_beta=FALSE,
                       verbose=FALSE, training_size=NULL,score = "max",
                       s_type = "st-dev",B=100,lambda=0,
                       tau = 1-(B+1)/(2*B)) {



  if(is.null(training_size) || length(training_size)!=B)
    training_size=rep(0.5,B)

  if (!is.null(seed)) set.seed(seed)

  n0=nrow(x0)
  p=ncol(x0)
  q=ncol(y)
  n=nrow(x)
  full=q*n0
  loB<-upB<-matrix(0,nrow=B,ncol=full)



  future::plan(future::multisession)
  options(future.rng.onMisuse="ignore")

  lo_up <- t(future.apply::future_sapply(1:B, function(bbb) {


    out<-conformal.multidim.split(x,y, x0, train.fun, predict.fun,
                               alpha*(1-tau) + (alpha*lambda)/B,
                             split, seed+bbb, randomized,seed_beta,
                             verbose, training_size[bbb] ,score, s_type)

    return(cbind(t(out$lo),t(out$up)))

      }))


  Y = rbind(lo_up[,1:full],lo_up[,-(1:full)])
  tr <- tau*B + .001
  finalInt = t(future.apply::future_sapply(1:full, function(kk) interval.build(Y[,kk],B,tr)))

  lo<-matrix(finalInt[,1], nrow = n0, ncol = q, byrow = TRUE)
  up<-matrix(finalInt[,2], nrow = n0, ncol = q, byrow = TRUE)


  ## To avoid CRAN check errors
  ## R CMD check: make sure any open connections are closed afterward
  future::plan(future::sequential)



  return(list(lo=lo,up=up,x0=x0))
}


