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
#' @param training_size Vector of length B with the proportion of data used for
#'   training in each split. Default NULL, i.e. 0.5 for every split.
#' @param s_type Type of modulation function: "identity", "st-dev", or "alpha-max". Default "st-dev".
#' @param B Number of repetitions. Default 100.
#' @param lambda Smoothing parameter. Default 0.
#' @param tau Parameter of the multi split procedure. Each split conformal
#'   region is computed at miscoverage level alpha*(1-tau+lambda/B). With
#'   aggregation = "componentwise", a value is retained in the final interval of
#'   a component when it is contained in more than tau*B of the B split
#'   intervals of that component:
#'   \describe{
#'     \item{tau = 1 - 1/B}{Bonferroni intersection method.}
#'     \item{tau = 0}{Unadjusted intersection.}
#'   }
#'   With aggregation = "depth", the floor(2*tau*B) deepest of the 2B vertices
#'   are retained. Default NULL, i.e. 1 - (B + 1)/(2 * B) for "componentwise"
#'   and 0.1 for "depth".
#' @param seed_beta Seed for the randomized version. Default FALSE.
#' @param score Nonconformity measure to use for the split conformal function.
#' @param aggregation How the B split conformal regions are combined, either
#'   "componentwise" (default) or "depth". See details.
#' @param mad.train.fun,mad.predict.fun Optional functions for local scaling of
#'   the residuals, passed to \code{\link{conformal.multidim.split}}. Default NULL.
#'
#' @return A list with components x0, lo, and up. lo and up are matrices of dimension n0 x q.
#'
#' @details This function extends the univariate Multi Split conformal approach to the multivariate case.
#'   The split conformal method is run B times, each run returning an
#'   axis-aligned box with lower vertex lo^b and upper vertex up^b.
#'   With aggregation = "componentwise" the univariate multi split procedure of
#'   Solari and Djordjilovic (2022) is applied to each component separately.
#'   The resulting box contains every y that belongs to more than tau*B of the
#'   B regions, so that its coverage is at least 1-alpha. If, for a component,
#'   the set of retained values is not an interval, its convex hull is returned.
#'   With aggregation = "depth" the 2B vertices lo^b, up^b are ranked according
#'   to the depth -max_j |z_j|, where z is the vertex standardised
#'   component-wise, the floor(2*tau*B) deepest vertices are retained and their
#'   bounding box is returned. The latter is a heuristic: no finite sample
#'   coverage guarantee is available for it.
#'   Parallelization is performed via the \code{\link[future.apply]{future_sapply}} function.
#'
#' @references Solari, Djordjilovic (2022), "Multi Split Conformal Prediction",
#'   Statistics and Probability Letters 184, 109395 (baseline for univariate case).
#'
#' @example inst/examples/ex.msplit.R
#' @export conformal.multidim.msplit




conformal.multidim.msplit = function(x,y, x0, train.fun, predict.fun, alpha=0.1,
                       split=NULL, seed=FALSE, randomized=FALSE,seed_beta=FALSE,
                       verbose=FALSE, training_size=NULL,score = "max",
                       s_type = "st-dev",B=100,lambda=0,
                       tau = NULL, aggregation = c("componentwise","depth"),
                       mad.train.fun = NULL, mad.predict.fun = NULL) {

  aggregation = match.arg(aggregation)
  check.pos.int(B)
  check.pos.num(lambda)

  if(is.null(tau))
    tau = switch(aggregation, "componentwise" = 1-(B+1)/(2*B), "depth" = 0.1)
  check.num.01(tau)

  if(is.null(training_size) || length(training_size)!=B)
    training_size=rep(0.5,B)

  if (!is.null(seed)) set.seed(seed)

  x0=as.matrix(x0)
  n0=nrow(x0)
  p=ncol(x0)
  q=ncol(y)
  n=nrow(x)
  full=q*n0


  oplan = future::plan(future::multisession)
  on.exit(future::plan(oplan), add = TRUE)
  options(future.rng.onMisuse="ignore")


  lo_up <- t(future.apply::future_sapply(1:B, function(bbb) {


    out<-conformal.multidim.split(x=x,y=y, x0=x0, train.fun=train.fun,
                                  predict.fun=predict.fun,
                                  alpha=alpha*(1-tau) + (alpha*lambda)/B,
                                  split=split, seed=seed+bbb,
                                  randomized=randomized, seed_tau=seed_beta,
                                  verbose=verbose, training_size=training_size[bbb],
                                  score=score, s_type=s_type,
                                  mad.train.fun=mad.train.fun,
                                  mad.predict.fun=mad.predict.fun)

    return(cbind(t(out$lo),t(out$up)))

      }))

  # lo_up is a B x (2 q n0) matrix: the first q*n0 columns contain the lower
  # vertices, the remaining ones the upper vertices (component index running
  # fastest)

  if(aggregation=="componentwise"){

    Y = rbind(lo_up[,1:full,drop=FALSE],lo_up[,-(1:full),drop=FALSE])
    tr <- tau*B + .001
    finalInt = t(future.apply::future_sapply(1:full, function(kk) interval.build(Y[,kk],B,tr)))

    lo<-matrix(finalInt[,1], nrow = n0, ncol = q, byrow = TRUE)
    up<-matrix(finalInt[,2], nrow = n0, ncol = q, byrow = TRUE)

  } else {

    tr = 2*tau*B + .001
    a = max(1,floor(tr))

    box = lapply(1:n0, function(i){
      cols = ((i-1)*q+1):(i*q)
      vertices = rbind(lo_up[,cols,drop=FALSE],lo_up[,full+cols,drop=FALSE])
      dep = depth.max(vertices)
      kept = vertices[order(dep,decreasing = TRUE)[1:a],,drop=FALSE]
      rbind(apply(kept,2,min),apply(kept,2,max))
    })

    lo = t(vapply(box, function(b) b[1,], numeric(q)))
    up = t(vapply(box, function(b) b[2,], numeric(q)))
    if(q==1){
      lo = matrix(lo, ncol=1)
      up = matrix(up, ncol=1)
    }

  }

  return(list(lo=lo,up=up,x0=x0))
}
