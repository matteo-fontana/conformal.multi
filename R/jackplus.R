#' Multivariate Response Jackknife + Prediction Regions
#'
#' Compute prediction regions using multivariate Jackknife + inference.
#'
#' @param x The feature variables, a matrix n x p.
#' @param y The matrix of multivariate responses (dimension n x q)
#' @param x0 The new points to evaluate, a matrix of dimension n0 x p.
#' @param train.fun A function to perform model training, i.e., to produce an
#'   estimator of E(Y|X), the conditional expectation of the response variable
#'   Y given features X. Its input arguments should be x: matrix of features,
#'   and y: matrix of responses.
#' @param predict.fun A function to perform prediction for the (mean of the)
#'   responses at new feature values. Its input arguments should be out: output
#'   produced by train.fun, and newx: feature values at which we want to make
#'   predictions.
#' @param alpha Miscoverage level. Default for alpha is 0.1.
#' @param method The multivariate extension of the jackknife+: "depth"
#'   (default) or "max". See details.
#' @param scale Only for method = "max": vector of q positive constants that
#'   scale the components of the residuals. They must not depend on the data
#'   for the coverage guarantee to hold. Default NULL, i.e. 1 for every
#'   component.
#' @return A list with components x0, lo and up. lo and up are matrices of
#'   dimension n0 x q containing the lower and upper bounds of the prediction
#'   region for each test point.
#'
#' @details The work is an extension of the univariate approach to jackknife +
#'  inference to a multivariate context. Let R_i be the vector of absolute
#'  leave-one-out residuals of observation i and mu_{-i} the leave-one-out
#'  regression function.
#'
#'  With method = "depth", the 2n vectors mu_{-i}(x0) -/+ R_i are standardised
#'  component-wise and ranked according to the depth -max_j |z_j|; the
#'  ceiling(2n(1-alpha/2)) deepest vectors are retained and their axis-aligned
#'  bounding box is returned. This is a heuristic: the finite sample guarantee
#'  of the univariate jackknife+ is not established for it, and the coverage
#'  of the region can be below 1-alpha, in particular when q is large relative
#'  to n.
#'
#'  With method = "max", the jackknife+ of Barber et al. (2021) is applied to
#'  the scalar score r_i = max_j R_ij / scale_j: for each component j, the
#'  lower bound is the floor(alpha(n+1))-th smallest value of
#'  mu_{-i,j}(x0) - scale_j r_i and the upper bound is the
#'  ceiling((1-alpha)(n+1))-th smallest value of mu_{-i,j}(x0) + scale_j r_i
#'  (-Inf and Inf when the index is out of range). The box contains the
#'  jackknife+ prediction set based on the score max_j |y_j - mu_j(x)| /
#'  scale_j, so that, under exchangeability and if the regression algorithm
#'  treats the observations symmetrically, its coverage is at least 1-2alpha
#'  (Barber et al., 2021, Theorem 1). For q = 1 it is the jackknife+ interval.
#' @details This function is based on the package future.apply to
#'  perform parallelisation.
#'
#' @references Barber, Candes, Ramdas, Tibshirani (2021), "Predictive inference
#'   with the jackknife+", The Annals of Statistics 49(1), 486-507.
#'
#' @example inst/examples/ex.jackplus.R
#' @export conformal.multidim.jackplus

conformal.multidim.jackplus = function(x,y,x0, train.fun, predict.fun, alpha=0.1,
                                       method = c("depth","max"), scale = NULL) {

  method = match.arg(method)

  ## Check Data
  check.split(x=x,y=y,x0=x0,train.fun=train.fun,
              predict.fun=predict.fun, alpha=alpha, seed=FALSE,
              training_size=0.5, seed.tau=FALSE, randomized=FALSE,
              score="max")

  y=as.matrix(y)
  n=dim(x)[1]
  q=dim(y)[2]
  n0=nrow(x0)

  ### Parallel sessions
  oplan = future::plan(future::multisession)
  on.exit(future::plan(oplan), add = TRUE)

  ## Fit the n leave-one-out models; each worker returns the absolute LOO
  ## residual of the left-out observation and the predictions at x0
  loo = future.apply::future_lapply(1:n,function(jj){
    mod_jj = train.fun(x[-jj,,drop=FALSE],y[-jj,,drop=FALSE])
    list(res = abs(c(as.matrix(predict.fun(mod_jj, x[jj,,drop=FALSE]))) - y[jj,]),
         fit0 = matrix(as.matrix(predict.fun(mod_jj, x0)), nrow=n0))
  }, future.seed = TRUE)

  ## LOO residuals (n x q)
  Loo = matrix(t(vapply(loo, function(l) l$res, numeric(q))), nrow=n, ncol=q)

  ## Fitted values of the LOO models at the test points (n x q for each x0)
  fitted = lapply(1:n0, function(i)
    matrix(t(vapply(loo, function(l) l$fit0[i,], numeric(q))), nrow=n, ncol=q))


  if(method=="max"){

    ## Jackknife+ based on the scalar score max_j |R_ij| / scale_j
    if(is.null(scale)) scale = rep(1,q)
    if(length(scale)!=q || any(!is.finite(scale)) || any(scale<=0))
      stop("scale must be a vector of q positive numbers")
    r = apply(sweep(Loo,2,scale,"/"),1,max)
    k_lo = floor(alpha*(n+1))
    k_up = ceiling((1-alpha)*(n+1))

    lo = t(vapply(1:n0, function(i) vapply(1:q, function(j)
      if(k_lo<1) -Inf else sort(fitted[[i]][,j]-scale[j]*r)[k_lo], numeric(1)),
      numeric(q)))
    up = t(vapply(1:n0, function(i) vapply(1:q, function(j)
      if(k_up>n) Inf else sort(fitted[[i]][,j]+scale[j]*r)[k_up], numeric(1)),
      numeric(q)))
    lo = matrix(lo, nrow=n0, ncol=q)
    up = matrix(up, nrow=n0, ncol=q)

    return(list(lo=lo,up=up,x0=x0))
  }

  ## Depth-based extension: number of retained vectors
  a = min(2*n, ceiling(2*n*(1-alpha/2)))

  box = lapply(1:n0, function(i){
    joint = rbind(fitted[[i]]-Loo, fitted[[i]]+Loo)
    dep = depth.max(joint)
    kept = joint[order(dep,decreasing = TRUE)[1:a],,drop=FALSE]
    rbind(apply(kept,2,min),apply(kept,2,max))
  })

  lo = matrix(t(vapply(box, function(b) b[1,], numeric(q))), nrow=n0, ncol=q)
  up = matrix(t(vapply(box, function(b) b[2,], numeric(q))), nrow=n0, ncol=q)

  return(list(lo=lo,up=up,x0=x0))

}


#' Depth of a set of vectors
#'
#' Standardises the columns of inp and returns, for every row, minus the
#' maximum absolute standardised component. Larger values correspond to
#' deeper (more central) vectors. Missing values are ignored in the
#' standardisation (as in base::scale) and rows containing them get depth -Inf.
#' Columns with numerically zero variance do not contribute.
#'
#' @param inp A matrix, each row being a vector.
#' @return A vector of length nrow(inp).
#' @noRd

depth.max=function(inp){
  inp=as.matrix(inp)
  mus=colMeans(inp,na.rm=TRUE)
  sds=apply(inp,2,stats::sd,na.rm=TRUE)
  # columns with (numerically) zero variance do not contribute
  const=!is.finite(sds) | sds==0 | sds <= sqrt(.Machine$double.eps)*abs(mus)
  sds[const]=1
  z=sweep(sweep(inp,2,mus),2,sds,"/")
  z[,const]=0
  d=-apply(abs(z),1,max)
  d[is.na(d)]=-Inf
  return(d)
}
