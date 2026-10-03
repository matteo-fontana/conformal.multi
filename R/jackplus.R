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
#' @param alpha Miscoverage level for the prediction intervals, i.e., intervals
#'   with coverage 1-alpha are formed. Default for alpha is 0.1.
#' @return A list with components x0, lo and up. lo and up are matrices of
#'   dimension n0 x q containing the lower and upper bounds of the prediction
#'   region for each test point.
#'
#' @details The work is an extension of the univariate approach to jackknife +
#'  inference to a multivariate context, exploiting the concept of depth measures.
#'  For each test point x0, the n leave-one-out models give the 2n vectors
#'  mu_{-i}(x0) -/+ R_i, where R_i is the vector of absolute leave-one-out
#'  residuals of observation i. The vectors are standardised component-wise and
#'  ranked according to the depth -max_j |z_j|; the ceiling(2n(1-alpha/2))
#'  deepest vectors are retained and their axis-aligned bounding box is
#'  returned. The finite sample guarantee of the univariate jackknife+ is not
#'  established for this extension.
#' @details This function is based on the package future.apply to
#'  perform parallelisation.
#'
#' @references Barber, Candes, Ramdas, Tibshirani (2021), "Predictive inference
#'   with the jackknife+", The Annals of Statistics 49(1), 486-507.
#'
#' @example inst/examples/ex.jackplus.R
#' @export conformal.multidim.jackplus

conformal.multidim.jackplus = function(x,y,x0, train.fun, predict.fun, alpha=0.1) {

  ## Check Data
  check.split(x=x,y=y,x0=x0,train.fun=train.fun,
              predict.fun=predict.fun, alpha=alpha, seed=FALSE,
              training_size=0.5, seed.tau=FALSE, randomized=FALSE,
              score="max")

  x=as.matrix(x)
  y=as.matrix(y)
  x0=as.matrix(x0)
  n=dim(x)[1]
  p=dim(x)[2]
  q=dim(y)[2]
  n0=nrow(x0)

  ### Parallel sessions
  oplan = future::plan(future::multisession)
  on.exit(future::plan(oplan), add = TRUE)
  options(future.rng.onMisuse="ignore")

  ## Compute models without each observation
  updated_models = future.apply::future_lapply(1:n,function(jj){
    mod_jj = train.fun(x[-jj,,drop=FALSE],y[-jj,,drop=FALSE])
    return(mod_jj)})

  ## Compute LOO residuals (n x q)
  Loo = t(vapply(1:n, function(jj)
    abs(c(predict.fun(updated_models[[jj]], x[jj,,drop=FALSE])) - y[jj,]),
    numeric(q)))
  if(q==1) Loo = matrix(Loo, ncol=1)


  ## Fitted values of the LOO models at the test points (n x q for each x0)
  fitted = lapply(1:n0, function(i)
    t(vapply(1:n, function(k)
      c(predict.fun(updated_models[[k]], x0[i,,drop=FALSE])), numeric(q))))


  ## Number of retained vectors
  a = min(2*n, ceiling(2*n*(1-alpha/2)))

  box = lapply(1:n0, function(i){
    fi = matrix(fitted[[i]], ncol=q)
    joint = rbind(fi-Loo, fi+Loo)
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
#'
#' @param inp A matrix, each row being a vector.
#' @return A vector of length nrow(inp).
#' @noRd

depth.max=function(inp){
  inp=as.matrix(inp)
  sds=apply(inp,2,stats::sd,na.rm=TRUE)
  sds[!is.finite(sds) | sds==0]=1
  z=sweep(sweep(inp,2,colMeans(inp,na.rm=TRUE)),2,sds,"/")
  d=-apply(abs(z),1,max)
  d[is.na(d)]=-Inf
  return(d)
}
