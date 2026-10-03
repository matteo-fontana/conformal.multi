#' Split conformal prediction intervals with Multivariate Response
#'
#' Compute prediction intervals using split conformal inference with multivariate
#' response.
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
#' @param split Indices that define the data-split to be used (i.e., the indices
#'   define the first half of the data-split, on which the model is trained).
#'   Default is NULL, in which case the split is chosen randomly.
#' @param seed Integer to be passed to set.seed before defining the random
#'   data-split to be used. Default is FALSE, which effectively sets no seed.
#'   If both split and seed are passed, the former takes priority and the latter
#'   is ignored.
#' @param randomized Should the randomized approach be used? Default is FALSE.
#' @param seed_tau The seed for the randomized version. Default is FALSE.
#' @param verbose Should intermediate progress be printed out? Default is FALSE.
#' @param training_size Split proportion between training and calibration set.
#' Default is 0.5.
#' @param score The non-conformity measure. It can either be "max", "l2", "mahalanobis".
#' The default is "l2".
#' @param s_type The type of modulation function.
#'  Currently we have 3 options: "identity","st-dev","alpha-max". Default is "st-dev"
#' @param mad.train.fun A function to perform training on the absolute residuals
#'   i.e., to produce an estimator of E(R|X) where R is the absolute residual
#'   R = |Y - m(X)|, and m denotes the estimator produced by train.fun.
#'   This is used to scale the conformal score, to produce a prediction interval
#'   with varying local width. The input arguments to mad.train.fun should be
#'   x: matrix of features, y: matrix of absolute residuals, and out: the output
#'   produced by a previous call to mad.train.fun, at the \emph{same} features
#'   x. The function mad.train.fun may (optionally) leverage this returned
#'   output for efficiency purposes. See details below. The default for
#'   mad.train.fun is NULL, which means that no training is done on the absolute
#'   residuals, and the usual (unscaled) conformal score is used. Note that if
#'   mad.train.fun is non-NULL, then so must be mad.predict.fun (next).
#' @param mad.predict.fun A function to perform prediction for the (mean of the)
#'   absolute residuals at new feature values. Its input arguments should be
#'   out: output produced by mad.train.fun, and newx: feature values at which we
#'   want to make predictions. The default for mad.predict.fun is NULL, which
#'   means that no local scaling is done for the conformal score, i.e., the
#'   usual (unscaled) conformal score is used.
#'
#' @return A list with the following components: x0,pred,k_s,s_type,s,alpha,randomized,tau,
#'   average_width,lo,up. In particular pred, lo, up are the matrices of
#'   dimension n0 x q, k_s is a scalar, s_type is a string, s is a vector of length q,
#'   alpha is a scalar between 0 and 1, randomized is a logical value,
#'   tau is a scalar between 0 and 1,and average_width is a positive scalar.
#'
#' @details The residuals of the calibration set are divided component-wise by
#'   the modulation function s (computed on the training set) and summarised by
#'   the chosen score: the max norm ("max"), the Euclidean norm ("l2") or the
#'   Mahalanobis norm with respect to the covariance matrix of the scaled
#'   training residuals ("mahalanobis"). k_s is the empirical quantile of the
#'   calibration scores. The prediction region is the set of y whose score does
#'   not exceed k_s; lo and up are the bounds of its smallest axis-aligned
#'   bounding box, i.e. pred -/+ k_s*s for "max" and "l2", and
#'   pred -/+ k_s*s*sqrt(diag(Sigma)) for "mahalanobis".
#' @details If the two mad functions are provided they take precedence over the s_type parameter,
#' and they force a local scoring via the mad function predicted values.
#'
#' @importFrom stats mad mahalanobis
#'
#' @seealso \code{\link{conformal.multidim.full}}
#'
#' @references The s_regression and the "max" score are taken from "Conformal Prediction Bands
#' for Multivariate Functional Data" by Diquigiovanni, Fontana, Vantini (2022).
#'
#' @example inst/examples/ex.split.R
#' @export conformal.multidim.split



conformal.multidim.split = function(x,y, x0, train.fun, predict.fun, alpha=0.1,
                                 split=NULL, seed=FALSE, randomized=FALSE,seed_tau=FALSE,
                                 verbose=FALSE, training_size=0.5,
                                 score="l2",s_type="st-dev", mad.train.fun = NULL,
                                 mad.predict.fun = NULL) {


  ## Check Data and Splits
  check.split(x=x,y=y,x0=x0,train.fun=train.fun,
              predict.fun=predict.fun, alpha=alpha, seed=seed, training_size
              =training_size, seed.tau=seed_tau, randomized=randomized,
              mad.train.fun = mad.train.fun, mad.predict.fun = mad.predict.fun, score = score)

  x=as.matrix(x)
  y=as.matrix(y)
  x0=as.matrix(x0)
  n=dim(x)[1]
  p=dim(x)[2]
  q=dim(y)[2]
  n0=nrow(x0)
  flag = FALSE #in general mad funs are not present

  # If mad funs exist, they take precedence over the modulation function
  if(is.function(mad.train.fun) && is.function(mad.predict.fun)){
    s_type = "identity"
    flag = TRUE
  }


  if (verbose == TRUE) txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }


  if(is.null(split)){

    if(ceiling(n*training_size) !=n )
      m=ceiling(n*training_size)
    else
      m=ceiling(n*training_size)-1

    if(seed!=FALSE){set.seed(seed)}

    training=sample(1:n,m)

  }

  else{
    training=split
  }

  calibration=setdiff(1:n,training)
  l=length(calibration)

  if(randomized==FALSE) {tau=1} else{
    if(seed_tau!=FALSE){set.seed(seed_tau)}
    tau=stats::runif(n=1,min=0,max=1)
  }

  check.tau(alpha=alpha,tau=tau,l=l)


  ###### TRAINING & RESIDUALS COMPUTATION
  if (verbose) {
    cat(sprintf("%sComputing models on first part ...\n",txt))
  }


  out = train.fun(x[training,,drop=FALSE],y[training,,drop=FALSE])
  fit = matrix(predict.fun(out,x),nrow=n)
  pred = matrix(predict.fun(out,x0),nrow=n0)

  if (verbose) {
    cat(sprintf("%sComputing residuals and quantiles on second part ...\n",txt))
  }


  res = y - fit
  s=computing_s_regression(mat_residual=res[training,,drop=FALSE],type=s_type,
                                 alpha=alpha,tau=tau)

  # Scale of each residual: the modulation function s, or the local scale
  # predicted by the mad functions (trained on the training set only)
  if(flag){ # with mad
    mad.out = mad.train.fun(x[training,,drop=FALSE],abs(res[training,,drop=FALSE]))
    scale.train = matrix(mad.predict.fun(mad.out,x[training,,drop=FALSE]),ncol=q)
    scale.cal = matrix(mad.predict.fun(mad.out,x[calibration,,drop=FALSE]),ncol=q)
    scale.x0 = matrix(mad.predict.fun(mad.out,x0),ncol=q)
  } else {
    scale.train = matrix(s,nrow=length(training),ncol=q,byrow=TRUE)
    scale.cal = matrix(s,nrow=l,ncol=q,byrow=TRUE)
    scale.x0 = matrix(s,nrow=n0,ncol=q,byrow=TRUE)
  }

  resc = res[calibration,,drop=FALSE] / scale.cal

  # The covariance matrix of the Mahalanobis score is estimated on the
  # training set, so that the score function does not depend on the
  # calibration data
  if(score=="mahalanobis"){
    sigma = stats::cov(res[training,,drop=FALSE] / scale.train)
    rho = sqrt(mahalanobis(resc,rep(0,q),sigma,tol=1e-12))
  } else if(score=="max"){
    rho = apply(abs(resc),1,max)
  } else { # "l2"
    rho = sqrt(rowSums(resc^2))
  }


  k_s=sort(rho,decreasing=FALSE)[ceiling(l+tau-(l+1)*alpha)]

  # Half-widths of the bounding box of {y : score(y) <= k_s}
  if(score=="mahalanobis")
    band = k_s * scale.x0 * matrix(sqrt(diag(sigma)),nrow=n0,ncol=q,byrow=TRUE)
  else
    band = k_s * scale.x0

  up=pred+band
  lo=pred-band
  average_width = mean(up-lo)


  return(structure(.Data=list(x0,pred,k_s,s_type,s,alpha,randomized,tau,
                              average_width,
                              lo,up),
                   names=c("x0","pred","k_s","s_type","s","alpha","randomized","tau"
                           ,"average_width", "lo", "up")))

}

utils::globalVariables(c("pval"))
