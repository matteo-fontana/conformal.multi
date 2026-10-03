#' Checks for the prediction methods
#'
#' Contains all the check functions used internally in the package.
#' All arguments match those of \code{\link{conformal.multidim.split}} in `split.multi.R`.
#'
#' @param x Feature matrix of dimension n x p.
#' @param y Response matrix of dimension n x q.
#' @param x0 New feature points to evaluate, matrix of dimension n0 x p.
#' @param train.fun Function to train the model, producing an estimator of E(Y|X).
#'   Input arguments: x (features), y (responses).
#' @param predict.fun Function to predict responses at new feature values.
#'   Input arguments: out (output from train.fun), newx (new features).
#' @param alpha Miscoverage level for prediction intervals (coverage = 1 - alpha). Default 0.1.
#' @param seed Integer seed for random split. Default FALSE (no seed). Ignored if split is provided.
#' @param randomized Logical. Should the randomized approach be used? Default FALSE.
#' @param seed.tau Seed for the randomized version. Default FALSE.
#' @param training_size Proportion of data used for training vs calibration. Default 0.5.
#' @param score Nonconformity measure: "max", "l2", or "mahalanobis". Default "l2".
#' @param mad.train.fun Optional function to train on absolute residuals (E(R|X)), where
#'   R = |Y - m(X)| and m is the estimator from train.fun. Overrides s_type if provided.
#'   Input arguments: x (features), y (absolute residuals), out (previous output). Default NULL.
#' @param mad.predict.fun Optional function to predict mean absolute residuals at new features.
#'   Input arguments: out (output from mad.train.fun), newx (features). Default NULL.
#'
#' @noRd


check.split=function(x,y,x0,train.fun,
               predict.fun, alpha, seed, training_size, seed.tau, randomized,mad.train.fun=NULL, mad.predict.fun=NULL, score){


  check.null.data(y)
  check.null.data(x)
  check.null.data(x0)

  if ( is.data.frame(x)==FALSE && is.matrix(x)==FALSE){

    stop("x must be a matrix or a df of dimension n x p")}



  if ( is.data.frame(y)==FALSE && is.matrix(y)==FALSE){

    stop("y must be a matrix or a df of dimension n x q")}


  if ( is.data.frame(x0)==FALSE && is.matrix(x0)==FALSE){

    stop("x0 must be a matrix or a df of dimension n0 x p")}


 if (dim(x)[1] != dim(y)[1]) stop("x and y must have the same number of rows")

  if (dim(x)[2] != dim(x0)[2]) stop("x and x0 must have the same number of columns")


  if (is.null(train.fun) || !is.function(train.fun))
    stop("train.fun must be a function")


  if (is.null(predict.fun) || !is.function(predict.fun))
    stop("predict.fun must be a function")

  if (!is.null(mad.train.fun) && !is.function(mad.train.fun))
    stop("mad.train.fun must be a function")
  if (!is.null(mad.predict.fun) && !is.function(mad.predict.fun))
    stop("mad.predict.fun must be a function")
  if ((!is.null(mad.train.fun) && is.null(mad.predict.fun)) ||
      (is.null(mad.train.fun) && !is.null(mad.predict.fun)))
    stop("mad.train.fun and mad.predict.fun must both be provided")


  check.num.01(alpha)


  if (is.null(seed)==TRUE || (seed!=FALSE & is.numeric(seed)==FALSE))
    stop("Argument 'seed' must be either FALSE or an integer.")



  if (is.null(training_size)==TRUE || (training_size!=FALSE & is.numeric(training_size)==FALSE)) stop("Argument 'training_size' must be either FALSE or an integer.")


  check.num.01(training_size)


  if (is.null(randomized)==TRUE || randomized %in% c("TRUE","FALSE")==FALSE)
    stop("Argument 'randomized' must be either TRUE or FALSE")

  if(score == "scaled.max")
    stop("Argument 'score' cannot take value 'scaled.max' in the multivariate conformal split")

  check.score(score)

}

check.full=function(x,y,x0,train.fun,predict.fun, alpha, num.grid.pts.dim,grid.factor,score,mad.train.fun=NULL, mad.predict.fun=NULL){

  check.null.data(y)
  check.null.data(x)
  check.null.data(x0)

  if ( is.data.frame(x)==FALSE && is.matrix(x)==FALSE){

    stop("x must be a matrix or a df of dimension n x p")}



  if ( is.data.frame(y)==FALSE && is.matrix(y)==FALSE){

    stop("y must be a matrix or a df of dimension n x q")}

  if(ncol(y)>2){
    stop("Too many dimensions for the response. The full conformal method is too computationally extensive. Try with the split conformal method.")
  }


  if ( is.data.frame(x0)==FALSE && is.matrix(x0)==FALSE){

    stop("x0 must be a matrix or a df of dimension n0 x p")}


  if (dim(x)[1] != dim(y)[1]) stop("x and y must have the same number of rows")

  if (dim(x)[2] != dim(x0)[2]) stop("x and x0 must have the same number of columns")


  if (is.null(train.fun) || !is.function(train.fun))
    stop("train.fun must be a function")


  if (is.null(predict.fun) || !is.function(predict.fun))
    stop("predict.fun must be a function")

  if (!is.null(mad.train.fun) && !is.function(mad.train.fun))
    stop("mad.train.fun must be a function")
  if (!is.null(mad.predict.fun) && !is.function(mad.predict.fun))
    stop("mad.predict.fun must be a function")
  if ((!is.null(mad.train.fun) && is.null(mad.predict.fun)) ||
      (is.null(mad.train.fun) && !is.null(mad.predict.fun)))
    stop("mad.train.fun and mad.predict.fun must both be provided")


  check.num.01(alpha)

  if (length(num.grid.pts.dim) != 1 || !is.numeric(num.grid.pts.dim)
      || num.grid.pts.dim <= 1 || num.grid.pts.dim >= 1000
      || round(num.grid.pts.dim) != num.grid.pts.dim) {
    stop("num.grid.pts must be an integer between 1 and 1000")
  }
  check.pos.num(grid.factor)


    possible_score_functions=c('l2','mahalanobis','max')
    if (is.null(score) || score %in% possible_score_functions==FALSE) {
      stop(c("The 'score' argument is not correct. Please select one of the following:",paste(possible_score_functions,collapse=", "),"."))
    }
}





check.null.data=function(y){


  if (is.null(y)) stop("y must be either be a matrix or a df")
}



check.bool = function(b) {
  if (is.null(b) || length(b)!=1 || !is.logical(b))
    stop(paste(deparse(substitute(b)),"must be a Boolean"))
}

check.num = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a))
    stop(paste(deparse(substitute(a)),"must be a number"))
}

check.int = function(i) {
  if (is.null(i) || length(i)!= 1 || !is.numeric(i) || round(i) != i)
    stop(paste(deparse(substitute(i)),"must be an integer"))
}

check.pos.num = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a) || a<0)
    stop(paste(deparse(substitute(a)),"must be a positive number"))
}

check.pos.int = function(i) {
  if (is.null(i) || length(i)!= 1 || !is.numeric(i) || round(i) != i || i<1)
    stop(paste(deparse(substitute(i)),"must be a positive integer"))
}

check.num.01 = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a) || a<0 || a>1)
    stop(paste(deparse(substitute(a)),"must be a number between 0 and 1"))
}

check.tau=function(alpha,tau,l){

  if (alpha<tau/(l+1) & alpha>0)
    stop ("The prediction band obtained with such a small value of alpha is the entire space.
                                       If you are using the non randomized version of the algorithm, non-trivial prediction bands can be obtained by using an alpha greater or equal than 1/(l+1) and less than 1.
                                       If you are using the randomized version of the algorithm, non-trivial prediction bands can be obtained by using an alpha greater or equal than tau/(l+1) and less than (tau+l)/(l+1).")


  if (alpha>=(l+tau)/(l+1) || alpha<=0)
    stop("The alpha value is not admissible.
                                                   If you are using the non randomized version of the algorithm, non-trivial prediction bands can be obtained by using an alpha greater or equal than 1/(l+1) and less than 1.
                                                   If you are using the randomized version of the algorithm, non-trivial prediction bands can be obtained by using an alpha greater or equal than tau/(l+1) and less than (tau+l)/(l+1).")

}

check.s_regression=function(mat_residual,type){


    if( is.matrix(mat_residual)==FALSE & is.data.frame(mat_residual)==FALSE & (is.atomic(mat_residual)==FALSE || is.vector(mat_residual)==FALSE)) stop("vec_residual must be either a matrix, a dataframe or an atomic vector (naive case).")



  #check on 'type' argument
  possible_s_functions=c("identity","st-dev","alpha-max")
  if (is.null(type) || type %in% possible_s_functions==FALSE) {
    stop(c("The 'type' argument is not correct. Please select one of the following:",paste(possible_s_functions,collapse=", "),"."))
  }

}

check.score = function(score){

  #check on 'score' argument
  possible_scores=c("l2","mahalanobis","max","scaled.max")
  if (is.null(score) || score %in% possible_scores==FALSE) {
    stop(c("The 'score' argument is not correct. Please select one of the following:",paste(possible_scores,collapse=", "),"."))
  }
}
