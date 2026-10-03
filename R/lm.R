#' Linear Modeling of Multivariate Response
#'
#' This model can be used with conformal prediction functions.
#' It returns a training function and a prediction function.
#'
#' @return A list with two components:
#' \item{train.fun}{Function to train the model. Fits a separate linear model for each dimension of the response.}
#' \item{predict.fun}{Function to make predictions on new data.}
#'
#' @details
#' The training function takes as input:
#' \describe{
#'   \item{x}{Feature matrix of dimension n x p.}
#'   \item{y}{Response matrix of dimension n x q.}
#' }
#'
#' The prediction function takes as input:
#' \describe{
#'   \item{out}{Output of a previous call to \code{train.fun}.}
#'   \item{newx}{New feature matrix to evaluate, dimension n0 x p.}
#' }
#'
#' @importFrom stats lm
#' @seealso \code{\link{conformal.multidim.split}}
#' @export lm_multi




lm_multi = function() {

  # Training function
  train.fun = function(x,y,out=NULL) {

    q=dim(y)[2]
    p=dim(x)[2]

    df1=data.frame(cbind(x,y[,1]))

    coeff=vapply(1:q, function(i) lm(formula = y[,i] ~  x)$coefficients,numeric(p+1))
    # dim (p+1) x q

    return(list(coeff=coeff))

  }

  # Prediction function
  predict.fun = function(out,newx) {

    c=out$coeff
    n0=dim(newx)[1]
    newxx=cbind(rep(1,n0),newx)
    sol=newxx%*%c


    return(sol)
  }



  return(list(train.fun=train.fun, predict.fun=predict.fun))
}
