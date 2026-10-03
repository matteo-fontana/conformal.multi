#' Mean of Multivariate Response
#'
#' This model can be used with conformal prediction functions.
#' It returns a training function and a prediction function.
#'
#' @return A list with two components:
#' \item{train.fun}{Function to train the model.}
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
#' @seealso \code{\link{conformal.multidim.split}}
#'
#' @export mean_multi


mean_multi = function() {

  # Training function
  train.fun = function(x,y,out=NULL) {

    m=colMeans(y)

    return(list(m=m))

    }

  # Prediction function
  predict.fun = function(out,newx) {

    temp=out$m
    n0=dim(newx)[1]
    len=length(temp)

    sol= t(matrix(rep(temp,n0),nrow=len))
    return(sol)
  }



  return(list(train.fun=train.fun, predict.fun=predict.fun))
}
