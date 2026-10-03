#' Plot Prediction Regions for a Multivariate Response
#'
#' Generate plots for the prediction regions produced by the multivariate
#' conformal prediction functions.
#'
#' @param split Output of \code{\link{conformal.multidim.split}},
#'   \code{\link{conformal.multidim.msplit}},
#'   \code{\link{conformal.multidim.jackplus}} or
#'   \code{\link{conformal.multidim.full}}.
#' @param same.scale Logical. Should all plots use the same y-axis scale? Default is FALSE.
#'
#' @return For the output of the split, multi split and jackknife+ functions,
#'   a list (of length p) of lists (of length q) of \code{ggplot} objects: the
#'   plot in position [[i]][[j]] displays the prediction interval of the j-th
#'   component of the response against the i-th covariate, at each test point.
#'   The plots are also drawn on a grid. For the output of
#'   \code{\link{conformal.multidim.full}}, the output of
#'   \code{\link{plot_multidim_full}}.
#'
#' @details This function uses the \code{\link[ggplot2]{ggplot2}} and
#'   \code{\link[gridExtra]{gridExtra}} packages for visualization. When the
#'   input does not contain point predictions (multi split and jackknife+), the
#'   midpoint of each interval is displayed.
#'
#' @example inst/examples/ex.split.R
#' @export plot_multidim



plot_multidim=function(split, same.scale = FALSE){

  if(!is.null(split$valid_points))
    return(plot_multidim_full(split))

  #Get Data
  x0 = as.matrix(split$x0)
  lo = as.matrix(split$lo)
  up = as.matrix(split$up)
  pred = split$pred
  if(is.null(pred))
    pred = (lo+up)/2
  pred = as.matrix(pred)

  # Find bounds for the plots

  if(same.scale){

    y_up = max(up) +0.01 * sd(up)
    y_lo = min(lo) -0.01 * sd(lo)

  }

  # Define dimensions
  p<-ncol(x0)
  q<-ncol(lo)



  gl <- lapply(1:p,function(ii) lapply(1:q, function(jj){

    df=data.frame(xd=x0[,ii],yg=pred[,jj],y_min=lo[,jj], y_max=up[,jj])

    ggg<- ggplot2::ggplot(df, ggplot2::aes(x=xd,y = yg)) + ggplot2::geom_pointrange(ggplot2::aes(ymin = y_min, ymax = y_max), color ="red") + ggplot2::xlab(paste("x ",ii)) + ggplot2::ylab(paste("y ",jj))


    if(same.scale)
      ggg = ggg + ggplot2::ylim(y_lo,y_up)

    return(ggg)

  }))



  glist <- do.call(c, gl)
  do.call(gridExtra::"grid.arrange", c(glist, ncol=q,top="Prediction Intervals"))
  return(gl)

}

utils::globalVariables(c( "xd", "y_max", "y_min", "yg"))
