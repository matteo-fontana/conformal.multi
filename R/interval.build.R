#' Build Multi-Split Interval for Each X0
#'
#' Constructs a multi-split prediction interval for each new point X0 by combining
#' results from multiple split conformal predictions using the multi-split algorithm.
#'
#' @param yyy Column vector containing B lower bounds and B upper bounds.
#' @param B Number of replications used in the multi-split procedure.
#' @param tr Truncation threshold for the algorithm.
#'
#' @return A vector with the lower and upper end of the convex hull of the set
#'   of values contained in more than tr of the B intervals. If the set is
#'   empty, c(NA, NA) is returned.
#'
#' @importFrom utils flush.console
#' @noRd




interval.build=function(yyy,B,tr){



  h=rep(1:0,each=B)

  o = order(yyy,2-h)

  ys <- yyy[o]
  hs <- h[o]

  count <- 0
  lo <- up <- NA


  for (j in 1:(2*B) ){
    if ( hs[j]==1 ) {
      count <- count + 1

      # first value at which the count exceeds the threshold
      if ( count > tr && (count - 1) <= tr && is.na(lo)) {
        lo <- ys[j]
      }

    }

    else {
      # last value at which the count drops below the threshold
      if ( count > tr && (count - 1) <= tr) {
        up <- ys[j]
      }

      count <- count - 1
    }
  }




  return(c(lo,up))
}
