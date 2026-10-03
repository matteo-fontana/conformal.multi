#' Plot Prediction Regions obtained from Full Conformal
#'
#' @param full It's the output of the multivariate full conformal prediction function
#' @return A list of ggplots (output[[i]] is the i-th observation prediction region).
#'   Each plot displays the trial values retained in the prediction region,
#'   coloured according to their conformal p-value, and the point prediction.
#'   For a univariate response, the p-values of the retained trial values are
#'   plotted against the trial values.
#' @details It exploits the package \code{\link[ggplot2]{ggplot2}}
#' to better visualize the results.
#'
#' @example inst/examples/ex.full.R
#'
#' @export plot_multidim_full
plot_multidim_full <- function(full) {

  # Get Data
  valid_points <- full$valid_points
  pred <- as.matrix(full$pred)
  n0 <- length(valid_points)

  if (ncol(pred) == 1) {
    # univariate response: p-values of the retained trial values
    plots <- lapply(seq_len(n0), function(k) {
      df <- valid_points[[k]]
      colnames(df) <- c("Var1", "pval")
      ggplot2::ggplot(data = df, ggplot2::aes(Var1, pval)) +
        ggplot2::geom_point(size = 0.8) +
        ggplot2::geom_vline(xintercept = pred[k, 1], linetype = 2) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("y") +
        ggplot2::ggtitle(paste("Test Observation", k))
    })
    return(plots)
  }

  plots <- lapply(seq_len(n0), function(k) {
    df <- valid_points[[k]]
    colnames(df) <- c("Var1", "Var2", "pval")
    df_pred <- data.frame(X1 = pred[k, 1], X2 = pred[k, 2])
    g_plot <- ggplot2::ggplot(data = df, ggplot2::aes(Var1, Var2)) +
      ggplot2::geom_raster(ggplot2::aes(fill = pval)) +
      ggplot2::scale_fill_distiller(palette = "RdPu", direction = 1) +
      ggplot2::theme_minimal() +
      ggplot2::xlab("y1") +
      ggplot2::ylab("y2") +
      ggplot2::ggtitle(paste("Test Observation", k)) +
      ggplot2::geom_point(
        data = df_pred,
        ggplot2::aes(x = X1, y = X2),
        shape = 8, size = 10
      )

    return(g_plot)
  })

  return(plots)
}

utils::globalVariables(c("Var1", "Var2", "X1", "X2", "pval"))
