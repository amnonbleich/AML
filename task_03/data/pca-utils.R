#' Plot a scree plot with percent explained variance on the y-axis
#'
#' @param x a prcomp object
#' @param npcs the number of principal components to plot
#' @return
#' @examples
#'
#' library("ALL")
#' data("ALL")
#' ex <- exprs(ALL)
#' p <- prcomp(t(ex), center = TRUE, scale. = TRUE)
#' screeplot_percent(p, npcs = 20)
#'
#' @author Matt Huska, \email{huska@@molgen.mpg.de}
screeplot_percent <- function(x, npcs = min(10, length(x$sdev)), ...) {
  idx <- seq_len(npcs)
  sum_var <- sum(x$sdev ^ 2)
  vars <- 100 * (x$sdev[idx] ^ 2 / sum_var)
  cumvar <- cumsum(vars)

  barplot(vars, width = 0.9, space = 0.1, names.arg = idx, ylim = c(0, 100),
       xlab = "Principal Component", ylab = "Percent Variance",
       xaxp = c(1, npcs, npcs - 1), las = 1)
  lines(x = idx - 0.5, y = cumvar, type = "b", lty = 2)
  legend("topleft", legend = c("Proportion", "Cumulative"), lty = c(1, 2),
         pch = c(19, 1))
}

#' Plot a scatter plot matrix of serveral components from a PCA
#'
#' @param x a prcomp object
#' @param npcs the number of principal components to plot
#' @param col the color of each point in the plot
#' @param pch the plotting character for each point in the plot
#' @return
#' @examples
#'
#' library("ALL")
#' data("ALL")
#' ex <- exprs(ALL)
#' pdata <- pData(ALL)
#' p <- prcomp(t(ex), center = TRUE, scale. = TRUE)
#' splom_pca(p, npcs = 6, col = ifelse(pdata[,"relapse"], "blue", "red"))
#'
#' @author Matt Huska, \email{huska@@molgen.mpg.de}
splom_pca <- function(x, npcs = min(5, ncol(x$x)), col = rgb(0, 0, 0, 0.5),
                      pch = 19, ...) {
  if(!require(lattice))
    stop("The lattice library is required for this function. Aborting.")
  idx <- seq_len(npcs)
  splom(x$x[,idx], col=col, pch=pch, ...)
}
