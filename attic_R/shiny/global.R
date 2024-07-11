# library(smallRNA)
load("object.RData")

plPCA <- function(x, normalize.genes = F, plot.components = c(1, 2),
                  plot.labels = T, tsplit = NULL, doclusplot = 0, points.color = NULL,
                  labels.color = "black", ...) {
  if (!(length(plot.components) %in% c(2, 3))) {
    stop("Plot components should be of size 2 or 3")
  }
  if (is.null(points.color)) {
    points.color <- ifelse(plot.labels, "white", "black")
  }
  x <- t(x[which(!(apply(x, 1, var) == 0)), ])
  if (normalize.genes) {
    x <- scale(x)
  }
  pca <- prcomp(x)
  xlab <- paste("PC ", plot.components[1], " (", round(pca$sdev[plot.components[1]] / sum(pca$sdev) *
    100, 0), "%)", sep = "")
  ylab <- paste("PC ", plot.components[2], " (", round(pca$sdev[plot.components[2]] / sum(pca$sdev) *
    100, 0), "%)", sep = "")
  if (length(plot.components) == 2) {
    if (doclusplot > 1) {
      library(cluster)
      clusplot(pam(-1 * pca$x[, plot.components], doclusplot),
        shade = TRUE, sub = "", xlab = xlab, ylab = ylab,
        ...
      )
    } else {
      plot(pca$x[, plot.components[1]], pca$x[, plot.components[2]],
        xlab = xlab, ylab = ylab, col = points.color,
        ...
      )
    }
    if (plot.labels) {
      text(pca$x[, plot.components[1]], pca$x[, plot.components[2]],
        labels = row.names(pca$x), cex = 0.8, font = 2,
        col = labels.color
      )
    }
    if (!is.null(tsplit)) {
      if (length(tsplit) == nrow(x) & length(unique(tsplit)) ==
        2) {
        divline <- supervisedPlotDivision(
          pca$x[, plot.components],
          tsplit
        )
        abline(
          a = divline$intercept, b = divline$slope,
          lty = "dashed", lwd = 2, col = "grey"
        )
      } else {
        warning("tsplit is not a vector of length=ncol(x) factorizable to two levels, and will therefore be ignored.")
      }
    }
  }
}

byheatmap <- function(x, ...) {
  x <- x[which(apply(x, 1, FUN = function(y) {
    !all(is.na(y))
  })), which(apply(x, 2, FUN = function(y) {
    !all(is.na(y))
  }))]
  pheatmap::pheatmap(x, color = colorRampPalette(c("blue", "black", "yellow"))(29), border_color = NA, ...)
}
