#' plotFeatureCoverage
#'
#' Plots the coverage across a transcript/feature
#'
#' @param o An object of class shortRNAexp
#' @param feature The name of the feature
#' @param exact Logical; whether the feature name should be interpretated as exact (default TRUE), otherwise as a regular expression.
#' @param includeAmbiguous Logical; whether to include ambiguous sequences (default TRUE; plotted a different color)
#' @param logCounts Logical; whether to log-transform read counts for plotting (defualt FALSE) color)
#'
#' @export
plotFeatureCoverage <- function(o, feature, exact = TRUE, includeAmbiguous = TRUE, logCounts = FALSE, main = "") {
  us <- getFeatureSeqs(o, feature, exact, allowAmbiguous = F)
  as <- setdiff(getFeatureSeqs(o, feature, exact, allowAmbiguous = T), us)
  u <- rowSums(o@seqcounts[us, , drop = F])
  a <- rowSums(o@seqcounts[as, , drop = F])
  if (logCounts) {
    u <- log(u + 1)
    a <- log(a + 1)
  }
  pmax <- max(o@sources[us, "posInFeature"] + o@sources[us, "length"] - 1)
  if (includeAmbiguous) {
    plot(u + a, type = "l", lwd = 2, pch = 20, xlab = "Position across transcript (nt)", ylab = ifelse(logCounts, "log(total read count +1)",
      "total read count"
    ), main = main)
    polygon(c(1, 1:length(u), length(u)), c(0, u + a, 0), col = "red")
    polygon(c(1, 1:length(u), length(u)), c(0, u, 0), col = "blue")
    legend("topright", fill = c("blue", "red"), legend = c("Uniquely assigned", "Ambiguous"))
  } else {
    plot(u, type = "l", lwd = 2, pch = 20, xlab = "Position across transcript (nt)", ylab = ifelse(logCounts, "log(total read count +1)",
      "total read count"
    ), main = main)
    polygon(c(1, 1:length(u), length(u)), c(0, u, 0), col = "blue")
  }
}



#' plotMAs
#'
#' Plots MAs relative to the mean
#'
#' @param o An object of class shortRNAexp
#' @param normalized Logical; whether normalized counts should be used (default TRUE)
#' @param type Types of sequences to be fetched (e.g. miRNA, tRNA, etc.), fetches all by default.
#' @param status Character vector of seq statuses to be selected (any combination of 'unmapped','unknown','ambiguous', and 'unique'). Fetches all mapped by default.
#' @param absfc Logical; whether to plot absolute log2-foldchange instead of log2-foldchange (default FALSE)
#' @param samples Optional numeric vector indicating which sample(s) to compare. If NULL (default), all samples are used.
#' @param ... any further arguments passed to the scatter.smooth funciton.
#'
#' @return Nothing, but plots multiple frames.
#'
#' @export
plotMAs <- function(o, normalized = TRUE, type = NULL, status = c("unknown", "ambiguous", "unique"), absfc = FALSE, samples = NULL,
                    average = FALSE, ...) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  if (!is.null(samples) & !all(samples %in% 1:ncol(o@seqcounts))) {
    stop("`samples` should be either NULL or an integer vector indicating which sample(s) to compare.")
  }
  m <- getSeqCounts(o, type = type, status = status, normalized = normalized)
  rm <- rowMeans(m)

  if (!is.null(samples)) {
    m <- m[, samples]
  }

  if (average) {
    fc <- m
    for (i in 1:ncol(m)) {
      if (absfc) {
        fc[, i] <- abs(log2((m[, i] + 0.1) / (rm + 0.1)))
      } else {
        fc[, i] <- log2((m[, i] + 0.1) / (rm + 0.1))
      }
    }
    fc <- rowMeans(fc)
    scatter.smooth(log10(rm + 1), fc,
      xlab = "log10(count)", ylab = ifelse(absfc, "absolute log2FC over mean", "log2(foldchange) over mean"),
      pch = 20, col = maketrans("black"), main = "Averaged M", lpars = list(lwd = 2, col = "blue"), ...
    )
    abline(h = 0, lwd = 2)
  } else {
    nf <- autoLayout(ncol(m))
    for (i in 1:ncol(m)) {
      if (absfc) {
        fc <- abs(log2((m[, i] + 1) / (rm + 1)))
      } else {
        fc <- log2((m[, i] + 1) / (rm + 1))
      }
      scatter.smooth(log10(rm + 1), fc,
        xlab = "log10(count)", ylab = ifelse(absfc, "absolute log2FC over mean", "log2(foldchange) over mean"),
        pch = 20, col = maketrans("black"), main = colnames(m)[i], lpars = list(lwd = 2, col = "blue"), ...
      )
      abline(h = 0, lwd = 2)
    }
    if (nf > ncol(m)) {
      for (i in 1:(nf - ncol(m))) plot.new()
    }
  }
}



#' plotComposition
#'
#' Plots the relative abundance of different types of RNAs populating a shortRNAexp dataset.
#'
#' @param o An object of class shortRNAexp
#' @param exclude Character vector of RNA classes to exclude (by default, excludes only ambiguous and unknown). Trumps `include`.
#' @param include Character vector of RNA classes to include (by default, all except those specified in `exclude`).
#' @param normalized Logical; whether to use normalized abundances (default FALSE).
#' @param abridged Logical; whether to return the abridged RNA types (default TRUE)
#' @param scale Logical; whether to linearly scale the abundances so that they sum to 1 (default TRUE)
#' @param ... any further arguments passed to the plotting functions.
#'
#' @return Nothing, but plots multiple frames.
#'
#' @export
plotComposition <- function(o, exclude = NULL, include = NULL, normalized = FALSE, abridged = TRUE, scale = TRUE, ...) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  if (abridged) {
    af <- row.names(o@composition$raw$abridged)
  } else {
    af <- row.names(o@composition$raw$all)
  }
  if (!is.null(include)) {
    if ("tRNA" %in% include & !abridged) {
      include <- c(include, tRNAtype())
    }
    if (!is.null(exclude)) {
      exclude <- unique(c(exclude, setdiff(af, include)))
    } else {
      exclude <- setdiff(af, include)
    }
  } else {
    if (is.null(exclude)) {
      exclude <- c("unknown", "ambiguous")
    }
  }
  su2 <- getComposition(o, abridged = abridged, normalized = normalized, scale = scale, exclude = exclude)

  cols <- getQualitativePalette(nrow(su2))
  od <- order(rowSums(su2), decreasing = T)
  su2 <- su2[od, ]
  cols <- cols[od]

  labs <- paste0(row.names(su2), " (", round(100 * rowSums(su2) / sum(rowSums(su2)), 1), "%)")

  layout(matrix(c(1, 2, 3, 2, 3, 2), nrow = 2))

  plot.new()
  legend("bottomleft", fill = rev(cols), legend = rev(labs), bty = "n", cex = 1.6)

  barplot(as.matrix(su2), col = cols, las = 3, ylab = "proportion of the reads", xpd = T, ...)

  pca <- prcomp(t(su2))
  pp <- pca$x[, 1:2]
  polcomp <- psych::polar(pp, sort = F)[, 2]
  xlab <- paste0("PC 1 (", round(pca$sdev[1] / sum(pca$sdev) * 100, 0), "%)")
  ylab <- paste0("PC 2 (", round(pca$sdev[2] / sum(pca$sdev) * 100, 0), "%)")

  plot(pp,
    pch = 20, col = .colorMap(polcomp), xlim = range(pp[, 1]) * 1.2, cex = 2.5, xlab = xlab, ylab = ylab, main = "Composition PCA",
    ...
  )
  text(pca$x[, 1:2], pos = 1, labels = row.names(pca$x))
  pd <- o@phenoData
  if ("group" %in% colnames(pd)) {
    group <- pd$group
    legend("topright", bty = "n", legend = paste0(c("PC1~group=", "polar~group="), c(format(summary(aov(pca$x[, 1] ~ group))[[1]][
      "group",
      5
    ], digits = 3), format(summary(aov(polcomp ~ group))[[1]]["group", 5], digits = 3))), ...)
  }
}

#' plotFeature
#'
#' Plots the read counts of sequences associated to a given feature (or potentially feature(s) if regular expressions are used), or the aggregated counts if `aggregated`=TRUE.
#'
#' @param o An object of class shortRNAexp
#' @param feature The name of the feature (e.g. gene name) or, if !`exact`), a regular expression to be evaluated against all feature names. One and only one of `sequences` and `feature` should be given.
#' @param sequences A character vector of the sequences to plot. One and only one of `sequences` and `feature` should be given.
#' @param exact Logical; indicates whether `feature` should be interpreted as an exact name (default) rather than a regular expression.
#' @param normalized Logical; whether to return the normalized counts (default TRUE)
#' @param aggregated Logical; whether to return aggregated counts rather than per-sequence counts (default FALSE)
#' @param plotType Either 'barplot' (default), 'heatmap', or 'zheatmap' (heatmap of feature-wise z-scores').
#' @param allowAmbiguous Logical; indicates whether to allow selection of ambiguous sequences when using exact=FALSE (default FALSE).
#' @param proportion Logical; whether to plot proportions (i.e. sequences summing to 1) rather than actual/normalized counts (default FALSE).
#' @param plotLegend Logical; whether to plot the legend (default TRUE); will be plotted on a distinct frame.
#' @param dolog Logical; whether to plot log-transformed counts (default FALSE if plotType='barplot', TRUE otherwise).
#' @param main Plot title, defaults to the name of the feature queried.
#' @param alignSeqs Logical; whether to show sequences as a multiple alignment (default TRUE).
#' @param barplotLayout a matrix populated by 1s and 2s indicating the layout for the barplot (2) and legend (1). Defaults to matrix(c(1,2,2),nrow=1). Ignored if plotType!='barplot'.
#' @param ... arguments passed to the plotting function.
#'
#' @return Nothings but generates plots.
#'
#' @export
plotFeature <- function(o, feature = NULL, sequences = NULL, exact = TRUE, normalized = TRUE, aggregated = FALSE, plotType = "barplot",
                        proportion = FALSE, plotLegend = TRUE, dolog = plotType != "barplot", main = feature, ylab = NULL, allowAmbiguous = FALSE,
                        alignSeqs = !aggregated, barplotLayout = matrix(c(1, 2, 2), nrow = 1), ...) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  if (is.null(feature) & is.null(sequences) | (!is.null(feature) & !is.null(sequences))) {
    stop("Exactly one of `feature` or `sequences` must be given.")
  }
  if (!is.null(sequences) & aggregated) {
    aggregated <- FALSE
    message("Setting `aggregated` to FALSE because `sequences` is given.")
  }
  if (!is.null(sequences) & !all(sequences %in% row.names(o@seqcounts))) {
    stop(paste("Some of the requested sequences are not found in this object:", paste(head(setdiff(sequences, row.names(o@seqcounts))),
      collapse = ", "
    )))
  }
  plotType <- match.arg(plotType, c("barplot", "heatmap", "zheatmap"))
  if (aggregated) {
    if (exact) {
      m <- o@agcounts[which(row.names(o@agcounts) == feature), , drop = F]
      if (nrow(m) == 0) {
        m <- o@agcounts_ambiguous[which(row.names(o@agcounts_ambiguous) == feature), , drop = F]
      }
    } else {
      m <- o@agcounts[grep(feature, row.names(o@agcounts), ignore.case = TRUE), , drop = F]
      if (allowAmbiguous) {
        m <- rbind(m, o@agcounts_ambiguous[grep(feature, row.names(o@agcounts_ambiguous), ignore.case = TRUE), , drop = F])
      }
    }
  } else {
    if (is.null(sequences)) {
      m <- getFeatureCounts(o, feature, exact = exact, normalized = normalized)
    } else {
      if (normalized) {
        m <- normalizeCounts(o@seqcounts[sequences, ], o@norm)
      } else {
        m <- o@seqcounts[sequences, ]
      }
    }
  }
  if (nrow(m) == 0) {
    stop("Feature not found!")
  }
  if (nrow(m) == 1 & plotType != "barplot") {
    message("Single element... reverting to barplot.")
    plotType <- "barplot"
  }
  if (nrow(m) > 1 & alignSeqs & !aggregated) {
    row.names(m) <- msaWrapper(row.names(m))
  } else {
    alignSeqs <- FALSE
  }
  if (!aggregated) {
    row.names(m) <- apply(cbind(row.names(m), o@sources[row.names(m), "cigar"]), 1, FUN = function(x) {
      capitalizeRead(x[1], x[2])
    })
  }
  if (plotType == "barplot") {
    m <- as.matrix(m[order(rowSums(m), decreasing = T), , drop = FALSE])
    if (nrow(m) > 22) {
      m <- as.matrix(rbind(m[1:21, ], as.data.frame(t(colSums(m[22:nrow(m), , drop = F])), row.names = paste0(
        "other (",
        nrow(m) - 21, ")"
      ))))
    }
    if (proportion) {
      if (is.null(ylab)) {
        ylab <- "Proportion of reads assigned to feature"
      }
      m <- t(t(m) / colSums(m))
    } else {
      if (dolog) {
        m <- log(m + 1)
        if (is.null(ylab)) {
          ylab <- paste0("log(", ifelse(normalized, "normalized read count", "read count"), ")")
        }
      } else {
        if (is.null(ylab)) {
          ylab <- ifelse(normalized, "normalized read count", "read count")
        }
      }
    }
    cols <- getQualitativePalette(nrow(m))
    if (plotLegend) {
      layout(barplotLayout)
      par(cex = 1)
      plot.new()

      pp <- par("family")
      if (alignSeqs) {
        par(family = "mono")
      }
      legend("topleft", fill = rev(cols), legend = rev(row.names(m)), cex = 1, bty = "n", xpd = T)
      par(family = pp)
    }
    barplot(m, col = cols, main = main, ylab = ylab, ...)
  } else {
    ra <- data.frame(row.names = row.names(m), abundance = log10(rowSums(m) + 1))
    if (dolog) {
      m <- log(m + 1)
    }
    if (!aggregated) {
      row.names(m) <- apply(cbind(row.names(m), o@sources[row.names(m), "cigar"]), 1, FUN = function(x) {
        capitalizeRead(x[1], x[2])
      })
    }
    .byheatmap(m,
      scale = ifelse(plotType == "zheatmap", "row", "none"), annotation_row = ra, annotation_col = phenoData(o),
      family = ifelse(alignSeqs, "mono", par("family")), ...
    )
  }
}


plotType <- function(o, types, aggregated = FALSE, ambiguous = FALSE, scaleFeatures = TRUE, dolog = TRUE, showNames = NULL, ...) {
  if (aggregated) {
    m <- getAggCounts(o, type = types, normalized = TRUE)
  } else {
    status <- "unique"
    if (ambiguous) {
      status <- c("unique", "ambiguous")
    }
    m <- getSeqCounts(o, type = types, status = status, normalized = TRUE)
  }
  if (nrow(m) > 1000) {
    warning("Only the top 1000 features will be plotted...")
    m <- m[order(rowSums(m), decreasing = T)[1:1000], ]
  }
  ra <- data.frame(row.names = row.names(m), abundance = log10(rowSums(m) + 1))
  if (dolog) {
    m <- log(m + 1)
  }
  if (is.null(showNames)) {
    showNames <- nrow(m) < 30
  }
  .byheatmap(m, scale = ifelse(scaleFeatures, "row", "none"), annotation_row = ra, annotation_col = phenoData(o), show_rownames = showNames)
}



.byheatmap <- function(x, ...) {
  x <- x[which(apply(x, 1, FUN = function(y) {
    !all(is.na(y))
  })), which(apply(x, 2, FUN = function(y) {
    !all(is.na(y))
  }))]
  pheatmap::pheatmap(x, color = colorRampPalette(c("blue", "black", "yellow"))(29), border_color = NA, ...)
}



#' checkSizeSelection
#'
#' Plots different measures of abundance according to sequence length across samples. This essentially calls `plotSizeAbundance` several times on different aggregations.
#'
#' @param o An object of class shortRNAexp
#' @param type Types of sequences to be fetched (e.g. miRNA, tRNA, etc.), fetches all by default.
#' @param status Character vector of seq statuses to be selected (any combination of 'unmapped','invalid','unknown','ambiguous', and 'unique'). Defaults to all mapped sequences.
#' @param trim Logical; whether to trim the 5\% most extreme values on each side before calculating mean (default TRUE).
#'
#' @export
checkSizeSelection <- function(o, type = NULL, status = c("unknown", "invalid", "ambiguous", "unique"), trim = TRUE, onlyDetectedInAll = TRUE) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  layout(matrix(1:4, nrow = 2))
  plotSizeAbundance(o, type, status, normalized = F, dolog10 = F, bty = "n", main = "sum of raw counts", plotLegend = TRUE)
  plotSizeAbundance(o, type, status, normalized = T, dolog10 = F, bty = "n", main = "sum of normalized counts")
  af <- ifelse(trim, "trimmed mean of ", "mean of ")
  cf <- ifelse(onlyDetectedInAll, " (common sequences)", "")
  plotSizeAbundance(o, type, status, normalized = F, dolog10 = F, agFun = ifelse(trim, function(x) {
    mean(x, trim = 0.05)
  }, mean), onlyDetectedInAll = onlyDetectedInAll, bty = "n", main = paste0(af, "raw counts", cf))
  plotSizeAbundance(o, type, status, normalized = T, dolog10 = F, agFun = ifelse(trim, function(x) {
    mean(x, trim = 0.05)
  }, mean), onlyDetectedInAll = onlyDetectedInAll, bty = "n", main = paste0(af, "normalized counts", cf))
}

#' plotSizeAbundance
#'
#' Plots the abundance according to sequence length across samples
#'
#' @param o An object of class shortRNAexp
#' @param type Types of sequences to be fetched (e.g. miRNA, tRNA, etc.), fetches all by default.
#' @param status Character vector of seq statuses to be selected (any combination of 'unmapped','unknown','ambiguous', and 'unique'). Defaults to all mapped sequences.
#' @param normalized Logical; whether to return the normalized counts (default FALSE)
#' @param plotLengthFrequencies Logical; whether to plot the number of sequences with each length (default TRUE)
#' @param onlyDetectedInAll Logical; whether use counts only from sequences detected in all samples (default FALSE)
#' @param dolog10 Logical; whether to use log10-counts (default FALSE).
#' @param ... arguments passed to the plot function.
#'
#' @return Nothing, but produces a plot.
#'
#' @export
plotSizeAbundance <- function(o, type = NULL, status = c("unknown", "ambiguous", "unique"), normalized = TRUE, plotLengthFrequencies = TRUE,
                              onlyDetectedInAll = FALSE, dolog10 = FALSE, agFun = sum, plotLegend = NULL, ...) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  m <- getSeqCounts(o, type, status, normalized)
  if (onlyDetectedInAll) {
    m <- m[which(apply(m, 1, FUN = function(x) {
      all(x > 0)
    })), ]
  }
  sl <- sapply(row.names(m), nchar)
  m <- aggregate(m, by = list(length = sl), FUN = agFun)
  x <- as.numeric(as.character(m[, 1]))
  o <- order(x)
  y <- t(m[o, -1, drop = F])
  if (dolog10) {
    y <- log10(y + 1)
    ylab <- ifelse(normalized, "log10(normalized count+1)", "log10(count+1)")
  } else {
    ylab <- ifelse(normalized, "normalized count", "raw count")
  }
  x <- x[o]
  if (nrow(y) > 22) {
    cols <- rep(maketrans("blue", 50), nrow(y))
    pch <- "."
    lwd <- 1
  } else {
    cols <- getQualitativePalette(nrow(y))
    pch <- 16
    lwd <- ifelse(nrow(y) > 12, 1, 2)
  }
  if (plotLengthFrequencies) {
    ln <- aggregate(sl, by = list(sl), FUN = length)
    ln[, 2] <- ln[, 2] * max(y) / max(ln[, 2])
    plot(ln[, 1], ln[, 2],
      type = "l", lty = "dashed", lwd = 3, col = "lightgrey", xlab = "length", ylab = ylab, ylim = range(y),
      ...
    )
  } else {
    plot(0, 0, col = "white", xlab = "length", ylab = ylab, ylim = range(y), ...)
  }
  for (i in 1:nrow(y)) lines(x, y[i, ], type = "b", pch = pch, col = cols[i], lwd = lwd)
  if (is.null(plotLegend)) {
    plotLegend <- nrow(y) <= 5
  }
  if (plotLegend) {
    legend("topleft", bty = "n", fill = cols, legend = row.names(y))
  }
}




#' getQualitativePalette
#'
#' Returns a qualitative color palette of given size
#'
#' @param nbcolors number of colors (from 1 to 22)
#'
#' @return A vector of colors
#'
# @import colorspace
#'
#' @export
getQualitativePalette <- function(nbcolors) {
  # based on Paul Tol's colors
  if (nbcolors > 22) {
    return(rainbow_hcl(nbcolors))
  }
  switch(as.character(nbcolors),
    `1` = c("#4477AA"),
    `2` = c("#4477AA", "#CC6677"),
    `3` = c("#4477AA", "#DDCC77", "#CC6677"),
    `4` = c("#4477AA", "#117733", "#DDCC77", "#CC6677"),
    `5` = c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677"),
    `6` = c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677", "#AA4499"),
    `7` = c(
      "#332288", "#88CCEE", "#44AA99", "#117733",
      "#DDCC77", "#CC6677", "#AA4499"
    ),
    `8` = c(
      "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677",
      "#AA4499"
    ),
    `9` = c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"),
    `10` = c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"),
    `11` = c(
      "#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255",
      "#AA4499"
    ),
    `12` = c(
      "#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677",
      "#AA4466", "#882255", "#AA4499"
    ),
    `13` = c(
      "#882E72", "#B178A6", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987",
      "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C"
    ),
    `14` = c(
      "#882E72", "#B178A6", "#D6C1DE", "#1965B0",
      "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C"
    ),
    `15` = c(
      "#114477",
      "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777",
      "#771144", "#AA4477", "#DD77AA"
    ),
    `16` = c(
      "#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711",
      "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA", "black"
    ),
    `17` = c(
      "#771155",
      "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#777711", "#AAAA44", "#DDDD77", "#774411",
      "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"
    ),
    `18` = c(
      "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA",
      "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122",
      "#AA4455", "#DD7788"
    ),
    `19` = c(
      "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA",
      "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788", "black"
    ),
    `20` = c(
      "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#117744", "#44AA77",
      "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"
    ),
    `21` = c(
      "#771155",
      "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA",
      "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"
    ),
    `22` = c(
      "#771155",
      "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA",
      "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788", "black"
    ),
    stop("Max 22 colors")
  )
}


#' maketrans
#'
#' Makes a given color transparent.
#'
#' @param tcol the color to be made transparent, specified by any R color specification (see \code{\link[grDevices]{col2rgb}}
#' @param alpha an integer from 0 to 255 indicating the degree of transparency (0 is totally transparent, 255 is totally opaque). Default 100.
#'
#' @return the transparent color code
#'
#' @examples
#' maketrans("black")
#' @export
maketrans <- function(tcol, alpha = 100) {
  c <- col2rgb(tcol)
  rgb(c["red", 1][[1]], c["green", 1][[1]], c["blue", 1][[1]], alpha, maxColorValue = 255)
}



#' autoLayout
#'
#' Creates a tiling layout for the given number of panels
#'
#' @param nb minimum number of panels
#'
#' @return The total number of frames in the layout
#'
#' @export
autoLayout <- function(nb, byrow = FALSE) {
  nc <- ceiling(sqrt(nb))
  nr <- ceiling(nb / nc)
  layout(matrix(1:(nr * nc), nrow = nr, byrow = byrow))
  return(nr * nc)
}

.colorMap <- function(x) {
  if (length(x) == 1) {
    return("blue")
  }
  pal <- colorRampPalette(c("blue", "red"), 0.5)(100)
  xmin <- min(x, na.rm = T)
  xmax <- max(x, na.rm = T) - xmin
  sapply(x, FUN = function(y) {
    pal[1 + floor((99) * (y - xmin) / xmax)]
  })
}


#' volcano
#'
#' Plots a volcano plot from a DEA result data.frame
#'
#' @param res The results data.frame
#' @param writeTop Number of top features to write on the plot (default 0).
#' @param useUncorrected Logical; whether to use uncorrected p-value instead of FDR (default FALSE).
#'
#' @export
volcano <- function(res, writeTop = 0, useUncorrected = FALSE) {
  if (class(res) != "data.frame") {
    stop("`res` should be a data.frame.")
  }
  fdr <- ifelse(useUncorrected, .getPvalField(res), .getFDRfield(res))
  LSD::heatscatter(res$logFC, -log10(res[[fdr]]), xlab = "log2(foldchange)", ylab = ifelse(useUncorrected, "-log10(p-value)",
    "-log10(FDR)"
  ), bty = "n", main = "")
  abline(h = -log10(0.05), lty = "dashed")
  if (writeTop > 0) {
    text(res$logFC[1:writeTop], -log10(res[[fdr]][1:writeTop]), labels = row.names(res)[1:writeTop], xpd = TRUE)
  }
}



#' plotGC
#'
#' Plots the expression across different GC proportions.
#'
#' @param o An object of class shortRNAexp
#' @param aggregateSamples Logical; whether to aggregated samples instead of plotting them individually (default).
#' @param ... Any filtering argument passed to the getSeqCounts function.
#'
#' @export
plotGC <- function(o, aggregateSamples = FALSE, ...) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  m <- getSeqCounts(o, ...)
  gc <- .gcContents(row.names(m))
  if (aggregateSamples) {
    return(LSD::heatscatter(gc, log10(rowSums(m) / sum(rowSums(m))), xlab = "GC content", ylab = "Log10 read proportion"))
  }
  autoLayout(ncol(m))
  for (i in 1:ncol(m)) {
    LSD::heatscatter(gc, log10(m[, i] / sum(m[, i])), xlab = "GC content", ylab = "Log10 read proportion", main = colnames(m)[i])
  }
}


#' msaWrapper
#'
#' Returns multiple alignment of the input sequences
#'
#' @param seqs A character vector containing DNA sequences
#'
#' @return A character vector of length=length(seqs), with each sequence offsetted so that they align
#'
# @import msa
#'
#' @export
msaWrapper <- function(seqs) {
  as.character(msa(toupper(seqs), type = "dna", order = "input")@unmasked)
}


#' plotClipping
#'
#' Plots the distribution of number of bases clipped.
#'
#' @param o An object of class shortRNAexp
#'
#' @export
plotClipping <- function(o) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  layout(matrix(1:2, nrow = 1))
  sc <- sapply(o@sources$cigar, FUN = function(x) {
    x <- .splitCigar(x)
    sum(as.numeric(x[which(x[, 1] %in% c("H", "S")), 2]))
  })
  if (all(sc == 0)) {
    stop("There does not seem to have been any soft-clipping allowed in this alignment.")
  }
  d2 <- aggregate(o@seqcounts[row.names(o@sources), ], by = list(clipping = sc), FUN = sum)
  row.names(d2) <- as.character(d2[, 1])
  d2[, 1] <- NULL
  d2 <- t(t(d2) / colSums(d2))
  scrang <- min(sc[which(sc > 0)]):max(sc)
  cols <- c("grey", .colorMap(scrang))
  scrang <- c("noClipping", paste0(scrang, "nt"))
  hist(sc, xlab = "Nucleotides clipped", ylab = "# unique sequences", col = cols, main = "Clipping summary")
  barplot(d2, col = cols, ylab = "Proportion of the reads", border = NA, las = 3)
  if (length(scrang) > 1) {
    if (length(scrang) > 2) {
      if (length(scrang) > 3) {
        ci <- c(1, 2, ceiling(length(scrang) / 2), length(scrang))
      } else {
        ci <- c(1, round(length(scrang) / 2), length(scrang))
      }
      cols <- cols[ci]
      scrang <- scrang[ci]
    }
    legend("top", inset = -0.1, legend = scrang, fill = cols, bty = "n", xpd = T, horiz = T)
  }
}


#' plotAlignStats
#'
#' Plots alignment statistics
#'
#' @param o An object of class shortRNAexp
#'
#' @export
plotAlignStats <- function(o) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  if (is.null(nrow(o@alignStats$uniqueSeqs)) || nrow(o@alignStats$uniqueSeqs) == 1) {
    m <- cbind(o@alignStats$uniqueSeqs, o@alignStats$reads)
    colnames(m)[1] <- "unique\nseqs"
  } else {
    nds <- nrow(o@alignStats$uniqueSeqs)
    m <- cbind(t(o@alignStats$uniqueSeqs), o@alignStats$reads)
    colnames(m)[1:nds] <- paste("unique\nseqs", 1:nds)
  }
  m <- t(t(m) / colSums(m))
  m <- m[rev(c("unmapped", "ambiguous", "uniqueSoftClipped", "uniqueFullMatch")), ]
  cols <- c("#4477AA", "#117733", "#CC6677", "#DDCC77")
  barplot(m, col = cols, las = 3, ylab = "Proportion", space = c(0, 1, rep(0.1, ncol(m) - 2)))
  legend("top", inset = -0.1, legend = row.names(m), col = cols, pch = 15, bty = "n", xpd = T, horiz = T, pt.cex = 2.2)
}
