#' runDEA
#'
#' Performs differential expression analysis (DEA) on a shortRNAexp object. This is a wrapper around common RNAseq DEA packages.
#'
#' @param o An object of class shortRNAexp
#' @param formula The formula to test. Either this or `groups` should be given.
#' @param groups Either NULL (if `formula` is given), or the name of a column in phenoData, or a vector of length equal to the number of samples.
#' @param method DEA method, either 'edgeR', 'voom', or 'DESeq2'.
#' @param types The types of RNA to test (default all). A separate analysis will anyway be performed for miRNAs and piRNAs, but this parameter dictates the RNAs that are included in the global analyses.
#' @param normalizeByType Logical; whether each RNA type should be analyzed on its own (default FALSE). If TRUE, the default normalization method is used.
#' @param replicates An optional vector indicating (technical or 'incomplete') replicates (ignored unless method='voom')
#' @param forceGLM Whether to force using GLM even when doing a simple comparison between two groups.
#' @param coef Coefficient to be tested (defaults to the last one of the formula).
#' @param ... Further arguments passed to the removeCorrelatedFeatures function.
#'
#' @return A list containing the DEA results on different sets of features.
# @import edgeR
#' @export
runDEA <- function(o, formula = NULL, groups = NULL, method = "edgeR", types = NULL, normalizeByType = FALSE, replicates = NULL,
                   forceGLM = FALSE, coef = NULL, filterFun = function(x) {
                     sum(x > 20) > 2
                   }, ...) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  if (any(is.na(o@norm$norm.factors))) {
    stop("Run `calcNormFactors` first.")
  }
  method <- match.arg(method, c("voom", "edgeR", "DESeq2"))
  if (method %in% c("voom", "edgeR")) {
  }

  if (!is.null(formula)) {
    if (class(formula) == "character") {
      formula <- as.formula(formula)
    }
  }

  res <- list(seqLevel = list(), aggregated = list(), info = list(
    types = types, norm = o@norm[c("method", "onType", "onStatus")],
    DEA.call = match.call(), assignmentRules = o@rules, phenoData = phenoData(o)
  ))
  if (is.null(res$info$DEA.call$method)) {
    res$info$DEA.call$method <- "edgeR"
  }
  if (is.null(res$info$DEA.call$formula)) {
    res$info$DEA.call$formula <- NA
  }


  if (o@norm$method == "cqn") {
    if (is.null(o@norm$agglmoffset)) {
      message("Calculating aggregated glm offset...")
      o@norm$agglmoffset <- getAggOffset(o, o@norm$glm.offset)
    }
  }

  if (normalizeByType) {
    message("normalizeByType=T; only miRNA and tRNA fragments will be tested.")

    res$seqLevel$miRNA <- .deaWrapper(method, getSeqCounts(o, type = "miRNA", status = "unique", formatCase = T), calcNormFactors.shortRNAexp(o,
      normalizeOnType = "miRNA"
    ), filterFun, formula, groups, forceGLM, coef, replicates)
    res$seqLevel$tRNA <- .deaWrapper(method, getSeqCounts(o, type = "tRNA", status = "unique", formatCase = T), calcNormFactors.shortRNAexp(o,
      normalizeOnType = "tRNA"
    ), filterFun, formula, groups, forceGLM, coef, replicates)
    res$seqLevel$tRNA$name <- o@sources[toupper(row.names(res$seqLevel$tRNA)), "gene_id"]
    res$aggregated$miRNA <- .deaWrapper(method, removeCorrelatedFeatures(
      getAggCounts(o, type = "miRNA", ambiguous = FALSE),
      o, ...
    ), calcNormFactors.shortRNAexp(o, normalizeOnType = "miRNA"), filterFun, formula, groups, forceGLM, coef, replicates)
    res$aggregated$tRNA <- .deaWrapper(method, removeCorrelatedFeatures(
      getAggCounts(o, type = "tRNA", ambiguous = FALSE),
      o, ...
    ), calcNormFactors.shortRNAexp(o, normalizeOnType = "tRNA"), filterFun, formula, groups, forceGLM, coef, replicates)
  } else {
    res$seqLevel$miRNA <- .deaWrapper(
      method, getSeqCounts(o, type = "miRNA", status = "unique", formatCase = T), o, filterFun,
      formula, groups, forceGLM, coef, replicates
    )
    res$seqLevel$tRNA <- .deaWrapper(
      method, getSeqCounts(o, type = "tRNA", status = "unique", formatCase = T), o, filterFun,
      formula, groups, forceGLM, coef, replicates
    )
    res$seqLevel$tRNA$name <- o@sources[toupper(row.names(res$seqLevel$tRNA)), "gene_id"]
    res$seqLevel$primary_piRNA <- NULL
    if ("primary_piRNA" %in% row.names(o@composition$raw$all)) {
      res$seqLevel$primary_piRNA <- .deaWrapper(
        method, getSeqCounts(o, type = "primary_piRNA", status = "unique", formatCase = T),
        o, filterFun, formula, groups, forceGLM, coef, replicates
      )
      if (!is.null(res$seqLevel$primary_piRNA)) {
        res$seqLevel$primary_piRNA$name <- o@sources[toupper(row.names(res$seqLevel$primary_piRNA)), "gene_id"]
      }
    }
    res$seqLevel$allKnownUnique <- .deaWrapper(
      method, getSeqCounts(o, type = types, status = "unique", formatCase = T),
      o, filterFun, formula, groups, forceGLM, coef, replicates
    )
    res$seqLevel$allKnownUnique$name <- o@sources[toupper(row.names(res$seqLevel$allKnownUnique)), "gene_id"]
    res$seqLevel$ambiguous <- .deaWrapper(
      method, getSeqCounts(o, type = types, status = "ambiguous", formatCase = T), o,
      filterFun, formula, groups, forceGLM, coef, replicates
    )
    res$seqLevel$ambiguous$name <- o@sources[toupper(row.names(res$seqLevel$ambiguous)), "gene_id"]
    res$seqLevel$unknown <- .deaWrapper(
      method, getSeqCounts(o, status = "unknown", formatCase = T), o, filterFun, formula,
      groups, forceGLM, coef, replicates
    )

    res$aggregated$miRNA <- .deaWrapper(method, removeCorrelatedFeatures(
      getAggCounts(o, type = "miRNA", ambiguous = FALSE),
      o, ...
    ), o, filterFun, formula, groups, forceGLM, coef, replicates)
    res$aggregated$tRNA <- .deaWrapper(method, removeCorrelatedFeatures(
      getAggCounts(o, type = "tRNA", ambiguous = FALSE),
      o, ...
    ), o, filterFun, formula, groups, forceGLM, coef, replicates)
    res$aggregated$primary_piRNA <- NULL
    res$aggregated$piRNA_clusters <- NULL

    if ("primary_piRNA" %in% row.names(o@composition$raw$all)) {
      a <- getSeqCounts(o, type = "primary_piRNA")
      if (nrow(a) > 1) {
        ag <- aggregate(a, by = list(feature = o@sources[row.names(a), "gene_id"]), FUN = sum)
        row.names(ag) <- ag[, 1]
        ag[, 1] <- NULL
        res$aggregated$primary_piRNA <- .deaWrapper(method, ag, o, filterFun, formula, groups, forceGLM, coef, replicates)
        ag <- aggregate(a, by = list(feature = sapply(o@sources[row.names(a), "gene_id"], FUN = function(x) {
          x <- strsplit(x, "-", fixed = T)[[1]]
          paste(x[-length(x)], collapse = "-")
        })), FUN = sum)
        row.names(ag) <- ag[, 1]
        ag[, 1] <- NULL
        res$aggregated$piRNA_clusters <- .deaWrapper(method, ag, o, filterFun, formula, groups, forceGLM, coef, replicates)
      }
    }
    res$aggregated$allUnique <- .deaWrapper(
      method, getAggCounts(o, type = types, ambiguous = FALSE), o, filterFun, formula,
      groups, forceGLM, coef, replicates
    )
    e <- getAggCounts(o, type = types, ambiguous = TRUE)
    e <- e[setdiff(row.names(e), row.names(res$aggregated$allUnique)), ]
    res$aggregated$ambiguous <- try(.deaWrapper(method, e, o, filterFun, formula, groups, forceGLM, coef, replicates), silent = T)
  }
  return(res)
}

.deaWrapper <- function(method, e, o, filterFun, formula = NULL, groups = NULL, forceGLM = FALSE, coef = NULL, replicates = NULL) {
  e <- e[which(apply(e, 1, FUN = filterFun)), , drop = F]
  if (nrow(e) == 0) {
    return(NULL)
  }
  return(switch(method,
    edgeR = .edgeRwrapper(e, o, formula, groups, forceGLM, coef),
    voom = .voomWrapper(
      e, o, formula, groups,
      coef, replicates
    ),
    DESeq2 = .deseq2Wrapper(e, o, formula = formula, groups = groups, coef = coef),
    stop("Unknown DEA method requested.")
  ))
}

#' .edgeRwrapper
#'
#' @param e 
#' @param o 
#' @param formula 
#' @param groups 
#' @param forceGLM 
#' @param coef 
#'
#' @return
#' @export
#' 
# @import edgeR
#'
#' @examples
.edgeRwrapper <- function(e, o, formula = NULL, groups = NULL, forceGLM = FALSE, coef = NULL) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  if (any(is.na(o@norm$norm.factors))) {
    stop("Run `calcNormFactors` first.")
  }
  if ((!is.null(formula) & !is.null(groups)) | (is.null(formula) & is.null(groups))) {
    stop("Exactly one of `formula` or `groups` must be given.")
  }
  if (is.null(groups)) {
    mm <- model.matrix(formula, data = phenoData(o))
  } else {
    mm <- NULL
    if (length(groups) == 1) {
      if (groups %in% colnames(phenoData(o))) {
        groups <- phenoData(o)[, groups]
      } else {
        stop("`groups` should either be the name of a column in phenoData, or a vector of length equal to the number of samples.")
      }
    } else {
      if (length(groups) != ncol(e)) {
        stop("`groups` should either be the name of a column in phenoData, or a vector of length equal to the number of samples.")
      }
    }
    if (forceGLM | o@norm$method == "cqn") {
      mm <- model.matrix(~groups)
    }
  }
  dds <- DGEList(e, lib.size = o@norm$lib.sizes, norm.factors = o@norm$norm.factors, group = groups)
  if (o@norm$method == "cqn") {
    dds$offset <- .getOffsets(row.names(e), o)
    if (is.null(dds$offset)) {
      warning("Could not find cqn offsets corresponding to features... aborting.")
      return(NULL)
    }
  }
  if (is.null(mm)) {
    dds <- estimateDisp(dds)
    return(as.data.frame(topTags(exactTest(dds), nrow(e))))
  } else {
    if (is.null(coef)) {
      coef <- colnames(mm)[ncol(mm)]
      message(paste("Coefficient tested:", coef))
    }
    dds <- estimateDisp(dds, mm, robust = T)
    return(as.data.frame(topTags(glmLRT(glmFit(dds, mm), coef = coef), nrow(e))))
  }
}

#' .deseq2wrapper
#'
#' @param e 
#' @param o 
#' @param formula 
#' @param groups 
#' @param coef 
#'
#' @return
#' @export
#' 
# @import DESeq2
#'
#' @examples
.deseq2wrapper <- function(e, o, formula = NULL, groups = NULL, coef = NULL) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  if (any(is.na(o@norm$norm.factors))) {
    stop("Run `calcNormFactors` first.")
  }
  if ((!is.null(formula) & !is.null(groups)) | (is.null(formula) & is.null(groups))) {
    stop("Exactly one of `formula` or `groups` must be given.")
  }
  pd <- as.data.frame(phenoData(o))
  if (!is.null(groups)) {
    if (length(groups) == 1) {
      if (!(groups %in% colnames(phenoData(o)))) {
        stop("`groups` should either be the name of a column in phenoData, or a vector of length equal to the number of samples.")
      }
    } else {
      if (length(groups) != ncol(e)) {
        stop("`groups` should either be the name of a column in phenoData, or a vector of length equal to the number of samples.")
      }
      pd$groups <- groups
    }
    formula <- ~groups
  }
  if (o@norm$method == "cqn") {
    stop("cqn normalization is not currently supported with DESeq")
  }
  dds <- DESeqDataSetFromMatrix(floor(as.matrix(e)), colData = pd, design = formula)
  sizeFactors(dds) <- 1 / (o@norm$norm.factors / o@norm$lib.sizes * mean(o@norm$lib.sizes))
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  if (!is.null(coef)) {
    return(as.data.frame(results(dds, name = coef)))
  }
  return(as.data.frame(results(dds)))
}


#' Title
#'
#' @param e 
#' @param o 
#' @param formula 
#' @param groups 
#' @param coef 
#' @param replicates 
#'
#' @return
#' @export
#' 
# @import limma
#'
#' @examples
.voomWrapper <- function(e, o, formula = NULL, groups = NULL, coef = NULL, replicates = NULL) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  if (any(is.na(o@norm$norm.factors))) {
    stop("Run `calcNormFactors` first.")
  }
  if ((!is.null(formula) & !is.null(groups)) | (is.null(formula) & is.null(groups))) {
    stop("Exactly one of `formula` or `groups` must be given.")
  }
  if (is.null(groups)) {
    mm <- model.matrix(formula, data = phenoData(o))
  } else {
    if (length(groups) == 1) {
      if (groups %in% colnames(phenoData(o))) {
        groups <- phenoData(o)[, groups]
      } else {
        stop("`groups` should either be the name of a column in phenoData, or a vector of length equal to the number of samples.")
      }
    } else {
      if (length(groups) != ncol(e)) {
        stop("`groups` should either be the name of a column in phenoData, or a vector of length equal to the number of samples.")
      }
    }
    mm <- model.matrix(~groups)
  }
  dds <- DGEList(e, lib.size = o@norm$lib.sizes, norm.factors = o@norm$norm.factors, group = groups)
  if (o@norm$method == "cqn") {
    dds$offset <- .getOffsets(row.names(e), o)
    if (is.null(dds$offset)) {
      warning("Could not find cqn offsets corresponding to features... aborting.")
      return(NULL)
    }
  }
  v <- voom(dds, mm)
  if (!is.null(replicates)) {
    if (length(replicates) == 1) {
      if (replicates %in% colnames(phenoData(o))) {
        replicates <- phenoData(o)[, replicates]
      } else {
        stop("If given, `replicates` should either be the name of a column in phenoData, or a vector of length equal to the number of samples.")
      }
    } else {
      if (length(replicates) != ncol(e)) {
        stop("If given, `replicates` should either be the name of a column in phenoData, or a vector of length equal to the number of samples.")
      }
    }
    dc <- duplicateCorrelation(v, mm, block = replicates)
    v <- voom(dds, mm, block = replicates, correlation = dc$consensus.correlation)
    fit <- lmFit(v, mm, block = replicates, correlation = dc$consensus.correlation)
  } else {
    fit <- lmFit(v, mm)
  }
  fit <- eBayes(fit)
  if (is.null(coef)) {
    coef <- colnames(mm)[ncol(mm)]
  }
  return(as.data.frame(topTable(eBayes(fit), coef = coef, number = nrow(e))))
}

.getOffsets <- function(features, o) {
  if (all(features %in% row.names(o@norm$glm.offset))) {
    return(o@norm$glm.offset[features, ])
  }
  if (!is.null(o@norm$agglmoffset) & all(features %in% row.names(o@norm$agglmoffset))) {
    return(o@norm$agglmoffset[features, ])
  }
  return(NULL)
}

#' writeDEA
#'
#' Writes the results of a differential expression analysis on a shortRNAexp object in an excel file.
#'
#' @param res A list such as produced by the `runDEA` function.
#' @param file Path to the file where to save the results. Defaults to DEA.xlsx in the working directory.
#' @param fdr.threshold Either NULL (default) or a numeric value between 0 and 1 indicating the maximum FDR value for a feature to be included in the file.
#' 
# @import xlsx
#' 
#' @export
writeDEA <- function(res, file = "DEA.xlsx", fdr.threshold = NULL) {
  if (!all(c("info", "seqLevel", "aggregated") %in% names(res))) {
    stop("`res` does not seem to be a list such as produced by the `runDEA` function.")
  }
  if (!is.null(fdr.threshold)) {
    ff <- .getFDRfield(res$seqLevel[[1]])
    if (is.null(ff)) {
      stop("Could not find FDR field!")
    }
  }
  write.xlsx(res$info$phenoData, file = file, sheetName = "design", row.names = TRUE)
  for (i in c("seqLevel", "aggregated")) {
    for (j in names(res[[i]])) {
      if (!is.null(res[[i]][[j]])) {
        e <- res[[i]][[j]]
        if (!is.null(fdr.threshold)) {
          e[which(e[, ff] < fdr.threshold), ]
        }
        message(paste("Creating sheet:", paste(i, j, sep = "_")))
        write.xlsx(e, file = file, sheetName = paste(i, j, sep = "_"), append = TRUE, row.names = TRUE)
      }
    }
  }
  pp <- function(x) {
    paste(x, collapse = ", ")
  }
  d <- data.frame(param = c(
    "Assignment rules", "Normalization on types", "Normalization on status", "Normalization method",
    "RNA types tested", "DEA method", "DEA formula", "DEA group"
  ), value = c(
    NA, pp(res$info$norm$onType), pp(res$info$norm$onStatus),
    pp(res$info$norm$method), pp(res$info$types), res$info$DEA.call$method, as.character(res$info$DEA.call$formula), pp(res$info$DEA.call$group)
  ))
  write.xlsx(d, file = file, sheetName = "info", append = TRUE, row.names = FALSE)
}

.getFDRfield <- function(e) {
  w <- which(tolower(colnames(e)) %in% c("adj.p.val", "fdr", "padj"))
  if (length(w) == 0) {
    return(NULL)
  }
  colnames(e)[w[1]]
}
.getPvalField <- function(e) {
  w <- which(tolower(colnames(e)) %in% c("pvalue", "p.value", "pval"))
  if (length(w) == 0) {
    return(NULL)
  }
  colnames(e)[w[1]]
}


#' removeCorrelatedFeatures
#'
#' Writes the results of a differential expression analysis on a shortRNAexp object in an excel file.
#'
#' @param counts A count matrix from which to remove correlated subfeatures.
#' @param o An object of class shortRNAexp
#' @param cort Max feature correlation (default 0.6) for splitting
#' @param cvt Max coefficient of variation (default 0.1) for splitting
#' @param minmc Min mean count (default 5) for splitting
#' @param minmdiff Min mean difference (default 1) for splitting

#' @return A filtered version of `counts`
#'
#' @export
removeCorrelatedFeatures <- function(counts, o, cort = 0.6, cvt = 0.1, minmc = 5, minmdiff = 1) {
  e <- normalizeCounts(counts, o@norm)
  a <- o@agdef[row.names(e), , drop = F]
  ll <- strsplit(a$sequences, ";")
  s <- unique(unlist(ll))
  m <- matrix(FALSE, nrow = nrow(a), ncol = length(s))
  for (i in 1:length(ll)) {
    m[i, which(s %in% ll[[i]])] <- TRUE
  }
  row.names(m) <- row.names(a)
  colnames(m) <- s
  gr <- lapply(1:nrow(m), m = m, FUN = function(x, m) {
    sort(row.names(m)[which(apply(m[, which(as.logical(m[x, ])), drop = F], 1, FUN = any))])
  })
  rm(m)
  gr <- gr[which(sapply(gr, FUN = length) > 1)]
  gr <- gr[which(!duplicated(sapply(gr, collapse = ",", FUN = paste)))]

  toRemove <- unique(unlist(lapply(gr, e = e, cort = cort, cvt = cvt, minmc = minmc, minmdiff = minmdiff, FUN = function(x,
                                                                                                                         e, cort, cvt, minmc, minmdiff) {
    a <- e[x, ]
    o <- matrix(0, nrow = length(x), ncol = length(x))
    for (i in 1:length(x)) {
      for (j in 1:length(x)) {
        if (i != j) {
          o[i, j] <- .compareCounts(a[i, ], a[j, ], cort, cvt, minmc, minmdiff)
        }
      }
    }
    x[which(apply(o, 1, FUN = function(i) {
      any(i == 1)
    }))]
  })))
  counts <- counts[setdiff(row.names(counts), toRemove), , drop = F]
  counts
}

.compareCounts <- function(x, y, cort = 0.6, cvt = 0.1, minmc = 5, minmdiff = 1) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (cor(x, y) > cort || all(c(sd(x) / mean(x), sd(y) / mean(y)) < cvt) || any(c(mean(x), mean(y)) < minmc) || mean(abs(x - y)) <
    minmdiff) {
    if (sum(x) < sum(y)) {
      return(1)
    }
    return(2)
  }
  return(0)
}
