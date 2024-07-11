setClass("shortRNAexp", slots = representation(
  seqcounts = "matrix", description = "character", sources = "data.frame", phenoData = "data.frame",
  composition = "list", norm = "list", alignStats = "list", agcounts = "matrix", agcounts_ambiguous = "matrix", agdef = "data.frame",
  creationDate = "Date", rules = "list", allsrcs = "data.frame", objvers = "numeric", annotationInfo = "list"
), prototype = prototype(
  description = NA_character_,
  alignStats = list(), agcounts = matrix(), agcounts_ambiguous = matrix(), agdef = data.frame(), norm = list(
    norm.factors = NA_real_,
    lib.sizes = NA_real_, cnq = NULL
  ), allsrcs = data.frame(), composition = list(), rules = list(), creationDate = Sys.Date(),
  objvers = 1, annotationInfo = list()
), validity = function(object) {
  errors <- c()
  if (!all(colnames(object@seqcounts) == row.names(object@phenoData))) {
    errors <- c(errors, "The colnames of `seqcounts` should correspond to row.names of `phenoData`.")
  }
  if (!is.numeric(object@seqcounts[1, 1])) {
    errors <- c(errors, "`seqcounts` should be a numeric matrix.")
  }
  if (length(errors) == 0) {
    return(TRUE)
  }
  errors
})


.recastSources <- function(sources) {
  return(sources)
}

setGeneric("phenoData", function(object, ...) {
  standardGeneric("phenoData", ...)
})
setMethod("phenoData", "shortRNAexp", function(object) {
  object@phenoData
})

setMethod("show", "shortRNAexp", function(object) {
  message(paste0("A `shortRNAexp` object of ", nrow(object@seqcounts), " small RNA sequences (", sum(object@sources$status ==
    "unique"), " of which can be uniquely assigned) across ", ncol(object@seqcounts), " samples."))
})


setMethod("[", "shortRNAexp", function(x, i) {
  x@phenoData <- x@phenoData[i, , drop = F]
  x@seqcounts <- x@seqcounts[, i, drop = F]
  x@agcounts <- x@agcounts[, i, drop = F]
  x@agcounts_ambiguous <- x@agcounts_ambiguous[, i, drop = F]
  x@composition$raw$all <- x@composition$raw$all[, i, drop = F]
  x@composition$raw$abridged <- x@composition$raw$abridged[, i, drop = F]
  x@alignStats$reads <- x@alignStats$reads[, i, drop = F]
  if (all(!is.na(x@norm$norm.factors))) {
    x@composition$normalized$all <- x@composition$normalized$all[, i, drop = F]
    x@composition$normalized$abridged <- x@composition$normalized$abridged[, i, drop = F]
    x@norm$norm.factors <- x@norm$norm.factors[i]
    x@norm$lib.sizes <- x@norm$lib.sizes[i]
    for (s in c("offset", "glm.offset", "agoffset", "agglmoffset")) {
      if (!is.null(x@norm[[s]])) {
        x@norm[[s]] <- x@norm[[s]][, i, drop = F]
      }
    }
  }
  return(x)
})

setMethod("names", "shortRNAexp", function(x) {
  colnames(x@seqcounts)
})

setMethod("names<-", "shortRNAexp", function(x, value) {
  if (length(value) != length(unique(value))) {
    stop("Some names are in duplicate!")
  }
  if (length(value) != ncol(x@seqcounts)) {
    stop("The number of names given does not match the number of samples.")
  }
  row.names(x@phenoData) <- value
  colnames(x@seqcounts) <- value
  colnames(x@agcounts) <- value
  colnames(x@agcounts_ambiguous) <- value
  colnames(x@composition$raw$all) <- value
  colnames(x@composition$raw$abridged) <- value
  colnames(o@alignStats$reads) <- value
  if (!is.null(x@composition$normalized)) {
    colnames(x@composition$normalized$all) <- value
    colnames(x@composition$normalized$abridged) <- value
  }
  return(x)
})
