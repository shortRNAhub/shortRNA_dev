#' calcNormFactors.shortRNAexp
#'
#' Calculates the normalization factors for a shortRNAexp dataset.
#'
#' @param object An object of class shortRNAexp
#' @param method Normalization method (any accepted by edgeR). Defaults to 'TMM'.
#' @param normalizeOnType Types of sequences (e.g. miRNA, tRNA, etc.) to consider for calculating normalization factors (default all).
#' @param normalizeOnStatus Sequence statuses to consider for calculating normalization factors (default all mapped)
#' @param ... Further arguments passed on to edgeR's generic calcNormFactors.
#'
#' @return The shortRNAexp object with updated norm.factors and lib.sizes slots.
#'
# @import edgeR
#'
#' @export
calcNormFactors.shortRNAexp <- function(object, method = "TMM", normalizeOnType = NULL, normalizeOnStatus = c(
                                          "unknown", "ambiguous",
                                          "unique"
                                        ), ...) {
  if (!is(object, "shortRNAexp")) {
    stop("`object` should be an object of class `shortRNAexp`.")
  }
  ss <- getSeqsByType(object, normalizeOnType, normalizeOnStatus)
  if (length(ss) == 0) {
    stop("No sequence matched the request!")
  }
  if (length(ss) == 1) {
    stop("A single sequence matched the request!")
  }
  if (length(ss) < 200) {
    warning(paste("Only", length(ss), "matched the request!"))
  }
  m <- object@seqcounts[ss, ]
  object@norm$onType <- normalizeOnType
  object@norm$onStatus <- normalizeOnStatus
  object@norm$lib.sizes <- colSums(m)
  if (tolower(method) == "cqn") {
    object <- .cqnNormWrapper(object, m)
  } else {
    dds <- calcNormFactors(DGEList(m), method = method, ...)
    object@norm$lib.sizes <- dds$samples$lib.size
    object@norm$norm.factors <- dds$samples$norm.factors
    object@norm$offset <- NULL
    object@norm$glm.offset <- NULL
    object@norm$agoffset <- NULL
    object@norm$agglmoffset <- NULL
    object@norm$method <- method
  }
  object@composition$normalized <- estimateComposition(object, normalizeCounts(object@seqcounts, object@norm))
  object
}

#' .cqnNormWrapper
#'
#' @param o 
#' @param m 
#' @param minc 
#'
#' @return
#' @export
#' 
# @import cqn
#'
#' @examples
.cqnNormWrapper <- function(o, m, minc = 10) {
  subs <- row.names(m)[which(rowMeans(m) > minc)]
  subs <- which(row.names(o@seqcounts) %in% subs)
  qn <- cqn(o@seqcounts,
    x = .gcContents(row.names(o@seqcounts)), lengths = sapply(row.names(o@seqcounts), nchar), sizeFactors = o@norm$lib.sizes,
    subindex = subs
  )
  o@norm$norm.factors <- rep(1, ncol(o@seqcounts))
  o@norm$offset <- qn$offset
  o@norm$glm.offset <- qn$glm.offset
  o@norm$agoffset <- getAggOffset(o, qn$offset)
  o@norm$agglmoffset <- getAggOffset(o, qn$glm.offset)
  o@norm$method <- "cqn"
  return(o)
}

#' getAggOffset
#'
#' Aggregates offsets for cqn normalization using weighted means
#'
#' @param o An object of class shortRNAexp
#' @param offset The offset matrix to aggregate
#'
#' @return A numeric matrix.
#'
#' @export
getAggOffset <- function(o, offset) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  o1 <- t(sapply(row.names(o@agcounts), ag = o@agdef, offset = offset, FUN = function(x, ag, offset) {
    s <- strsplit(as.character(ag[x, "sequences"]), ";", fixed = T)[[1]]
    offset <- offset[s, , drop = F]
    w <- o@seqcounts[s, , drop = F]
    m <- rep(0, ncol(offset))
    for (i in 1:length(m)) {
      w[which(!(w >= 0))] <- 0
      if (any(w[, i] > 0)) {
        m[i] <- weighted.mean(offset[, i], w[, i], na.rm = T)
      } else {
        m[i] <- mean(offset[, i], na.rm = T)
      }
    }
    return(m)
  }))
  o2 <- t(sapply(row.names(o@agcounts_ambiguous), o = o, offset = offset, FUN = function(x, o, offset) {
    s <- row.names(o@sources)[which(o@sources$gene_id == x)]
    offset <- offset[s, , drop = F]
    w <- o@seqcounts[s, , drop = F]
    m <- rep(0, ncol(offset))
    for (i in 1:length(m)) {
      w[which(!(w >= 0))] <- 0
      if (any(w[, i] > 0)) {
        m[i] <- weighted.mean(offset[, i], w[, i], na.rm = T)
      } else {
        m[i] <- mean(offset[, i], na.rm = T)
      }
    }
    return(m)
  }))
  return(rbind(o1, o2))
}


#' normalizeCounts
#'
#' Normlalizes a count matrix using the given norm.factors
#'
#' @param counts A numeric count matrix (with samples as columns)
#' @param normParams A list of normalization parameters (e.g. the slot `norm` of a `shortRNAexp` object).
#'
#' @return A normalized matrix.
#'
#' @export
normalizeCounts <- function(counts, normParams) {
  if (any(is.na(normParams$norm.factors))) {
    stop("Run `calcNormFactors` first.")
  }
  if (is.null(normParams$lib.sizes)) {
    normParams$lib.sizes <- colSums(counts)
  }
  if (is.null(normParams$offset)) {
    return(cpm(DGEList(counts = counts, norm.factors = normParams$norm.factors, lib.size = normParams$lib.sizes), normalized.lib.sizes = TRUE) *
      mean(normParams$lib.sizes) / 1e+06)
  } else {
    if (all(row.names(counts) %in% row.names(normParams$offset))) {
      o <- normParams$offset[row.names(counts), ]
    } else {
      if (is.null(normParams$agoffset)) {
        stop("Could not find cqn offset for features")
      }
      if (all(row.names(counts) %in% row.names(normParams$agoffset))) {
        o <- normParams$agoffset[row.names(counts), ]
      } else {
        stop("It appears that the features do not have a defined cqn offset!")
      }
    }
    x <- 2^(log2(t(1e+06 * t(counts) / normParams$lib.sizes)) + o)
    x <- x - min(min(x), 1)
    x <- t(t(x) * mean(normParams$lib.sizes) / 1e+06)
    return(x)
  }
}



#' Conditional Quantile Normalization (CQN)
#'
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
# @import cqn cqn
#'
#' @seealso cqn cqn
#'
#' @param IDs A GEO sample accession or SRA run accession, `GSM...` or `GSE...` or `SRX...` or `SRR...` or `SRP....`. Multiple IDs input is possible.
#' @param sqlfile A SQL database of SRAdb. If not provides, it will be downloaded from internet.
#' @param nParallel No. of processes for dumping fastq files. Default: 4
#' @param suppGEO If supplementary GEO files to be downloaded. Default: TRUE
#'
#' @return A Directory/ directories with the same name as ID/ IDs and downloaded data in it (or log file if data is not availabls via SRAdb).
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
# cqnSR <- function()
