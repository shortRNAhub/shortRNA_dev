#' shortRNAexp
#'
#' Creates a shortRNAexp object
#'
#' @param seqcounts Either a numeric matrix (with samples as columns and short sequences as rows) or the path to such a file.
#' @param alignments An alignment file (on the custom genome created along with the `annotation`) for each unique sequence, in BAM format.
#' @param annotation A feature annotation object of type `GRanges`, as created by the `prepareAnnotation` function.
#' @param phenoData Either a data.frame of phenotypic data (with samples as rows and features as columns) or the path to such a file.
#' @param description Character vector describing the experiment/dataset.
#' @param rules Assignment rules (see `?defaultAssignRules`).
#' @param keepAllSrcs Logical; whether to keep information on all sequences' potential sources (before assignment).
#' @param do.aggregate Logical; whether to run aggregation by feature name and type (default TRUE).
#'
#' @return a shortRNAexp object.
#'
#' @export
shortRNAexp <- function(seqcounts, alignments, annotation, phenoData = NULL, description = NA_character_, rules = defaultAssignRules(),
                        keepAllSrcs = TRUE, do.aggregate = TRUE) {
  if (is.character(seqcounts)) {
    seqcounts <- read.delim(seqcounts, header = T, row.names = 1)
  }
  if (is.null(phenoData)) {
    phenoData <- data.frame(row.names = colnames(seqcounts))
  } else {
    if (is.character(phenoData)) {
      phenoData <- read.delim(phenoData, header = T, row.names = 1, stringsAsFactors = F)
    }
  }

  message("
        Finding overlaps between alignments and the annotation...
    ")

  srcs <- findoverlaps.bam.featureAnnotation(alignments, annotation, overlapBy = rules$overlapBy)
  srcs <- srcs[which(srcs$seq %in% row.names(seqcounts)), ]
  srcs$transcript_type <- as.character(srcs$transcript_type)

  message("
        Assigning reads to features...
    ")
  sources <- assignReads(srcs, rules)

  message("
        Building the shortRNAexp object...
    ")
  o <- new("shortRNAexp",
    seqcounts = as.matrix(seqcounts), description = description, sources = sources, phenoData = phenoData,
    rules = rules
  )
  if (!is.null(attr(annotation, "info"))) {
    o@annotationInfo <- attr(annotation, "info")
  }

  o@alignStats <- .getAlignStats(o, srcs)

  if (keepAllSrcs) {
    o@allsrcs <- srcs
  }
  rm(srcs)

  if (do.aggregate) {
    message("
          Aggregating counts by features...
      ")

    o <- aggregateSequences(o)

    o@composition$raw <- estimateComposition(o)
  }

  return(o)
}

.getDefaultAbridgedTypes <- function(types = NULL) {
  newcats <- c(
    "Mt_rRNA", "Mt_tRNA", "rRNA", "protein_coding", "tRNA", "miRNA", "transposable_elements", "transposable_elements",
    "piRNA_precursor", "transposable_elements", "other", "snRNA", "other", "other", "other", "snRNA", "ambiguous", "unknown",
    "primary_piRNA", rep("tRNA", length(tRNAtype()))
  )
  names(newcats) <- c(
    "Mt_rRNA", "Mt_tRNA", "rRNA", "protein_coding", "tRNA", "miRNA", "SINE", "LTR", "piRNA_precursor", "LINE",
    "Simple_repeat", "snoRNA", "Low_complexity", "other", "unprocessed_pseudogene", "snRNA", "ambiguous", "unknown", "primary_piRNA",
    tRNAtype()
  )
  if (!is.null(types)) {
    for (i in setdiff(types, names(newcats))) newcats[i] <- "other"
  }
  return(newcats)
}


#' estimateComposition
#'
#' Estimates the relative abundance of different types of RNAs populating a dataset.
#'
#' @param o An object of class shortRNAexp
#' @param m An optional seqcount table. If ommited, o@seqcounts will be used.
#' @param abridgedTypes An optional named vector of RNA types to be aggregated.
#'
#' @return A list containing two data.frames.
#'
#' @export
estimateComposition <- function(o, m = NULL, abridgedTypes = NULL) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  if (is.null(m)) {
    m <- o@seqcounts
  }

  m <- m[intersect(row.names(m), row.names(o@sources)), ]
  tm <- aggregate(m, by = list(type = o@sources[row.names(m), "transcript_type"]), na.rm = T, FUN = sum)
  row.names(tm) <- tm[, 1]
  tm[, 1] <- NULL

  if (is.null(abridgedTypes)) {
    abridgedTypes <- .getDefaultAbridgedTypes(row.names(tm))
  }

  su <- aggregate(tm, by = list(type = abridgedTypes[row.names(tm)]), FUN = sum)
  row.names(su) <- su$type
  su$type <- NULL

  return(list(all = tm, abridged = su))
}

#' getComposition
#'
#' Get the relative abundance of different types of RNAs populating a shortRNAexp dataset.
#'
#' @param o An object of class shortRNAexp
#' @param abridged Logical; whether to return the abridged RNA types (default FALSE)
#' @param exclude Character vector of RNA classes to exclude (by default, excludes only ambiguous and unknown).
#' @param normalized Logical; whether to use normalized counts rather than raw counts to calculate abundances (default FALSE)
#' @param scale Logical; whether to linearly scale the abundances so that they sum to 1 (default FALSE)
#'
#' @return A matrix of relative abundances.
#'
#' @export
getComposition <- function(o, abridged = FALSE, exclude = c("unknown", "ambiguous"), normalized = FALSE, scale = TRUE) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  if (any(is.null(exclude) || is.na(exclude))) {
    exclude <- c()
  }
  if (normalized) {
    co <- o@composition$normalized
  } else {
    co <- o@composition$raw
  }
  if (abridged) {
    co <- co$abridged
  } else {
    co <- co$all
  }
  if ("unknown" %in% exclude) {
    exclude <- c(exclude, "NA")
  }
  if ("tRNA" %in% exclude & !abridged) {
    exclude <- c(exclude, tRNAtype())
  }
  for (i in exclude) co <- co[which(row.names(co) != i), , drop = F]
  if (scale) {
    co <- t(t(co) / colSums(co))
  }
  return(as.matrix(co))
}

#' getCompositionComponent
#'
#' Retrieves component variables explaining the variation in RNA type composition.
#'
#' @param o An object of class shortRNAexp
#' @param summarytype The type of dimensionality reduction to use. Either 'pca' (reporting the first two principal components) or 'mds' (using a single variable). Default is 'mds'.
#' @param ... Any argument passed to the `getComposition` function.
#'
#' @return A matrix with the composition varible(s).
#'
#' @export
getCompositionComponent <- function(o, summarytype = "MDS", ...) {
  cc <- getComposition(o, ...)
  summarytype <- match.arg(tolower(summarytype), c("pca", "mds"))
  a <- as.matrix(switch(tolower(summarytype),
    mds = cmdscale(dist(t(cc)), 1),
    pca = prcomp(t(cc))$x[, 1:2]
  ))
  colnames(a) <- paste0("compositionV", 1:ncol(a))
  a / max(abs(a))
}


#' getSeqsByType
#'
#' Returns sequences matching filters
#'
#' @param o An object of class shortRNAexp
#' @param type Types of sequences to be fetched (e.g. miRNA, tRNA, etc.), fetches all by default.
#' @param status Character vector of seq statuses to be selected (any combination of 'unmapped','invalid','unknown','ambiguous', and 'unique'). Fetches all by default.
#'
#' @return A character vector of sequences.
#'
#' @export
getSeqsByType <- function(o, type = NULL, status = NULL) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  type <- .parseTypes(o, type)
  if (!is.null(status)) {
    status <- match.arg(status, c("unmapped", "invalid", "unknown", "ambiguous", "unique"), several.ok = TRUE)
  }
  if (is.null(type) & is.null(status)) {
    return(row.names(o@sources))
  }
  if (!is.null(type) & !is.null(status)) {
    return(row.names(o@sources)[which(o@sources$status %in% status & o@sources$transcript_type %in% type)])
  }
  if (is.null(type)) {
    return(row.names(o@sources)[which(o@sources$status %in% status)])
  }
  return(row.names(o@sources)[which(o@sources$transcript_type %in% type)])
}

.parseTypes <- function(o, type) {
  if (is.null(type)) {
    return(NULL)
  }
  if ("tRNA" %in% type) {
    type <- c(setdiff(type, "tRNA"), intersect(row.names(o@composition$raw$all), tRNAtype()))
  }
  type <- match.arg(type, sort(row.names(o@composition$raw$all)), several.ok = TRUE)
  return(type)
}

#' tRNAtype
#'
#' Returns the list of tRNA fragment types
#'
#' @export
tRNAtype <- function() {
  c("tRNA_3p_fragment", "tRNA_5p_fragment", "tRNA_3p_half", "tRNA_5p_half", "tRNA_internal_fragment", "tRNA_fragment")
}

#' getSeqCounts
#'
#' Returns a count matrix for requested sequences
#'
#' @param o An object of class shortRNAexp
#' @param type Types of sequences to be fetched (e.g. miRNA, tRNA, etc.), fetches all by default.
#' @param status Character vector of seq statuses to be selected (any combination of 'unmapped','invalid', 'unknown', 'ambiguous', and 'unique'). Defaults to all mapped sequences.
#' @param normalized Logical; whether to return the normalized counts (default FALSE).
#' @param formatCase Logical; whether to format sequences' case on the basis of clipping (default FALSE).
#'
#' @return A count matrix.
#'
#' @export
getSeqCounts <- function(o, type = NULL, status = c("unknown", "invalid", "ambiguous", "unique"), normalized = FALSE, formatCase = FALSE) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  m <- o@seqcounts[getSeqsByType(o, type, status), , drop = FALSE]
  if (formatCase) {
    row.names(m) <- apply(cbind(row.names(m), o@sources[row.names(m), "cigar"]), 1, FUN = function(x) {
      capitalizeRead(x[1], x[2])
    })
  }
  if (normalized) {
    m <- normalizeCounts(m, o@norm)
  }
  return(m)
}

#' getFeatureSeqs
#'
#' Returns the sequences associated to a feature
#'
#' @param o An object of class shortRNAexp
#' @param feature The name of the feature (e.g. gene name) or, if !`exact`), a regular expression to be evaluated against all feature names.
#' @param exact Logical; indicates whether `feature` should be interpreted as an exact name (default) rather than a regular expression.
#' @param allowAmbiguous Logical; indicates whether to allow selection of ambiguous sequences when using exact=FALSE (default FALSE).
#' @param ignore.case Logical; indicates whether to ignore case when using exact=FALSE (default TRUE)
#'
#' @return A vector of sequences
#'
#' @export
getFeatureSeqs <- function(o, feature, exact = TRUE, allowAmbiguous = FALSE, ignore.case = TRUE) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  if (exact) {
    w <- which(o@sources$transcript_id == feature || o@sources$gene_id == feature)
    s <- row.names(o@sources)[w]
  } else {
    w <- union(grep(feature, o@sources$transcript_id, ignore.case = ignore.case), grep(feature, o@sources$gene_id, ignore.case = ignore.case))
    if (allowAmbiguous) {
      s <- row.names(o@sources)[w]
    } else {
      s <- row.names(o@sources)[intersect(w, which(o@sources$status == "unique"))]
    }
  }
  return(s)
}

#' getFeatureCounts
#'
#' Returns the read counts of sequences associated to a feature
#'
#' @param o An object of class shortRNAexp
#' @param feature The name of the feature (e.g. gene name) or, if !`exact`), a regular expression to be evaluated against all feature names.
#' @param exact Logical; indicates whether `feature` should be interpreted as an exact name (default) rather than a regular expression.
#' @param normalized Logical; whether to return the normalized counts (default FALSE)
#' @param allowAmbiguous Logical; indicates whether to allow selection of ambiguous sequences when using exact=FALSE (default FALSE).
#'
#' @return A count matrix.
#'
#' @export
getFeatureCounts <- function(o, feature, exact = TRUE, normalized = FALSE, allowAmbiguous = FALSE) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  m <- o@seqcounts[getFeatureSeqs(o, feature, exact = exact), , drop = FALSE]
  if (normalized) {
    m <- normalizeCounts(m, o@norm)
  }
  return(m)
}

#' getAggCounts
#'
#' Returns the read counts of aggregated features
#'
#' @param o An object of class shortRNAexp
#' @param type Types of features to be fetched (e.g. miRNA, tRNA, etc.), fetches all by default.
#' @param ambiguous Logical; whether to include ambiguous features (default TRUE if `type` is NULL).
#' @param normalized Logical; whether to return the normalized counts (default FALSE)
#'
#' @return A count matrix.
#'
#' @export
getAggCounts <- function(o, type = NULL, ambiguous = is.null(type), normalized = FALSE) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  type <- .parseTypes(o, type)
  if (!is.null(type)) {
    w <- row.names(o@agdef)[which(as.character(o@agdef$type) %in% type)]
    m <- o@agcounts[w, ]
  } else {
    m <- o@agcounts
  }
  if (ambiguous) {
    m <- rbind(m, o@agcounts_ambiguous)
  }
  if (normalized) {
    if (!is.null(o@norm$offset)) {
      if (is.null(o@norm$agoffset)) {
        message("Aggregating offsets...")
        o@norm$agoffset <- getAggOffset(o, o@norm$offset)
      }
    }
    m <- normalizeCounts(m, o@norm)
  }
  return(m)
}


# computes per-sequence and per-read alignment stats; o is object, sr is seqs.srcs
.getAlignStats <- function(o, sr) {
  sr <- sr[which(!is.na(sr$cigar)), ]
  sr <- sr[which(sr$seq %in% row.names(o@seqcounts)), ]
  unmapped <- setdiff(row.names(o@seqcounts), sr$seq)
  uniq <- row.names(o@sources)[which(o@sources$status == "unique")]
  unknown <- row.names(o@sources)[which(o@sources$status == "unknown")]
  ambi <- row.names(o@sources)[which(o@sources$status == "ambiguous" & o@sources$chrPos == "ambiguous")]
  sr <- sr[order(sr$seq, sr$overlap, decreasing = T), 1:2]
  sr <- sr[!duplicated(sr$seq), ]
  sr <- sr[which(sr$seq %in% uniq), ]
  pmap <- sr$seq[grep("S", sr$cigar, invert = T)]
  scmap <- sr$seq[grep("S", sr$cigar)]
  res <- list(uniqueSeqs = c(
    unmapped = length(unmapped), uniqueFullMatch = length(pmap), uniqueSoftClipped = length(scmap),
    ambiguous = length(ambi), unknown = length(unknown)
  ))
  res$reads <- rbind(colSums(o@seqcounts[unmapped, ]), colSums(o@seqcounts[pmap, ]), colSums(o@seqcounts[scmap, ]), colSums(o@seqcounts[ambi, ]), colSums(o@seqcounts[unknown, ]))
  row.names(res$reads) <- names(res$uniqueSeqs)
  return(res)
}

.mergeAndCollapseColumns <- function(x, y, by = "row.names", checkCol = 1, ...) {
  o <- ncol(x)
  z <- .merge2(x, y, by = by, ...)
  z <- as.data.frame(t(apply(x, 1, o = o, checkCol = checkCol, FUN = function(x, o, checkCol) {
    i <- 1:o
    if (is.na(x[[checkCol]])) {
      i <- i + o
    }
    return(x[i])
  })), stringsAsFactors = FALSE)
  colnames(z) <- colnames(x)
  z
}

.merge2 <- function(x, y, by = "row.names", fill.missing = FALSE, ...) {
  x <- merge(x, y, by = by, ...)
  row.names(x) <- x[, 1]
  x[, 1] <- NULL
  if (fill.missing) {
    x <- as.matrix(x)
    x[which(is.na(x))] <- 0
  }
  return(x)
}



#' mergeshortRNADatasets
#'
#' Merges two shortRNAexp objects. NEEDS TO BE RE-WRITTEN!!!
#'
#' @param o1 An object of class shortRNAexp.
#' @param o2 An object of class shortRNAexp.
#' @param by Level at which to aggregate; either 'feature' (default) or 'sequence'.
#' @param ... Arguments passed to the merge function calls. This can be used to pass for instance `all.x` and `all.y` arguments, although the result object is likely not the be compatible with all the packages' functions afterwards.
#'
#' @return An object of class shortRNAexp.
#'
#' @export
mergeshortRNADatasets <- function(o1, o2, by = "feature", ...) {
  by <- match.arg(by, c("sequence", "feature"))
  if (o@objvers != o@objvers) {
    warning("Objects were built with different version of the class - some trouble might arise!")
  }

  o <- o1
  o@seqcounts <- as.matrix(.merge2(o1@seqcounts, o2@seqcounts, ...))
  if (any(is.na(o@seqcounts))) {
    message("Missing counts replaced by 0.")
    o@seqcounts[is.na(o@seqcounts)] <- 0
  }
  o@description <- paste(o1@description, o2@description, sep = " | ")
  o@phenoData <- as.data.frame(t(.merge2(t(o1@phenoData), t(o2@phenoData), ...)), stringsAsFactors = F)
  o@norm <- list(norm.factors = NA_real_, lib.sizes = NA_real_, cnq = NULL)
  o@allsrcs <- data.frame()
  o@alignStats$reads <- cbind(o1@alignStats$reads, o2@alignStats$reads)
  o@alignStats$uniqueSeqs <- rbind(o1@alignStats$uniqueSeqs, o2@alignStats$uniqueSeqs)
  o@composition$raw$all <- .merge2(o1@composition$raw$all, o2@composition$raw$all, all = TRUE, fill.missing = TRUE)
  o@composition$raw$abridged <- .merge2(o1@composition$raw$abridged, o2@composition$raw$abridged, all = TRUE, fill.missing = TRUE)
  o@composition$normalized <- list()


  sr <- .recastSources(.mergeAndCollapseColumns(o1@sources, o2@sources, checkCol = 6, ...))
  o@sources <- sr[intersect(row.names(sr), row.names(o@seqcounts)), ]

  if (by == "feature") {
    o@agcounts <- .merge2(o1@agcounts, o2@agcounts, by = "row.names", ...)
    o@agcounts_ambiguous <- .merge2(o1@agcounts_ambiguous, o2@agcounts_ambiguous, by = "row.names", ...)
    o@agdef <- .mergeAndCollapseColumns(o1@agdef, o2@agdef, ...)
  } else {
    message("Re-aggregating sequences into features...")
    o <- aggregateSequences(o)
  }

  message("Normalization factors must be recalculated; please run `calcNormFactors`.")

  return(o)
}
