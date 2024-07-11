#' aggregateSequences
#'
#' Aggregates individual sequences into features
#'
#' @param o An object of class shortRNAexp.
#' @param quiet Logical; indicates whether to suppresses messages (default FALSE).
#'
#' @return An updated object of class shortRNAexp.
#'
#' @export
aggregateSequences <- function(o, quiet = FALSE) {
  if (!quiet) {
    message("- sequences of unique origins")
  }
  o <- .aggUnambiguous(o)

  if (!quiet) {
    message("- sequences of ambiguous origins")
  }
  o <- .aggAmbiguous(o)

  return(o)
}

.aggUnambiguous <- function(o) {
  if (!is(o, "shortRNAexp")) {
    stop("`o` should be an object of class `shortRNAexp`.")
  }
  m <- o@seqcounts
  sources <- o@sources
  l1 <- .doag(m, sources, which(sources$status == "unique"))
  l2 <- .doag(m, sources, which(sources$transcript_type %in% c("tRNA", tRNAtype()) & sources$status == "unique"), nFUN = .tRNAbasename)
  l3a <- .doag(m, sources, which(sources$transcript_type == "miRNA" & sources$status == "unique" & !grepl(";", sources$transcript_id)),
    by = "transcript_id"
  )
  l3b <- .doag(m, sources, which(sources$transcript_type == "miRNA" & sources$status == "unique"), nFUN = .miRNAbasename)
  l4 <- .doag(m, sources, which(sources$transcript_type == "miRNA" & sources$status == "unique"), nFUN = function(x) {
    .miRNAbasename(x, T)
  })
  w <- which(sources$transcript_type %in% tRNAtype() & sources$status == "unique")
  l5 <- list(os = sources$gene_id[w], as = paste(.tRNAbasename(sources$gene_id[w]), gsub("^tRNA_", "", sources$transcript_type[w]),
    sep = "_"
  ), ss = row.names(sources)[w])
  l5$ag <- aggregate(m[l5$ss, , drop = F], by = list(l5$as), FUN = sum)
  l5$tt <- as.character(sources$transcript_type[w])

  ag <- aggregate(rbind(l1$ag, l2$ag, l3a$ag, l3b$ag, l4$ag, l5$ag)[, -1, drop = F], by = list(c(l1$ag[, 1], l2$ag[, 1], l3a$ag[
    ,
    1
  ], l3b$ag[, 1], l4$ag[, 1], l5$ag[, 1])), FUN = max)
  row.names(ag) <- ag[, 1]
  ag[, 1] <- NULL

  # remove duplicated rows if names are a subset
  tsig <- apply(ag, 1, collapse = ";", FUN = paste)
  toRemove <- unlist(sapply(unique(tsig[which(duplicated(tsig))]), tsig, ag, FUN = function(x, tsig, ag) {
    e <- row.names(ag[which(tsig == x), ])
    e <- e[order(nchar(e))]
    e[which(sapply(e[1:(length(e) - 1)], e = e, FUN = function(y, e) {
      any(grepl(y, e, fixed = T))
    }))]
  }))
  ag <- ag[which(!(row.names(ag) %in% toRemove)), ]

  o@agcounts <- as.matrix(ag)

  ll <- list(l1, l2, l3a, l3b, l4, l5)
  rm(l1, l2, l3a, l3b, l4, l5)
  os <- unlist(lapply(ll, FUN = function(x) {
    x$os
  }))
  as <- unlist(lapply(ll, FUN = function(x) {
    x$as
  }))
  tt <- unlist(lapply(ll, FUN = function(x) {
    x$tt
  }))
  ss <- unlist(lapply(ll, FUN = function(x) {
    x$ss
  }))

  ag <- aggregate(cbind(os, ss, tt), by = list(as), FUN = function(x) {
    paste(sort(unique(as.character(x))), collapse = ";")
  })
  row.names(ag) <- ag[, 1]
  ag <- ag[row.names(o@agcounts), -1]
  colnames(ag) <- c("originalNames", "sequences", "type")
  o@agdef <- ag

  o
}

.aggAmbiguous <- function(o) {
  u <- unique(o@sources$gene_id[which(o@sources$status == "ambiguous" & !is.na(o@sources$gene_id) & o@sources$gene_id != "ambiguous")])
  ag <- matrix(0, nrow = length(u), ncol = ncol(o@seqcounts))
  row.names(ag) <- u
  for (i in 1:length(u)) {
    g <- u[i]
    sr <- o@sources[which(o@sources$gene_id == g), , drop = F]
    ag[i, ] <- colSums(o@seqcounts[row.names(sr), , drop = F])
    if (all(sr$transcript_type %in% tRNAtype())) {
      type <- NULL
      if (length(unique(sr$transcript_type)) == 1) {
        type <- sr$transcript_type[1]
      } else {
        if (all(grepl("5p", sr$transcript_type))) {
          type <- "5p_fragment"
        } else {
          if (all(grepl("3p", sr$transcript_type))) {
            type <- "3p_fragment"
          }
        }
      }
      if (!is.null(type)) {
        row.names(ag)[i] <- paste(g, gsub("tRNA_", "", type, fixed = T))
      }
    }
  }
  o@agcounts_ambiguous <- ag

  return(o)
}

.doag <- function(m, sources, w, by = "gene_id", nFUN = function(x) {
                    x
                  }) {
  os <- sources[[by]][w]
  as <- nFUN(sources[[by]][w])
  ss <- row.names(sources)[w]
  ag <- aggregate(m[ss, , drop = F], by = list(as), FUN = sum)
  tt <- as.character(sources$transcript_type[w])
  return(list(os = os, as = as, ss = ss, ag = ag, tt = tt))
}
