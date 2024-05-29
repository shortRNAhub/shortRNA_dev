#' mergeAtomicLists
#'
#' Pairwise merging of two atomic lists of the same length. Does not perform
#' any matching of the names, but triggers an error if names are present and
#' not in the same order.
#'
#' @param x An `AtomicList` object
#' @param y An `AtomicList` object
#'
#' @return An `AtomicList` object of the same length as `x` and `y`
#' @export
#'
#' @importFrom IRanges CharacterList
#' @importFrom S4Vectors splitAsList
#' @examples
#' f1 <- FactorList(list(LETTERS[1:3], LETTERS[1:4]))
#' f2 <- FactorList(list(LETTERS[4:6], LETTERS[7:9]))
#' mergeAtomicLists(f1, f2)
mergeAtomicLists <- function(x, y) {
  stopifnot(length(y) == length(x))
  if (is.null(names(x)) || is.null(names(y))) {
    n <- seq_along(x)
    if (!is.null(names(x))) {
      n <- names(x)
    }
    if (!is.null(names(y))) {
      n <- names(y)
    }
  } else {
    stopifnot(identical(names(x), names(y)))
    n <- names(x)
  }
  if (!is.integer(n)) {
    n <- as.factor(n)
  }
  if (is(x, "FactorList") && is(y, "FactorList")) {
    lvls <- union(levels(x[[1]]), levels(y[[1]]))
    fx <- factor(unlist(x, use.names = FALSE), lvls)
    fy <- factor(unlist(y, use.names = FALSE), lvls)
    return(splitAsList(
      factor(c(fx, fy), levels = seq_along(lvls), lvls),
      c(rep(n, lengths(x)), rep(n, lengths(y)))
    ))
  } else {
    if (is(x, "FactorList")) {
      x <- CharacterList(x)
    }
    if (is(y, "FactorList")) {
      y <- CharacterList(y)
    }
  }
  splitAsList(
    c(unlist(x, use.names = FALSE), unlist(y, use.names = FALSE)),
    c(rep(n, lengths(x)), rep(n, lengths(y)))
  )
}

#' longestOrderedOverlap
#'
#' Returns the longest ordered overlap between elements of an atomic list
#'
#' @param x An AtomicList
#'
#' @return An atomic vector
#' @export
#' @importFrom IRanges IntegerList
#'
#' @examples
#' il <- IntegerList(list(1:4, 1:5, c(1:3, 6)))
#' longestOrderedOverlap(il)
longestOrderedOverlap <- function(x) {
  if (!is(x, "AtomicList")) {
    if (is(x, "list")) {
      # x <- AtomicList(x)
      x <- IntegerList(x)
    } else {
      stop("x should be an atomic list")
    }
  }
  x1 <- x[[1]]
  i <- rep(seq_along(x), lengths(x))
  n <- length(x)
  x <- unlist(x)
  d <- which(!duplicated(i))
  j <- 0
  while (length(x) > 0 && length(d) == n && length(u <- unique(x[d])) == 1) {
    j <- j + 1
    x <- x[-d]
    i <- i[-d]
    d <- which(!duplicated(i))
  }
  x1[seq_len(j)]
}

#' fList2tree
#'
#' Converts a `FactorList` to a `phylo` object.
#'
#' @param fL A `FactorList` object
#' @param addRoot Logical; whether to add a root node
#' @param collapseSingles Logical; wether to remove intermediate nodes of single
#' childs
#' @importFrom ape collapse.singles
#' 
#' @return A `phylo` object
#' @export
fList2tree <- function(fL, addRoot = TRUE,
                       collapseSingles = FALSE, root = "ROOT") {
  fL <- relist(droplevels(unlist(fL, use.names = FALSE)), fL)
  lvls <- levels(fL[[1]])
  # temporarily get rid of levels
  x <- IntegerList(fL)
  maxLen <- max(lengths(x))
  f <- rep(seq_along(x), lengths(x))
  x <- unlist(x)
  w <- which(!duplicated(f))
  xfrom <- integer(0)
  xto <- integer(0)
  while (length(w) > 0) {
    xfrom <- c(xfrom, x[w])
    x <- x[-w]
    f <- f[-w]
    w <- which(!duplicated(f))
    xto <- c(xto, x[w])
    wdup <- which(f %in% unique(f[-w]))
    x <- x[wdup]
    f <- f[wdup]
    w <- which(!duplicated(f))
  }
  m <- matrix(c(xfrom, xto), nrow = length(xfrom))
  rm(xfrom, xto)
  m <- m[!duplicated(m), , drop = FALSE]
  if (addRoot) {
    roots <- setdiff(unique(m[, 1]), unique(m[, 2]))
    m <- rbind(cbind(rep(length(lvls) + 1L, length(roots)), roots), m)
    colnames(m) <- NULL
    lvls <- c(lvls, root)
  }
  tips <- setdiff(unique(m[, 2]), unique(m[, 1]))
  ro <- unique(c(tips, max(m), unique(m[, 1])))
  conv <- seq_along(ro)[order(ro)]
  m <- matrix(conv[as.integer(m)], ncol = 2)
  lvls <- lvls[ro]
  tips <- setdiff(unique(m[, 2]), unique(m[, 1]))
  nodes <- setdiff(unique(m[, 1]), tips)
  o <- list()
  class(o) <- "phylo"
  o$edge <- m
  o$Nnode <- length(nodes)
  o$node.label <- lvls[nodes]
  o$tip.label <- lvls[tips]
  if (collapseSingles) {
    o <- ape::collapse.singles(tree = o, root.edge = TRUE)
  }
  o
}


.addNames2fL <- function(fL) {
  nlvls <- length(levels(fL[[1]]))
  lvls <- c(levels(fL[[1]]), names(fL))
  x <- c(as.integer(unlist(fL, use.names = FALSE)), nlvls + seq_len(length(fL)))
  x <- factor(x, seq_len(nlvls + length(fL)), lvls)
  splitAsList(x, c(rep(seq_along(fL), lengths(fL)), seq_along(fL)))
}



#' Add reads to the tree
#'
#' @param fL factorList of the annotation
#' @param mappedFeaturesDF `DFrame` with Sequences and transcript IDs
#' @param featuresCol Feature name column
#' @param readsCol Column with sequences, if omitted will use the row.names
#' @param unassigned If unassigned reads to be included in the tree or not
#' @param extraTreeBranch Additional sequences as `FactorList` to be added to
#'  the tree. One can obtain unaligned reads using `getReadsFromBam` to obtain
#'  reads from bam files.
#' @param ... Other parameters to be passed to `fList2tree`
#'
#' @importFrom IRanges IntegerList FactorList
#' @importFrom S4Vectors splitAsList
#'
#' @export
addReadsToTree <- function(fL,
                           mappedFeaturesDF,
                           featuresCol = "transcript_id",
                           readsCol = NULL,
                           unassigned = FALSE,
                           extraTreeBranch = NULL,
                           ...) {

  # Parallel processing of apply functions
  # plan(multisession)

  # Features per sequence
  features <- mappedFeaturesDF[, featuresCol]
  if(is.null(readsCol)){
    names(features) <- row.names(mappedFeaturesDF)
  }else{
    names(features) <- mappedFeaturesDF[, readsCol]
  }

  # # Filtering out non-existing features
  # featuresPerRead <- lapply(features, function(x) intersect(x = as.character(x), y = names(fl)))
  # sm <- sum(lengths(featuresPerRead) == 0)
  #
  # if (sm > 0) {
  #   warning(paste(sm, "features were removed as they were not found in the tree!"))
  # }
  #
  # featuresPerRead <- featuresPerRead[lengths(featuresPerRead) > 0]

  # Features per sequence factorList
  # fpr_fl <- FactorList(featuresPerRead)
  fpr_fl <- features

  # For sequences overlapping with 1 feature
  fpr_fl_s <- fpr_fl[lengths(fpr_fl) == 1]

  fpr_fl_s_o <- fl[unlist(fpr_fl_s)]
  names(fpr_fl_s_o) <- names(fpr_fl_s)
  fpr_fl_s_o <- .addNames2fL(fL = fpr_fl_s_o)

  # fpr_s_tree <- fList2tree(fL = fpr_fl_s_o, addRoot = T, collapseSingles = F)

  # For sequences overlapping with more than 1 feature
  fpr_fl_m <- fpr_fl[lengths(fpr_fl) > 1]

  if (length(fpr_fl_m) > 0) {
    fli <- IntegerList(fL)
    fpr_fl_m_o <- lapply(fpr_fl_m, FUN = function(x) {
      longestOrderedOverlap(fli[x])
    })
    rm(fli)

    fpr_fl_m_o1 <- splitAsList(
      x = factor(as.numeric(unlist(fpr_fl_m_o)),
        levels = seq_along(levels(fl[[1]])),
        labels = levels(fL[[1]])
      ),
      f = rep(names(fpr_fl_m_o), lengths(fpr_fl_m_o))
    )

    fpr_fl_m_o2 <-
      FactorList(as.list(as.character(paste(fpr_fl_m, collapse = "/"))))

    fpr_fl_m_of <- mergeAtomicLists(x = fpr_fl_m_o1, y = fpr_fl_m_o2)
    names(fpr_fl_m_of) <- names(fpr_fl_m_o1)

    fpr_fl_m_of <- .addNames2fL(fpr_fl_m_of)
    # Features for the tree
    fpr_tree <- c(fpr_fl_s_o, fpr_fl_m_of)
  } else {
    fpr_tree <- fpr_fl_s_o
  }

  # Features for the tree
  # fpr_tree <- c(fpr_fl_s_o, fpr_fl_m_of)

  # Reads not assigned to any features
  if (unassigned) {
    fpr_fl_n <- fpr_fl[lengths(fpr_fl) == 0]
    fpr_fl_n <- .addNames2fL(fpr_fl_n)
    fpr_fl_n <- as.character(unlist(fpr_fl_n))
    n <- length(fpr_fl_n)
    fpr_fl_n <- splitAsList(as.character(unlist(fpr_fl_n)), seq_len(n))

    x1 <- splitAsList(
      rep(factor(as.character(fpr_fl_m_of[[1]][1])), n),
      seq_len(n)
    )
    x2 <- splitAsList(rep(factor("unassigned"), n), seq_len(n))
    x <- mergeAtomicLists(x1, x2)

    fpr_fl_n_o <- mergeAtomicLists(x, fpr_fl_n)

    fpr_tree <- c(fpr_tree, fpr_fl_n_o)
  }


  # If additional branch is provided
  if (!is.null(extraTreeBranch)) fpr_tree <- c(fpr_tree, extraTreeBranch)

  tree <- fList2tree(fL = fpr_tree, addRoot = FALSE, collapseSingles = F, ...)

  # plan(sequential)

  return(tree)
}


#' Get reads from BAM file for a specific samFlag
#'
#' @param bam path to the BAM file.
#' @param flag Flag for which reads are to be obtained. Default: 4 (Unaligned)
#' @param bam A label to be given to reads
#'
#' @return A `FactorList` of reads with labels
#'
#' @importFrom Rsamtools BamFile scanBam
#' @importFrom IRanges FactorList
#'
#' @export
getReadsFromBam <- function(bam, flag = 4, label = "unaligned") {
  bamFile <- BamFile(bam)
  aln <- scanBam(bamFile)[[1]]
  ua <- aln$qname[aln$flag == flag]

  n <- length(ua)
  uafl <- FactorList(splitAsList(ua, seq_len(n)))

  x <- splitAsList(rep(factor(label), n), seq_len(n))
  res <- mergeAtomicLists(x, uafl)

  return(res)
}
