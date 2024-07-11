#' findoverlaps.bam.featureAnnotation
#'
#' Map a BAM file to the features file
#'
#' @param bamFile A character vector of path to BAM file.
#' @param featureAnnotation A GRanges object as produced by prepareAnnotation function.
#' @param ignoreStrand Logical; whether to ignore strand when searching for overlaps (default TRUE)
#'
#' @return A dataframe.
#' 
# @import GenomicAlignments GenomicRanges Rsamtools
#' 
#' @export
findoverlaps.bam.featureAnnotation <- function(bamFile, featureAnnotation, ignoreStrand = TRUE) {
  # convert bam to GRanges
  param <- ScanBamParam(what = c("cigar", "pos", "seq"))
  bam <- readGAlignments(bamFile, param = param)
  bam <- GRanges(bam)

  # Finding overlaps
  overlap <- suppressWarnings(findOverlaps(bam, featureAnnotation, ignore.strand = ignoreStrand))
  read <- data.frame(bam[overlap@from, ], stringsAsFactors = F)
  colnames(read) <- paste0(colnames(read), "Read")
  feature <- data.frame(featureAnnotation[overlap@to], stringsAsFactors = F)
  colnames(feature) <- paste0(colnames(feature), "Feature")

  m <- cbind(read, feature)

  # Position of the featureAnnotation
  m$posInFeature <- apply(m, 1, function(x) {
    ifelse(x[13] == "+", as.numeric(x["startRead"]) - as.numeric(x["startFeature"]),
      as.numeric(x["endFeature"]) - as.numeric(x["endRead"])
    )
  })

  m$overlap <- apply(m, 1, function(x) {
    min(as.numeric(x["endRead"]), as.numeric(x["endFeature"])) - max(
      as.numeric(x["startRead"]),
      as.numeric(x["startFeature"])
    ) + 1
  })

  m$percentOverlap <- apply(m, 1, function(x) {
    (as.numeric(x["overlap"]) / (as.numeric(x["endRead"]) - as.numeric(x["startRead"]) +
      1)) * 100
  })

  m$featureWidth <- apply(m, 1, function(x) as.numeric(x["endFeature"]) - as.numeric(x["startFeature"]) + 1)

  # Non overlapping reads
  nonOverlapRead <- data.frame(bam[-overlap@from, ], stringsAsFactors = F)
  colnames(nonOverlapRead) <- paste0(colnames(nonOverlapRead), "Read")
  colsAddition <- colnames(m)[!colnames(m) %in% colnames(nonOverlapRead)]
  newCols <- data.frame(matrix(nrow = nrow(nonOverlapRead), ncol = length(colsAddition)), stringsAsFactors = F)
  colnames(newCols) <- colsAddition
  newCols[, "overlap"] <- 0
  newCols[, "percentOverlap"] <- 0

  nonOverlapRead <- cbind(nonOverlapRead, newCols)

  # Final table
  m <- rbind(m, nonOverlapRead)
  m$chrPos <- paste(m[, 1], m[, 7], sep = ":")

  columns <- c(
    "seqRead", "cigarRead", "chrPos", "strandRead", "strandFeature", "posInFeature", "overlap", "percentOverlap",
    "transcript_idFeature", "gene_idFeature", "transcript_typeFeature", "featureWidth"
  )
  r <- m[, columns]
  names(r) <- c(
    "seq", "cigar", "chrPos", "strandRead", "strandFeature", "posInFeature", "overlap", "percentOverlap", "transcript_id",
    "gene_id", "transcript_type", "featureWidth"
  )
  r <- r[order(r[, "seq"]), ]
  return(r)
}
