#' Prepare annotation for short RNA-seq analysis
#'
#' This function prepares annotation data for short RNA-seq analysis by processing
#' transcript information, handling extra sequences and genomic ranges, and
#' building a genome index.
#'
#' @param ensdb An 'EnsDb' object, as obtained from AnnotationHub
#' @param genome A genome fasta, BSGenome or TwoBit object. If missing, will be
#'   fetched using `ensdb`
#' @param output_dir The path where to save the files
#' @param extra.gr A list of GRanges or GRangesList giving the coordinates of
#'   extra features
#' @param extra.seqs Extra transcript sequences to include
#' @param resolveSplicing The biotypes for which to resolve splicing
#' @param rules Rules for reads assignment
#' @param tRNAEnsembleRemove Remove tRNA annotations from Ensembl database
#' @param clusterMiRNA Whether to cluster miRNA annotations
#' @param ... Additional arguments passed to `Rsubread::buildindex`
#'
#' @return A GRangesList of features
#'
#' @import S4Vectors
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom Biostrings DNAStringSet writeXStringSet
#' @importFrom Rsubread buildindex
#' @importFrom dplyr bind_rows
#' @importFrom ensembldb genes exonsBy seqlevelsStyle getGenomeFaFile
#' @importFrom R.utils gzip
#' @importFrom rtracklayer import
#'
#' @export
prepareAnnotation <- function(ensdb, genome = NULL, output_dir = "",
                              extra.gr = list(), extra.seqs = NULL,
                              resolveSplicing = NULL,
                              rules = defaultAssignRules(),
                              tRNAEnsembleRemove = TRUE,
                              clusterMiRNA = TRUE,
                              ...) {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  tx <- process_transcripts(ensdb, resolveSplicing, rules, tRNAEnsembleRemove)
  tx <- process_extra_sequences(tx, extra.seqs)
  tx <- process_extra_granges(tx, extra.gr)
  
  tx <- finalize_transcripts(tx, clusterMiRNA)
  
  anno.out <- save_features(tx, output_dir)
  
  extra.seqs <- prepare_extra_seqs(extra.seqs)
  
  genome <- process_genome(genome, ensdb, extra.seqs, output_dir)
  
  build_index(genome, output_dir, ...)
  
  return(tx)
}

#' Process Transcripts
#'
#' This function processes transcript information from an EnsDb object.
#'
#' @param ensdb An 'EnsDb' object
#' @param resolveSplicing The biotypes for which to resolve splicing
#' @param rules Rules for reads assignment
#' @param tRNAEnsembleRemove Whether to remove tRNA annotations
#'
#' @return A GRanges object containing processed transcript information
#'
#' @importFrom ensembldb genes exonsBy
process_transcripts <- function(ensdb, resolveSplicing, rules, tRNAEnsembleRemove) {
  if (is.null(resolveSplicing)) {
    p <- rules$priorities
    resolveSplicing <- names(p)[p >= 0]
  }
  
  if (!is.null(resolveSplicing)) {
    gs <- get_genes(ensdb, resolveSplicing)
    anofilter <- ~ tx_biotype == resolveSplicing
  } else {
    gs <- NULL
    anofilter <- NULL
  }
  
  tx <- get_exons(ensdb, anofilter)
  
  if ("tRNA" %in% names(extra.seqs) & tRNAEnsembleRemove) {
    tx <- tx[grep("tRNA", tx$tx_biotype, invert = TRUE)]
  }
  
  tx$tx_biotype[tx$tx_biotype == "miRNA"] <- ifelse(
    test = width(tx[tx$tx_biotype == "miRNA"]) > 25,
    yes = "miRNA_precursor", no = "miRNA"
  )
  
  if (!is.null(resolveSplicing)) tx <- c(tx, gs)
  
  tx <- tx[, c("tx_id", "tx_biotype", "symbol")]
  colnames(mcols(tx))[colnames(mcols(tx)) == "tx_biotype"] <- "tx_type"
  
  return(tx)
}

#' Get Genes
#'
#' This function retrieves gene information from an EnsDb object.
#'
#' @param ensdb An 'EnsDb' object
#' @param resolveSplicing The biotypes for which to resolve splicing
#'
#' @return A GRanges object containing gene information
#'
#' @importFrom ensembldb genes
get_genes <- function(ensdb, resolveSplicing) {
  gs <- genes(ensdb,
              filter = ~ tx_biotype != resolveSplicing,
              columns = c("tx_id", "gene_id", "gene_biotype", "symbol")
  )
  gs$tx_id <- gs$symbol
  gs$tx_biotype <- gs$gene_biotype
  gs$gene_id <- gs$gene_biotype <- NULL
  gs <- unique(gs)
  return(gs)
}

#' Get Exons
#'
#' This function retrieves exon information from an EnsDb object.
#'
#' @param ensdb An 'EnsDb' object
#' @param anofilter A filter for annotations
#'
#' @return A GRanges object containing exon information
#'
#' @importFrom ensembldb exonsBy
get_exons <- function(ensdb, anofilter) {
  tx <- exonsBy(ensdb,
                filter = anofilter,
                columns = c("tx_id", "tx_biotype", "symbol")
  )
  tx <- unlist(tx)
  return(tx)
}

#' Process Extra Sequences
#'
#' This function processes extra sequences and adds them to the transcript information.
#'
#' @param tx A GRanges object containing transcript information
#' @param extra.seqs Extra transcript sequences to include
#'
#' @return An updated GRanges object with extra sequences included
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Biostrings DNAStringSet
process_extra_sequences <- function(tx, extra.seqs) {
  if (!is.null(extra.seqs)) {
    extra.seqs <- prepare_extra_seqs(extra.seqs)
    m <- mcols(extra.seqs)
    
    if (is.null(m$tx_id)) m$tx_id <- names(extra.seqs)
    names(extra.seqs) <- paste0("pseudoChr_", names(extra.seqs))
    stopifnot(all(c("tx_id", "tx_biotype") %in% colnames(mcols(extra.seqs))))
    
    if (is.null(m$symbol)) m$symbol <- m$tx_id
    gr <- GRanges(names(extra.seqs), IRanges(1L, width = nchar(extra.seqs)),
                  strand = "+", tx_id = m$tx_id, tx_biotype = m$tx_biotype,
                  symbol = m$symbol
    )
    tx <- c(tx, gr)
  }
  return(tx)
}

#' Prepare Extra Sequences
#'
#' This function prepares extra sequences for inclusion in the transcript information.
#'
#' @param extra.seqs Extra transcript sequences to include
#'
#' @return A DNAStringSet object containing prepared extra sequences
#'
#' @importFrom Biostrings DNAStringSet
#' @importFrom dplyr bind_rows
prepare_extra_seqs <- function(extra.seqs) {
  stopifnot(is.list(extra.seqs) || is(extra.seqs, "DNAStringSet"))
  
  if (is.list(extra.seqs)) {
    m <- bind_rows(
      lapply(extra.seqs, FUN = function(x) {
        data.frame(
          row.names = names(x), tx_id = names(x),
          seq = as.character(x)
        )
      }),
      .id = "tx_biotype"
    )
    
    extra.seqs <- DNAStringSet(m$seq)
    m$seq <- NULL
    names(extra.seqs) <- row.names(m)
    mcols(extra.seqs) <- m
  }
  
  return(extra.seqs)
}

#' Process Extra GRanges
#'
#' This function processes extra GRanges and adds them to the transcript information.
#'
#' @param tx A GRanges object containing transcript information
#' @param extra.gr A list of GRanges or GRangesList giving the coordinates of extra features
#'
#' @return An updated GRanges object with extra GRanges included
#'
#' @importFrom IRanges overlapsAny
process_extra_granges <- function(tx, extra.gr) {
  if (length(extra.gr) > 0) {
    extra.gr1 <- process_extra_gr_list(extra.gr)
    tx2 <- tx[!overlapsAny(tx, do.call(c, extra.gr1))]
    tx <- c(tx2, do.call(c, extra.gr1))
  }
  return(tx)
}

#' Process Extra GRanges List
#'
#' This function processes a list of extra GRanges.
#'
#' @param extra.gr A list of GRanges or GRangesList
#'
#' @return A list of processed GRanges
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom ensembldb seqlevelsStyle
process_extra_gr_list <- function(extra.gr) {
  extra.gr1 <- lapply(names(extra.gr), FUN = function(x) {
    gr <- extra.gr[[x]]
    seqlevelsStyle(gr) <- "ensembl"
    
    if (is(gr, "GRangesList")) {
      stopifnot(!is.null(names(gr)))
      if (!("tx_type" %in% colnames(gr@unlistData))) {
        gr@unlistData$tx_biotype <- x
      }
      if (!("symbol" %in% colnames(gr@unlistData))) {
        gr@unlistData$symbol <- gr@unlistData$tx_id
      }
      gr@unlistData <- gr@unlistData[, c("tx_id", "tx_type", "symbol")]
    } else {
      stopifnot("tx_id" %in% colnames(mcols(gr)))
      if (!("tx_type" %in% colnames(mcols(gr)))) mcols(gr)$tx_biotype <- x
      if (!("symbol" %in% colnames(mcols(gr)))) mcols(gr)$symbol <- x
      if (!is.null(gr$type) && any(gr$type == "exon")) gr <- gr$exon
      mcols(gr) <- mcols(gr)[, c("tx_id", "tx_type", "symbol")]
    }
    return(gr)
  })
  return(extra.gr1)
}

#' Finalize Transcripts
#'
#' This function finalizes the transcript information, including processing miRNA transcripts.
#'
#' @param tx A GRanges object containing transcript information
#' @param clusterMiRNA Whether to cluster miRNA annotations
#'
#' @return A finalized GRanges object containing transcript information
finalize_transcripts <- function(tx, clusterMiRNA) {
  colnames(mcols(tx))[colnames(mcols(tx)) == "tx_type"] <- "tx_biotype"
  names(tx) <- NULL
  
  tx2 <- tx[grep("miRNA", tx$tx_biotype, invert = TRUE)]
  tx3 <- tx[grep("miRNA", tx$tx_biotype)]
  tx4 <- process_mir_transcripts(tx3)
  
  tx <- Reduce(c, list(tx2, tx3, tx4))
  tx <- tx[, c("tx_id", "tx_biotype", "symbol")]
  
  if (clusterMiRNA) tx <- miRNAcluster(tx)
  
  return(tx)
}

#' Process miRNA Transcripts
#'
#' This function processes miRNA transcripts.
#'
#' @param tx3 A GRanges object containing miRNA transcript information
#'
#' @return A processed GRanges object containing miRNA transcript information
process_mir_transcripts <- function(tx3) {
  tx4 <- tx3[startsWith(tx3$symbol, "Mir")]
  tx3 <- tx3[!startsWith(tx3$symbol, "Mir")]
  
  tx4$symbol1 <- paste(tx4$symbol, tx4$tx_id, sep = "_")
  tx4$tx_id1 <- tx4$symbol
  tx4$tx_id1 <- gsub("Mir", "miR-", tx4$tx_id1)
  tx4$symbol <- tx4$symbol1
  tx4$tx_id <- tx4$tx_id1
  
  tx4 <- tx4[, c("tx_id", "tx_biotype", "symbol")]
  
  return(tx4)
}

#' Save Features
#'
#' This function saves the processed features to a file.
#'
#' @param tx A GRanges object containing transcript information
#' @param output_dir The directory to save the features file
#'
#' @return The path to the saved features file
save_features <- function(tx, output_dir) {
  anno.out <- file.path(output_dir, "features.rds")
  saveRDS(tx, file = anno.out)
  message("Features saved in \n", anno.out)
  return(anno.out)
}

#' Process Genome
#'
#' This function processes the genome data, including handling extra sequences.
#'
#' @param genome A genome fasta, BSGenome or TwoBit object
#' @param ensdb An 'EnsDb' object
#' @param extra.seqs Extra transcript sequences to include
#' @param output_dir The directory to save the processed genome file
#'
#' @return The path to the processed genome file
#'
#' @importFrom ensembldb getGenomeTwoBitFile
#' @importFrom Biostrings getSeq writeXStringSet
#' @importFrom rtracklayer import
#' @importFrom R.utils gzip
process_genome <- function(genome, ensdb, extra.seqs, output_dir) {
  if (is.null(genome)) genome <- getSeq(getGenomeTwoBitFile(ensdb))
  if (is.character(genome) && length(genome) != 1) {
    stop("`genome` should be the path to a fasta(.gz) file, or a BSgenome/TwoBitFile object")
  }
  
  if (is(genome, "TwoBitFile") || (length(genome) == 1 && grepl(pattern = "\\.2bit", x = genome))) {
    genome <- import(genome)
  }
  
  isCompressed <- !is.character(genome) || grepl("\\.gz$", genome)
  genome.out <- paste0("customGenome.fasta", ifelse(isCompressed, ".gz", ""))
  genome.out <- file.path(output_dir, genome.out)
  
  save_genome(genome, extra.seqs, genome.out, isCompressed)
  
  return(genome.out)
}

#' Save Genome
#'
#' This function saves the processed genome data to a file.
#'
#' @param genome A genome object
#' @param extra.seqs Extra transcript sequences to include
#' @param genome.out The output file path for the genome
#' @param isCompressed Whether the output should be compressed
#'
#' @importFrom Biostrings writeXStringSet
#' @importFrom R.utils gzip
save_genome <- function(genome, extra.seqs, genome.out, isCompressed) {
  if (is.null(extra.seqs)) {
    if (!is.character(genome)) {
      writeXStringSet(genome, genome.out, compress = TRUE)
    } else {
      file.copy(genome, genome.out)
      if (!isCompressed) gzip(genome.out)
    }
  } else {
    if (!is.character(genome)) {
      genome <- c(genome, extra.seqs)
      writeXStringSet(genome, genome.out, compress = TRUE)
    } else {
      file.copy(genome, genome.out)
      writeXStringSet(extra.seqs,
                      filepath = genome.out,
                      append = TRUE
      )
      if (!isCompressed) gzip(genome.out)
    }
  }
  message(
    "Genome including eventual extra chromosomes was saved in:\n",
    genome.out
  )
}

#' Build Index
#'
#' This function builds the genome index using Rsubread.
#'
#' @param genome The path to the genome file
#' @param output_dir The directory to save the index
#' @param ... Additional arguments passed to `Rsubread::buildindex`
#'
#' @importFrom Rsubread buildindex
build_index <- function(genome, output_dir, ...) {
  message("Now building the index...")
  buildindex(
    basename = paste0(output_dir, "/customGenome"),
    reference = genome,
    ...
  )
}

#' Cluster miRNA Annotations
#'
#' This function clusters miRNA annotations.
#'
#' @param tx A GRanges object containing transcript information
#'
#' @return A GRanges object with clustered miRNA annotations
#'
#' @importFrom GenomicRanges reduce
miRNAcluster <- function(tx) {
  miRNA <- tx[tx$tx_biotype %in% c("miRNA", "miRNA_precursor")]
  other <- tx[!tx$tx_biotype %in% c("miRNA", "miRNA_precursor")]
  
  miRNA_reduced <- reduce(miRNA, min.gapwidth = 0L)
  mcols(miRNA_reduced) <- DataFrame(
    tx_id = paste(miRNA_reduced$tx_id, collapse = ","),
    tx_biotype = "miRNA_cluster",
    symbol = paste(miRNA_reduced$symbol, collapse = ",")
  )
  
  tx <- c(other, miRNA_reduced)
  return(tx)
}

# Export only the main function
#' @export
# prepareAnnotation

# Do not export helper functions

# 
# db_mmu <- getDB(tRNA_includeMt = FALSE)
# 
# mm10_annoprep <- prepareAnnotation(
#   ensdb = db_mmu$ensdb,
#   output_dir = "./pr1/",
#   extra.gr = list(piRNA = db_mmu$piRNA_GR, miRNA = db_mmu$miRNA_GR),
#   extra.seqs = list(rRNA = db_mmu$rRNA_fa, tRNA = db_mmu$tRNA_fa),
#   resolveSplicing = NULL,
#   rules = defaultAssignRules(),
#   tRNAEnsembleRemove = TRUE,
#   clusterMiRNA = TRUE
# )
# 
# 
# mm10_annoprep2 <- prepareAnnotation2(
#   ensdb = db_mmu$ensdb,
#   output_dir = "./pr2/",
#   extra.gr = list(piRNA = db_mmu$piRNA_GR, miRNA = db_mmu$miRNA_GR),
#   extra.seqs = list(rRNA = db_mmu$rRNA_fa, tRNA = db_mmu$tRNA_fa),
#   resolveSplicing = NULL,
#   rules = defaultAssignRules(),
#   tRNAEnsembleRemove = TRUE,
#   clusterMiRNA = TRUE
# )
