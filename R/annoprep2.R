#' #' prepareAnnotation
#' #'
#' #' @param ensdb An 'EnsDb' object, as obtained from AnnotationHub
#' #' @param genome A genome fasta, BSGenome or TwoBit object. If missing, will be
#' #' fetched using `ensdb`
#' #' @param output_dir The path where to save the files
#' #' @param extra.gr A list of GRanges or GRangesList giving the coordinates of
#' #' extra features. If a GRanges containing a 'type' column, exons will be used
#' #' and split by transcript. The objects should minimally contain a `tx_id`
#' #' field; if `tx_biotype` is missing and `extra.gr` is named, the names will be
#' #' used as biotype. `symbol` can also be provided.
#' #' @param extra.seqs Extra transcript sequences to include. Can either be a
#' #' named list of named characters (in which case the names of the list are used
#' #' as biotype, and the names of the character vector as tx_id), or a
#' #' `DNAStringSet` object with `mcols` columns 'tx_id' and 'tx_biotype' (`symbol`
#' #' can also be provided).
#' #' @param resolveSplicing The biotypes (in the `ensdb`) for which to resolve
#' #' splicing. If NULL, will resolve all but the low-priority (i.e. long RNA)
#' #' transcripts.
#' #' @param rules Rules for reads assignment
#' #' @param tRNAEnsembleRemove removes the tRNA annotations from Ensembl database
#' #' only if names of extra.seqs contains tRNA
#' #' @param ... passed to `Rsubread::buildindex`
#' #'
#' #' @return Produces a Rsubread index and returns a GRanglesList of features
#' #'
#' #' @import S4Vectors
#' #' @importFrom GenomicRanges GRanges
#' #' @importFrom IRanges IRanges overlapsAny
#' #' @importFrom Biostrings DNAStringSet getSeq writeXStringSet readDNAStringSet readRNAStringSet
#' #' @importFrom Rsubread buildindex
#' #' @importFrom dplyr bind_rows
#' #' @importFrom ensembldb genes exonsBy seqlevelsStyle getGenomeFaFile
#' #' @importFrom R.utils gzip gunzip downloadFile
#' #' @importFrom rtracklayer import import.chain liftOver
#' #' @importFrom plyr ldply
#' #' @importFrom tRNAdbImport import.mttRNAdb
#' #' @importFrom data.table fread
#' #' @importFrom AnnotationHub AnnotationHub query
#' #' @import httr
#' #' @export
#' prepareAnnotation <- function(ensdb, genome = NULL, output_dir = "",
#'                               extra.gr = list(), extra.seqs = NULL,
#'                               resolveSplicing = NULL,
#'                               rules = defaultAssignRules(),
#'                               tRNAEnsembleRemove = TRUE,
#'                               clusterMiRNA = TRUE,
#'                               ...) {
#' 
#'   if (!dir.exists(output_dir)) dir.create(output_dir)
#' 
#'   # If splicing is not to be resolved
#'   if (is.null(resolveSplicing)) {
#'     p <- rules$priorities
#'     resolveSplicing <- names(p)[p >= 0]
#'   }
#' 
#'   # If splicing is to be resolved
#'   if (!is.null(resolveSplicing)) {
#'     gs <- genes(ensdb,
#'       filter = ~ tx_biotype != resolveSplicing,
#'       columns = c("tx_id", "gene_id", "gene_biotype", "symbol")
#'     )
#'     gs$tx_id <- gs$symbol
#'     gs$tx_biotype <- gs$gene_biotype
#'     gs$gene_id <- gs$gene_biotype <- NULL
#'     gs <- unique(gs)
#'     anofilter <- ~ tx_biotype == resolveSplicing
#'   } else {
#'     anofilter <- NULL
#'   }
#' 
#'   tx <- exonsBy(ensdb,
#'     filter = anofilter,
#'     columns = c("tx_id", "tx_biotype", "symbol")
#'   )
#' 
#'   tx <- unlist(tx)
#' 
#'   if ("tRNA" %in% names(extra.seqs) & tRNAEnsembleRemove) {
#'     tx <- tx[grep("tRNA", tx$tx_biotype, invert = T)]
#'   }
#' 
#'   tx$tx_biotype[tx$tx_biotype == "miRNA"] <- ifelse(
#'     test = width(tx[tx$tx_biotype == "miRNA"]) > 25,
#'     yes = "miRNA_precursor", no = "miRNA"
#'   )
#' 
#' 
#'   # If additional sequence files are provided
#'   if (!is.null(extra.seqs)) {
#'     stopifnot(is.list(extra.seqs) || is(extra.seqs, "DNAStringSet"))
#' 
#' 
#'     # If additional sequence files are provided as a list
#'     if (is.list(extra.seqs)) {
#'       m <- bind_rows(
#'         lapply(extra.seqs, FUN = function(x) {
#'           data.frame(
#'             row.names = names(x), tx_id = names(x),
#'             seq = as.character(x)
#'           )
#'         }),
#'         .id = "tx_biotype"
#'       )
#' 
#'       extra.seqs <- DNAStringSet(m$seq)
#'       m$seq <- NULL
#'       names(extra.seqs) <- row.names(m)
#'       mcols(extra.seqs) <- m
#'     }
#'     m <- mcols(extra.seqs)
#' 
#'     if (is.null(m$tx_id)) m$tx_id <- names(extra.seqs)
#'     names(extra.seqs) <- paste0("pseudoChr_", names(extra.seqs))
#'     stopifnot(all(c("tx_id", "tx_biotype") %in% colnames(mcols(extra.seqs))))
#' 
#'     if (is.null(m$symbol)) m$symbol <- m$tx_id
#'     gr <- GRanges(names(extra.seqs), IRanges(1L, width = nchar(extra.seqs)),
#'       strand = "+", tx_id = m$tx_id, tx_biotype = m$tx_biotype,
#'       symbol = m$symbol
#'     )
#'   }
#' 
#' 
#'   if (!is.null(resolveSplicing)) tx <- c(tx, gs)
#' 
#'   if (!is.null(extra.seqs)) tx <- c(tx, gr)
#' 
#'   tx <- tx[, c("tx_id", "tx_biotype", "symbol")]
#'   colnames(mcols(tx))[colnames(mcols(tx)) == "tx_biotype"] <- "tx_type"
#' 
#'   if (length(extra.gr) > 0) {
#'     extra.gr1 <- lapply(names(extra.gr), FUN = function(x) {
#'       gr <- extra.gr[[x]]
#'       seqlevelsStyle(gr) <- "ensembl"
#' 
#'       if (is(gr, "GRangesList")) {
#'         stopifnot(!is.null(names(gr)))
#'         if (!("tx_type" %in% colnames(gr@unlistData))) {
#'           gr@unlistData$tx_biotype <- x
#'         }
#'         if (!("symbol" %in% colnames(gr@unlistData))) {
#'           gr@unlistData$symbol <- gr@unlistData$tx_id
#'         }
#'         gr@unlistData <- gr@unlistData[, c("tx_id", "tx_type", "symbol")]
#'       } else {
#'         stopifnot("tx_id" %in% colnames(mcols(gr)))
#'         if (!("tx_type" %in% colnames(mcols(gr)))) mcols(gr)$tx_biotype <- x
#'         if (!("symbol" %in% colnames(mcols(gr)))) mcols(gr)$symbol <- x
#'         if (!is.null(gr$type) && any(gr$type == "exon")) gr <- gr$exon
#'         mcols(gr) <- mcols(gr)[, c("tx_id", "tx_type", "symbol")]
#'         # gr <- split(gr, gr$tx_id)
#'       }
#'       return(gr)
#'     })
#' 
#'     tx2 <- tx[!overlapsAny(tx, do.call(c, extra.gr1))]
#'     tx <- c(tx2, do.call(c, extra.gr1))
#'     anno.out <- file.path(output_dir, "features.rds")
#'   }
#' 
#'   # tx[grep(pattern = "let-7f-1|Mirlet7f-1", x = tx$symbol)]
#' 
#'   colnames(mcols(tx))[colnames(mcols(tx)) == "tx_type"] <- "tx_biotype"
#' 
#'   names(tx) <- NULL
#' 
#'   tx2 <- tx[grep("miRNA", tx$tx_biotype, invert = TRUE)]
#' 
#'   tx3 <- tx[grep("miRNA", tx$tx_biotype)]
#'   tx4 <- tx3[startsWith(tx3$symbol, "Mir")]
#'   tx3 <- tx3[!startsWith(tx3$symbol, "Mir")]
#' 
#'   tx4$symbol1 <- paste(tx4$symbol, tx4$tx_id, sep = "_")
#'   tx4$tx_id1 <- tx4$symbol
#'   tx4$tx_id1 <- gsub("Mir", "miR-", tx4$tx_id1)
#'   tx4$symbol <- tx4$symbol1
#'   tx4$tx_id <- tx4$tx_id1
#' 
#'   tx4 <- tx4[, c("tx_id", "tx_biotype", "symbol")]
#' 
#'   tx <- Reduce(c, list(tx2, tx3, tx4))
#' 
#'   tx <- tx[, c("tx_id", "tx_biotype", "symbol")]
#' 
#'   if (clusterMiRNA) tx <- miRNAcluster(tx)
#' 
#'   # tx <- as(tx, "GRangesList")
#' 
#'   saveRDS(tx, file = anno.out)
#'   message("Features saved in \n", anno.out)
#' 
#'   # If genome is not provided
#'   if (is.null(genome)) genome <- getSeq(getGenomeTwoBitFile(ensdb))
#'   if (is.character(genome) && length(genome) != 1) {
#'     stop(
#'       "`genome` should be the path to a fasta(.gz) file, or a ",
#'       "BSgenome/TwoBitFile object"
#'     )
#'   }
#' 
#'   if (is(genome, "TwoBitFile") || (length(genome) == 1 && grepl(pattern = "\\.2bit", x = genome))) {
#'     genome <- import(genome)
#'   }
#' 
#'   isCompressed <- !is.character(genome) || grepl("\\.gz$", genome)
#'   genome.out <- paste0("customGenome.fasta", ifelse(isCompressed, ".gz", ""))
#'   genome.out <- file.path(output_dir, genome.out)
#' 
#' 
#'   if (is.null(extra.seqs)) { # in case of no extra seq
#'     if (!is.character(genome)) { # if genome is not a file
#'       writeXStringSet(genome, genome.out, compress = TRUE)
#'     } else { # if genome is a file
#'       file.copy(genome, genome.out)
#'       if (!isCompressed) gzip(genome.out)
#'     }
#'   } else { # in case of extra seq
#'     if (!is.character(genome)) { # if genome is not a file
#'       genome <- c(genome, extra.seqs)
#'       writeXStringSet(genome, genome.out, compress = TRUE)
#'     } else { # If genome is a file
#'       file.copy(genome, genome.out)
#'       writeXStringSet(extra.seqs,
#'         filepath = genome.out,
#'         append = TRUE
#'       )
#'       if (!isCompressed) gzip(genome.out)
#'     }
#'   }
#' 
#'   rm(genome)
#'   message(
#'     "Genome including eventual extra chromosomes was saved in:\n",
#'     genome.out, "\nNow building the index..."
#'   )
#' 
#'   buildindex(
#'     basename = paste0(output_dir, "/customGenome"),
#'     reference = ifelse(test = isCompressed,
#'       yes = genome.out,
#'       no = paste0(genome.out, ".gz")
#'     ), ...
#'   )
#'   return(tx)
#' }



#' Obtain miRNA GRanges from mirBase
#'
#' @param sp 3 alphabet species code. Please check `data(species)`
#'
#' @return a list of GRanges annotation of miRNAs, miRNAs DNAStringSet,
#' and miRNA hairpins miRNAs DNAStringSet
#' 
#' @importFrom rtracklayer import
#' @importFrom plyr ldply
#' 
#' @export
#'
#' @examples
#' mir_mouse <- getmiRNA(sp = "mmu")
#' mir_human <- getmiRNA(sp = "hsa")
getmiRNA <- function(sp = "mmu") {
  sp <- tolower(sp)
  link <- paste0("https://mirbase.org/download_version_genome_files/22.1/", sp, ".gff3")
  gr <- rtracklayer::import(link)

  gr$transcript_id <- gr$Name
  idm <- gr$Name
  names(idm) <- gr$ID
  gr$gene_id <- idm[gr$Derives_from]
  w <- which(is.na(gr$gene_id))
  gr$gene_id[w] <- gr$Name[w]
  gr$transcript_type <- gr$type
  gr$transcript_type <- as.character(gr$transcript_type)
  gr$transcript_type[gr$transcript_type == "miRNA_primary_transcript"] <-
    "miRNA_precursor"

  gr$transcript_id <- gsub(
    pattern = paste0(sp, "-"),
    replacement = "",
    x = gr$transcript_id
  )
  gr$gene_id <- gsub(
    pattern = paste0(sp, "-"),
    replacement = "",
    x = gr$gene_id
  )


  df <- data.frame(gr[, c("transcript_id", "gene_id", "ID", "transcript_type")])
  df_sp <- split(df, df$gene_id)

  gtf <- GRanges(
    ldply(
      lapply(df_sp, function(x) {
        id <- x$ID[x$transcript_type == "miRNA_precursor"]
        x$symbol <- paste(id, x$gene_id, sep = "_")
        return(x)
      }),
      .id = NULL
    )
  )


  gtf$transcript_id[gtf$transcript_type == "miRNA_precursor"] <-
    gtf$ID[gtf$transcript_type == "miRNA_precursor"]

  gtf$tx_id <- gtf$transcript_id
  gtf$tx_type <- gtf$transcript_type

  gtf <- gtf[, c("tx_id", "symbol", "tx_type")]

  miRNA <- rtracklayer::import("https://mirbase.org/download/mature.fa",
    format = "fasta"
  )
  miRNA <- miRNA[grep(pattern = paste0("^", sp), x = names(miRNA))]
  names(miRNA) <- sapply(names(miRNA), function(x) {
    gsub(
      pattern = paste0(sp, "-"),
      replacement = "",
      x = strsplit(x, " ")[[1]][1]
    )
  })

  miRNA_h <- rtracklayer::import("https://mirbase.org/download/hairpin.fa",
    format = "fasta"
  )
  miRNA_h <- miRNA_h[grep(pattern = paste0("^", sp), x = names(miRNA_h))]
  names(miRNA_h) <- sapply(names(miRNA_h), function(x) strsplit(x, " ")[[1]][2])

  return(list(gtf = gtf, mature = miRNA, hairpin = miRNA_h))
}


#' Cluster miRNAs located close to each other
#'
#' @param gr GRanges of miRNAs
#' @param minGap minimum gap between the miRNAs for grouping
#'
#' @return A GRanges object
#' 
#' @importFrom GenomicRanges reduce findOverlaps GRanges
#' @importFrom S4Vectors queryHits subjectHits
#' 
#' @export
miRNAcluster <- function(gr, minGap = 10000) {
  gr_miRNA <- gr[grep("miRNA", gr$tx_biotype)]
  gr_o <- gr[grep("miRNA", gr$tx_biotype, invert = T)]

  if (length(gr_o) > 0) gr_o$miRNAcluster <- NA

  rd <- reduce(gr_miRNA, min.gapwidth = minGap)
  hits <- findOverlaps(rd, gr_miRNA)
  idx1 <- subjectHits(hits)
  idx2 <- queryHits(hits)
  r1 <- data.frame(rd[idx2])
  r1$miRNAcluster <- paste0(
    "miRNAcluster_", r1$seqnames, ":",
    r1$start, "-", r1$end, r1$strand
  )
  values <- data.frame(gr_miRNA[idx1])

  df <- GRanges(data.frame(values, miRNAcluster = r1$miRNAcluster))
  grn <- c(df, gr_o)
  return(grn)
}



#' Obtain mitiochondrial tRNAs from mitotRNA batabase
#'
#' @param sp Scientific name of species
#' @param addCCA Whether to add CCA modifications to sequences or not.
#' Default: TRUE
#'
#' @importFrom tRNAdbImport import.mttRNAdb
#' @importFrom plyr ldply
#' @importFrom Biostrings DNAStringSet
#'
#' @export
getMttRNA <- function(sp = "Mus musculus", addCCA = TRUE) {
  tab <- data.frame(
    import.mttRNAdb(organism = sp)
  )[, c("tRNA_length", "tRNA_type", "tRNA_anticodon", "tRNA_seq")]

  tab <- tab[!duplicated(tab), ]

  tab$tRNA_type <- gsub("[0-9]", "", tab$tRNA_type)

  tab_sp <- split(tab, tab$tRNA_type)

  tab <- ldply(
    lapply(tab_sp, function(x) {
      sp_x <- split(x, x$tRNA_anticodon)
      a <- ldply(lapply(sp_x, function(y) {
        if (nrow(y) > 1) {
          y$n <- 1:table(x$tRNA_anticodon)
        } else {
          y$n <- 1
        }
        y
      }), .id = NULL)

      a$name <- paste("mt_tRNA",
        x$tRNA_type,
        x$tRNA_anticodon,
        a$n,
        1,
        sep = "-"
      )
      a[, !colnames(a) %in% "n"]
    }),
    .id = NULL
  )

  if ("His" %in% tab$tRNA_type) {
    his <- tab$tRNA_seq[tab$tRNA_type == "His"]
    his <- sapply(his, function(x) paste0("G", x))
    tab$tRNA_seq[tab$tRNA_type == "His"] <- his
  }

  mt <- tab$tRNA_seq
  names(mt) <- tab$name

  if (addCCA) {
    mt_CCA <- paste0(mt, "CCA")
    # names(mt_CCA) <- paste(names(mt), "CCA", sep = "_")
    names(mt_CCA) <- names(mt)
    mt <- mt_CCA
  }

  return(DNAStringSet(mt))
}



#' Obtain tRNA sequences from the GttRNA database
#'
#' @param sp Species for which tRNA sequences to be obtained. Currently
#' available for mm10, mm39, hg19, and hg38
#' @param mt Whether to also download and include mitochondrial tRNAs from
#' mitotRNA database. Please see also `getMttRNA`
#' @param addCCA Whether to add CCA modifications to sequences or not.
#'
#' @importFrom Biostrings readDNAStringSet DNAStringSet
#' 
#' @export
gettRNA <- function(sp = "mm10", mt = TRUE, addCCA = TRUE) {
  match.arg(sp, c("hg19", "hg38", "mm10", "mm39"))

  url <- "http://gtrnadb.ucsc.edu/genomes/eukaryota/"
  mt_sp <- NULL

  if (sp == "hg19") {
    url <- paste0(url, "Hsapi19/hg19-tRNAs.fa")
    mt_sp <- "Homo sapiens"
  } else if (sp == "hg38") {
    url <- paste0(url, "Hsapi38/hg38-tRNAs.fa")
    mt_sp <- "Homo sapiens"
  } else if (sp == "mm10") {
    url <- paste0(url, "Mmusc10/mm10-tRNAs.fa")
    mt_sp <- "Mus musculus"
  } else if (sp == "mm39") {
    url <- paste0(url, "Mmusc39/mm39-tRNAs.fa")
    mt_sp <- "Mus musculus"
  }
  
  trna <- Biostrings::readDNAStringSet(url)
  names(trna) <- trimws(gsub(
    pattern = ".*_|\\(.*", replacement = "",
    x = names(trna)
  ))

  if (any(grepl("His", names(trna)))) {
    trna <- as.character(trna)

    # https://github.com/junchaoshi/sports1.1#trna_mappingpl

    his <- trna[grepl("His", names(trna))]
    his <- sapply(his, function(x) paste0("G", x))
    trna[names(his)] <- his
    trna <- Biostrings::DNAStringSet(trna)
  }

  if (addCCA) {
    trna_CCA <- paste0(as.character(trna), "CCA")
    # names(trna_CCA) <- paste(names(trna), "CCA", sep = "_")
    names(trna_CCA) <- names(trna)
    trna <- Biostrings::DNAStringSet(trna_CCA)
  }

  if (mt) {
    mt <- getMttRNA(sp = mt_sp, addCCA = addCCA)
    trna <- c(trna, mt)
  }
  return(trna)
}


#' Obtain rRNAs from SILVA database
#'
#' @param sp Scientific names of the species
#' @param release Which release should be used from 128, 132 and 138.1
#'
#' @importFrom data.table fread
#' @importFrom Biostrings readRNAStringSet DNAStringSet
#'
#' @return DNAStringSet of rRNAs
#' @export
getrRNA <- function(sp = "Mus musculus", release = "138.1") {
  match.arg(release, c("128", "132", "138.1"))

  url_files <- paste0("https://ftp.arb-silva.de/release_", release, "/Exports/")
  files <- listFilesFTP(url_files)

  ss_file <- grep(
    pattern = "SSURef_N[R|r]99_tax_silva.fasta.gz$",
    x = files,
    value = TRUE
  )

  if (length(ss_file) == 0) {
    ss_file <- grep(
      pattern = "SSURef_tax_silva.fasta.gz$",
      x = files,
      value = TRUE
    )
  }

  ss <- readRNAStringSet(ss_file)
  ss <- ss[grepl(pattern = sp, x = names(ss))]
  ss <- DNAStringSet(ss)
  names(ss) <- gsub(pattern = "\\..*", replacement = "", x = names(ss))

  ls_file <- grep(
    pattern = "LSURef_N[R|r]99_tax_silva.fasta.gz$",
    x = files,
    value = TRUE
  )
  if (length(ls_file) == 0) {
    ls_file <- grep(
      pattern = "LSURef_tax_silva.fasta.gz$",
      x = files,
      value = TRUE
    )
  }

  ls <- readRNAStringSet(ls_file)
  ls <- ls[grepl(pattern = sp, x = names(ls))]
  ls <- DNAStringSet(ls)
  names(ls) <- gsub(pattern = "\\..*", replacement = "", x = names(ls))

  if (release != "128") {
    url_meta <- paste0(url_files, "full_metadata/")
    meta <- listFilesFTP(url_meta)

    ssm <- grep(
      pattern = "SSURef_N[R|r]99.full_metadata.gz$",
      x = meta,
      value = TRUE
    )

    if (length(ssm) == 0) {
      ssm <- grep(
        pattern = "SSURef.full_metadata.gz$",
        x = meta,
        value = TRUE
      )
    }

    ssm <- fread(input = ssm)
    ssm <- data.frame(ssm[ssm$acc %in% names(ss), ])
    ssm <- ssm[!duplicated(ssm$acc), ]
    rownames(ssm) <- ssm$acc


    lsm <- grep(
      pattern = "LSURef_N[R|r]99.full_metadata.gz$",
      x = meta,
      value = TRUE
    )

    if (length(lsm) == 0) {
      lsm <- grep(
        pattern = "LSURef.full_metadata.gz$",
        x = meta,
        value = TRUE
      )
    }

    lsm <- fread(input = lsm)
    lsm <- data.frame(lsm[lsm$acc %in% names(ls), ])
    lsm <- lsm[!duplicated(lsm$acc), ]
    rownames(lsm) <- lsm$acc

    names(ss) <- paste0("SSU-", names(ss), "(", ssm[names(ss), "product"], ")")
    names(ls) <- paste0("LSU-", names(ls), "(", lsm[names(ls), "product"], ")")
  } else {
    names(ss) <- paste0("SSU-", names(ss))
    names(ls) <- paste0("LSU-", names(ls))
  }

  res <- c(ss, ls)
  names(res) <- gsub(pattern = "\\()", replacement = "", x = names(res))
  names(res) <- make.unique(names(res), sep = "_")
  names(res) <- gsub(pattern = " ", replacement = "_", x = names(res))
  return(res)
}


#' Obtain databases for species (currently supported for mouse only)
#'
#' @param species 3 alphabet code for the species. Default: "mmu"
#' @param genomeVersion genome version to be used for the species.
#' Default: "GRCm38"
#' @param ensemblVer Ensemble version to be used for genome. Default: 102
#' @param tRNA_addCCA Wether to add "CCA" modifications to the tRNA sequences.
#' Default: TRUE
#' @param tRNA_includeMt Whether to include mitochondrial tRNAs. Default: TRUE
#' @param rRNA_release the version of rRNA database to be used. Default: "138.1"
#'
#' @import httr
#' @importFrom AnnotationHub AnnotationHub query
#' @importFrom R.utils gunzip downloadFile
#' @importFrom rtracklayer import.chain liftOver
#' 
#' @return A list of databases (DNAstringSet or GRanges)
#'
#' @export
getDB <- function(species = "mmu", genomeVersion = "GRCm38",
                  ensemblVer = "102",
                  tRNA_addCCA = TRUE, tRNA_includeMt = TRUE,
                  rRNA_release = "138.1") {
  
  # config <- httr::config(ssl_verifypeer = 0, ssl_verifyhost = 0)
  # tt <- options()$timeout
  # options(timeout = 6000)
  options(download.file.method="curl", download.file.extra="-k -L")
  
  if (species == "mmu") {

    # EnsDb
    ah <- AnnotationHub()
    # ensdb <- rev(query(ah, genomeVersion, "Ensdb"))[[1]]
    ensdb <- query(ah, c(genomeVersion, "EnsDb", ensemblVer))[[1]]

    # miRNA
    miRNA <- getmiRNA(sp = species)$gtf

    # tRNA
    if (genomeVersion == "GRCm38") {
      tRNA <- gettRNA(sp = "mm10", addCCA = tRNA_addCCA, mt = tRNA_includeMt)
    } else if (genomeVersion == "GRCm39") {
      tRNA <- gettRNA(sp = "mm39", addCCA = tRNA_addCCA, mt = tRNA_includeMt)
    }

    # piRNA precursors
    data("piRNA_mmu")
    piRNA <- piRNA_mmu
    piRNA <- piRNA[piRNA$type == "transcript"]
    colnames(mcols(piRNA))[5:7] <- c("symbol", "tx_id", "tx_type")
    piRNA <- piRNA[, c("symbol", "tx_id", "tx_type")]

    if (genomeVersion == "GRCm39") {
      R.utils::gunzip(
        R.utils::downloadFile(paste0(
          "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/",
          "liftOver/mm10ToMm39.over.chain.gz"
        ))
      )
      chain <- rtracklayer::import.chain("mm10ToMm39.over.chain")
      miRNA <- unlist(rtracklayer::liftOver(miRNA, chain))
      piRNA <- unlist(rtracklayer::liftOver(piRNA, chain))
      unlink("mm10ToMm39.over.chain")
    }

    # rRNA
    rRNA <- getrRNA(sp = "Mus musculus", release = rRNA_release)

    db <- list(
      ensdb = ensdb,
      miRNA_GR = miRNA,
      tRNA_fa = tRNA,
      piRNA_GR = piRNA,
      rRNA_fa = rRNA
    )
    return(db)
  } else if (species == "hsa") {

    # EnsDb
    ah <- AnnotationHub()
    # ensdb <- rev(query(ah, genomeVersion, "Ensdb"))[[1]]
    ensdb <- query(ah, c(genomeVersion, "EnsDb", ensemblVer))[[1]]

    # miRNA
    miRNA <- getmiRNA(sp = species)$gtf

    # tRNA
    if (genomeVersion == "GRCh38") {
      tRNA <- gettRNA(sp = "hg38", addCCA = tRNA_addCCA, mt = tRNA_includeMt)
    } else if (genomeVersion == "GRCh19") {
      tRNA <- gettRNA(sp = "hg19", addCCA = tRNA_addCCA, mt = tRNA_includeMt)
    }

    # rRNA
    rRNA <- getrRNA(sp = "Homo sapiens", release = rRNA_release)

    db <- list(
      ensdb = ensdb,
      miRNA_GR = miRNA,
      tRNA_fa = tRNA,
      rRNA_fa = rRNA
    )
    
    # config <- httr::config(ssl_verifypeer = 1, ssl_verifyhost = 1)
    # options(timeout = tt)
    
    return(db)
  }
}




# devtools::load_all("./")

# db_hg38 <- getDB(species = "hsa", genomeVersion = "GRCh38", ensemblVer = "102")
#
# hg38_annoprep <- prepareAnnotation(
#   ensdb = db_hg38$ensdb,
#   genome = "/mnt/IM/reference/genome/gencode/hg38/fasta/GRCh38.p13.genome.fa",
#   output_dir = "../../shortRNA_reports/schratt_human/shortRNA/genome/",
#   extra.gr = list(miRNA = db_hg38$miRNA_GR),
#   extra.seqs = list(rRNA = db_hg38$rRNA_fa, tRNA = db_hg38$tRNA_fa),
#   resolveSplicing = NULL,
#   rules = defaultAssignRules(),
#   tRNAEnsembleRemove = TRUE,
#   clusterMiRNA = TRUE
# )
#
#
# db_mmu <- getDB()
#
# ensdb <- db_mmu$ensdb
# genome <- "/mnt/IM/reference/genome/gencode/fasta/GRCm38.p5.genome.fa"
# output_dir <- "/mnt/IM/projects/software/shortRNA/genome"
# extra.gr <- list(piRNA = db_mmu$piRNA_GR, miRNA = db_mmu$miRNA_GR)
# extra.seqs <- list(rRNA = db_mmu$rRNA_fa, tRNA = db_mmu$tRNA_fa)
# resolveSplicing <- NULL
# rules <- defaultAssignRules()
# tRNAEnsembleRemove <- TRUE
# clusterMiRNA <- TRUE




# mm10_annoprep <- prepareAnnotation(
#   ensdb = db_mmu$ensdb,
#   genome = "/mnt/IM/reference/genome/gencode/fasta/GRCm38.p5.genome.fa",
#   output_dir = "./",
#   extra.gr = list(piRNA = db_mmu$piRNA_GR, miRNA = db_mmu$miRNA_GR),
#   extra.seqs = list(rRNA = db_mmu$rRNA_fa, tRNA = db_mmu$tRNA_fa),
#   resolveSplicing = NULL,
#   rules = defaultAssignRules(),
#   tRNAEnsembleRemove = TRUE,
#   clusterMiRNA = TRUE
# )
# 
# 
# # devtools::load_all("./")
# db_mmu <- getDB()
# mm10_annoprep <- prepareAnnotation(
#   ensdb = db_mmu$ensdb,
#   # genome = "/mnt/IM/reference/genome/gencode/fasta/GRCm38.p5.genome.fa",
#   output_dir = "~/Desktop",
#   extra.gr = list(piRNA = db_mmu$piRNA_GR, miRNA = db_mmu$miRNA_GR),
#   extra.seqs = list(rRNA = db_mmu$rRNA_fa, tRNA = db_mmu$tRNA_fa),
#   resolveSplicing = NULL,
#   rules = defaultAssignRules(),
#   tRNAEnsembleRemove = TRUE,
#   clusterMiRNA = TRUE
# )
