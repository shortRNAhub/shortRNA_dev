#' getAnnotationFiles
#'
#' Download (and untar) the files necessary to build a short RNA annotation, and returns a named list indicating the different feature files. Is called by `prepareAnnotation`.
#'
#' @param org An organism assembly to use. Currently supported are: mm10, hg38.
#' @param destination The destination folder for the files (defaults to current working directory)
#'
#' @return A list.
#' 
# @import Biostrings rtracklayer R.utils
#'
#' @export
getAnnotationFiles <- function(org = "mm10", destination = getwd()) {
  if (!dir.exists(destination)) dir.create(destination)
  
  # Mouse
  if (org == "mm10") {
    ## list of files to be doanloaded
    ff <- list(
      c(
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/gencode.vM17.tRNAs.gtf.gz",
        file.path(destination, "gencode.tRNAs.gtf.gz")
      ),
      c(
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf.gz",
        file.path(destination, "gencode.features.gtf.gz")
      ),
      c(
        "ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz",
        file.path(destination, "mirbase.fa.gz")
      ),
      c(
        "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/mmu.gff3",
        file.path(destination, "mirbase.gtf")
      ),
      c(
        "http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz",
        file.path(destination, "gtRNAdb.tar.gz")
      )
    )

    srcs <- sapply(ff, FUN = function(x) {
      download.file(x[[1]], ifelse(length(x) > 1, x[[2]], NULL))
      return(x[[1]])
    })
    untar(file.path(destination, "gtRNAdb.tar.gz"), exdir = destination)

    fa <- readRNAStringSet(file.path(destination, "mirbase.fa.gz"), "fasta")
    fa <- fa[grep("^mmu-", names(fa))]
    writeXStringSet(fa, file.path(destination, "mm10.mature.miRNAs.fa"), format = "fasta")
    unlink(file.path(destination, "mature.miRNA.fa.gz"))

    return(list(
      features.gtf = file.path(destination, "gencode.features.gtf.gz"),
      mirbase.fa = file.path(destination, "mm10.mature.miRNAs.fa"),
      mirbase.gtf = file.path(destination, "mirbase.gtf"),
      pi_precursors.gtf = file.path(path.package("shortRNA"), "extdata", "mm10.piRNA.precursors.gtf"),
      gtRNAdb.fa = file.path(destination, "mm10-tRNAs.fa"),
      gtRNAdb.bed = file.path(destination, "mm10-tRNAs.bed"),
      tRNA.gtf = file.path(destination, "gencode.tRNAs.gtf.gz"),
      equivalences = NA_character_,
      srcs = srcs
    ))
  }

  # Rat
  if (org == "rn6") {
    ff <- list(
      c(
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/gencode.vM17.tRNAs.gtf.gz",
        file.path(destination, "gencode.tRNAs.gtf.gz")
      ),
      c(
        "ftp://ftp.ensembl.org/pub/release-94/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.94.gtf.gz",
        file.path(destination, "features.gtf.gz")
      ),
      c(
        "ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz",
        file.path(destination, "mirbase.fa.gz")
      ),
      c(
        "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/rno.gff3",
        file.path(destination, "mirbase.gtf")
      ),
      c(
        "http://gtrnadb.ucsc.edu/genomes/eukaryota/Rnorv5/rn5-tRNAs.tar.gz",
        file.path(destination, "rn5.gtRNAdb.tar.gz")
      ),
      c(
        "http://hgdownload.soe.ucsc.edu/goldenPath/rn5/liftOver/rn5ToRn6.over.chain.gz",
        file.path(destination, "rn5ToRn6.chain.gz")
      )
    )

    srcs <- sapply(ff, FUN = function(x) {
      download.file(x[[1]], ifelse(length(x) > 1, x[[2]], NULL))
      return(x[[1]])
    })

    untar(file.path(destination, "rn5.gtRNAdb.tar.gz"), exdir = destination)
    gunzip(file.path(destination, "rn5ToRn6.chain.gz"))
    chain <- import.chain(file.path(destination, "rn5ToRn6.chain"))
    b <- import.bed("rn5-tRNAs.bed")
    b <- unlist(liftOver(b, chain))
    b <- b[which(width(b) == width(b$thick)), ]
    b$thick <- b@ranges
    export.bed(b, file.path(destination, "rn6.gtRNAdb.bed"))
    b2 <- bed2gtf(b)
    b2$gene_id <- strsplit(b2$transcript_id, "-", fixed = T)[[1]][[2]]
    b2$gene_name <- paste0("tRNA_", b$gene_id)
    b2$transcript_type <- "tRNA"
    export.gff(b2, file.path(destination, "rn6.tRNA.gtf"))

    fa <- readDNAStringSet(file.path(destination, "rn5-tRNAs.fa"), "fasta")
    names(fa) <- sapply(names(fa), FUN = function(x) {
      strsplit(gsub("Rattus_norvegicus_", "", x), " ", fixed = T)[[1]][[1]]
    })
    fa <- fa[which(names(fa) %in% as.character(b$name))]
    writeXStringSet(fa, file.path(destination, "rn6-tRNAs.fa"), format = "fasta")

    fa <- readRNAStringSet(file.path(destination, "mirbase.fa.gz"), "fasta")
    fa <- fa[grep("^rno-", names(fa))]
    writeXStringSet(fa, file.path(destination, "rno.mature.miRNAs.fa"), format = "fasta")
    unlink(file.path(destination, "mature.miRNA.fa.gz"))
    unlink(file.path(destination, "rn5ToRn6.chain"))
    return(list(
      features.gtf = file.path(destination, "features.gtf.gz"),
      mirbase.fa = file.path(destination, "rno.mature.miRNAs.fa"),
      mirbase.gtf = file.path(destination, "mirbase.gtf"),
      gtRNAdb.fa = file.path(destination, "rn6-tRNAs.fa"),
      gtRNAdb.bed = file.path(destination, "rn6.gtRNAdb.bed"),
      tRNA.gtf = file.path(destination, "gencode.tRNAs.gtf.gz"),
      equivalences = NA_character_,
      srcs = srcs
    ))
  }

  # Human
  if (org == "hg38") {
    ff <- list( # c("http://www.smallrnagroup.uni-mainz.de/piRNAclusterDB/homo_sapiens/proTRAC_normal_ovary_generic/BED_files.zip", file.path(destination,"piRNAs.clusters.zip")),
      c(
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.tRNAs.gtf.gz",
        file.path(destination, "gencode.tRNAs.gtf.gz")
      ),
      c(
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.chr_patch_hapl_scaff.annotation.gtf.gz",
        file.path(destination, "gencode.features.gtf.gz")
      ),
      c(
        "ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz",
        file.path(destination, "mirbase.fa.gz")
      ),
      c(
        "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3",
        file.path(destination, "mirbase.gtf")
      ),
      c(
        "http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz",
        file.path(destination, "gtRNAdb.tar.gz")
      )
    )

    srcs <- sapply(ff, FUN = function(x) {
      download.file(x[[1]], ifelse(length(x) > 1, x[[2]], NULL))
      return(x[[1]])
    })

    # unzip(file.path(destination,"piRNAs.clusters.zip"),exdir=destination)
    # e <- read.delim(file.path(destination,"homo_sapiens","proTRAC_normal_ovary_generic","piRNA_clusters.BED"),header=F,skip=1,stringsAsFactors=F)
    # e$name <- sapply(as.character(e$V4),FUN=function(x){ x <- strsplit(x,"|",fixed=T)[[1]]; paste0("locus",x[length(x)])})
    # e$score <- 0
    # e$strand <- "*"
    # e$strand[grep("mono:minus",e$V4)] <- "-"
    # e$strand[grep("mono:plus",e$V4)] <- "+"
    # write.table(e[,-4],file.path(destination,"piRNA.clusters.bed"),col.names=F,row.names=F,sep="\t",quote=F)
    # unlink(file.path(destination,"homo_sapiens/proTRAC_normal_ovary_generic"),T)

    untar(file.path(destination, "gtRNAdb.tar.gz"), exdir = destination)
    fa <- readRNAStringSet(file.path(destination, "mirbase.fa.gz"), "fasta")
    fa <- fa[grep("^hsa-", names(fa))]
    writeXStringSet(fa, file.path(destination, "hg38.mature.miRNAs.fa"),
      format = "fasta"
    )
    unlink(file.path(destination, "mature.miRNA.fa.gz"))
    return(list(
      features.gtf = file.path(destination, "gencode.features.gtf.gz"),
      mirbase.fa = file.path(destination, "hg38.mature.miRNAs.fa"),
      mirbase.gtf = file.path(destination, "mirbase.gtf"),
      #  pi_precursors=file.path(destination,"piRNA.clusters.bed"),
      gtRNAdb.fa = file.path(destination, "hg38-tRNAs.fa"),
      gtRNAdb.bed = file.path(destination, "hg38-tRNAs.bed"),
      tRNA.gtf = file.path(destination, "gencode.tRNAs.gtf.gz"),
      equivalences = NA_character_,
      srcs = srcs
    ))
  }
  stop("Unknown organism!")
}

#' prepareAnnotation
#'
#' Prepares an annotation GRanges and modified miRNA/tRNA sequences for a given organism.
#' Either `filelist` or `species` should be provided.
#'
#' @param species The species for automatic fetching and creation of the annotation.
#' See `?getAnnotationFiles` for the supported species/assemblies.
#' @param filelist A list of named feature files, as produced by 'getAnnotationFiles', indicating the necessary files.
#' The following slots are required:
#' 'features.gtf','features.gtf','mirbase.gtf','gtRNAdb.fa','gtRNAdb.bed','tRNA.gtf'; the 'pi_precursors' is optional.
#' @param destination The destination folder for the files (defaults to current working directory)
#' @param description An optional character vector containing a description of the annotation (default empty)
#' @param resolveSplicing A character vector indicating for which RNA types to
#' consider splicing when calculating overlaps and position within transcript,
#' or a logical indicating whether to do so for all transcripts. By default,
#' this is done only for short RNAs (see `getDefaultShortTypes()`), since it is
#' often unimportant whether reads coming from long RNAs stem from precursors or
#' from processed RNAs. Doing it for all transcripts will increase computing
#' time, while disabling it entirely will yield imperfect overlaps and positions
#' within spliced transcripts.
#'
#' @return A GRangesList object, and produces modified fasta files.
#' 
# @import Biostrings rtracklayer GenomicRanges tools data.table
#'
#' @export
prepareAnnotation <- function(species = NULL,
                              filelist = NULL,
                              destination = getwd(),
                              description = "",
                              includeIsomirs = TRUE,
                              resolveSplicing = defaultAssignRules()$highPriorityTypes,
                              unresolvedLevel = "gene") {
  if ((is.null(species) && is.null(filelist)) ||
    (!is.null(species) && !is.null(filelist))
  ) {
    stop("Exactly one of `species` or `filelist` must be given.")
  }
  dir.create(file.path(destination), showWarnings = FALSE, recursive = TRUE)
  if (!is.null(species)) {
    filelist <- getAnnotationFiles(species, destination = destination)
  }
  ff <- c(
    "features.gtf", "mirbase.gtf", "gtRNAdb.fa", "gtRNAdb.bed",
    "tRNA.gtf", "pi_precursors.gtf"
  )
  if (!is.list(filelist) || !all(ff %in% names(filelist))) {
    stop(paste(
      "`filelist` should be a list with the following elements:",
      paste(ff, collapse = ", ")
    ))
  }
  
  # miRNAs
  if (includeIsomirs) {
    mirs <- .buildIsomirs(filelist$mirbase.fa)
  } else {
    mirs <- DNAStringSet(complement(readRNAStringSet(filelist$mirbase.fa)))
  }
  names(mirs) <- paste0("pseudoChr_", gsub("-", "_", names(mirs)))
  mn <- sapply(names(mirs), FUN = function(x) {
    strsplit(x, " ", fixed = T)[[1]][[1]]
  })
  writeXStringSet(mirs, file.path(destination, "miRNAs.modified.fa"),
    format = "fasta"
  )
  mn2 <- sapply(mn, FUN = function(x) {
    strsplit(x, ".", fixed = T)[[1]][[1]]
  })
  mn3 <- sapply(mn2, FUN = function(x) {
    x <- strsplit(x, "-", fixed = T)[[1]]
    paste(x[-length(x)], collapse = "-")
  })
  le <- length(mn)
  mirs <- GRanges(mn, IRanges(rep(1, le), width(mirs)), rep("+", le),
    transcript_id = mn2, gene_id = mn3, transcript_type = "miRNA"
  )
  mirbase.gtf <- .prepare.mirbase.gtf(filelist$mirbase.gtf)
  mirs <- suppressWarnings(c(mirbase.gtf, mirs))
  levels(mirs$transcript_type)[which(
    levels(mirs$transcript_type) ==
      "miRNA_primary_transcript"
  )] <- "miRNA_hairpin"

  # tRNAs
  trnas <- .preparetRNAsequences(filelist$gtRNAdb.fa, filelist$gtRNAdb.bed)
  names(trnas) <- paste0("pseudoChr_", gsub("-", "_", names(trnas)))
  writeXStringSet(trnas, file.path(destination, "tRNAs.modified.fa"),
    format = "fasta"
  )
  le <- length(trnas)
  mn <- sapply(names(trnas), FUN = function(x) {
    x <- strsplit(x, " ", fixed = T)[[1]][[1]]
    x <- strsplit(x, "_", fixed = T)[[1]]
    x[length(x)]
  })
  trnas <- GRanges(names(trnas), IRanges(rep(1, le), width(trnas)),
    rep("+", le),
    transcript_id = names(trnas), gene_id = mn,
    transcript_type = "tRNA"
  )
  tr <- import.bed(filelist$gtRNAdb.bed)
  tr <- bed2gtf(tr)
  tr <- tr[which(tr$type == "exon"), ]
  tr$gene_id <- sapply(tr$transcript_id, FUN = function(x) {
    x <- strsplit(x, "-", fixed = T)[[1]]
    paste(x[1:3], collapse = "-")
  })
  tr$transcript_type <- "tRNA"
  g <- import.gff(filelist$tRNA.gtf)
  # add only tRNAs not in gtRNAdb and mark them as putative
  ov <- findOverlaps(g, tr)
  g$transcript_type <- "putative_tRNA"
  g$transcript_type[which(g$gene_type == "Pseudo_tRNA")] <- "pseudo_tRNA"
  g$gene_id <- g$transcript_name
  if (length(ov@from) > 0) g <- g[-unique(ov@from), ]
  trnas <- suppressWarnings(c(
    tr[, c("transcript_id", "gene_id", "transcript_type")],
    trnas,
    g[, c("transcript_id", "gene_id", "transcript_type")]
  ))

  # piRNAs
  if ("pi_precursors.gtf" %in% names(filelist) &&
    !is.na(filelist$pi_precursors.gtf)) {
    pirnas <- import(filelist$pi_precursors.gtf)
    if ("type" %in% names(pirnas@elementMetadata)) {
      pirnas <- pirnas[
        which(pirnas$type == "transcript"),
        c("transcript_id", "gene_id", "transcript_type")
      ]
    } else {
      pirnas$transcript_id <- pirnas$name
      pirnas$gene_id <- pirnas$name
      pirnas$transcript_type <- "piRNA_precursor"
      pirnas <- pirnas[, c("transcript_id", "gene_id", "transcript_type")]
    }
  }

  # GENES
  meta <- c("transcript_id", "gene_id", "transcript_type")
  g <- import.gff(filelist$features.gtf)
  g <- g[which(g$type %in% c("gene", "transcript", "exon")), ]
  if (!("transcript_type" %in% colnames(mcols(g))) &&
    ("transcript_biotype" %in% colnames(mcols(g)))) {
    g$transcript_type <- g$transcript_biotype
    g$transcript_biotype <- NULL
  }
  mcols(g)$transcript_type <- as.character(mcols(g)$transcript_type)

  w <- which(g$type == "transcript" & g$transcript_type == "miRNA")
  if (length(w) > 0) {
    # remove miRNAs that are already in mirbase
    gmirs <- g[w, ]
    ov <- findOverlaps(gmirs, mirbase.gtf)
    gmirs <- gmirs[unique(ov@from), ]
    g <- g[which(!(g$transcript_id %in% gmirs$transcript_id | g$gene_id %in% gmirs$gene_id)), ]
  }

  w <- which(g$type == "transcript" & g$transcript_type == "tRNA")
  if (length(w) > 0) {
    # remove tRNAs that are in previous annotations
    gtr <- g[w, ]
    ov <- findOverlaps(gtr, trnas)
    gtr2 <- gtr[unique(ov@from), ]
    if (length(ov) > 0) gtr <- gtr[-unique(ov@from), ]
    g <- g[which(!(g$transcript_id %in% gtr2$transcript_id | g$gene_id %in% gtr2$gene_id)), ]
    g[which(g$transcript_id %in% gtr$transcript_id), "transcript_type"] <- "putative_tRNA"
  }

  if (all(resolveSplicing == FALSE)) {
    g <- .toGRangesList(g[which(g$type == "transcript"), meta])
  } else {
    if (!all(resolveSplicing == TRUE)) {
      if (unresolvedLevel == "gene") {
        g1 <- g[which(!(g$transcript_type %in% resolveSplicing) &
          g$type == "gene"), meta]
        g1$transcript_id <- g1$gene_id
        g1$transcript_type <- "long_RNA"
      } else {
        g1 <- g[which(!(g$transcript_type %in% resolveSplicing) &
          g$type == "transcript"), meta]
      }
      g1 <- .toGRangesList(g1)
      g2 <- g[which(g$transcript_type %in% resolveSplicing &
        g$type == "exon"), meta]
      g2 <- split(g2, g2$transcript_id)
      g <- c(g1, g2)
      rm(g1, g2)
    } else {
      g <- g[which(g$type == "exon"), meta]
      g <- split(g, g$transcript_id)
    }
  }

  if ("equivalences" %in% names(filelist) && !is.na(filelist$equivalences)) {
    gid <- sapply(as.character(g$gene_id), FUN = function(x) {
      x <- strsplit(x, ".", fixed = T)[[1]]
      if (length(x) == 1) {
        return(x)
      }
      paste(x[-length(x)], collapse = ".")
    })
    eqv <- read.delim(filelist$equivalences, header = F, row.names = 1)
    w <- which(gid %in% row.names(eqv))
    g$gene_id[w] <- eqv[gid[w], 1]
  }

  g <- g[, meta]
  if ("pi_precursors.gtf" %in% names(filelist) &&
    !is.na(filelist$pi_precursors.gtf)) {
    g2 <- suppressWarnings(c(mirs, trnas, pirnas))
  } else {
    g2 <- suppressWarnings(c(mirs, trnas))
  }
  g2 <- .toGRangesList(g2[, meta])
  g <- suppressWarnings(c(g, g2))
  rm(g2)

  # set transcript-level metadata from exon-level one where missing:
  ma <- as.data.frame(mcols(g, use.names = F), stringsAsFactors = F)
  w <- which(is.na(ma$gene_id))
  if (length(w) > 0) {
    ma[w, ] <- t(sapply(g[w], FUN = function(x) {
      as.matrix(mcols(x))[1, ]
    }))
  }
  ma$transcript_type <- factor(ma$transcript_type)

  # remove exon-level metadata:
  tx <- factor(rep(1:length(g), elementNROWS(g)))
  g <- GenomicRanges::split(unlist(g)[, c()], tx)
  names(g) <- NULL
  mcols(g) <- ma

  if ("srcs" %in% names(filelist)) {
    srcs <- filelist$srcs
  } else {
    srcs <- NA
  }
  attr(g, "description") <- description
  attr(g, "srcs") <- srcs
  attr(g, "md5") <- md5sum(
    as.character(filelist[which(names(filelist) != "srcs" & !is.na(filelist))])
  )
  message(paste0("Saved modified sequences in:
", file.path(destination, "tRNAs.modified.fa"), "
", file.path(destination, "miRNAs.modified.fa")))
  g
}

.prepare.mirbase.gtf <- function(g) {
  if (is.character(g)) g <- import.gff(g)
  g$transcript_id <- g$Name
  idm <- g$Name
  names(idm) <- g$ID
  g$gene_id <- idm[g$Derives_from]
  w <- which(is.na(g$gene_id))
  g$gene_id[w] <- g$Name[w]
  g$transcript_type <- g$type
  g[, c("transcript_id", "gene_id", "transcript_type")]
}

# where mature.fa is either a RNAStringSet or a path to a fasta file of RNA sequences
#' .buildIsomirs
#'
#' @param mature.fa 
#'
#' @return
#' @export
#' 
# @import Biostrings
#'
#' @examples
.buildIsomirs <- function(mature.fa) {
  if (is.character(mature.fa)) mature.fa <- readRNAStringSet(mature.fa)
  fa <- DNAStringSet(complement(mature.fa))
  a <- DNAStringSet(unlist(lapply(fa, FUN = function(x) {
    sapply(c("A", "C", "G", "T"), y = x, FUN = function(y, x) {
      xscat(y, x)
    })
  })))
  names(a) <- paste(sapply(names(a), FUN = function(x) {
    strsplit(x, " ", fixed = T)[[1]][[1]]
  }), rep(c("A", "C", "G", "T"), length(a) / 4), sep = ".")
  c(fa, a)
}

#' .preparetRNAsequences
#'
#' @param cDNA 
#' @param bed 
#'
#' @return
#' @export
#' 
# @import Biostrings
#'
#' @examples
.preparetRNAsequences <- function(cDNA, bed = NULL) {
  if (is.character(cDNA)) cDNA <- readDNAStringSet(cDNA)
  if (is.null(bed)) {
    warning(
      "No tRNA bed file provided -
      assuming that the fasta is stranded and spliced."
    )
  } else {
    # a bed file is provided, so we will splice the sequences
    # and invert if on negative strand
    if (is.character(bed)) bed <- import.bed(bed)
    names(cDNA) <- sapply(names(cDNA), FUN = function(x) {
      x <- strsplit(x, " ", fixed = T)[[1]][[1]]
      x <- strsplit(x, "_", fixed = T)[[1]]
      x[length(x)]
    })
    names(bed) <- bed$name
    i <- intersect(names(cDNA), names(bed))
    bed <- bed[i, ]
    cDNA <- cDNA[i]
    # splicing
    for (i in which(sapply(bed$blocks, length) > 1)) {
      cDNA[[i]] <- DNAString(paste(apply(as.data.frame(bed$blocks[[i]]), 1,
        fa = as.character(cDNA[[i]]),
        FUN = function(x, fa) {
          substr(fa, x[[1]], x[[2]])
        }
      ), collapse = ""))
    }
    # inverting sequence
    for (i in which(as.character(strand(bed)) == "-")) {
      cDNA[[i]] <- reverse(cDNA[[i]])
    }
  }
  # we add the CCA tail
  x <- xscat(cDNA, "CCA")
  names(x) <- names(cDNA)
  x
}




#' bed2gtf
#'
#' Converts a bed-format GRanges object to a GTF format GRanges object.
#'
#' @param b A GRanges object in bed format
#' (as produced by `rtracklayer::import.bed`)
#' @param fieldsToInclude The additional columns of `b` to include
#' in the results. By default, only the transcript_id is taken
#'  (from the `name` column)
#'
#' @return A GRanges object that can be saved with `rtracklayer::export.gff`
#' @export
bed2gtf <- function(b, fieldsToInclude = c("transcript_id")) {
  b$type <- "transcript"
  if (!("transcript_id" %in% names(b@elementMetadata))) b$transcript_id <- b$name
  nbl <- sapply(b$blocks, FUN = length)
  bl <- unlist(b$blocks)
  b2 <- GRanges(
    seqnames = rep(seqnames(b), nbl),
    ranges = IRanges(
      start(bl) + rep(start(b), nbl) - 1,
      end(bl) + rep(start(b), nbl) - 1
    ),
    strand = rep(as.character(strand(b)), nbl),
    type = rep("exon", length(bl))
  )
  for (f in fieldsToInclude) {
    b2@elementMetadata[[f]] <- rep(b@elementMetadata[[f]], nbl)
  }
  b2 <- sort(c(b[, c("type", fieldsToInclude)], b2))
}



#' buildSrcTree
#'
#' @param a An annotation GRanges object
#' @param metatypes An optional named vector linking different values of
#' `transcript_type` in `a` to meta-types. If omitted, will be built
#' using pre-built grouping scheme.
#' @param clusters An optional named vector linking different values of
#' `transcript_id` or `gene_id` in `a` to clusters.
#'
#' @return a `data.`
#' 
# @import data.tree
#' 
#' @export
buildSrcTree <- function(a, metatypes = NULL, clusters = NULL) {
  a <- a@elementMetadata
  a <- a[, intersect(colnames(a), c(
    "transcript_id", "gene_id",
    "transcript_type", "metatype", "cluster"
  ))]
  a <- as.data.frame(a[!duplicated(a), ])
  if (!("metatype" %in% colnames(a))) {
    a$metatype <- as.character(a$transcript_type)
    if (!is.null(metatypes)) {
      a$metatype[which(a$metatype %in% names(metatypes))] <-
        metatypes[a$metatype]
    } else {
      a$metatype[which(a$metatype == "miRNA_primary_transcript")] <- "miRNA"
      a$metatype[which(a$metatype == "pseudo_tRNA")] <- "tRNA"
      a$metatype[which(a$metatype == "Mt_tRNA")] <- "tRNA"
      stypes <- c(
        "miRNA", "snoRNA", "snRNA", "tRNA", "pseudo_tRNA",
        "putative_tRNA", "Mt_tRNA", "scaRNA", "sRNA"
      )
      a$metatype[which(a$metatype %in% stypes)] <- "other_short_RNA"
      a$metatype[which(!(a$metatype %in% stypes))] <- "other"
    }
    a$metatype <- as.factor(a$metatype)
  }
  if (!("cluster" %in% colnames(a)) &&
    !is.null(clusters) && length(clusters) > 0) {
    a$cluster <- clusters[a$gene_id]
    a$cluster[which(a$transcript_id %in% names(clusters))] <-
      clusters[a$transcript_id]
  }
  fields <- c(
    "metatype", "cluster", "gene_id",
    "transcript_type", "transcript_id"
  )
  fields <- fields[which(fields %in% colnames(a))]
  message(paste(
    "Building source tree with fields",
    paste(fields, collapse = ", ")
  ))
  a$pathString <- apply(a[, fields], 1, FUN = function(x) {
    x <- as.character(x)
    x <- unique(x[which(!is.na(x))])
    paste("RNA", paste(x, collapse = "/"), sep = "/")
  })
  as.Node(a)
}

.toGRangesList <- function(x, name.field = NULL) {
  if ("transcript_type" %in% names(x@elementMetadata)) {
    x@elementMetadata$transcript_type <- as.character(
      x@elementMetadata$transcript_type
    )
  }
  y <- as(x, "GRangesList")
  if (!is.null(name.field)) {
    names(y) <- mcols(x, use.names = F)[[name.field]]
  }
  return(y)
}
