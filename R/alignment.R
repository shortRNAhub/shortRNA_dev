#' Build index for Rsubread
#' @author Deepak Tanwar
#' @param fastaGenome path to a gzipped fasta file
#' (optional if want to build genome other than mm9, mm10, hg19 or hg38).
#' @param refGenome mm9/mm10/hg19/hg38.
#' `fasta` file is downloaded from gencode (https://gencodegenes.org)
#' @param gencodeRelease Release of the genome to be downloaded.
#' @param outPath Output direcory where genome is downloaded and index is build.
#' @param basename Basename for the Rsubread index.
#' @param extraFasta `vector`. Additional fasta files location.
#' @param gzExtra if additional fasta files are gzip. Default: FALSE.
#' @param ... Other options for \code{\link[Rsubread]{buildindex}}.
#' @return `index` for aligning fastq files with Rsubread.
#' @importFrom Rsubread buildindex
#' @importFrom R.utils gunzip gzip
#'
#' @examples
#' indexRsubread(
#'   fastaGenome = "../test/reference.fa.gz",
#'   basename = "../test/reference"
#' )
#' @export
indexRsubread <- function(fastaGenome = NULL,
                          refGenome = "mm10",
                          gencodeRelease = "M18",
                          outPath = getwd(),
                          basename = "refRsubread",
                          extraFasta = NULL,
                          gzExtra = FALSE,
                          ...) {
  # If fasta file is provided
  if (!is.null(fastaGenome)) {
    message("...building index for Rsubread...\n")
    Rsubread::buildindex(
      basename = basename,
      reference = R.utils::gunzip(fastaGenome),
      ...
    )
    R.utils::gzip(gsub(pattern = "\\.gz$", "", fastaGenome))
  } else {
    # If not provided
    message(paste0(
      "As the `fasta` file is not provided via fastaGenome, ",
      refGenome,
      " will be dowloaded and index will be build from GENCODE version ",
      gencodeRelease, "."
    ))

    # Various checks
    if (!grepl(pattern = "hg19|hg38|mm9|mm10", x = refGenome)) {
      stop(
        "Please provide path to the `fasta` file!
      This function can only build indexes without a `fasta` file for:
      mm9, mm10, hg19 or hg38, or from a `fasta file`!"
      )
    }

    if (grepl(pattern = "hg[0-9][0-9]", x = refGenome) &
      !grepl(pattern = "[0-9][0-9]", x = gencodeRelease)) {
      stop("Please see: https://www.gencodegenes.org/human/releases.html")
    }

    if (grepl(pattern = "mm[0-9]|mm[0-9][0-9]", x = refGenome) &
      !grepl(pattern = "M[0-9]|M[0-9][0-9]", x = gencodeRelease)) {
      stop("Please see: https://www.gencodegenes.org/mouse/releases.html")
    }

    # Obtaining link to download fasta file
    url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_"
    ftp <- NULL
    if (refGenome == "mm9") {
      url <- paste0(url, "mouse/release_M1/NCBIM37.genome.fa.gz")
      ftp <- url
    } else if (refGenome == "mm10") {
      u <- paste0(url, "mouse/release_", gencodeRelease, "/")
      files <- listFilesFTP(url = u)
      ftp <- files[grep(
        pattern = "GRCm38.p[0-9].genome.fa.gz|GRCm38.p[0-9][0-9].genome.fa.gz",
        x = files
      )]
    } else if (refGenome == "hg19") {
      url <- paste0(url, "human/release_19/GRCh37.p13.genome.fa.gz")
      ftp <- url
    } else if (refGenome == "hg38") {
      u <- paste0(url, "human/release_", gencodeRelease, "/")
      files <- listFilesFTP(url = u)
      ftp <- files[grep(
        pattern = "GRCh38.p[0-9].genome.fa.gz|GRCh38.p[0-9][0-9].genome.fa.gz",
        x = files
      )]
    }

    # Getting name to save fasta file
    str <- strsplit(x = ftp, split = "/")[[1]]
    n <- str[length(str)]
    fa <- paste0(outPath, "/", n)
    message("...downloading file...\n")
    system(paste("mkdir -p", outPath))
    download.file(url = ftp, destfile = fa)
    message("...file downloaded...\n\n...building index for Rsubread...\n")

    genomeFasta <- NULL
    # Extra fasta
    if (!is.null(extraFasta)) {
      if (gzExtra) {
        system(paste0(
          "cat ", fa, " ", paste(extraFasta, collapse = " "),
          " > ", outPath, "/genome.fa.gz"
        ))
        genomeFasta <- paste(outPath, "genome.fa.gz", sep = "/")
      } else {
        tmp <- as.character(sapply(extraFasta, gzip))
        system(paste0(
          "cat ", fa, " ", paste(tmp, collapse = " "),
          " > ", outPath, "/genome.fa.gz"
        ))
        genomeFasta <- paste(outPath, "genome.fa.gz", sep = "/")
      }
    } else {
      genomeFasta <- fa
    }

    # Index building
    Rsubread::buildindex(
      basename = paste(outPath, basename, sep = "/"),
      reference = genomeFasta, ...
    )
  }
}




#' Alignment using Rsubread
#' @author Deepak Tanwar
#'
#' @param fastq character vector of file names.
#' @param fastq2 in case the sequencing is paired-end.
#' @param nThreads numeric value of how many cores to be used for alignment.
#' @param index path to the Rsubread index
#' @param numeric value specifying the maximal number of equally-best mapping
#'  locations that will be reported for a multi-mapping read.
#' @param GTF path to the GTF file.
#' @param outDir output directory. Default: current working directory
#' @param outFile basename of the output file. Default: basename of input file
#' @param mode function to be used from Rsubread package for
#' performing alignment. Default: align
#' @param ... other parameters specific to \code{\link[Rsubread]{align}} or
#' \code{\link[Rsubread]{subjunc}}.
#' @return Stores a `BAM` file.
#'
#' @importFrom Rsubread align subjunc
#' @importFrom parallel detectCores
#'
#' @export
alignShortRNA <- function(fastq,
                          fastq2 = NULL,
                          index,
                          nBestLocations = 100,
                          nThreads = parallel::detectCores(),
                          GTF = NULL,
                          outDir = getwd(),
                          outFile = NULL,
                          mode = "align",
                          ...) {
  gtfOption <- NULL
  if (!is.null(GTF)) {
    gtfOption <- TRUE
  } else {
    gtfOption <- FALSE
  }

  out <- NULL
  if (!is.null(outFile)) {
    out <- paste0(outDir, "/", outFile, ".bam")
  } else {
    o <- gsub(pattern = ".f.*", replacement = "", x = basename(fastq))
    out <- paste0(outDir, "/", o, ".bam")
  }

  align_func <- if (mode == "align") {
    Rsubread::align
  } else if (mode == "subjunc") {
    Rsubread::subjunc
  } else {
    stop("Please provide one option: align or subjunc!")
  }

  align_func(
    index = index, readfile1 = fastq, readfile2 = fastq2,
    nBestLocations = nBestLocations, output_file = out,
    nthreads = nThreads, sortReadsByCoordinates = TRUE,
    annot.ext = GTF, isGTF = gtfOption, ...
  )
}
