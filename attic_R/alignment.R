#' shortRNAexp_align
#'
#' A wrapper for the 2-steps alignment method.
#'
#' @param fasta Path to a fasta file containing the unique reads.
#' @param outputfile Path to the output bam file.
#' @param bowtie1index Path to the base of the bowtie1 index.
#' @param starindex Path to the STAR index folder.
#' @param bowtie1 Executable for bowtie1 (default 'bowtie').
#' @param star Executable for STAR (default 'STAR').
#' @param samtools Executable for samtools (default 'samtools').
#' @param m Maximum alignments for a read to be considered (default 1000).
#' @param nthreads Number of threads for alignment (default 4).
#'
#' @export
shortRNAexp_align <- function(fasta, outputfile, bowtie1index,
                              starindex, bowtie1 = "bowtie",
                              star = "STAR", samtools = "samtools",
                              m = 1000, nthreads = 4) {
  of2 <- gsub(".bam", "", outputfile, fixed = T)
  cmd <- paste0(
    bowtie1, " -p ", nthreads, " -v 0 -S -a --best --strata -m ",
    m, " -f --un ", of2, ".unmapped.fasta ", bowtie1index,
    " ", fasta, " | ", samtools, " view -bh > ", of2, ".unsorted.bam &&
    ", samtools, " sort -@ ", nthreads, " -m 2G ",
    of2, ".unsorted.bam > ", of2, ".perfectMatch.bam &&
    rm ", of2, ".unsorted.bam"
  )
  print(cmd)
  system(cmd)
  cmd <- paste0(
    star, " --genomeDir ", starindex, " --runThreadN ", nthreads,
    " --readFilesIn ", of2,
    ".unmapped.fasta --alignIntronMax 1 --outFilterMultimapNmax ", m,
    " --outSAMattributes NH HI NM --outSAMtype BAM SortedByCoordinate",
    "--outFileNamePrefix ", of2,
    ".imperfect --outSAMprimaryFlag AllBestScore &&",
    samtools, " merge -f ", outputfile, " ", of2, ".perfectMatch.bam ",
    of2, ".imperfect*bam && ", samtools, " index ", outputfile
  )
  print(cmd)
  system(cmd)
}


#' shortRNAexp_parseBam
#'
#' Extracts reads' possible sources from an alignment file.
#'
#' @param bam Path to the sorted bam file.
#' @param elements Path to the genomic elements bed file.
#' @param outputfile Path for the output file.
#' @param shell Shell executable, default 'bash'.
#' @param samtools Executable for samtools (default 'samtools').
#' @param bedtools Executable for bedtools (default 'bedtools').
#'
#' @export
shortRNAexp_parseBam <- function(bam,
                                 elements,
                                 outputfile = NULL,
                                 shell = "bash",
                                 samtools = "samtools",
                                 bedtools = "bedtools") {
  if (is.null(outputfile)) {
    f <- tempfile("srcs")
  } else {
    f <- outputfile
  }
  p <- paste0(path.package("shortRNA"), "/parseBam.sh")
  cmd <- paste0(shell, p, bam, elements, f, samtools, bedtools)
  system(cmd)
  if (is.null(outputfile)) {
    ret <- read.delim(f, header = F, stringsAsFactors = F)
    unlink(f)
    return(ret)
  }
}


#' Build index for Rsubread
#' @author Deepak Tanwar (tanward@ethz.ch)
# @import Rsubread R.utils
#' @param fastaGenome path to a gzipped fasta file
#' (optional if want to build genome other than mm9, mm10, hg19 or hg38).
#' @param refGenome mm9/mm10/hg19/hg38.
#' `fasta` file is downloaded from gencode (https://gencodegenes.org)
#' @param gencodeRelease Release of the genome to be downloaded.
#' @param outPath Output direcory where genome is downloaded and index is build.
#' @param basename Basename for the Rsubread index.
#' @param extraFasta `vector`. Additional fasta files location.
#' @param gzExtra if additional fasta files are gzip. Default: FALSE.
#' @param ... Other options for \code{\link[Rsubread:buildindex]}.
#'
#' @return **index** for aligning fastq files with Rsubread.
#' 
# @import Rsubread R.utils
#' 
#' \dontrun{
#' @sexamples
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
    buildindex(basename = basename, reference = gunzip(fastaGenome), ...)
    gzip(gsub(pattern = "\\.gz$", "", fastaGenome))
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
    buildindex(
      basename = paste(outPath, basename, sep = "/"),
      reference = genomeFasta, ...
    )
  }
}




#' Alignment using Rsubread
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
# @import Rsubread parallel
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
#' @param ... other parameters specific to \code{\link[Rsubread:align]} or
#' \code{\link[Rsubread:subjunc]}.
#' @return Stores a `BAM` file.
#' 
# @import Rsubread parallel
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

  if (mode == "align") {
    align(
      index = index, readfile1 = fastq, readfile2 = fastq2,
      nBestLocations = nBestLocations, output_file = out,
      nthreads = nThreads, sortReadsByCoordinates = TRUE,
      annot.ext = GTF, isGTF = gtfOption, ...
    )
  } else if (mode == "subjunc") {
    subjunc(
      index = index, readfile1 = fastq, readfile2 = fastq2,
      nBestLocations = nBestLocations, output_file = out,
      nthreads = nThreads, sortReadsByCoordinates = TRUE,
      annot.ext = GTF, isGTF = gtfOption, ...
    )
  } else {
    stop("Please provide one option: align or subjunc!")
  }
}
