#' Get Genome Versions
#'
#' This function returns a list of genome versions for different species.
#'
#' @return A list containing genome versions for mouse (mmu), human (hsa), and rat (rno).
#' @export
#'
#' @examples
#' genome_versions <- get_genome_versions()
#' mouse_versions <- genome_versions$mmu
get_genome_versions <- function() {
  genome_versions <- list(
    mmu = c("GRCm38", "mm10", "GRCm39", "mm39"),
    hsa = c("GRCh38", "hg38"),
    rno = c("mRatBN7.2", "rn7")
  )
  return(genome_versions)
}

#' Get miRNA Data
#'
#' This function retrieves miRNA data for a specified species, genome version, and miRBase version.
#'
#' @param sp A character string specifying the species. Must be one of "mmu" (mouse), "hsa" (human), or "rno" (rat).
#' @param genome_version A character string specifying the genome version. Must be a valid version for the specified species.
#' @param mirbase_version A character string specifying the miRBase version. If NULL, it will be inferred from the genome version.
#'
#' @return A list containing:
#'   \item{gtf}{A GRanges object with miRNA genomic features}
#'   \item{mature}{A DNAStringSet object with mature miRNA sequences}
#'   \item{hairpin}{A DNAStringSet object with miRNA hairpin sequences}
#'   \item{genome_version}{The specified genome version}
#'   \item{mirbase_version}{The used miRBase version}
#'
#' @import rtracklayer
#' @importFrom plyr ldply
#'
#' @export
#'
#' @examples
#' mouse_mirna <- getmiRNA("mmu", "GRCm39")
#' human_mirna <- getmiRNA("hsa", "GRCh38", mirbase_version = "22.1")
#' rat_mirna <- getmiRNA("rno", "rn7", mirbase_version = "21")
getmiRNA <- function(sp = "mmu", mirbase_version = "CURRENT") {
  # Validate input
  valid_species <- c("mmu", "hsa", "rno")
  if (!sp %in% valid_species) {
    stop("Invalid species code. Use 'mmu' for mouse, 'hsa' for human, or 'rno' for rat.")
  }
  
  # If mirbase_version is not provided, infer it from the genome version
  if (is.null(mirbase_version)) {
    mirbase_version <- switch(genome_version,
                              "GRCm38" = "21",
                              "mm10" = "21",
                              "GRCm39" = "22.1",
                              "mm39" = "22.1",
                              "GRCh38" = "22.1",
                              "hg38" = "22.1",
                              "mRatBN7.2" = "22.1",
                              "rn7" = "22.1",
                              "22.1")  # Default to latest version if no match
  }
  
  # Validate mirbase_version (this is a simplified check, you might want to expand it)
  valid_mirbase_versions <- c("21", "22", "22.1")
  if (!mirbase_version %in% valid_mirbase_versions) {
    stop(paste("Invalid miRBase version. Valid versions are:", paste(valid_mirbase_versions, collapse = ", ")))
  }
  
  # Construct URLs based on miRBase version
  base_url <- paste0("https://mirbase.org/download_version_genome_files/", mirbase_version, "/")
  gff3_url <- paste0(base_url, sp, ".gff3")
  mature_url <- paste0("https://mirbase.org/download_version/", mirbase_version, "/mature.fa")
  hairpin_url <- paste0("https://mirbase.org/download_version/", mirbase_version, "/hairpin.fa")
  
  # Import GFF3 data
  gr <- rtracklayer::import(gff3_url)
  gr$transcript_id <- gr$Name
  idm <- gr$Name
  names(idm) <- gr$ID
  gr$gene_id <- idm[gr$Derives_from]
  w <- which(is.na(gr$gene_id))
  gr$gene_id[w] <- gr$Name[w]
  gr$transcript_type <- gr$type
  gr$transcript_type <- as.character(gr$transcript_type)
  gr$transcript_type[gr$transcript_type == "miRNA_primary_transcript"] <- "miRNA_precursor"
  gr$transcript_id <- gsub(pattern = paste0(sp, "-"), replacement = "", x = gr$transcript_id)
  gr$gene_id <- gsub(pattern = paste0(sp, "-"), replacement = "", x = gr$gene_id)
  
  df <- data.frame(gr[, c("transcript_id", "gene_id", "ID", "transcript_type")])
  df_sp <- split(df, df$gene_id)
  gtf <- GRanges(
    plyr::ldply(
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
  
  # Import mature miRNA sequences
  miRNA <- rtracklayer::import(mature_url, format = "fasta")
  miRNA <- miRNA[grep(pattern = paste0("^", sp), x = names(miRNA))]
  names(miRNA) <- sapply(names(miRNA), function(x) {
    gsub(
      pattern = paste0(sp, "-"),
      replacement = "",
      x = strsplit(x, " ")[[1]][1]
    )
  })
  
  # Import hairpin sequences
  miRNA_h <- rtracklayer::import(hairpin_url, format = "fasta")
  miRNA_h <- miRNA_h[grep(pattern = paste0("^", sp), x = names(miRNA_h))]
  names(miRNA_h) <- sapply(names(miRNA_h), function(x) strsplit(x, " ")[[1]][2])
  
  return(list(gtf = gtf, mature = miRNA, hairpin = miRNA_h, genome_version = genome_version, mirbase_version = mirbase_version))
}