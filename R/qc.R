#' QC of fastq files
#' 
#' @author Deepak Tanwar
#' 
#' @importFrom Rfastp rfastp
#' @importFrom parallel detectCores
#' 
#' @param file character vector of file names
#' @param outDir output directory. Default: current working directory
#' @param name character vector of output file names.
#' @param adapterSeq Adapter sequence for read1
#' @param adaptersFa Adapter sequence fasta file
#' @inheritParams Rfastp::rfastp
#' @param ... Additional arguments passed to \code{\link[Rfastp]{rfastp}}
#' 
#' @return A list of results and location to the JSON results file
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#' # Analysis
#' qc_res <- qcFastq(file = fq1, outDir = outDir)
#' # Output
#' qc_res
#' 
#' @export
qcFastq <- function(file,
                    outDir = ".",
                    name = gsub(
                      pattern = "\\.fq|\\.fq.gz|\\.fastq|\\.fastq.gz",
                      replacement = "", x = basename(file)
                    ),
                    adapterSeq = "auto",
                    adaptersFa = "",
                    minLength = 15,
                    nThread = parallel::detectCores(), ...) {
  if (adapterSeq == "auto") {
    message(
      "The adapter sequence is automatically recognised.
    If the adapter sequence is already known, it is better to specify it.\n\n"
    )
  }
  res <- rfastp(
    read1 = file,
    outputFastq = paste(outDir, name, sep = "/"),
    thread = nThread,
    cutMeanQual = 20,
    cutFrontMeanQual = 20,
    cutTailMeanQual = 20,
    cutLowQualTail = T,
    cutLowQualFront = T,
    qualityFiltering = T,
    maxNfilter = 5,
    qualityFilterPhred = 20,
    qualityFilterPercent = 40,
    lengthFiltering = T,
    minReadLength = minLength,
    maxIndexMismatch = 1,
    adapterTrimming = TRUE,
    adapterSequenceRead1 = adapterSeq,
    adapterFasta = adaptersFa,
    ...
  )

  json <- list.files(path = outDir, pattern = name, full.names = T)
  json <- grep(pattern = "json", x = json, value = TRUE)

  ret <- list(
    result = res,
    json = json
  )
  return(ret)
}



#' Make general information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @importFrom rjson fromJSON
#'
#' @param json A JSON file with QC information from \code{\link{qcFastq}}
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qcFastq(file = fq1, outDir = outDir)
#' df <- makeGeneralDF(json = qc_res$json)
#'
#' # Output
#' df
#' @export
makeGeneralDF <- function(json) {
  js <- fromJSON(file = json)

  pkg <- as.character(packageVersion("Rfastp"))
  seq <- js$read1_before_filtering$total_cycles
  seq_type <- ifelse(test = any(grepl(pattern = "read2", x = names(js))),
    yes = "paired-end", no = "single-end"
  )
  cycle <- paste0(seq_type, " (", seq, " cycles)")

  mean_len_bf <- js$summary$before_filtering$read1_mean_length
  mean_len_af <- js$summary$after_filtering$read1_mean_length
  dup_rate <- ifelse(test = seq_type == "paired-end",
    yes = paste0(round(x = (js$duplication$rate * 100), digits = 2), "%"),
    no = paste0(
      round(x = (js$duplication$rate * 100), digits = 2),
      "% (may be overestimated since this is SE data)"
    )
  )

  df <- data.frame(
    test = c(
      "Rfastp version",
      "Sequencing",
      "Mean length before filtering",
      "Mean length after filtering",
      "Duplication rate"
    ),
    value = c(
      pkg,
      cycle,
      mean_len_bf,
      mean_len_af,
      dup_rate
    )
  )
  return(df)
}



#' Make before QC information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @importFrom rjson fromJSON
#'
#' @param json A JSON file with QC information from \code{\link{qcFastq}}
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qcFastq(file = fq1, outDir = outDir)
#' df <- makeBeforeFiltDF(json = qc_res$json)
#'
#' # Output
#' df
#' @export
makeBeforeFiltDF <- function(json) {
  js <- fromJSON(file = json)
  tot_reads <- js$read1_before_filtering$total_reads
  tot_bases <- js$read1_before_filtering$total_bases
  q_20_bases <- js$read1_before_filtering$q20_bases
  q_20_rate <- js$summary$before_filtering$q20_rate * 100
  q_30_bases <- js$read1_before_filtering$q30_bases
  q_30_rate <- js$summary$before_filtering$q30_rate * 100
  gc <- paste0(round(
    x = js$summary$before_filtering$gc_content * 100,
    digits = 2
  ), "%")

  df <- data.frame(
    test = c(
      "Total reads",
      "Total bases",
      "Q20 bases",
      "Q30 bases",
      "GC content"
    ),
    value = c(
      tot_reads,
      tot_bases,
      paste0(q_20_bases, " (", q_20_rate, "%)"),
      paste0(q_30_bases, " (", q_30_rate, "%)"),
      gc
    )
  )
  return(df)
}



#' Make after QC information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @importFrom rjson fromJSON
#'
#' @param json A JSON file with QC information from \code{\link{qcFastq}}
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qcFastq(file = fq1, outDir = outDir)
#' df <- makeAfterFiltDF(json = qc_res$json)
#'
#' # Output
#' df
#' @export
makeAfterFiltDF <- function(json) {
  js <- fromJSON(file = json)
  tot_reads <- js$read1_after_filtering$total_reads
  tot_bases <- js$read1_after_filtering$total_bases
  q_20_bases <- js$read1_after_filtering$q20_bases
  q_20_rate <- js$summary$after_filtering$q20_rate * 100
  q_30_bases <- js$read1_after_filtering$q30_bases
  q_30_rate <- js$summary$after_filtering$q30_rate * 100
  gc <- paste0(
    round(x = js$summary$after_filtering$gc_content * 100, digits = 2),
    "%"
  )

  df <- data.frame(
    test = c(
      "Total reads",
      "Total bases",
      "Q20 bases",
      "Q30 bases",
      "GC content"
    ),
    value = c(
      tot_reads,
      tot_bases,
      paste0(q_20_bases, " (", q_20_rate, "%)"),
      paste0(q_30_bases, " (", q_30_rate, "%)"),
      gc
    )
  )
  return(df)
}



#' Make filtering information information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @importFrom rjson fromJSON
#'
#' @param json A JSON file with QC information from \code{\link{qcFastq}}
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qcFastq(file = fq1, outDir = outDir)
#' df <- makeFiltDF(json = qc_res$json)
#'
#' # Output
#' df
#' @export
makeFiltDF <- function(json) {
  js <- fromJSON(file = json)

  reads <- js$summary$before_filtering$total_reads
  read_pass <- js$filtering_result$passed_filter_reads
  read_pass_pct <- round((read_pass / reads) * 100, 2)
  read_low_Q <- js$filtering_result$low_quality_reads
  read_low_Q_pct <- round((read_low_Q / reads) * 100, 2)
  read_with_N <- js$filtering_result$too_many_N_reads
  read_with_N_pct <- round((read_with_N / reads) * 100, 2)
  read_short <- js$filtering_result$too_short_reads
  read_short_pct <- round((read_short / reads) * 100, 2)
  read_long <- js$filtering_result$too_long_reads
  read_long_pct <- round((read_pass / reads) * 100, 2)

  df <- data.frame(
    test = c(
      "Reads passed filters",
      "Reads with low quality",
      "Reads with too many N",
      "Reads too short",
      "Reads too long"
    ),
    value = c(
      paste0(read_pass, " (", read_pass_pct, "%)"),
      paste0(read_low_Q, " (", read_low_Q_pct, "%)"),
      paste0(read_with_N, " (", read_with_N_pct, "%)"),
      paste0(read_short, " (", read_short_pct, "%)"),
      paste0(read_long, " (", read_long_pct, "%)")
    )
  )
  return(df)
}



#' Reads duplication plot
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @importFrom dplyr `%>%`
#' @importFrom plotly plot_ly add_trace add_lines layout
#'
#' @param json A JSON file with QC information from \code{\link{qcFastq}}
#' @param df A data frame
#' @param dupRate Column with duplication rate in df
#' @param GCpercent Column with GC percentage in df
#' @param dupLevel Column with duplication levels in df
#' 
#' @return An interactive plot
#'
#' @examples
#' # Input
#' df <- data.frame(
#'   dupRate = rev(sort(c(70, 10, 1:5, rep(0, 13)))),
#'   dupLevel = 1:20
#' )
#' df$GCpercent <- c(runif(n = 5, min = 51, max = 55), rep(0, 15))
#'
#' # Analysis
#' ## Plot with a data frame
#' p <- readDuplicationPlot(
#'   df = df,
#'   dupRateCol = "dupRate",
#'   GCpercentCol = NULL,
#'   dupLevelCol = "dupLevel"
#' )
#'
#' ## Plot with a data frame with GC
#' p_GC <- readDuplicationPlot(
#'   df = df,
#'   dupRateCol = "dupRate",
#'   GCpercentCol = "GCpercent",
#'   dupLevelCol = "dupLevel"
#' )
#'
#' ## Plot with a json file
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#' qc_res <- qc_fastq(file = fq1, outDir = outDir)
#' p_js <- readDuplicationPlot(json = qc_res$json)
#'
#' # Output
#' p
#' p_GC
#' p_js
#' @export
readDuplicationPlot <- function(json = NULL,
                                df = NULL,
                                dupRate = "dupRate",
                                GCpercent = NULL,
                                dupLevel = "dupLevel") {
  if (!is.null(json)) {
    df <- dupPlotDF(json)
    GCpercent <- "GCpercent"
  }

  p <- plot_ly(df, x = df[, dupLevel]) %>%
    add_trace(
      y = df[, dupRate],
      type = "bar",
      name = "Read percent",
      hoverinfo = "text+x",
      hovertext = paste(
        "<b>", df[, dupRate], "</b>"
      )
    )

  if (!is.null(GCpercent)) {
    df[, GCpercent][df[, GCpercent] == 0] <- ""
    p <- p %>%
      add_lines(
        y = df[, GCpercent],
        name = "Mean GC ratio (%)",
        line = list(color = "red", width = 3),
        hoverinfo = "text+x",
        hovertext = paste(
          "<b>", df[, GCpercent], "</b>"
        )
      ) %>%
      layout(
        title = paste0(
          "Mean duplication rate: ",
          round(mean(df[, dupRate]), 2),
          "%"
        ),
        xaxis = list(title = "Duplication level"),
        yaxis = list(title = "Reads percent & GC ratio (%)"),
        hovermode = "x closest"
      )
  } else {
    p <- p %>%
      layout(
        title = paste0(
          "Mean duplication rate: ",
          round(mean(df[, dupRate]), 2),
          "%"
        ),
        xaxis = list(title = "Duplication level"),
        yaxis = list(title = "Reads percent"),
        hovermode = "x closest"
      )
  }

  return(p)
}

# df <- data.frame(
#   dupRate = rev(sort(c(70, 10, 1:5, rep(0, 13)))),
#   dupLevel = 1:20
# )
# df$GCpercent <- c(runif(n = 5, min = 51, max = 55), rep(0, 15))
# 
# p <- readDuplicationPlot(
#   df = df,
#   dupRate = "dupRate",
#   GCpercent = NULL,
#   dupLevel = "dupLevel"
# )
# 
# p_GC <- readDuplicationPlot(
#   df = df,
#   dupRate = "dupRate",
#   GCpercent = "GCpercent",
#   dupLevel = "dupLevel"
# )
# 


#' Reads duplication data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @importFrom rjson fromJSON
#'
#' @param json A JSON file with QC information from \code{\link{qcFastq}}
#'
#' @return An interactive plot
#'
#' #' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qcFastq(file = fq1, outDir = outDir)
#' df <- dupPlotDF(qc_res$json)
#'
#' # Output
#' df
dupPlotDF <- function(json) {

  js <- fromJSON(file = json)

  bars <- js$duplication$histogram
  gc <- js$duplication$mean_gc
  xaxis <- 1:length(bars)

  data <- data.frame(
    dupRate = round(bars / sum(bars) * 100, 2),
    GCpercent = round(gc * 100, 2),
    dupLevel = xaxis
  )

  return(data)
}



#' Read quality score plot
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @importFrom dplyr `%>%`
#' @importFrom plotly plot_ly add_trace add_lines layout
#' @importFrom rcartocolor carto_pal
#'
#' @param json A JSON file with QC information from \code{\link{qcFastq}}
#' @param df A data frame
#' @param bpPos Column name with base-pair position
#' @param Aa Column with Adenine quality score
#' @param Tt Column with Thymine quality score
#' @param Gg Column with Guanine quality score
#' @param Cc Column with Cytosine quality score
#' @param meanQual Column with `mean` quality score
#'
#' @return An interactive plot
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qc_fastq(file = fq1, outDir = outDir)
#' p <- readsQualityPlot(json = qc_res$json)
#'
#' # Output
#' p
#' @export
readsQualityPlot <- function(json = NULL,
                             df = NULL,
                             bpPos = "bpPos",
                             Aa = "Aa",
                             Tt = "Tt",
                             Gg = "Gg",
                             Cc = "Cc",
                             meanQual = NULL) {

  if (!is.null(json)) {
    df <- readQualDF(json)
    meanQual <- "meanQual"
  }

  if (is.null(meanQual)) {
    df$meanQual <- rowMeans(df[, c(Aa, Tt, Gg, Cc)])
    meanQual <- "meanQual"
  }

  df$Q20 <- 20
  df$Q28 <- 8
  df$Q50 <- 18

  cols <- carto_pal(n = 10, name = "Safe")[c(2, 1, 3, 4, 10)]

  plot_ly(df, x = df[, bpPos]) %>%
    add_trace(
      y = df[, "Q20"], name = "Not good", fillcolor = "rgba(255,0,0, 0.2)",
      type = "scatter", mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_trace(
      y = df[, "Q28"], name = "OK", fillcolor = "rgba(255,255,0,0.2)",
      type = "scatter", mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_trace(
      y = df[, "Q50"], name = "Good", fillcolor = "rgba(0,128,0,0.2)",
      type = "scatter", mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_lines(
      y = df[, Aa],
      name = "A",
      line = list(color = cols[1], width = 2),
      hoverinfo = "text",
      hovertext = paste("<b>", df[, Aa], "</b>")
    ) %>%
    add_lines(
      y = df[, Tt],
      name = "T",
      line = list(color = cols[2], width = 2),
      hoverinfo = "text",
      hovertext = paste("<b>", df[, Tt], "</b>")
    ) %>%
    add_lines(
      y = df[, Gg],
      name = "G",
      line = list(color = cols[3], width = 2),
      hoverinfo = "text",
      hovertext = paste("<b>", df[, Gg], "</b>")
    ) %>%
    add_lines(
      y = df[, Cc],
      name = "C",
      line = list(color = cols[4], width = 2),
      hoverinfo = "text",
      hovertext = paste("<b>", df[, Cc], "</b>")
    ) %>%
    add_lines(
      y = df[, meanQual],
      name = "Mean",
      line = list(color = "black", width = 3, dash = "dot"),
      hoverinfo = "text",
      hovertext = paste("<b>", df[, meanQual], "</b>")
    ) %>%
    layout(
      title = "Per base sequence quality",
      xaxis = list(title = "Position in read (bp)"),
      yaxis = list(title = "Quality score"),
      hovermode = "x closest"
    )
}



#' Reads qualoty data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @importFrom rjson fromJSON
#'
#' @param json A JSON file with QC information from \code{\link{qcFastq}}
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qc_fastq(file = fq1, outDir = outDir)
#' df <- readsQualDF(json = qc_res$json)
#'
#' # Output
#' df
readQualDF <- function(json, which = "before") {
  js <- fromJSON(file = json)

  qual <- NULL
  if (which == "before") {
    qual <- js$read1_before_filtering$quality_curves
  } else {
    qual <- js$read1_after_filtering$quality_curves
  }

  data <- data.frame(
    bpPos = 1:length(qual$mean),
    Aa = qual$A,
    Tt = qual$T,
    Gg = qual$G,
    Cc = qual$C,
    meanQual = qual$mean
  )

  return(data)
}



#' Base content line plot
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @importFrom dplyr `%>%`
#' @importFrom plotly plot_ly add_trace add_lines layout
#' @importFrom rcartocolor carto_pal
#'
#' @param json A JSON file with QC information from \code{\link{qcFastq}}
#' @param df A data frame
#' @param bpPos Column name with base-pair position
#' @param Aa Column with Adenine quality score
#' @param Tt Column with Thymine quality score
#' @param Gg Column with Guanine quality score
#' @param Cc Column with Cytosine quality score
#' @param GC Column with Mean Guanine & Cytosine quality score
#' @param Nn Column with N (could be any of A,T,G, or C) score
#' @param isRatio If the Quality scores are ratios or percentages
#'
#' @return An interactive plot
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qc_fastq(file = fq1, outDir = outDir)
#' p <- baseRatioLinePlot(json = qc_res$json)
#'
#' # Output
#' p
#' @export
baseRatioLinePlot <- function(json = NULL,
                              df = NULL,
                              bpPos = "bpPos",
                              Aa = "Aa",
                              Tt = "Tt",
                              Gg = "Gg",
                              Cc = "Cc",
                              GC = NULL,
                              Nn = "Nn",
                              isRatio = FALSE) {


  if (!is.null(json)) {
    df <- baseRatioDF(json)
    GC <- "GC"
  }

  if (isRatio) {
    df <- df[, c(bpPos, Aa, Tt, Gg, Cc, Nn)] * 100
  }

  if (is.null(GC)) {
    df$GC <- rowMeans(df[, c(Gg, Cc)])
    GC <- "GC"
  }

  cols <- carto_pal(n = 11, name = "Safe")[c(2, 1, 3, 4, 10)]

  plot_ly(df, x = df[, bpPos]) %>%
    add_lines(
      y = df[, Aa],
      name = "A",
      line = list(color = cols[1], width = 2),
      hoverinfo = "text+x",
      hovertext = paste("<b>", df[, Aa], "</b>")
    ) %>%
    add_lines(
      y = df[, Tt],
      name = "T",
      line = list(color = cols[2], width = 2),
      hoverinfo = "text+x",
      hovertext = paste("<b>", df[, Tt], "</b>")
    ) %>%
    add_lines(
      y = df[, Gg],
      name = "G",
      line = list(color = cols[3], width = 2),
      hoverinfo = "text+x",
      hovertext = paste("<b>", df[, Gg], "</b>")
    ) %>%
    add_lines(
      y = df[, Cc],
      name = "C",
      line = list(color = cols[4], width = 2),
      hoverinfo = "text+x",
      hovertext = paste("<b>", df[, Cc], "</b>")
    ) %>%
    add_lines(
      y = df[, GC],
      name = "GC",
      line = list(color = cols[5], width = 3, dash = "dot"),
      hoverinfo = "text+x",
      hovertext = paste("<b>", df[, GC], "</b>")
    ) %>%
    add_lines(
      y = df[, Nn],
      name = "N",
      line = list(color = "black", width = 3),
      hoverinfo = "text+x",
      hovertext = paste("<b>", df[, Nn], "</b>")
    ) %>%
    layout(
      title = "Per base sequence content",
      xaxis = list(title = "Position in read (bp)/ Cycle"),
      yaxis = list(title = "Percentage (%)"),
      hovermode = "x closest"
    )
}



#' Base content data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @importFrom rjson fromJSON
#'
#' @param json A JSON file with QC information from \code{\link{qcFastq}}
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qc_fastq(file = fq1, outDir = outDir)
#' df <- baseRatioDF(json = qc_res$json)
#'
#' # Output
#' df
baseRatioDF <- function(json, which = "before") {
  js <- fromJSON(file = json)

  qual <- NULL
  if (which == "before") {
    qual <- js$read1_before_filtering$content_curves
  } else {
    qual <- js$read1_after_filtering$content_curves
  }

  data <- data.frame(
    bpPos = 1:length(qual$A),
    Aa = qual$A * 100,
    Tt = qual$T * 100,
    Gg = qual$G * 100,
    Cc = qual$C * 100,
    GC = qual$GC * 100,
    Nn = qual$N * 100
  )

  return(data)
}



#' Base content ratio proportion plot
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @importFrom dplyr `%>%`
#' @importFrom plotly plot_ly add_trace add_lines layout
#' @importFrom rcartocolor carto_pal
#'
#' @param json A JSON file with QC information from \code{\link{qcFastq}}
#' @param df A data frame
#' @param bpPos Column name with base-pair position
#' @param Aa Column with Adenine quality score
#' @param Tt Column with Thymine quality score
#' @param Gg Column with Guanine quality score
#' @param Cc Column with Cytosine quality score
#' @param GC Column with Mean Guanine & Cytosine quality score
#' @param Nn Column with N (could be any of A,T,G, or C) score
#' @param isRatio If the Quality scores are ratios or percentages
#'
#' @return An interactive plot
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qc_fastq(file = fq1, outDir = outDir)
#' p <- baseRatioProportionPlot(json = qc_res$json)
#'
#' # Output
#' p
#' @export
baseRatioProportionPlot <- function(json = NULL,
                                    df = NULL,
                                    bpPos = "bpPos",
                                    Aa = "Aa",
                                    Tt = "Tt",
                                    Gg = "Gg",
                                    Cc = "Cc",
                                    GC = NULL,
                                    Nn = "Nn",
                                    isRatio = FALSE) {


  if (!is.null(json)) {
    df <- baseRatioDF(json)
    GC <- "GC"
  }

  if (isRatio) {
    df <- df[, c(bpPos, Aa, Tt, Gg, Cc, Nn)] * 100
  }

  if (is.null(GC)) {
    df$GC <- rowMeans(df[, c(Gg, Cc)])
    GC <- "GC"
  }

  cols <- carto_pal(n = 11, name = "Safe")[c(2, 1, 3, 4, 10)]

  plot_ly(df, x = df[, bpPos]) %>%
    add_trace(
      y = df[, Nn], name = "N", fillcolor = "black", type = "scatter",
      mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_trace(
      y = df[, Cc], name = "C", fillcolor = cols[4], type = "scatter",
      mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_trace(
      y = df[, Gg], name = "G", fillcolor = cols[3], type = "scatter",
      mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_trace(
      y = df[, Tt], name = "T", fillcolor = cols[2], type = "scatter",
      mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_trace(
      y = df[, Aa], name = "A", fillcolor = cols[1], type = "scatter",
      mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_lines(
      y = df[, GC],
      name = "GC",
      line = list(color = cols[5], width = 3, dash = "dot"),
      hoverinfo = "text",
      hovertext = paste("<b>", df[, GC], "</b>")
    ) %>%
    layout(
      title = "Per base sequence content",
      xaxis = list(title = "Position in read (bp)/ Cycle"),
      yaxis = list(title = "Percentage (%)"),
      hovermode = "x closest"
    )
}



#' Extract reads length information from a fastq file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @importFrom seqTools fastqq
#' @importFrom qckitfastq read_length
#'
#' @param fq A fastq file
#'
#' @return A data frame with reads length (read_length)
#' and number of reads (num_reads)
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#'
#' # Analysis
#' df <- readLengthFastq(fq = fq1)
#'
#' # Output
#' df
readLengthFastq <- function(fq) {

  out <- list()

  infile <- fq
  fseq <- fastqq(infile)
  read_len <- read_length(fseq)

  return(read_len)
}



#' Read length distribution (barplot)
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @importFrom dplyr `%>%`
#' @importFrom plotly plot_ly add_trace add_lines layout
#'
#' @param fq A fastq file
#' @param df A data frame with reads length and number of reads
#' @param read_length Column with reads length
#' @param num_reads Column with number of reads of a particulat length
#' @param col Color of the bars
#' @param name Either the fastq file is trimmed or not-trimmed
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#'
#' # Analysis
#' p <- readLengthHist(fq = fq1)
#'
#' # Output
#' p
#' @export
readLengthHist <- function(fq = NULL,
                           df = NULL,
                           read_length = "read_length",
                           num_reads = "num_reads",
                           col = "red",
                           name = "Not trimmed") {


  if (!is.null(fq)) {
    df <- readLengthFastq(fq)
    df[, name] <- name
  }

  plot_ly(hoverinfo = "text", type = "bar", textposition = "auto") %>%
    add_trace(
      data = df,
      legendgroup = df[, name],
      x = df[, read_length],
      y = df[, num_reads],
      name = df[, name],
      marker = list(color = col),
      hovertext = paste(
        "<b>Length:</b> ", df[, read_length],
        "<br><b>nReads:</b> ", df[, num_reads]
      )
    ) %>%
    layout(
      title = "Histogram of read length distribution",
      xaxis = list(title = "Reads length"),
      yaxis = list(title = "No. of reads"),
      hovermode = "x closest"
    )
}
