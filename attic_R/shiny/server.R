# library(shiny)
# library(edgeR)

server <- shinyServer(function(input, output) {
  getCols <- function(x) {
    ff <- tolower(c("name", "logCPM", "AvgExp", "logFC", "PValue", "P.Value", "adj.p.val", "fdr", "padj"))
    x$feature <- row.names(x)
    row.names(x) <- NULL
    x[, c("feature", colnames(x)[which(tolower(colnames(x)) %in% ff)])]
  }

  output$title2 <- renderText({
    o@description
  })

  output$composition <- renderPlot({
    plotComposition(o)
  })

  output$pca <- renderPlot({
    ft <- strsplit(input$featureType, ".", fixed = T)[[1]]
    if (ft[[1]] == "seq") {
      if (ft[[2]] == "all") {
        type <- NULL
      } else {
        type <- ft[[2]]
      }
      e <- getSeqCounts(o, type = type, status = c("unknown", "ambiguous", "unique"), normalized = TRUE)
    } else {
      if (ft[[2]] == "all") {
        type <- NULL
      } else {
        type <- ft[[2]]
      }
      e <- getAggCounts(o, type = type, normalized = TRUE)
    }
    if (input$pcaType == "pca") {
      plPCA(log(e + 1))
    } else {
      e <- cor(log(e + 1))
      for (i in 1:ncol(e)) e[i, i] <- NA
      byheatmap(e, annotation_col = o@phenoData)
    }
  })

  s <- reactive({
    if (nchar(input$findFeature) > 3) {
      return(getFeatureSeqs(o, input$findFeature, exact = input$exactFeature))
    } else {
      return(c())
    }
  })
  output$seqc <- renderPlot({
    if (length(s()) == 0) {
      return("Nothing found!")
    } else {
      if (length(s()) > 100) {
        plot.new()
        legend("topleft", legend = "Too many matching sequences to plot...")
      } else {
        plotFeature(o, input$findFeature, exact = input$exactFeature, normalized = input$normFeature, main = input$findFeature, plotType = input$fplotType, aggregated = input$fagg, las = 3)
      }
    }
  })
  output$seqtable <- renderDataTable({
    if (length(s()) == 0) {
      return("No data.")
    } else {
      getCols(res$seqLevel$allKnownUnique[which(row.names(res$seqLevel$allKnownUnique) %in% s()), ])
    }
  })

  output$volcano <- renderPlot({
    tt <- strsplit(input$deaTypes, ".", fixed = T)[[1]]
    if (tt[[1]] == "agg") {
      rr <- res$aggregated[[tt[[2]]]]
    } else {
      rr <- res$seqLevel[[tt[[2]]]]
    }
    if (nrow(rr) > 0) {
      volcano(rr)
    } else {
      return("No data.")
    }
  })

  output$deatable <- renderDataTable(
    {
      tt <- strsplit(input$deaTypes, ".", fixed = T)[[1]]
      if (tt[[1]] == "agg") {
        rr <- res$aggregated[[tt[[2]]]]
      } else {
        rr <- res$seqLevel[[tt[[2]]]]
      }
      if (is.null(rr)) {
        return("No data.")
      } else {
        getCols(rr)
      }
    },
    options = list(lengthMenu = c(20, 50, 100), pageLength = 20)
  )

  output$sizeplot <- renderPlot({
    checkSizeSelection(o)
  })

  output$info <- renderPrint({
    res$info
  })
})
