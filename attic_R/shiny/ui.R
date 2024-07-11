# library(shiny)


shinyUI(fluidPage(
  titlePanel(textOutput("title2"), windowTitle = "shiny smallRNAseq explorer"),
  mainPanel(
    tabsetPanel(
      tabPanel("Overview", fluidRow(
        column(
          3, selectInput("featureType", label = "Feature type", choices = c(
            "seq-level: all known" = "seq.all", "seq-level: miRNA" = "seq.miRNA", "seq-level: tRNA" = "seq.tRNA",
            "aggregated: all known" = "agg.all", "aggregated: miRNA" = "agg.miRNA", "aggregated: tRNA" = "agg.tRNA"
          )),
          selectInput("pcaType", label = h3("Plot type"), choices = c("PCA" = "pca", "Correlation heatmap" = "corheatmap"))
        ),
        column(7, plotOutput("pca", height = 500, width = 600))
      )),
      tabPanel("Composition", plotOutput("composition", height = 800)),
      tabPanel(
        "Sequences", fluidRow(
          column(
            2, textInput("findFeature", label = h3("Find feature(s)"), value = ""),
            checkboxInput("exactFeature", label = "Exact", value = FALSE),
            checkboxInput("normFeature", label = "Normalized", value = TRUE),
            checkboxInput("fagg", label = "Aggregated", value = FALSE),
            selectInput("fplotType", label = h3("Plot type"), choices = c("barplot" = "barplot", "heatmap" = "heatmap", "heatmap of z-scores" = "zheatmap"))
          ),
          column(10, plotOutput("seqc", height = 700))
        ),
        dataTableOutput("seqtable")
      ),
      tabPanel(
        "Differential expression", fluidRow(
          column(3, selectInput("deaTypes", label = "Feature type", choices = c(
            "seq-level: all known" = "seq.allKnownUnique", "seq-level: miRNA" = "seq.miRNA", "seq-level: tRNA" = "seq.tRNA",
            "seq-level: primary piRNA" = "seq.primary_piRNA", "seq-level: ambiguous" = "seq.ambiguous", "seq-level: unknown" = "seq.unknown",
            "aggregated: all non-ambiguous" = "agg.allUnique", "aggregated: miRNA" = "agg.miRNA", "aggregated: ambiguous" = "agg.ambiguous"
          ))),
          column(7, plotOutput("volcano", height = 300, width = 400))
        ),
        dataTableOutput("deatable")
      ),
      tabPanel("Abundances by size", plotOutput("sizeplot", height = 800, width = 1200)),
      tabPanel("info", verbatimTextOutput("info"))
    )
  )
))
