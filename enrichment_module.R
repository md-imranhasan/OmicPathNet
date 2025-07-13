# modules/enrichment_module.R

enrichmentUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "enrichment",
          fluidPage(
            selectInput(ns("organism"), "Organism:",
                        choices = list("Human" = "hsa", "Mouse" = "mmu", "Rat" = "rno"),
                        selected = "hsa"),
            numericInput(ns("topn"), "Top N pathways:", value = 10, min = 5, max = 50, step = 5),
            actionButton(ns("run"), "Run Dynamic Enrichment"),
            br(), br(),
            uiOutput(ns("dynamic_ui"))
          )
  )
}

enrichmentServer <- function(id, user_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    library(clusterProfiler)
    library(ggplot2)
    library(enrichplot)
    library(circlize)
    library(openxlsx)
    library(DT)
    library(stringr)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    library(org.Rn.eg.db)
    
    observeEvent(input$run, {
      req(user_data$seq_lists())
      
      OrgDb <- switch(input$organism,
                      "hsa" = org.Hs.eg.db,
                      "mmu" = org.Mm.eg.db,
                      "rno" = org.Rn.eg.db)
      
      n_lists <- length(user_data$seq_lists())
      showNotification(paste("Running enrichment for", n_lists, "lists"), type = "message")
      
      results <- vector("list", n_lists)
      
      withProgress(message = "Running enrichment", value = 0, {
        for (i in seq_len(n_lists)) {
          incProgress(1 / n_lists, detail = paste("List", i))
          
          mapped <- tryCatch({
            bitr(user_data$seq_lists()[[i]],
                 fromType = user_data$id_type(),
                 toType = "ENTREZID",
                 OrgDb = OrgDb)
          }, error = function(e) NULL)
          
          if (is.null(mapped) || nrow(mapped) == 0) {
            showNotification(paste("List", i, "failed to map."), type = "error")
            next
          }
          
          ekegg <- enrichKEGG(gene = unique(mapped$ENTREZID), organism = input$organism)
          if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
            showNotification(paste("List", i, "has no significant pathways."), type = "warning")
            next
          }
          
          ekegg <- setReadable(ekegg, OrgDb, keyType = "ENTREZID")
          ekegg <- pairwise_termsim(ekegg)
          
          results[[i]] <- list(
            ekegg = ekegg,
            df = as.data.frame(ekegg)
          )
        }
      })
      
      ## Dynamic UI
      output$dynamic_ui <- renderUI({
        all_ui <- lapply(seq_len(n_lists), function(i) {
          if (is.null(results[[i]])) return(NULL)
          tagList(
            tags$hr(),
            h3(paste("Gene List", i)),
            h4("Enrichment Table:"),
            DTOutput(ns(paste0("table_", i))),
            downloadButton(ns(paste0("dl_table_csv_", i)), "Download CSV"),
            downloadButton(ns(paste0("dl_table_xlsx_", i)), "Download XLSX"),
            h4("Bar Plot:"),
            plotOutput(ns(paste0("barplot_", i))),
            downloadButton(ns(paste0("dl_barplot_", i)), "Download Barplot"),
            h4("Dot Plot:"),
            plotOutput(ns(paste0("dotplot_", i))),
            downloadButton(ns(paste0("dl_dotplot_", i)), "Download Dotplot"),
            h4("Enrichment Map:"),
            plotOutput(ns(paste0("emap_", i))),
            downloadButton(ns(paste0("dl_emap_", i)), "Download Enrichment Map"),
            h4("Chord Diagram:"),
            plotOutput(ns(paste0("chord_", i))),
            downloadButton(ns(paste0("dl_chord_", i)), "Download Chord"),
            h4("Gene–Pathway Network (Cnetplot):"),
            plotOutput(ns(paste0("cnet_", i))),
            downloadButton(ns(paste0("dl_cnet_", i)), "Download Cnetplot")
          )
        })
        do.call(tagList, all_ui)
      })
      
      ## Outputs for each list
      lapply(seq_len(n_lists), function(i) {
        if (is.null(results[[i]])) return(NULL)
        df <- results[[i]]$df
        ekegg <- results[[i]]$ekegg
        
        output[[paste0("table_", i)]] <- renderDT({
          datatable(df, options = list(scrollX = TRUE))
        })
        output[[paste0("dl_table_csv_", i)]] <- downloadHandler(
          filename = function() paste0("list_", i, ".csv"),
          content = function(file) write.csv(df, file, row.names = FALSE)
        )
        output[[paste0("dl_table_xlsx_", i)]] <- downloadHandler(
          filename = function() paste0("list_", i, ".xlsx"),
          content = function(file) openxlsx::write.xlsx(df, file)
        )
        
        ## Barplot
        bar <- reactive({
          df_bar <- ekegg %>% as.data.frame()
          df_bar$logP <- -log10(df_bar$p.adjust)
          df_bar <- df_bar[order(df_bar$logP, decreasing = TRUE), ]
          df_bar <- head(df_bar, input$topn)
          df_bar$Description <- str_wrap(df_bar$Description, width = 40)
          max_len <- max(nchar(df_bar$Description))
          font_size <- ifelse(max_len > 50, 8, ifelse(max_len > 30, 10, 12))
          ggplot(df_bar, aes(x = reorder(Description, logP), y = logP)) +
            geom_bar(stat = "identity", fill = "steelblue") +
            coord_flip() +
            labs(x = "Pathway", y = "-log10(p.adjust)") +
            theme_minimal(base_size = 16) +
            theme(axis.text.y = element_text(size = font_size, hjust = 1))
        })
        output[[paste0("barplot_", i)]] <- renderPlot({ bar() })
        output[[paste0("dl_barplot_", i)]] <- downloadHandler(
          filename = function() paste0("barplot_", i, ".png"),
          content = function(file) ggsave(file, plot = bar(), width = 12, height = 7, dpi = 300)
        )
        
        ## Dotplot
        dot <- reactive({ dotplot(ekegg, showCategory = input$topn) })
        output[[paste0("dotplot_", i)]] <- renderPlot({ dot() })
        output[[paste0("dl_dotplot_", i)]] <- downloadHandler(
          filename = function() paste0("dotplot_", i, ".png"),
          content = function(file) ggsave(file, plot = dot(), width = 10, height = 6, dpi = 300)
        )
        
        ## Enrichment map
        emap <- reactive({
          if (!is.null(ekegg@termsim)) {
            emapplot(ekegg, showCategory = input$topn)
          } else {
            plot.new(); text(0.5, 0.5, "No enrichment map available")
          }
        })
        output[[paste0("emap_", i)]] <- renderPlot({ emap() })
        output[[paste0("dl_emap_", i)]] <- downloadHandler(
          filename = function() paste0("emap_", i, ".png"),
          content = function(file) {
            png(file, width = 1000, height = 800, res = 300)
            print(emap())
            dev.off()
          }
        )
        
        ## Chord diagram
        chord <- reactive({
          circos.clear()
          topn <- min(input$topn, nrow(df))
          if (topn == 0) return(NULL)
          pathgenes <- strsplit(df$geneID[1:topn], "/")
          pathways <- df$Description[1:topn]
          links <- data.frame(
            Pathway = rep(pathways, sapply(pathgenes, length)),
            Gene = unlist(pathgenes)
          )
          mat <- table(links)
          gap <- ifelse(length(unique(pathways)) > 20, 1, 4)
          circos.par(gap.degree = gap, track.height = 0.08)
          tryCatch({
            chordDiagram(
              mat,
              annotationTrack = "grid",
              preAllocateTracks = list(track.height = 0.1)
            )
            circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
              sector.name <- get.cell.meta.data("sector.index")
              xlim <- get.cell.meta.data("xlim")
              ylim <- get.cell.meta.data("ylim")
              circos.text(
                mean(xlim), ylim[1] + .1,
                sector.name,
                facing = "clockwise",
                niceFacing = TRUE,
                adj = c(0, 0.5),
                cex = 0.7
              )
            }, bg.border = NA)
          }, error = function(e) {
            showNotification(paste("Chord diagram failed:", e$message), type = "error")
          })
          circos.clear()
        })
        output[[paste0("chord_", i)]] <- renderPlot({ chord() })
        output[[paste0("dl_chord_", i)]] <- downloadHandler(
          filename = function() paste0("chord_", i, ".png"),
          content = function(file) {
            png(file, width = 1600, height = 1600, res = 300)
            chord(); dev.off()
          }
        )
        
        ## Cnetplot
        cnet <- reactive({ cnetplot(ekegg, showCategory = input$topn) })
        output[[paste0("cnet_", i)]] <- renderPlot({ cnet() })
        output[[paste0("dl_cnet_", i)]] <- downloadHandler(
          filename = function() paste0("cnet_", i, ".png"),
          content = function(file) {
            png(file, width = 1200, height = 1000, res = 300)
            print(cnet())
            dev.off()
          }
        )
      })
    })
  })
}
