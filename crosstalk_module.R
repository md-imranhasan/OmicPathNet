# modules/crosstalk_module.R

crosstalkUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "crosstalk",
          fluidPage(
            numericInput(ns("topn"), "Top N shared pathways:", value = 10, min = 5, max = 50, step = 5),
            selectInput(ns("list1"), "Select List 1:", choices = NULL),
            selectInput(ns("list2"), "Select List 2:", choices = NULL),
            actionButton(ns("analyze"), "Analyze Crosstalk"),
            br(), br(),
            
            h4("Shared Pathways Table"),
            DT::DTOutput(ns("overlap_table")),
            downloadButton(ns("dl_table_csv"), "Download (.csv)"),
            downloadButton(ns("dl_table_xlsx"), "Download (.xlsx)"),
            
            br(), br(),
            h4("Overlap Significance Barplot"),
            plotOutput(ns("barplot"), height = "600px"),
            downloadButton(ns("dl_barplot"), "Download Barplot (.png)"),
            
            br(), br(),
            h4("UpSet Plot"),
            plotOutput(ns("upset"), height = "600px"),
            downloadButton(ns("dl_upset"), "Download UpSet (.png)"),
            
            br(), br(),
            h4("Pathway–Gene Chord Diagram"),
            plotOutput(ns("chord"), height = "800px"),
            downloadButton(ns("dl_chord"), "Download Chord (.png)")
          )
  )
}

crosstalkServer <- function(id, user_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observe({
      req(user_data$seq_lists())
      n_lists <- length(user_data$seq_lists())
      choices <- as.character(seq_len(n_lists))
      updateSelectInput(session, "list1", choices = choices, selected = "1")
      updateSelectInput(session, "list2", choices = choices, selected = ifelse(n_lists >= 2, "2", "1"))
    })
    
    observeEvent(input$analyze, {
      req(user_data$seq_lists())
      
      if (length(user_data$seq_lists()) < 2) {
        showNotification("Upload at least two lists for crosstalk!", type = "error")
        return(NULL)
      }
      if (input$list1 == input$list2) {
        showNotification("Please select two different lists.", type = "error")
        return(NULL)
      }
      
      genes1 <- user_data$seq_lists()[[as.integer(input$list1)]]
      genes2 <- user_data$seq_lists()[[as.integer(input$list2)]]
      
      OrgDb <- org.Hs.eg.db
      g1 <- bitr(genes1, fromType = user_data$id_type(), toType = "ENTREZID", OrgDb = OrgDb)
      g2 <- bitr(genes2, fromType = user_data$id_type(), toType = "ENTREZID", OrgDb = OrgDb)
      
      if (nrow(g1) == 0 || nrow(g2) == 0) {
        showNotification("One or both lists failed to map.", type = "error")
        return(NULL)
      }
      
      ekegg1 <- enrichKEGG(gene = unique(g1$ENTREZID), organism = "hsa")
      ekegg2 <- enrichKEGG(gene = unique(g2$ENTREZID), organism = "hsa")
      
      if (nrow(as.data.frame(ekegg1)) == 0 || nrow(as.data.frame(ekegg2)) == 0) {
        showNotification("One or both lists have no enriched pathways.", type = "warning")
        return(NULL)
      }
      
      kegg1 <- ekegg1@result
      kegg2 <- ekegg2@result
      shared <- intersect(kegg1$Description, kegg2$Description)
      
      if (length(shared) == 0) {
        showNotification("No shared pathways found.", type = "warning")
        return(NULL)
      }
      
      df <- data.frame(
        Pathway = shared,
        Pval1 = kegg1$pvalue[match(shared, kegg1$Description)],
        Pval2 = kegg2$pvalue[match(shared, kegg2$Description)],
        Genes1 = kegg1$geneID[match(shared, kegg1$Description)],
        Genes2 = kegg2$geneID[match(shared, kegg2$Description)]
      )
      
      ## Add gene symbols
      all_ids1 <- unique(unlist(strsplit(df$Genes1, "/")))
      map1 <- bitr(all_ids1, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = OrgDb)
      df$Genes1_Symbol <- sapply(df$Genes1, function(x) {
        ids <- unlist(strsplit(x, "/"))
        syms <- map1$SYMBOL[match(ids, map1$ENTREZID)]
        paste(na.omit(syms), collapse = "/")
      })
      
      all_ids2 <- unique(unlist(strsplit(df$Genes2, "/")))
      map2 <- bitr(all_ids2, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = OrgDb)
      df$Genes2_Symbol <- sapply(df$Genes2, function(x) {
        ids <- unlist(strsplit(x, "/"))
        syms <- map2$SYMBOL[match(ids, map2$ENTREZID)]
        paste(na.omit(syms), collapse = "/")
      })
      
      df <- df[order(df$Pval1), ][1:min(input$topn, nrow(df)), ]
      
      ## === Table ===
      output$overlap_table <- DT::renderDT(df, options = list(scrollX = TRUE))
      output$dl_table_csv <- downloadHandler("crosstalk.csv", function(f) write.csv(df, f, row.names = FALSE))
      output$dl_table_xlsx <- downloadHandler("crosstalk.xlsx", function(f) openxlsx::write.xlsx(df, f))
      
      ## === ggplot2 Barplot ===
      output$barplot <- renderPlot({
        library(ggplot2)
        df$logP <- -log10(df$Pval1)
        ggplot(df, aes(x = reorder(Pathway, logP), y = logP)) +
          geom_bar(stat = "identity", fill = "skyblue") +
          coord_flip() +
          labs(x = "Pathway", y = "-log10(Pval1)") +
          theme_minimal(base_size = 16) +
          theme(axis.text.y = element_text(size = 12, hjust = 1))
      })
      output$dl_barplot <- downloadHandler("crosstalk_barplot.png", function(f) {
        ggsave(f, plot = last_plot(), width = 12, height = 7, dpi = 300)
      })
      
      ## === UpSetR Plot ===
      output$upset <- renderPlot({
        pathlist <- list(List1 = kegg1$Description, List2 = kegg2$Description)
        UpSetR::upset(
          UpSetR::fromList(pathlist),
          order.by = "freq",
          main.bar.color = "steelblue",
          sets.bar.color = "skyblue",
          text.scale = 1.4
        )
      })
      output$dl_upset <- downloadHandler("crosstalk_upset.png", function(f) {
        png(f, width = 1000, height = 700, res = 150)
        pathlist <- list(List1 = kegg1$Description, List2 = kegg2$Description)
        UpSetR::upset(
          UpSetR::fromList(pathlist),
          order.by = "freq",
          main.bar.color = "steelblue",
          sets.bar.color = "skyblue",
          text.scale = 1.4
        )
        dev.off()
      })
      
      ## === Chord Diagram ===
      output$chord <- renderPlot({
        library(circlize)
        circos.clear()
        links <- data.frame(
          Pathway = rep(df$Pathway, sapply(strsplit(df$Genes1_Symbol, "/"), length)),
          Gene = unlist(strsplit(df$Genes1_Symbol, "/"))
        )
        n_pathways <- length(unique(df$Pathway))
        gap <- ifelse(n_pathways > 20, 1, 4)
        circos.par(gap.degree = gap, track.height = 0.08)
        
        chordDiagram(
          table(links),
          annotationTrack = "grid",
          preAllocateTracks = list(track.height = 0.1)
        )
        
        ## Add readable sector names
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
            cex = 0.8
          )
        }, bg.border = NA)
        
        circos.clear()
      })
      output$dl_chord <- downloadHandler("crosstalk_chord.png", function(f) {
        png(f, width = 1600, height = 1600, res = 300)
        library(circlize)
        circos.clear()
        links <- data.frame(
          Pathway = rep(df$Pathway, sapply(strsplit(df$Genes1_Symbol, "/"), length)),
          Gene = unlist(strsplit(df$Genes1_Symbol, "/"))
        )
        n_pathways <- length(unique(df$Pathway))
        gap <- ifelse(n_pathways > 20, 1, 4)
        circos.par(gap.degree = gap, track.height = 0.08)
        
        chordDiagram(
          table(links),
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
            cex = 0.8
          )
        }, bg.border = NA)
        
        circos.clear()
        dev.off()
      })
    })
  })
}
