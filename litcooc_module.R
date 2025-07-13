# modules/litcooc_module.R

# Load required libraries
library(shiny)
library(DT)
library(ggplot2)
library(httr)
library(jsonlite)
library(visNetwork)
library(igraph)
library(webshot2)
library(dplyr)
library(purrr)

litcoocUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "litcooc",
          fluidPage(
            h3("📚 Literature Co-Occurrence Semantic Network"),
            p("Analyze the co-occurrence of your PPI hub genes with a specified term in Europe PMC, visualizing the relationships as a network."),
            
            # Input for search term
            textInput(ns("term"), "Co-Occurrence Term (e.g., 'cancer', 'apoptosis'):", "cancer"),
            
            # Action button to trigger the search and network generation
            actionButton(ns("search"), "Generate Co-Occurrence Network"),
            br(), br(),
            
            # Output for the interactive network visualization
            h4("Co-Occurrence Network (Genes and Term)"),
            visNetwork::visNetworkOutput(ns("lit_net"), height = "700px"),
            downloadButton(ns("dl_net_png"), "Download Network (.png)"),
            
            br(), br(),
            
            # Output for the co-mention table
            h4("Co-Mention Count Table"),
            DT::DTOutput(ns("lit_table")),
            
            # Output for the top co-mentions bar plot
            br(), br(),
            h4("Top Co-Mentions Barplot"),
            plotOutput(ns("lit_plot"), height = "500px"),
            downloadButton(ns("dl_barplot"), "Download Barplot (.png)")
          )
  )
}

litcoocServer <- function(id, ppi_clusters) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive value to store the co-mention data frame
    cooc_data <- reactiveVal(NULL)
    # Reactive value to store the visNetwork data (nodes/edges)
    vis_data <- reactiveVal(NULL)
    
    observeEvent(input$search, {
      req(ppi_clusters())
      genes <- unique(ppi_clusters()$Gene)
      term <- input$term
      
      if (length(genes) == 0) {
        showNotification("No genes found in PPI clusters.", type = "error")
        return()
      }
      
      showNotification("Searching Europe PMC for co-mentions...", type = "message")
      
      # 1. Fetch co-mention counts
      cooc_counts <- withProgress(message = "Querying Europe PMC...", value = 0, {
        
        counts <- sapply(genes, function(gene) {
          incProgress(1/length(genes), detail = paste("Processing:", gene))
          query <- paste0(gene, " AND ", term)
          # Note: Europe PMC REST API is generally preferred for systematic searches
          url <- paste0("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=", 
                        URLencode(query), "&format=json&resulttype=lite")
          
          res <- tryCatch(GET(url), error = function(e) NULL)
          
          if (!is.null(res) && res$status_code == 200) {
            json <- fromJSON(content(res, "text"))
            if (!is.null(json$hitCount)) return(as.numeric(json$hitCount))
          }
          return(0)
        })
        names(counts) <- genes
        return(counts)
      })
      
      # 2. Create the data frame for the table and barplot
      df <- data.frame(
        Gene = names(cooc_counts),
        CoMentionCount = as.numeric(cooc_counts)
      )
      df <- df[order(-df$CoMentionCount), ]
      df <- df[df$CoMentionCount > 0, ] # Only keep genes with at least one co-mention
      
      if (nrow(df) == 0) {
        showNotification(paste("No co-mentions found for term '", term, "' and the selected genes."), type = "warning")
        cooc_data(NULL)
        vis_data(NULL)
        return()
      }
      
      cooc_data(df)
      
      # 3. Build the Network data (for visNetwork)
      
      # Nodes: Genes and the search term
      # Add the search term as a central node
      nodes_df <- data.frame(
        id = c(df$Gene, input$term),
        label = c(df$Gene, input$term),
        group = c(rep("Gene", nrow(df)), "Term"),
        value = c(df$CoMentionCount, max(df$CoMentionCount) * 1.5), # Scale term node for emphasis
        shape = c(rep("ellipse", nrow(df)), "square")
      )
      
      # Add links for UniProt and gene description if possible
      nodes_df$title <- paste0("<b>", nodes_df$label, "</b><br>Group: ", nodes_df$group)
      
      # Edges: Co-mention relationships between genes and the term
      edges_df <- data.frame(
        from = df$Gene,
        to = input$term,
        value = df$CoMentionCount, # Use co-mention count as edge thickness
        title = paste("Co-mentions:", df$CoMentionCount)
      )
      
      # Store network data
      vis_data(list(nodes = nodes_df, edges = edges_df))
      
      # 4. Render outputs
      
      # Render the network plot
      output$lit_net <- visNetwork::renderVisNetwork({
        visNetwork(nodes = vis_data()$nodes, edges = vis_data()$edges) %>%
          visIgraphLayout(layout = "layout_with_fr") %>%
          visNodes(
            color = list(background = "skyblue", border = "darkblue", highlight = "red"),
            font = list(size = 20)
          ) %>%
          visEdges(smooth = TRUE, color = list(color = "gray", highlight = "red"), scaling = list(min = 1, max = 15)) %>%
          visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
          visLegend()
      })
      
      # Render the table
      output$lit_table <- DT::renderDT(df)
      
      # Render the bar plot
      output$lit_plot <- renderPlot({
        top_20 <- head(df, 20)
        ggplot(top_20, aes(x = reorder(Gene, CoMentionCount), y = CoMentionCount)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          coord_flip() +
          labs(x = "Gene", y = paste("Co-Mentions with:", term)) +
          theme_minimal(base_size = 14)
      })
    })
    
    # --- Download Handlers ---
    
    # Download Barplot
    output$dl_barplot <- downloadHandler(
      filename = function() "cooc_barplot.png",
      content = function(file) {
        ggsave(file, plot = last_plot(), width = 10, height = 8, dpi = 300)
      }
    )
    
    # Download Network Snapshot
    output$dl_net_png <- downloadHandler(
      filename = function() { "cooc_network.png" },
      content = function(file) {
        if (is.null(vis_data())) return(NULL)
        tmp_html <- tempfile(fileext = ".html")
        
        # Save the network visualization to a temporary HTML file
        visNetwork::visSave(
          visNetwork::visNetwork(vis_data()$nodes, vis_data()$edges) %>%
            visIgraphLayout(layout = "layout_with_fr") %>%
            visNodes(
              color = list(background = "skyblue", border = "darkblue"),
              font = list(size = 20)
            ) %>%
            visEdges(smooth = TRUE, color = list(color = "gray"), scaling = list(min = 1, max = 15)) %>%
            visOptions(highlightNearest = TRUE),
          file = tmp_html
        )
        
        # Add a short delay to ensure the HTML is ready before taking the snapshot
        Sys.sleep(1)
        
        # Use webshot2 to capture the high-resolution image
        webshot2::webshot(url = tmp_html, file = file, vwidth = 1800, vheight = 1400)
      }
    )
  })
}