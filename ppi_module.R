# modules/ppi_module.R

# Load required libraries at the start if this file is sourced directly, 
# although in a typical Shiny app structure, they might be loaded globally or within the server function.
# We keep the UI definition as provided.

ppiUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "ppi",
          fluidPage(
            selectInput(ns("which_list"), "Which List(s):",
                        choices = NULL, selected = "Combine All"),
            actionButton(ns("run"), "Build PPI Network (STRINGdb)"),
            br(), br(),
            h4("Interactive PPI Network (High-Res + Gene Labels)"),
            visNetwork::visNetworkOutput(ns("ppi_net"), height = "800px"),
            br(),
            downloadButton(ns("dl_ppi_png"), "Download PPI Network (.png)"),
            br(), br(),
            h4("Cluster Assignments"),
            DT::DTOutput(ns("cluster_table")),
            downloadButton(ns("dl_cluster_csv"), "Download Cluster Table (.csv)"),
            downloadButton(ns("dl_cluster_xlsx"), "Download Cluster Table (.xlsx)"),
            br(), br(),
            h4("PPI Interactions (with Score)"),
            DT::DTOutput(ns("ppi_table")),
            downloadButton(ns("dl_ppi_table_csv"), "Download PPI Table (.csv)"),
            downloadButton(ns("dl_ppi_table_xlsx"), "Download PPI Table (.xlsx)")
          )
  )
}

ppiServer <- function(id, user_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values to store results
    ppi_vis <- reactiveVal()
    cluster_data <- reactiveVal()
    ppi_table_data <- reactiveVal()
    # New reactive value to store hub data (including centrality measures)
    cluster_hubs <- reactiveVal(NULL)
    
    # Update the 'Which List(s)' dropdown based on user_data
    observe({
      req(user_data$seq_lists())
      n_lists <- length(user_data$seq_lists())
      choices <- c("Combine All", paste0("List ", seq_len(n_lists)))
      updateSelectInput(session, "which_list", choices = choices, selected = "Combine All")
    })
    
    # Main logic triggered when the 'Run' button is clicked
    observeEvent(input$run, {
      # Ensure necessary libraries are loaded (good practice within observeEvent for localized dependencies)
      library(STRINGdb)
      library(visNetwork)
      library(webshot2)
      library(igraph)
      library(openxlsx)
      
      # 1. Determine which gene list to use
      if (input$which_list == "Combine All") {
        genes <- unique(unlist(user_data$seq_lists()))
      } else {
        list_num <- as.integer(gsub("List ", "", input$which_list))
        genes <- unique(user_data$seq_lists()[[list_num]])
      }
      
      # Validate gene list size
      if (length(genes) < 2) {
        showNotification("Need at least 2 valid genes!", type = "error")
        return(NULL)
      }
      
      # 2. Map genes to STRING IDs
      string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)
      mapped <- string_db$map(data.frame(Gene = genes), "Gene", removeUnmappedRows = TRUE)
      
      if (nrow(mapped) < 2) {
        showNotification("Not enough mapped genes.", type = "error")
        return(NULL)
      }
      
      # 3. Get PPI interactions from STRINGdb
      ppi_raw <- string_db$get_interactions(mapped$STRING_id)
      if (!is.data.frame(ppi_raw) || nrow(ppi_raw) == 0) {
        showNotification("No interactions found.", type = "error")
        return(NULL)
      }
      
      # 4. Map STRING IDs back to gene symbols and clean data
      node_map <- mapped[, c("Gene", "STRING_id")]
      ppi_raw$from <- node_map$Gene[match(ppi_raw$from, node_map$STRING_id)]
      ppi_raw$to <- node_map$Gene[match(ppi_raw$to, node_map$STRING_id)]
      
      ppi_clean <- na.omit(ppi_raw[, c("from", "to", "combined_score")])
      ppi_edges <- unique(ppi_clean)
      
      if (nrow(ppi_edges) == 0) {
        showNotification("No valid edges after filtering.", type = "error")
        return(NULL)
      }
      
      # 5. Build igraph object
      g <- igraph::graph_from_data_frame(ppi_edges, directed = FALSE)
      
      # 6. Perform Community detection (Louvain)
      comm <- cluster_louvain(g)
      membership <- membership(comm)
      
      # 7. Calculate centrality metrics (Degree and Betweenness)
      V(g)$degree <- igraph::degree(g)
      V(g)$betweenness <- igraph::betweenness(g)
      
      # 8. Create and store cluster/hub data
      cluster_df <- data.frame(
        Gene = V(g)$name,
        Cluster = membership,
        Degree = V(g)$degree,
        Betweenness = V(g)$betweenness
      )
      
      # Sort by cluster and then by degree (to identify hubs)
      cluster_df <- cluster_df[order(cluster_df$Cluster, -cluster_df$Degree), ]
      
      # Store data in reactive values
      cluster_data(cluster_df)
      ppi_table_data(ppi_edges)
      cluster_hubs(cluster_df) # Store the hub data specifically for export/downstream use
      
      # 9. Prepare data for visNetwork
      nodes <- data.frame(
        id = V(g)$name,
        label = V(g)$name,
        group = as.factor(membership),
        # Enhanced title for tooltips, including cluster and degree
        title = paste0(
          "<b>", V(g)$name, "</b><br>",
          "Cluster: ", membership, "<br>",
          "Degree: ", V(g)$degree, "<br>",
          "<a href='https://www.uniprot.org/uniprot/?query=", V(g)$name, "' target='_blank'>UniProt</a>"
        ),
        shape = "ellipse"
      )
      
      edges <- data.frame(
        from = ppi_edges$from,
        to = ppi_edges$to,
        value = ppi_edges$combined_score / 1000
      )
      
      ppi_vis(list(nodes = nodes, edges = edges))
      
      # 10. Render the PPI network
      output$ppi_net <- visNetwork::renderVisNetwork({
        visNetwork(nodes, edges) %>%
          visIgraphLayout(layout = "layout_with_fr") %>%
          visNodes(font = list(size = 30, face = "bold")) %>%
          visEdges(smooth = TRUE) %>%
          visOptions(highlightNearest = TRUE) %>%
          visLegend()
      })
      
      # 11. Render tables
      output$cluster_table <- DT::renderDT({ cluster_df }, options = list(scrollX = TRUE))
      output$ppi_table <- DT::renderDT({ ppi_edges }, options = list(scrollX = TRUE))
    })
    
    # --- Download Handlers ---
    
    output$dl_cluster_csv <- downloadHandler(
      filename = function() { "ppi_clusters.csv" },
      content = function(file) {
        write.csv(cluster_data(), file, row.names = FALSE)
      }
    )
    
    output$dl_cluster_xlsx <- downloadHandler(
      filename = function() { "ppi_clusters.xlsx" },
      content = function(file) {
        openxlsx::write.xlsx(cluster_data(), file)
      }
    )
    
    output$dl_ppi_table_csv <- downloadHandler(
      filename = function() { "ppi_interactions.csv" },
      content = function(file) {
        write.csv(ppi_table_data(), file, row.names = FALSE)
      }
    )
    
    output$dl_ppi_table_xlsx <- downloadHandler(
      filename = function() { "ppi_interactions.xlsx" },
      content = function(file) {
        openxlsx::write.xlsx(ppi_table_data(), file)
      }
    )
    
    # Download handler for high-resolution PNG snapshot
    output$dl_ppi_png <- downloadHandler(
      filename = function() { "ppi_network.png" },
      content = function(file) {
        if (is.null(ppi_vis())) return(NULL)
        tmp_html <- tempfile(fileext = ".html")
        
        # Save the network visualization to a temporary HTML file
        visNetwork::visSave(
          visNetwork::visNetwork(ppi_vis()$nodes, ppi_vis()$edges) %>%
            visIgraphLayout(layout = "layout_with_fr") %>%
            visNodes(font = list(size = 30, face = "bold")) %>%
            visEdges(smooth = TRUE) %>%
            visOptions(highlightNearest = TRUE) %>%
            visLegend(),
          file = tmp_html
        )
        
        # Add a short delay to ensure the HTML is ready before taking the snapshot
        Sys.sleep(1)
        
        # Use webshot2 to capture the high-resolution image
        webshot2::webshot(url = tmp_html, file = file, vwidth = 2200, vheight = 1800)
      }
    )
    
    # Return reactive values for use in other modules
    return(list(
      ppi_vis = ppi_vis,
      cluster_data = cluster_data,
      ppi_table_data = ppi_table_data,
      cluster_hubs = cluster_hubs # Added 'cluster_hubs' to the return list
    ))
  })
}