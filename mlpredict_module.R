# modules/mlpredict_module.R

# Define the UI for the ML Hub Ranking module
mlpredictUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "mlpredict",
          fluidPage(
            h3("🧬 ML Hub Ranking (Topological & FDA Link)"),
            # Change the button ID to 'rank' to match the provided integration code
            actionButton(ns("rank"), "Run Hub Ranking"), 
            br(), br(),
            
            # Add a plot output for visualization
            h4("Top 10 Hubs by ML Score"),
            plotOutput(ns("rank_plot"), height = "400px"),
            br(),
            
            # Change the DT output ID to 'rank_table'
            h4("Hub Ranking Table"),
            DT::DTOutput(ns("rank_table")), 
            downloadButton(ns("dl_rank_csv"), "Download (.csv)"),
            downloadButton(ns("dl_rank_xlsx"), "Download (.xlsx)")
          )
  )
}

# Define the server logic for the ML Hub Ranking module
# ppi_clusters: Reactive containing PPI cluster data (now includes Degree and Betweenness from ppi_module.R)
# openfda_genes: Reactive containing the list of genes identified in OpenFDA (if applicable)
mlpredictServer <- function(id, ppi_clusters, openfda_genes) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Load necessary libraries (DT and ggplot2 were requested in the integration code)
    library(DT)
    library(ggplot2)
    library(openxlsx)
    
    # Reactive value to store the ranking results
    rank_data <- reactiveVal(NULL)
    
    # Observe the action button to trigger the ranking process
    # The action button ID is now 'rank'
    observeEvent(input$rank, {
      
      # Get the cluster data (which now includes Degree and Betweenness)
      clusters <- ppi_clusters()
      
      # Ensure cluster data is available and has the necessary columns
      if (is.null(clusters) || !all(c("Gene", "Cluster", "Degree", "Betweenness") %in% names(clusters))) {
        showNotification("No PPI cluster data with topological features found! Please run PPI analysis first.", type = "error")
        return(NULL)
      }
      
      # Start with the cluster data frame
      df <- clusters
      
      # Get the list of genes found in OpenFDA
      fda_genes_list <- openfda_genes()
      
      # Add 'FDA_Link' feature: TRUE if the gene is in the OpenFDA list, FALSE otherwise
      # This serves as a drug-target indicator feature for hub ranking.
      df$FDA_Link <- df$Gene %in% fda_genes_list
      
      # Calculate a custom "HubScore" using topological features (Degree, Betweenness) and the FDA link status.
      # This implements the ML-like ranking logic without an external XGBoost model file, 
      # as requested in the integration instructions.
      df$HubScore <- df$Degree + 
        log(df$Betweenness + 1) + 
        ifelse(df$FDA_Link, 2, 0) # Add a boost if the gene is linked to an FDA drug
      
      # Sort the results by HubScore in descending order
      df <- df[order(-df$HubScore), ]
      
      # Store the results in the reactive value
      rank_data(df)
      
      # Render the Hub Ranking table
      output$rank_table <- DT::renderDT(df, options = list(scrollX = TRUE))
      
      # Render the Hub Ranking plot (Top 10)
      output$rank_plot <- renderPlot({
        top <- head(df, 10)
        ggplot(top, aes(x = reorder(Gene, HubScore), y = HubScore)) +
          geom_bar(stat = "identity", fill = "steelblue") +
          coord_flip() +
          labs(x = "Gene", y = "Hub Score") +
          theme_minimal()
      })
    })
    
    # --- Download Handlers ---
    
    # Download CSV
    output$dl_rank_csv <- downloadHandler(
      filename = function() { "ml_hub_ranking.csv" },
      content = function(file) {
        write.csv(rank_data(), file, row.names = FALSE)
      }
    )
    
    # Download XLSX
    output$dl_rank_xlsx <- downloadHandler(
      filename = function() { "ml_hub_ranking.xlsx" },
      content = function(file) {
        openxlsx::write.xlsx(rank_data(), file)
      }
    )
  })
}