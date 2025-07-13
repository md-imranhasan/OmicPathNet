# modules/report_module.R
library(shiny)
library(DT)
library(openxlsx)
library(rmarkdown)
library(htmltools)

reportUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "report",
          fluidPage(
            h3("📄 Comprehensive Report Generation"),
            
            # Button for generating the comprehensive report
            downloadButton(ns("generate_report"), "Generate Comprehensive Report (HTML/PDF)"),
            p("Note: This report will include summaries of Enrichment and PPI analysis."),
            
            hr(),
            
            # Interactive tables for quick review (similar to previous)
            h4("Interactive Data Tables:"),
            
            h5("Enrichment Summary (Example List)"),
            DT::DTOutput(ns("enrich_table")),
            br(),
            
            h5("PPI Cluster Assignments"),
            DT::DTOutput(ns("ppi_clusters")),
            downloadButton(ns("dl_ppi_clusters"), "Download Clusters (.csv)"),
            br(),
            
            h5("PPI Interaction Table"),
            DT::DTOutput(ns("ppi_interactions")),
            downloadButton(ns("dl_ppi_table"), "Download PPI Interactions (.csv)")
          )
  )
}

reportServer <- function(id, enrich_results, ppi_clusters, ppi_table) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # --- Interactive Tables (Quick View) ---
    
    ## Enrichment summary (Example: List 1)
    output$enrich_table <- DT::renderDT({
      # Assuming enrich_results contains a list of data frames from the enrichment server
      # We extract results for the first list (List 1) if available.
      
      # The enrich_results object passed from app.R is complex, checking for a key like 'ekegg1' or 'df'
      # If enrich_results is a list of results as returned by enrichmentServer, we can access them.
      
      # Fallback mechanism if enrich_results is structured differently than expected
      if (is.null(enrich_results) || length(enrich_results) == 0) {
        return(data.frame(Message = "No enrichment results available."))
      }
      
      # Since we only have the 'enrich_results' object from app.R which seems to rely on the enrichmentServer
      # We rely on 'enrich_results$kegg1' as defined in the original reportServer.
      
      if (!is.null(enrich_results$kegg1)) {
        datatable(as.data.frame(enrich_results$kegg1), options = list(scrollX = TRUE))
      } else {
        # Display a placeholder if no data is available
        datatable(data.frame(Message = "Run enrichment analysis first."), options = list(scrollX = TRUE))
      }
    })
    
    ## PPI Clusters
    output$ppi_clusters <- DT::renderDT({
      req(ppi_clusters())
      datatable(ppi_clusters(), options = list(scrollX = TRUE))
    })
    output$dl_ppi_clusters <- downloadHandler(
      filename = function() { "ppi_clusters.csv" },
      content = function(file) {
        write.csv(ppi_clusters(), file, row.names = FALSE)
      }
    )
    
    ## PPI Interactions
    output$ppi_interactions <- DT::renderDT({
      req(ppi_table())
      datatable(ppi_table(), options = list(scrollX = TRUE))
    })
    output$dl_ppi_table <- downloadHandler(
      filename = function() { "ppi_interactions.csv" },
      content = function(file) {
        write.csv(ppi_table(), file, row.names = FALSE)
      }
    )
    
    # --- Comprehensive Report Generation (RMarkdown) ---
    
    output$generate_report <- downloadHandler(
      filename = function() {
        paste("CoPathNet_Bioinformatics_Report-", Sys.Date(), ".html", sep="")
      },
      content = function(file) {
        
        # Ensure required data is available for the report
        req(ppi_clusters(), ppi_table(), !is.null(enrich_results$kegg1))
        
        # Path to the R Markdown template file
        # We assume the user has a report_template.Rmd file in the Shiny app directory.
        # This Rmd file will contain the layout and code to render the data.
        src <- normalizePath('report_template.Rmd')
        
        # Create a temporary directory for the generated report
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        
        # Parameters to pass to the Rmd document
        params <- list(
          enrichment_data = as.data.frame(enrich_results$kegg1),
          ppi_clusters = ppi_clusters(),
          ppi_interactions = ppi_table(),
          report_date = Sys.Date()
        )
        
        # Render the RMarkdown report
        # This will knit the Rmd template using the data parameters
        rmarkdown::render(
          src,
          output_file = file,
          params = params,
          envir = new.env(parent = globalenv())
        )
      }
    )
  })
}

# ⚠️ Note: For this code to work, you must create a file named `report_template.Rmd` 
# in your Shiny app's main directory. 
# This Rmd file will utilize the 'params' passed by the server function to generate the report.