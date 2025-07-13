# modules/openfda_module.R

openfdaUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "openfda",
          fluidPage(
            h3("🔎 OpenFDA Connected to PPI Clusters"),
            checkboxGroupInput(ns("selected_genes"), "Select Genes to Query:", choices = NULL),
            selectInput(ns("endpoint"), "Data Type:",
                        choices = c("Drug Events" = "drug/event",
                                    "Drug Labeling" = "drug/label")),
            actionButton(ns("search"), "Search OpenFDA"),
            br(), br(),
            h4("OpenFDA Results"),
            DT::DTOutput(ns("fda_table")),
            downloadButton(ns("dl_fda_csv"), "Download (.csv)"),
            downloadButton(ns("dl_fda_xlsx"), "Download (.xlsx)")
          )
  )
}

openfdaServer <- function(id, ppi_hubs) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    library(httr)
    library(jsonlite)
    library(openxlsx)
    library(dplyr)  # ✅ Use modern bind_rows
    
    ## Update gene choices when PPI hubs change
    observe({
      req(ppi_hubs())
      updateCheckboxGroupInput(session, "selected_genes",
                               choices = ppi_hubs(),
                               selected = ppi_hubs())
    })
    
    results <- reactiveVal(NULL)
    
    observeEvent(input$search, {
      req(input$selected_genes)
      all_results <- list()
      
      for (gene in input$selected_genes) {
        base_url <- "https://api.fda.gov/"
        endpoint <- input$endpoint
        search_term <- URLencode(gene)
        
        url <- paste0(base_url, endpoint, ".json?search=",
                      if (endpoint == "drug/event") {
                        "patient.drug.medicinalproduct:"
                      } else {
                        "openfda.generic_name:"
                      },
                      search_term, "&limit=5")
        
        resp <- tryCatch({ GET(url) }, error = function(e) NULL)
        if (is.null(resp) || resp$status_code != 200) next
        
        json_data <- tryCatch({ fromJSON(content(resp, "text", encoding = "UTF-8")) }, error = function(e) NULL)
        
        if (!is.null(json_data$results)) {
          df <- as.data.frame(json_data$results)
          df$QueriedGene <- gene
          all_results[[gene]] <- df
        }
      }
      
      if (length(all_results) == 0) {
        showNotification("No FDA results found for selected genes.", type = "warning")
        results(NULL)
      } else {
        combined <- bind_rows(all_results)  # ✅ Modern robust bind
        results(combined)
      }
    })
    
    output$fda_table <- DT::renderDT({
      req(results())
      datatable(results(), options = list(scrollX = TRUE))
    })
    
    output$dl_fda_csv <- downloadHandler(
      filename = function() { "openfda_results.csv" },
      content = function(file) {
        write.csv(results(), file, row.names = FALSE)
      }
    )
    
    output$dl_fda_xlsx <- downloadHandler(
      filename = function() { "openfda_results.xlsx" },
      content = function(file) {
        openxlsx::write.xlsx(results(), file)
      }
    )
  })
}
