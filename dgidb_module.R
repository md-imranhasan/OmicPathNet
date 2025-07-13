# modules/dgidb_module.R

dgidbUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "dgidb",
          fluidPage(
            h3("💊 DGIdb Drug–Gene Interactions"),
            actionButton(ns("run"), "Query DGIdb"),
            br(), br(),
            DT::DTOutput(ns("dgidb_table")),
            downloadButton(ns("dl_dgidb_csv"), "Download (.csv)"),
            downloadButton(ns("dl_dgidb_xlsx"), "Download (.xlsx)"),
            br(), br(),
            h4("Drug–Gene Interaction Network"),
            visNetwork::visNetworkOutput(ns("dgidb_net"), height = "700px"),
            downloadButton(ns("dl_net_png"), "Download Network (.png)")
          )
  )
}

dgidbServer <- function(id, ppi_clusters) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    library(httr)
    library(jsonlite)
    library(DT)
    library(openxlsx)
    library(visNetwork)
    library(webshot2)
    
    dgidb_results <- reactiveVal(NULL)
    vis_data <- reactiveVal(NULL)
    
    observeEvent(input$run, {
      clusters <- ppi_clusters()
      if (is.null(clusters)) {
        showNotification("No clusters found!", type = "error")
        return(NULL)
      }
      
      genes <- unique(na.omit(clusters$Gene))
      if (length(genes) == 0) {
        showNotification("No genes found in clusters.", type = "error")
        return(NULL)
      }
      
      if (length(genes) > 50) {
        genes <- genes[1:50]
        showNotification("Limiting to first 50 genes for DGIdb.", type = "message")
      }
      
      query <- '
      query($genes: [String!]) {
        genes(names: $genes) {
          nodes {
            name
            interactions {
              drug { name }
              interactionScore
              interactionTypes { type }
              sources { sourceDbName }
            }
          }
        }
      }'
      
      body <- list(query = query, variables = list(genes = genes))
      
      resp <- httr::POST(
        "https://dgidb.org/api/graphql",
        encode = "json",
        body = body,
        httr::add_headers("Content-Type" = "application/json")
      )
      
      if (httr::http_error(resp)) {
        showNotification("DGIdb GraphQL call failed!", type = "error")
        return(NULL)
      }
      
      raw <- httr::content(resp, as = "text", encoding = "UTF-8")
      data <- jsonlite::fromJSON(raw, simplifyVector = FALSE)
      
      if (is.null(data$data) || is.null(data$data$genes$nodes)) {
        showNotification("No DGIdb interactions found.", type = "warning")
        dgidb_results(NULL)
        return(NULL)
      }
      
      nodes <- data$data$genes$nodes
      if (!is.list(nodes) || length(nodes) == 0) {
        showNotification("No DGIdb interactions found.", type = "warning")
        dgidb_results(NULL)
        return(NULL)
      }
      
      all_hits <- lapply(nodes, function(node) {
        gene <- node$name %||% NA
        interactions <- node$interactions
        if (is.null(interactions) || length(interactions) == 0) return(NULL)
        
        do.call(rbind, lapply(interactions, function(x) {
          data.frame(
            Gene = gene,
            Drug = x$drug$name %||% NA,
            InteractionTypes = if (!is.null(x$interactionTypes) && length(x$interactionTypes) > 0) {
              paste(sapply(x$interactionTypes, function(t) t$type), collapse = "; ")
            } else { NA },
            Sources = if (!is.null(x$sources) && length(x$sources) > 0) {
              paste(sapply(x$sources, function(s) s$sourceDbName), collapse = "; ")
            } else { NA },
            Score = x$interactionScore %||% NA,
            stringsAsFactors = FALSE
          )
        }))
      })
      
      final_df <- do.call(rbind, all_hits)
      if (is.null(final_df) || nrow(final_df) == 0) {
        showNotification("No interactions parsed.", type = "warning")
        dgidb_results(NULL)
        return(NULL)
      }
      
      dgidb_results(final_df)
      
      ## === Build visNetwork nodes & edges ===
      genes <- unique(final_df$Gene)
      drugs <- unique(final_df$Drug)
      
      nodes <- data.frame(
        id = c(genes, drugs),
        label = c(genes, drugs),
        group = c(rep("Gene", length(genes)), rep("Drug", length(drugs))),
        shape = c(rep("ellipse", length(genes)), rep("box", length(drugs))),
        font.size = 28
      )
      
      edges <- data.frame(
        from = final_df$Gene,
        to = final_df$Drug,
        title = paste0("InteractionTypes: ", final_df$InteractionTypes, "<br>",
                       "Sources: ", final_df$Sources, "<br>",
                       "Score: ", round(as.numeric(final_df$Score), 2)),
        value = as.numeric(final_df$Score) / max(as.numeric(final_df$Score), na.rm = TRUE)
      )
      
      vis_data(list(nodes = nodes, edges = edges))
      
      output$dgidb_net <- visNetwork::renderVisNetwork({
        visNetwork::visNetwork(nodes, edges) %>%
          visIgraphLayout(layout = "layout_with_fr") %>%
          visNodes(font = list(size = 28)) %>%
          visEdges(smooth = TRUE) %>%
          visOptions(highlightNearest = TRUE) %>%
          visLegend()
      })
      
      output$dgidb_table <- renderDT({
        datatable(final_df, options = list(scrollX = TRUE))
      })
    })
    
    ## Downloads
    output$dl_dgidb_csv <- downloadHandler(
      filename = function() { "dgidb_interactions.csv" },
      content = function(file) {
        write.csv(dgidb_results(), file, row.names = FALSE)
      }
    )
    output$dl_dgidb_xlsx <- downloadHandler(
      filename = function() { "dgidb_interactions.xlsx" },
      content = function(file) {
        openxlsx::write.xlsx(dgidb_results(), file)
      }
    )
    
    ## Download network snapshot
    output$dl_net_png <- downloadHandler(
      filename = function() { "dgidb_network.png" },
      content = function(file) {
        if (is.null(vis_data())) return(NULL)
        tmp_html <- tempfile(fileext = ".html")
        visNetwork::visSave(
          visNetwork::visNetwork(vis_data()$nodes, vis_data()$edges) %>%
            visIgraphLayout(layout = "layout_with_fr") %>%
            visNodes(font = list(size = 28)) %>%
            visEdges(smooth = TRUE) %>%
            visOptions(highlightNearest = TRUE) %>%
            visLegend(),
          file = tmp_html
        )
        Sys.sleep(1)
        webshot2::webshot(url = tmp_html, file = file, vwidth = 1800, vheight = 1400)
      }
    )
  })
}

# Safe pipe helper
`%||%` <- function(x, y) if (is.null(x)) y else x
