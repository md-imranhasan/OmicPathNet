structureUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "structure",
          fluidPage(
            h3("🧬 UniProt & AlphaFold Structure Coverage"),
            actionButton(ns("run"), "Check Structures"),
            br(), br(),
            DT::DTOutput(ns("structure_table")),
            downloadButton(ns("dl_structure_csv"), "Download (.csv)"),
            downloadButton(ns("dl_structure_xlsx"), "Download (.xlsx)")
          )
  )
}

structureServer <- function(id, ppi_clusters) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    library(httr)
    library(jsonlite)
    library(openxlsx)
    library(DT)
    
    struct_data <- reactiveVal(NULL)
    
    observeEvent(input$run, {
      clusters <- ppi_clusters()
      if (is.null(clusters) || nrow(clusters) == 0) {
        showNotification("No clusters found!", type = "error")
        return(NULL)
      }
      
      genes <- unique(clusters$Gene)
      showNotification(paste("Checking structures for", length(genes), "genes..."), type = "message")
      
      withProgress(message = "Fetching UniProt data...", value = 0, {
        all_results <- lapply(seq_along(genes), function(i) {
          gene <- genes[i]
          incProgress(1/length(genes), detail = paste("Gene:", gene))
          
          af_link <- paste0("https://alphafold.ebi.ac.uk/entry/", gene)
          uniprot_id <- NA
          uniprot_name <- "NA"
          uniprot_length <- "NA"
          uniprot_reviewed <- "NA"
          
          tryCatch({
            q <- paste0("gene_exact:", gene, " AND organism_id:9606")
            base_url <- "https://rest.uniprot.org/uniprotkb/search?"
            params <- list(
              query = q,
              format = "json",
              fields = "accession,protein_name,length,reviewed",
              size = 1
            )
            full_url <- paste0(base_url, URLencode(paste(names(params), params, sep = "=", collapse = "&")))
            res <- GET(full_url)
            if (!http_error(res)) {
              js <- fromJSON(content(res, "text", encoding = "UTF-8"))
              if (!is.null(js$results) && nrow(js$results) > 0) {
                uniprot_id <- js$results$primaryAccession[1]
                if (!is.null(js$results$proteinDescription$recommendedName$fullName$value)) {
                  uniprot_name <- js$results$proteinDescription$recommendedName$fullName$value[1]
                }
                if (!is.null(js$results$sequence$length)) {
                  uniprot_length <- as.character(js$results$sequence$length[1])
                }
                if (!is.null(js$results$entryType)) {
                  uniprot_reviewed <- js$results$entryType[1]
                }
              }
            }
          }, error = function(e) { })
          
          data.frame(
            Gene = gene,
            UniProt_ID = ifelse(is.na(uniprot_id), "NA", uniprot_id),
            UniProt_Protein = uniprot_name,
            Protein_Length = uniprot_length,
            Reviewed = uniprot_reviewed,
            UniProt_Link = ifelse(is.na(uniprot_id), "NA", paste0("https://www.uniprot.org/uniprotkb/", uniprot_id)),
            AlphaFold_Link = af_link,
            stringsAsFactors = FALSE
          )
        })
      })
      
      df <- do.call(rbind, all_results)
      struct_data(df)
      
      output$structure_table <- renderDT({
        datatable(df, escape = FALSE, options = list(scrollX = TRUE),
                  rownames = FALSE,
                  callback = DT::JS(  # ✅ DT::JS!
                    "table.rows().every(function(){",
                    "  var data = this.data();",
                    "  if (data[1] !== 'NA') {",
                    "    data[1] = '<a href=\"' + data[5] + '\" target=\"_blank\">' + data[1] + '</a>';",
                    "  }",
                    "  if (data[5] !== 'NA') {",
                    "    data[5] = '<a href=\"' + data[5] + '\" target=\"_blank\">UniProt</a>';",
                    "  }",
                    "  if (data[6] !== 'NA') {",
                    "    data[6] = '<a href=\"' + data[6] + '\" target=\"_blank\">AlphaFold</a>';",
                    "  }",
                    "  this.invalidate();",
                    "});"
                  ))
      })
    })
    
    output$dl_structure_csv <- downloadHandler(
      filename = function() { "structure_table.csv" },
      content = function(file) {
        write.csv(struct_data(), file, row.names = FALSE)
      }
    )
    
    output$dl_structure_xlsx <- downloadHandler(
      filename = function() { "structure_table.xlsx" },
      content = function(file) {
        openxlsx::write.xlsx(struct_data(), file)
      }
    )
  })
}
