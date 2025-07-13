# modules/upload_module.R

uploadUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "upload",
          fluidPage(
            h3("Upload Gene Lists"),
            fileInput(ns("files"), "Upload one or more files (.txt or .csv):", multiple = TRUE),
            selectInput(ns("id_type"), "ID Type:",
                        choices = c("SYMBOL", "ENSEMBL", "UNIPROT", "ENTREZID", "CUSTOM")),
            actionButton(ns("load"), "Load Data"),
            br(),
            verbatimTextOutput(ns("preview"))
          )
  )
}

uploadServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    seq_lists <- reactiveVal(NULL)
    
    observeEvent(input$load, {
      req(input$files)
      
      lists <- lapply(seq_len(nrow(input$files)), function(i) {
        fpath <- input$files$datapath[i]
        ext <- tools::file_ext(fpath)
        
        if (ext == "csv") {
          # Assume 1st column
          read.csv(fpath, header = TRUE)[, 1]
        } else {
          # Plain text list
          readLines(fpath, warn = FALSE)
        }
      })
      
      seq_lists(lists)
      showNotification(paste("Loaded", length(lists), "file(s)."), type = "message")
    })
    
    output$preview <- renderPrint({
      req(seq_lists())
      lapply(seq_lists(), function(x) {
        list(
          First5 = head(x, 5),
          Total_N = length(x)
        )
      })
    })
    
    
    ## Return both the lists and the selected ID type
    return(list(
      seq_lists = seq_lists,
      id_type = reactive(input$id_type)
    ))
  })
}
