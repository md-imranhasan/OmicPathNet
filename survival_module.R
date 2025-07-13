# modules/survival_module.R

survivalUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "survival",
          fluidPage(
            h3("⏳ Survival Analysis (TCGA)"),
            
            # Cancer type selection (expanded options for flexibility)
            selectInput(ns("cancer"), "Select TCGA Cancer Type:",
                        choices = c("BRCA" = "TCGA-BRCA", "LUAD" = "TCGA-LUAD",
                                    "LAML" = "TCGA-LAML", "COAD" = "TCGA-COAD",
                                    "STAD" = "TCGA-STAD", "KIRC" = "TCGA-KIRC")),
            
            # Analysis type selection (Single Gene vs. Gene Module)
            radioButtons(ns("analysis_type"), "Analysis Type:",
                         choices = c("Single Gene" = "single", "Gene Module/Cluster" = "module"),
                         selected = "single", inline = TRUE),
            
            # Input for single gene (dynamically updated)
            conditionalPanel(
              condition = "input.analysis_type == 'single'",
              selectInput(ns("gene"), "Select Gene:", choices = NULL)
            ),
            
            # Input for gene module (shows all genes in the provided list)
            conditionalPanel(
              condition = "input.analysis_type == 'module'",
              div(class = "well", 
                  p("Module analysis uses the median expression of all genes in the current list/cluster."),
                  textOutput(ns("module_gene_count"))
              )
            ),
            
            actionButton(ns("run"), "Run Survival Analysis"),
            br(), br(),
            
            # Status messages for data loading
            textOutput(ns("status_message")), 
            
            plotOutput(ns("kmplot"), height = "600px"),
            DT::DTOutput(ns("surv_table")),
            downloadButton(ns("dl_table"), "Download (.csv)")
          )
  )
}

survivalServer <- function(id, gene_list) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    library(TCGAbiolinks)
    library(SummarizedExperiment)
    library(survival)
    library(survminer)
    library(DT)
    library(dplyr)
    library(tidyr)
    
    # Reactive value to cache the loaded TCGA data (SummarizedExperiment)
    tcga_data_cache <- reactiveVal(list())
    
    # Dynamic gene choices update and module gene count display
    observe({
      req(gene_list())
      genes <- gene_list()
      updateSelectInput(session, "gene", choices = genes)
      output$module_gene_count <- renderText({
        paste("Selected gene list/module contains", length(genes), "genes.")
      })
    })
    
    # Reactive expression to handle data loading and caching
    loaded_tcga_data <- reactive({
      cancer_type <- input$cancer
      
      # 1. Check cache
      if (!is.null(tcga_data_cache()[[cancer_type]])) {
        output$status_message <- renderText(paste("✅ Data for", cancer_type, "loaded from cache."))
        return(tcga_data_cache()[[cancer_type]])
      }
      
      # 2. If not cached, download and prepare data robustly
      output$status_message <- renderText(paste("⬇️ Downloading and preparing data for", cancer_type, "from GDC... This may take a few minutes."))
      
      data <- withProgress(message = paste("Loading TCGA data for", cancer_type), value = 0, {
        
        incProgress(0.1, detail = "Querying available workflows...")
        
        # Identify available workflows for gene expression quantification
        meta_query <- tryCatch({
          GDCquery(
            project = cancer_type,
            data.category = "Transcriptome Profiling",
            data.type = "Gene Expression Quantification"
          )
        }, error = function(e) {
          showNotification(paste("❌ Error querying GDC for project", cancer_type, ":", e$message), type = "error")
          return(NULL)
        })
        
        if (is.null(meta_query)) return(NULL)
        
        # Determine valid workflows and select the preferred one
        workflows <- unique(meta_query$results[[1]]$workflow.type)
        
        # Prioritize workflows known for reliable FPKM/counts
        preferred_workflows <- c("STAR - Counts", "HTSeq - Counts", "HTSeq - FPKM-UQ", "HTSeq - FPKM")
        
        selected_workflow <- intersect(preferred_workflows, workflows)
        
        if (length(selected_workflow) == 0) {
          output$status_message <- renderText("❌ No suitable gene expression workflow found for this project.")
          return(NULL)
        }
        
        selected_workflow <- selected_workflow[1]
        incProgress(0.3, detail = paste("✅ Using workflow:", selected_workflow))
        
        # Define the final query and download
        query <- GDCquery(
          project = cancer_type,
          data.category = "Transcriptome Profiling",
          data.type = "Gene Expression Quantification",
          workflow.type = selected_workflow
        )
        
        incProgress(0.5, detail = "Downloading data...")
        GDCdownload(query)
        
        incProgress(0.8, detail = "Preparing data...")
        prepared_data <- GDCprepare(query)
        
        incProgress(1, detail = "Done!")
        return(prepared_data)
      })
      
      # Store in cache
      if (!is.null(data)) {
        current_cache <- tcga_data_cache()
        current_cache[[cancer_type]] <- data
        tcga_data_cache(current_cache)
      }
      
      return(data)
    })
    
    # 4. Survival Analysis (Triggered by 'Run' button)
    surv_data <- reactiveVal(NULL)
    
    observeEvent(input$run, {
      req(loaded_tcga_data())
      
      # Determine the gene set based on analysis type
      if (input$analysis_type == "single") {
        req(input$gene)
        target_genes <- input$gene
        analysis_name <- input$gene
      } else {
        req(gene_list())
        target_genes <- gene_list()
        analysis_name <- "Gene Module"
      }
      
      withProgress(message = paste("Analyzing survival for", analysis_name), value = 0, {
        
        data <- loaded_tcga_data()
        expr <- assay(data)
        clinical <- colData(data)
        
        # Filter expression matrix to only target genes
        common_genes <- intersect(target_genes, rownames(expr))
        if (length(common_genes) == 0) {
          showNotification("Selected gene(s) not found in the expression data.", type = "error")
          return(NULL)
        }
        
        incProgress(0.4, detail = "Calculating signature score...")
        
        # Calculate expression score for the analysis
        if (input$analysis_type == "single") {
          gene_expr_score <- as.numeric(expr[common_genes, ])
        } else {
          # For a module, use the median expression across all genes in the module
          # Use t(expr[common_genes, ]) to get samples in rows and genes in columns
          gene_expr_score <- apply(t(expr[common_genes, , drop = FALSE]), 1, median)
        }
        
        # Prepare survival dataframe
        incProgress(0.6, detail = "Preparing survival data...")
        df <- data.frame(
          Sample = names(gene_expr_score),
          Expression = gene_expr_score,
          OS.time = clinical$days_to_death,
          OS.event = ifelse(clinical$vital_status == "Dead", 1, 0)
        )
        
        # Filter for valid survival times and dichotomize (median split)
        df <- df[!is.na(df$OS.time) & df$OS.time > 0, ]
        df$Group <- ifelse(df$Expression >= median(df$Expression, na.rm = TRUE), "High", "Low")
        
        incProgress(0.8, detail = "Fitting survival model...")
        
        # Run survival analysis (Kaplan-Meier and Cox model)
        fit <- survfit(Surv(OS.time, OS.event) ~ Group, data = df)
        
        output$kmplot <- renderPlot({
          ggsurvplot(fit, data = df,
                     pval = TRUE,
                     risk.table = TRUE,
                     ggtheme = theme_minimal(),
                     legend.title = analysis_name)
        })
        
        cox <- coxph(Surv(OS.time, OS.event) ~ Expression, data = df)
        summary_cox <- summary(cox)
        
        # Prepare results table
        out_tab <- data.frame(
          Analysis = analysis_name,
          HR = round(summary_cox$coefficients[,"exp(coef)"], 3),
          p.value = summary_cox$coefficients[,"Pr(>|z|)"]
        )
        surv_data(out_tab)
        output$surv_table <- renderDT(out_tab)
        
        incProgress(1, detail = "Done!")
      })
    })
    
    # Download handler for results table
    output$dl_table <- downloadHandler(
      filename = function() "survival_results.csv",
      content = function(file) {
        write.csv(surv_data(), file, row.names = FALSE)
      }
    )
  })
}