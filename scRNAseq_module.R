# modules/scRNAseq_module.R

scRNAseqUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "scrna",
          fluidPage(
            h3("🔬 scRNA-seq Explorer (Raw 10X / Preprocessed)"),
            
            radioButtons(ns("input_type"), "Data Type:",
                         choices = c("Preprocessed Seurat (.rds)" = "rds",
                                     "Raw 10X Folder (.zip)" = "zip")),
            
            conditionalPanel(
              condition = sprintf("input['%s'] == 'rds'", ns("input_type")),
              fileInput(ns("seurat_rds"), "Upload Seurat Object (.rds):")
            ),
            
            conditionalPanel(
              condition = sprintf("input['%s'] == 'zip'", ns("input_type")),
              fileInput(ns("raw_zip"), "Upload 10X Data Folder (.zip):")
            ),
            
            numericInput(ns("resolution"), 
                         "Clustering Resolution:", 
                         value = 0.5, min = 0.1, step = 0.1),
            
            numericInput(ns("dims"), 
                         "Number of PCA/UMAP Dimensions:", 
                         value = 10, min = 5, step = 1),
            
            actionButton(ns("load"), "Load & Process"),
            br(), br(),
            
            h4("📌 UMAP"),
            plotOutput(ns("umap"), height = "600px"),
            downloadButton(ns("dl_umap_png"), "Download UMAP (.png)"),
            br(), br(),
            
            h4("📌 Cluster Heatmap (Top 5 per cluster)"),
            plotOutput(ns("heatmap"), height = "600px"),
            downloadButton(ns("dl_heatmap_png"), "Download Heatmap (.png)"),
            br(), br(),
            
            h4("📌 FeaturePlot of Top Gene"),
            plotOutput(ns("featureplot"), height = "600px"),
            downloadButton(ns("dl_featureplot_png"), "Download FeaturePlot (.png)"),
            br(), br(),
            
            h4("📌 ViolinPlot of Selected Gene"),
            selectizeInput(ns("violin_gene"), "Select Gene:", choices = NULL, multiple = FALSE),
            plotOutput(ns("violinplot"), height = "400px"),
            downloadButton(ns("dl_violin_png"), "Download ViolinPlot (.png)"),
            br(), br(),
            
            h4("📌 Marker Genes"),
            DT::DTOutput(ns("marker_table")),
            downloadButton(ns("dl_markers_csv"), "Download Markers (.csv)"),
            downloadButton(ns("dl_markers_xlsx"), "Download Markers (.xlsx)")
          ))
}

scRNAseqServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    library(Seurat)
    library(DT)
    library(openxlsx)
    library(glue)
    library(tidyverse)
    
    sc_object <- reactiveVal(NULL)
    marker_data <- reactiveVal(NULL)
    top5_genes <- reactiveVal(NULL)
    
    observeEvent(input$load, {
      withProgress(message = "Processing scRNA-seq data", value = 0, {
        
        dims <- as.integer(input$dims)
        resolution <- input$resolution
        
        # ------------------------------
        # 1. Load data
        # ------------------------------
        if (input$input_type == "rds") {
          incProgress(0.1, detail = "Reading RDS...")
          req(input$seurat_rds)
          obj <- readRDS(input$seurat_rds$datapath)
          
        } else if (input$input_type == "zip") {
          incProgress(0.1, detail = "Unzipping...")
          req(input$raw_zip)
          temp_dir <- tempfile()
          unzip(input$raw_zip$datapath, exdir = temp_dir)
          data_dir <- list.dirs(temp_dir, recursive = FALSE)[1]
          
          incProgress(0.2, detail = "Reading 10X data...")
          counts <- Read10X(data.dir = data_dir)
          obj <- CreateSeuratObject(counts = counts)
        } else {
          showNotification("Invalid input type.", type = "error")
          return(NULL)
        }
        
        # ------------------------------
        # 2. Preprocessing
        # ------------------------------
        incProgress(0.4, detail = "QC & Preprocessing...")
        obj <- NormalizeData(obj)
        obj <- FindVariableFeatures(obj)
        obj <- ScaleData(obj)
        obj <- RunPCA(obj)
        obj <- RunUMAP(obj, dims = 1:dims)
        obj <- FindNeighbors(obj, dims = 1:dims)
        obj <- FindClusters(obj, resolution = resolution)
        
        # ✅ Safe: make cluster identities active
        Idents(obj) <- obj$seurat_clusters
        
        sc_object(obj)
        showNotification(glue("Seurat object ready! Resolution = {resolution}, Dims = {dims}"), type = "message")
        
        # ------------------------------
        # 3. Find DEGs
        # ------------------------------
        incProgress(0.7, detail = "Finding DEGs...")
        markers <- FindAllMarkers(obj, assay = "RNA", logfc.threshold = 0.25, only.pos = TRUE)
        markers <- markers %>% distinct(gene, .keep_all = TRUE)  # remove duplicates
        marker_data(markers)
        
        top5 <- markers %>%
          group_by(cluster) %>%
          top_n(n = 5, wt = avg_log2FC)
        top5_genes(top5)
        
        showNotification(glue("Found {nrow(markers)} unique DEGs"), type = "message")
        
        # ✅ Update violin gene input with server-side selectize
        updateSelectizeInput(session, "violin_gene",
                             choices = unique(markers$gene),
                             server = TRUE)
        
        # ------------------------------
        # 4. Render Plots
        # ------------------------------
        output$umap <- renderPlot({
          DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE)
        })
        
        output$heatmap <- renderPlot({
          DoHeatmap(obj, features = top5$gene) + NoLegend()
        })
        
        output$featureplot <- renderPlot({
          FeaturePlot(obj, features = top5$gene[1])
        })
        
        output$violinplot <- renderPlot({
          req(input$violin_gene)
          VlnPlot(obj, features = input$violin_gene)
        })
        
        output$marker_table <- renderDT({
          datatable(markers, options = list(scrollX = TRUE))
        })
        
        incProgress(1, detail = "Done!")
      })
    })
    
    # ------------------------------
    # 5. Downloads
    # ------------------------------
    output$dl_umap_png <- downloadHandler(
      filename = function() "scRNA_UMAP.png",
      content = function(file) {
        png(file, width = 1200, height = 900)
        DimPlot(sc_object(), reduction = "umap", label = TRUE, repel = TRUE)
        dev.off()
      }
    )
    output$dl_heatmap_png <- downloadHandler(
      filename = function() "scRNA_Heatmap.png",
      content = function(file) {
        png(file, width = 1200, height = 900)
        DoHeatmap(sc_object(), features = top5_genes()$gene) + NoLegend()
        dev.off()
      }
    )
    output$dl_featureplot_png <- downloadHandler(
      filename = function() "scRNA_FeaturePlot.png",
      content = function(file) {
        png(file, width = 1200, height = 900)
        FeaturePlot(sc_object(), features = top5_genes()$gene[1])
        dev.off()
      }
    )
    output$dl_violin_png <- downloadHandler(
      filename = function() "scRNA_ViolinPlot.png",
      content = function(file) {
        png(file, width = 1200, height = 900)
        VlnPlot(sc_object(), features = input$violin_gene)
        dev.off()
      }
    )
    output$dl_markers_csv <- downloadHandler(
      filename = function() "scRNA_markers.csv",
      content = function(file) write.csv(marker_data(), file, row.names = FALSE)
    )
    output$dl_markers_xlsx <- downloadHandler(
      filename = function() "scRNA_markers.xlsx",
      content = function(file) openxlsx::write.xlsx(marker_data(), file)
    )
    
    # ✅ Return unique DEGs for downstream
    return(reactive({
      req(marker_data())
      unique(marker_data()$gene)
    }))
  })
}
