# app.R
options(shiny.maxRequestSize = 1024 * 1024 * 1000)
# ✅ Load required libraries
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(UpSetR)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(STRINGdb)
library(purrr)
library(circlize)
library(stringr)
library(visNetwork)
library(igraph)
library(webshot2)
library(poweRlaw)
library(rmarkdown)
library(htmltools)
library(WGCNA)
library(Seurat)
library(DESeq2)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tidyverse)

# ✅ Load your connected modules
source("modules/upload_module.R")
source("modules/scrnaseq_module.R")
source("modules/enrichment_module.R")
source("modules/crosstalk_module.R")
source("modules/ppi_module.R")
source("modules/dgidb_module.R")
source("modules/structure_module.R")
source("modules/mlpredict_module.R")
source("modules/survival_module.R")
source("modules/openfda_module.R")
source("modules/litcooc_module.R")
source("modules/report_module.R")

# ✅ UI
ui <- dashboardPage(
  dashboardHeader(title = "OmicPathNet"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload Data", tabName = "upload"),
      menuItem("scRNA-seq", tabName = "scrna"),
      menuItem("Enrichment", tabName = "enrichment"),
      menuItem("Pathway Crosstalk", tabName = "crosstalk"),
      menuItem("PPI Network", tabName = "ppi"),
      menuItem("DGIdb Drug Targets", tabName = "dgidb"),
      menuItem("Structures", tabName = "structure"),
      menuItem("ML Hub Ranking", tabName = "mlpredict"),
      menuItem("Survival Analysis", tabName = "survival"),
      menuItem("OpenFDA", tabName = "openfda"),
      menuItem("Lit Co-Occurrence", tabName = "litcooc"),
      menuItem("Download Report", tabName = "report")
    )
  ),
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    tabItems(
      uploadUI("upload_ui"),
      scRNAseqUI("scrna_ui"),
      enrichmentUI("enrich_ui"),
      crosstalkUI("crosstalk_ui"),
      ppiUI("ppi_ui"),
      dgidbUI("dgidb_ui"),
      structureUI("structure_ui"),
      mlpredictUI("mlpredict_ui"),
      survivalUI("survival_ui"),
      openfdaUI("openfda_ui"),
      litcoocUI("litcooc_ui"),
      reportUI("report_ui")
    )
  )
)

# ✅ Server
server <- function(input, output, session) {
  ## 1️⃣ Upload: base module
  user_data <- uploadServer("upload_ui")
  
  ## 2️⃣ scRNA-seq: DEGs module
  scrna_genes <- scRNAseqServer("scrna_ui")
  
  ## 3️⃣ ✅ Master Gene Container — smart multi-list aware
  master_gene_container <- reactive({
    if (!is.null(user_data$seq_lists()) && length(user_data$seq_lists()) > 0) {
      # Multiple uploaded lists → use them AS-IS for enrichment & PPI
      list(seq_lists = user_data$seq_lists,
           id_type = user_data$id_type)
    } else if (!is.null(scrna_genes()) && length(scrna_genes()) > 0) {
      # Fallback: scRNA-seq DEGs as single list
      list(seq_lists = reactive(list(scrna_genes())),
           id_type = user_data$id_type)
    } else {
      # Fallback empty
      list(seq_lists = reactive(list()), id_type = user_data$id_type)
    }
  })
  
  ## 4️⃣ ✅ Enrichment: works with multi-lists!
  enrich_results <- enrichmentServer("enrich_ui", master_gene_container())
  
  ## 5️⃣ ✅ Crosstalk
  crosstalkServer("crosstalk_ui", master_gene_container())
  
  ## 6️⃣ ✅ PPI
  ppi_results <- ppiServer("ppi_ui", master_gene_container())
  
  ## 7️⃣ ✅ PPI hub cluster & fallback gene list
  ppi_hubs_data <- reactive({
    ppi_results$cluster_hubs()
  })
  ppi_gene_list <- reactive({
    if (!is.null(ppi_hubs_data()) && nrow(ppi_hubs_data()) > 0) {
      unique(ppi_hubs_data()$Gene)
    } else {
      # fallback to master merged if no clusters yet
      unlist(master_gene_container()$seq_lists())
    }
  })
  
  ## 8️⃣ ✅ DGIdb: uses PPI cluster data
  dgidbServer("dgidb_ui", ppi_results$cluster_data)
  
  ## 9️⃣ ✅ Structures: uses PPI clusters
  structureServer("structure_ui", ppi_results$cluster_data)
  
  ## 🔟 ✅ ML Hub Ranking
  mlpredictServer("mlpredict_ui", ppi_results$cluster_data, ppi_gene_list)
  
  ## 1️⃣1️⃣ ✅ Survival
  survivalServer("survival_ui", ppi_gene_list)
  
  ## 1️⃣2️⃣ ✅ OpenFDA
  openfdaServer("openfda_ui", ppi_gene_list)
  
  ## 1️⃣3️⃣ ✅ Literature Co-Occurrence
  litcoocServer("litcooc_ui", ppi_results$cluster_data)
  
  ## 1️⃣4️⃣ ✅ Report
  reportServer("report_ui",
               enrich_results,
               ppi_results$cluster_data,
               ppi_results$ppi_table_data)
}


# ✅ Run the app
shinyApp(ui, server)