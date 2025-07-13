# modules/timeline_module.R
timelineUI <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "timeline",
          fluidPage(
            textInput(ns("term"), "Enter gene/pathway:"),
            actionButton(ns("trend"), "Fetch Trends"),
            plotOutput(ns("trend_plot"))
          )
  )
}

timelineServer <- function(id, ...) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$trend, {
      # Dummy plot. In practice, call rentrez or easyPubMed.
      trend_years <- 2000:2024
      counts <- sample(10:100, length(trend_years))
      output$trend_plot <- renderPlot({
        plot(trend_years, counts, type = "l", col = "blue",
             main = paste("PubMed trends for:", input$term))
      })
    })
  })
}
