library(shiny)
library(shinydashboard)
library(DT)



ui <- dashboardPage(
  dashboardHeader(title = "Variable clustering methods"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data", tabName = "data", icon = icon("database")),
      menuItem("K-means", tabName = "kmeans", icon = icon("arrow-down-1-9")),
      menuItem("K-modes", tabName = "kmodes", icon = icon("arrow-down-a-z")),
      menuItem("K-prototypes", tabName = "kprototypes", icon = icon("shapes"))
    )
  ),
  dashboardBody(
    tabItems(

      # Data tab
      tabItem(
        tabName = "data",
        h2("Data loading"),
        fluidRow(
          # Column 1: file input and preview
          column(
            width = 6,
            box(
              title = "Load data", status = "primary", solidHeader = TRUE,
              width = 12,
              fileInput("file_kmeans", "Choose CSV File",
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv"))
            ),
            box(
              title = "Preview", status = "primary", solidHeader = TRUE,
              width = 12,
              DTOutput("table_csv")
            )
          ),
          # Column 2: active and descriptive variables selection
          column(
            width = 6,
            box(
              title = "Active variables",
              status = "info", solidHeader = TRUE, width = 12,
              selectInput(
                "active_vars",
                "Choose active variables:",
                choices = NULL,
                multiple = TRUE
              )
            ),
            box(
              title = "Descriptive variables",
              status = "info", solidHeader = TRUE, width = 12,
              selectInput(
                "desc_vars",
                "Choose descriptive variables:",
                choices = NULL,
                multiple = TRUE
              )
            )
          )
        )
      ),

      tabItem(
        tabName = "kmeans",
        h2("K-means Clustering"),
        actionButton("print_active", "An action button")
      ),
      tabItem(tabName = "kmodes", h2("K-modes Clustering")),
      tabItem(tabName = "kprototypes", h2("K-prototypes Clustering"))
    )
  )
)

server <- function(input, output, session) {

  # read data reactively
  data_csv <- reactive({
    req(input$file_kmeans)
    read.csv(input$file_kmeans$datapath)
  })

  # update select menus when data is available
  observeEvent(data_csv(), {
    vars <- names(data_csv())
    updateSelectInput(session, "active_vars", choices = vars)
    updateSelectInput(session, "desc_vars",  choices = vars)
  })

  # store selections
  active.variables      <- reactive(input$active_vars)
  descriptive.variables <- reactive(input$desc_vars)

  # update data table
  output$table_csv <- renderDT({
    data <- data_csv()

    dt <- datatable(
      data,
      options = list(
        pageLength = 10,
        dom = "rtip",
        scrollX = TRUE
      )
    )

    # hightlight active variables
    for (col in active.variables()) {
      if (col %in% names(data)) {
        dt <- dt %>%
          formatStyle(
            columns = col,
            backgroundColor = "#ffe8a1",
            fontWeight = "bold"
          )
      }
    }
    # hightlight descriptive variables
    for (col in descriptive.variables()) {
      if (col %in% names(data)) {
        dt <- dt %>%
          formatStyle(
            columns = col,
            backgroundColor = "#98eda4",
            fontWeight = "bold"
          )
      }
    }

    dt
  })

  observeEvent(input$print_active, {
    cat("Active:", active.variables(), "\n")
  })
}

shinyApp(ui = ui, server = server)
