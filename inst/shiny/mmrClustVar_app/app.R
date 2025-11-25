# inst/shiny/mmrClustVar_app/app.R

library(shiny)
library(mmrClustVar)

ui <- fluidPage(

  titlePanel("mmrClustVar - Clustering de variables"),

  sidebarLayout(

    sidebarPanel(

      h4("1. Charger un fichier CSV"),

      fileInput("file", "Sélectionner un fichier CSV"),
      checkboxInput("header", "En-tête présente", TRUE),
      selectInput("sep", "Séparateur",
                  choices = c("," = ",", ";" = ";", "\t" = "\t")),
      selectInput("encoding", "Encodage",
                  choices = c("UTF-8", "latin1", "CP1252"), selected = "UTF-8"),

      hr(),
      h4("2. Sélection des variables"),

      uiOutput("var_selector_active"),
      uiOutput("var_selector_descr"),

      hr(),
      h4("3. Paramètres du modèle"),

      selectInput(
        "method", "Méthode de clustering",
        choices = c(
          "k-means (numériques)"            = "kmeans",
          "k-modes (catégorielles)"         = "kmodes",
          "k-prototypes (mixtes)"           = "kprototypes"
        ),
        selected = "kmeans"
      ),

      sliderInput("K", "Nombre de clusters K", min = 2, max = 10, value = 3),
      checkboxInput("scale", "Standardiser les variables numériques", TRUE),
      numericInput("lambda", "Lambda (k-prototypes)", value = 1, min = .1, step = .1),

      hr(),
      h4("4. Actions"),

      actionButton("run_fit", "Lancer le clustering", class = "btn-primary"),
      br(), br(),
      actionButton("run_predict", "Rattacher les variables descriptives"),

      width = 3
    ),

    mainPanel(

      tabsetPanel(

        tabPanel("Résumé du modèle",
                 h4("Résumé"),
                 verbatimTextOutput("summary_out")
        ),

        tabPanel("Clusters",
                 h4("Variables actives et adhésion"),
                 dataTableOutput("clusters_table")
        ),

        tabPanel("Graphiques",
                 h4("Visualisation"),
                 selectInput(
                   "plot_type", "Type de graphique",
                   choices = c(
                     "Taille des clusters"               = "clusters",
                     "Adhésion des variables"            = "membership",
                     "Courbe d'inertie en fonction de K" = "inertia"
                   )
                 ),
                 conditionalPanel(
                   condition = "input.plot_type == 'inertia'",
                   sliderInput("Ks_range", "Plage de K", min = 2, max = 10, value = c(2,6))
                 ),
                 plotOutput("plot_out", height = "450px")
        ),

        tabPanel("Variables descriptives",
                 h4("Résultat du rattachement"),
                 dataTableOutput("descr_table"),
                 br(),
                 downloadButton("download_descr", "Exporter en CSV")
        )
      )
    )
  )
)

server <- function(input, output, session) {

  # --- Lecture CSV ----------------------------------------------------------------
  data_reactive <- reactive({
    req(input$file)
    read.csv(input$file$datapath,
             sep = input$sep,
             header = input$header,
             fileEncoding = input$encoding)
  })

  # --- Sélection variables ---------------------------------------------------------
  output$var_selector_active <- renderUI({
    req(data_reactive())
    selectInput("vars_active", "Variables actives",
                choices = names(data_reactive()),
                multiple = TRUE)
  })

  output$var_selector_descr <- renderUI({
    req(data_reactive())
    selectInput("vars_descr", "Variables descriptives",
                choices = names(data_reactive()),
                multiple = TRUE)
  })

  # objet du modèle
  model <- reactiveVal(NULL)

  # --- FIT ------------------------------------------------------------------------
  observeEvent(input$run_fit, {
    req(input$vars_active)

    X <- data_reactive()[, input$vars_active, drop = FALSE]

    obj <- mmrClustVar$new(
      method = input$method,
      K      = input$K,
      scale  = input$scale,
      lambda = input$lambda
    )

    obj$fit(X)
    model(obj)
  })

  # --- Résumé ---------------------------------------------------------------------
  output$summary_out <- renderPrint({
    req(model())
    model()$summary()
  })

  # --- Table clusters --------------------------------------------------------------
  output$clusters_table <- renderDataTable({
    req(model())
    model()$summary()
  })

  # --- Graphiques -----------------------------------------------------------------
  output$plot_out <- renderPlot({
    req(model())

    type <- input$plot_type

    if (type == "inertia") {
      Ks <- seq(input$Ks_range[1], input$Ks_range[2])
      model()$plot(type = "inertia", Ks = Ks)
    } else {
      model()$plot(type = type)
    }
  })

  # --- Predict (variables descriptives) -------------------------------------------
  descr_res <- reactiveVal(NULL)

  observeEvent(input$run_predict, {
    req(model(), input$vars_descr)

    X_descr <- data_reactive()[, input$vars_descr, drop = FALSE]

    # mmrClustVar$predict() attend un data.frame complet
    res <- model()$predict(X_descr)

    descr_res(res)
  })

  output$descr_table <- renderDataTable({
    req(descr_res())
    descr_res()
  })

  output$download_descr <- downloadHandler(
    filename = function() "variables_descriptives_clusters.csv",
    content = function(file) {
      write.csv(descr_res(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
