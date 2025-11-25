library(shiny)
library(mmrClustVar)

ui <- fluidPage(
    
    titlePanel("mmrClustVar — Clustering de variables"),
    
    sidebarLayout(
        sidebarPanel(
            
            h4("1. Données"),
            
            selectInput(
                "data_source", "Source des données",
                choices = c("Jeu interne", "Fichier utilisateur (CSV/XLSX)"),
                selected = "Jeu interne"
            ),
            
            conditionalPanel(
                "input.data_source == 'Jeu interne'",
                selectInput(
                    "builtin_dataset", "Jeu interne :",
                    choices = c("iris_num", "iris_mixed", "mtcars", "airquality"),
                    selected = "iris_num"
                )
            ),
            
            conditionalPanel(
                "input.data_source == 'Fichier utilisateur (CSV/XLSX)'",
                fileInput("file", "Sélectionner un fichier"),
                checkboxInput("header", "En-tête ? (CSV uniquement)", TRUE),
                selectInput("sep", "Séparateur (CSV uniquement)",
                            choices = c("," = ",", ";" = ";", "\t" = "\t")),
                selectInput("encoding", "Encodage (CSV uniquement)",
                            choices = c("UTF-8", "latin1", "CP1252"))
            ),
            
            hr(),
            
            h4("2. Sélection des variables"),
            uiOutput("var_select"),
            
            hr(),
            
            h4("3. Paramètres du modèle"),
            
            selectInput(
                "method", "Méthode",
                choices = list(
                    "Auto"          = "auto",
                    "k-means"       = "kmeans",
                    "k-modes"       = "kmodes",
                    "k-prototypes"  = "kprototypes",
                    "k-medoids"     = "kmedoids"
                ),
                selected = "auto"
            ),
            
            sliderInput("K", "Nombre de clusters K :", min = 2, max = 10, value = 3),
            checkboxInput("scale", "Standardiser les quantitatives", TRUE),
            
            numericInput("lambda", "λ (k-prototypes)", value = 1, min = 0.1, step = 0.1),
            
            hr(),
            
            actionButton("run_fit", "Lancer le clustering", class = "btn-primary"),
            actionButton("run_predict", "Rattacher les variables supplémentaires"),
            width = 3
        ),
        
        mainPanel(
            
            tabsetPanel(
                
                tabPanel("Résumé",
                         verbatimTextOutput("model_print"),
                         hr(),
                         verbatimTextOutput("model_summary")
                ),
                
                tabPanel("Clusters",
                         tableOutput("cluster_table"),
                         tableOutput("cluster_sizes")
                ),
                
                tabPanel("Graphiques",
                         selectInput("plot_type", "Type de graphique :",
                                     choices = c("Inertie" = "inertia",
                                                 "Clusters" = "clusters",
                                                 "Adhésion" = "membership",
                                                 "Profils"  = "profiles")),
                         plotOutput("plot")
                ),
                
                tabPanel("Variables supplémentaires",
                         tableOutput("predict_table")
                ),
                
                tabPanel("Diagnostics",
                         verbatimTextOutput("diag_infos")
                )
            )
        )
    )
)


server <- function(input, output, session) {
    
    # --- Chargement des données ---------------------------------------------
    
    dataset <- reactive({
        
        if (input$data_source == "Jeu interne") {
            
            env <- new.env()
            data(list = input$builtin_dataset,
                 package = "mmrClustVar",
                 envir = env)
            
            get(input$builtin_dataset, envir = env)
            
        } else {
            
            req(input$file)
            
            ext <- tools::file_ext(input$file$name)
            ext <- tolower(ext)
            
            if (ext %in% c("csv", "txt")) {
                read.csv(
                    input$file$datapath,
                    header      = input$header,
                    sep         = input$sep,
                    fileEncoding = input$encoding
                )
            } else if (ext %in% c("xlsx", "xls")) {
                if (!requireNamespace("readxl", quietly = TRUE)) {
                    stop("Le package 'readxl' est requis pour lire les fichiers Excel.")
                }
                as.data.frame(readxl::read_excel(input$file$datapath))
            } else {
                stop("Extension de fichier non supportée. Utilisez un CSV ou un XLSX.")
            }
        }
    })
    
    
    # --- Sélection variables -------------------------------------------------
    
    output$var_select <- renderUI({
        df <- dataset()
        choices <- names(df)
        
        tagList(
            selectInput("active_vars", "Variables actives",
                        choices = choices, multiple = TRUE),
            selectInput("suppl_vars", "Variables supplémentaires",
                        choices = choices, multiple = TRUE)
        )
    })
    
    
    # --- Clustering ----------------------------------------------------------
    
    model_obj <- reactiveVal(NULL)
    
    observeEvent(input$run_fit, {
        
        df <- dataset()
        
        req(input$active_vars)
        
        X <- df[, input$active_vars, drop = FALSE]
        
        obj <- mmrClustVar$new(
            method = input$method,
            K      = input$K,
            scale  = input$scale,
            lambda = input$lambda
        )
        
        obj$fit(X)
        model_obj(obj)
        
        showNotification("Clustering terminé.", type = "message")
    })
    
    
    # --- Rattachement des variables supplémentaires -------------------------
    
    observeEvent(input$run_predict, {
        
        obj <- model_obj()
        req(obj)
        
        df <- dataset()
        req(input$suppl_vars)
        
        X_suppl <- df[, input$suppl_vars, drop = FALSE]
        
        pred <- obj$predict(X_suppl)
        
        output$predict_table <- renderTable(pred)
        
        showNotification("Variables supplémentaires rattachées.", type = "message")
    })
    
    
    # --- Sorties texte : print et summary -----------------------------------
    
    output$model_print <- renderPrint({
        obj <- model_obj()
        req(obj)
        obj$print()
    })
    
    output$model_summary <- renderPrint({
        obj <- model_obj()
        req(obj)
        obj$summary()
    })
    
    
    # --- Table clusters ------------------------------------------------------
    
    output$cluster_table <- renderTable({
        obj <- model_obj()
        req(obj)
        
        df   <- dataset()
        vars <- input$active_vars
        clusters <- obj$get_clusters()
        
        if (length(vars) != length(clusters)) {
            return(data.frame())
        }
        
        data.frame(
            variable = vars,
            cluster  = clusters,
            stringsAsFactors = FALSE
        )
    })
    
    
    # --- Graphiques ----------------------------------------------------------
    
    output$plot <- renderPlot({
        obj <- model_obj()
        req(obj)
        
        if (input$plot_type == "inertia")
            obj$plot(type = "inertia")
        else if (input$plot_type == "clusters")
            obj$plot(type = "clusters")
        else if (input$plot_type == "membership")
            obj$plot(type = "membership")
        else if (input$plot_type == "profiles")
            obj$plot(type = "profiles")
    })
    
    
    # --- Diagnostics ---------------------------------------------------------
    
    output$diag_infos <- renderText({
        obj <- model_obj()
        req(obj)
        
        paste(
            "Méthode :", obj$get_method(), "\n",
            "Convergence :", obj$get_convergence(), "\n",
            "Inertie :", obj$get_inertia()
        )
    })
}

shinyApp(ui, server)
