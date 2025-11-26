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
            helpText("La standardisation est ignorée si toutes les variables actives sont qualitatives."),
            
            # λ uniquement pour k-prototypes
            conditionalPanel(
                "input.method == 'kprototypes'",
                numericInput("lambda", "λ (k-prototypes)", value = 1, min = 0.1, step = 0.1)
            ),
            
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
                         tableOutput("cluster_sizes"),
                         downloadButton("download_clusters", "Exporter les clusters (CSV)")
                ),
                
                tabPanel("Graphiques",
                         selectInput("plot_type", "Type de graphique :",
                                     choices = c("Inertie"   = "inertia",
                                                 "Clusters"  = "clusters",
                                                 "Adhésion"  = "membership",
                                                 "Profils"   = "profiles")),
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
                    header       = input$header,
                    sep          = input$sep,
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
        p_actives <- ncol(X)
        
        # 1) Vérifier K <= nb variables actives
        if (input$K > p_actives) {
            showNotification(
                "K ne peut pas dépasser le nombre de variables actives sélectionnées.",
                type = "error"
            )
            return()
        }
        
        # 2) Déterminer types des variables actives
        is_num <- vapply(X, is.numeric, logical(1L))
        has_num <- any(is_num)
        has_cat <- any(!is_num)
        
        # 3) Gérer scale : ignoré si aucune variable numérique
        scale_arg <- input$scale
        if (!has_num && isTRUE(scale_arg)) {
            scale_arg <- FALSE
            showNotification(
                "Standardisation ignorée : toutes les variables actives sont qualitatives.",
                type = "message"
            )
        }
        
        # 4) Création de l'objet façade
        lambda_arg <- if (!is.null(input$lambda)) input$lambda else 1
        
        obj <- mmrClustVar$new(
            method = input$method,
            K      = input$K,
            scale  = scale_arg,
            lambda = lambda_arg
        )
        
        # 5) Apprentissage avec gestion propre des erreurs
        tryCatch({
            obj$fit(X)
            model_obj(obj)
            showNotification("Clustering terminé.", type = "message")
        }, error = function(e) {
            showNotification(
                paste("Erreur pendant le clustering :", conditionMessage(e)),
                type = "error"
            )
        })
    })
    
    
    # --- Rattachement des variables supplémentaires -------------------------
    
    observeEvent(input$run_predict, {
        
        obj <- model_obj()
        req(obj)
        
        df <- dataset()
        req(input$suppl_vars)
        
        X_suppl <- df[, input$suppl_vars, drop = FALSE]
        
        pred <- tryCatch({
            obj$predict(X_suppl)
        }, error = function(e) {
            showNotification(
                paste("Erreur pendant le rattachement :", conditionMessage(e)),
                type = "error"
            )
            return(NULL)
        })
        
        if (!is.null(pred)) {
            output$predict_table <- renderTable(pred)
            showNotification("Variables supplémentaires rattachées.", type = "message")
        }
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
        
        # Cas où rien n'est sélectionné / appris
        if (is.null(vars) || length(vars) == 0L) {
            return(data.frame())
        }
        if (is.null(clusters) || length(clusters) == 0L) {
            return(data.frame())
        }
        
        # Cas incohérent (sécurité, doesn't occure normalement)
        if (length(vars) != length(clusters)) {
            return(data.frame(
                message = "Incohérence entre le nombre de variables actives et le nombre d'affectations de clusters."
            ))
        }
        
        data.frame(
            variable = vars,
            cluster  = clusters,
            stringsAsFactors = FALSE
        )
    })

    # --- Tailles de clusters -------------------------------------------------
    
    output$cluster_sizes <- renderTable({
        obj <- model_obj()
        req(obj)
        
        clusters <- obj$get_clusters()
        if (is.null(clusters) || length(clusters) == 0L) {
            return(data.frame())
        }
        
        as.data.frame(table(cluster = clusters))
    })

    # --- Export des clusters -------------------------------------------------
    
    output$download_clusters <- downloadHandler(
        filename = function() {
            paste0("mmrClustVar_clusters_", Sys.Date(), ".csv")
        },
        content = function(file) {
            obj <- model_obj()
            req(obj)
            
            df   <- dataset()
            vars <- input$active_vars
            clusters <- obj$get_clusters()
            
            if (is.null(vars) || length(vars) == 0L ||
                is.null(clusters) || length(clusters) == 0L ||
                length(vars) != length(clusters)) {
                utils::write.csv(
                    data.frame(message = "Aucun résultat exploitable (clustering non réalisé ou incohérent)."),
                    file,
                    row.names = FALSE
                )
            } else {
                out <- data.frame(
                    variable = vars,
                    cluster  = clusters,
                    stringsAsFactors = FALSE
                )
                utils::write.csv(out, file, row.names = FALSE)
            }
        }
    )

    
    
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
            "Méthode    :", obj$get_method(),      "\n",
            "Convergence:", obj$get_convergence(), "\n",
            "Inertie    :", obj$get_inertia()
        )
    })
}

shinyApp(ui, server)
