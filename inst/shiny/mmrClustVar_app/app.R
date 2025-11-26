library(shiny)
library(mmrClustVar)

ui <- fluidPage(
    
    titlePanel("mmrClustVar — Clustering de variables"),
    
    sidebarLayout(
        sidebarPanel(
            
            h4("1. Données"),
            
            selectInput(
                "data_source", "Source des données",
                choices  = c("Jeu interne", "Fichier utilisateur (CSV/XLSX)"),
                selected = "Jeu interne"
            ),
            
            conditionalPanel(
                "input.data_source == 'Jeu interne'",
                selectInput(
                    "builtin_dataset", "Jeu interne :",
                    choices  = c(
                        "iris_num",
                        "iris_mixed",
                        "mtcars_num",
                        "airquality_num",
                        "arthritis_cat",
                        "titanic_cat",
                        "housevotes_cat",
                        "credit_mix",
                        "adult_small_mix"
                    ),
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
            
            numericInput(
                "K", "Nombre de clusters K", value = 3, min = 2, step = 1
            ),
            
            checkboxInput(
                "scale", "Standardiser les variables numériques ?", TRUE
            ),
            
            conditionalPanel(
                "input.method == 'kprototypes' || input.method == 'kmedoids'",
                numericInput("lambda", "λ (pondération parties catégorielles)",
                             value = 1, min = 0.1, step = 0.1)
            ),
            
            hr(),
            
            actionButton("run_fit",     "Lancer le clustering", class = "btn-primary"),
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
                                     choices = c("Inertie (méthode du coude)" = "inertia",
                                                 "Clusters"                    = "clusters",
                                                 "Adhésion"                    = "membership",
                                                 "Profils / heatmaps"          = "profiles")),
                         plotOutput("plot")
                ),
                
                tabPanel("Variables supplémentaires",
                         tableOutput("predict_table")
                ),
                
                tabPanel("Diagnostics & Export",
                         verbatimTextOutput("diag_infos"),
                         br(),
                         downloadButton("download_report",  "Exporter le rapport texte"),
                         br(), br(),
                         downloadButton("download_bundle",  "Exporter tous les résultats (.zip)")
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
            data(list    = input$builtin_dataset,
                 package = "mmrClustVar",
                 envir   = env)
            
            get(input$builtin_dataset, envir = env)
            
        } else {
            
            req(input$file)
            ext <- tools::file_ext(input$file$name)
            
            if (tolower(ext) == "csv") {
                utils::read.csv(
                    input$file$datapath,
                    header   = isTRUE(input$header),
                    sep      = input$sep,
                    fileEncoding = input$encoding,
                    stringsAsFactors = FALSE,
                    check.names = TRUE
                )
            } else if (tolower(ext) %in% c("xls", "xlsx")) {
                if (!requireNamespace("readxl", quietly = TRUE)) {
                    stop("Le package 'readxl' est requis pour lire les fichiers Excel.")
                }
                as.data.frame(readxl::read_excel(input$file$datapath))
            } else {
                stop("Extension de fichier non supportée. Utilisez un CSV ou un XLS/XLSX.")
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
        
        # 2) Types des variables actives
        is_num  <- vapply(X, is.numeric, logical(1L))
        has_num <- any(is_num)
        
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
        
        # 5) Ajustement
        res <- tryCatch({
            obj$fit(X)
            obj
        }, error = function(e) {
            showNotification(
                paste("Erreur pendant le clustering :", conditionMessage(e)),
                type = "error"
            )
            NULL
        })
        
        if (!is.null(res)) {
            model_obj(res)
            showNotification("Clustering terminé.", type = "message")
        }
    })
    
    
    # --- Affichage résumé ----------------------------------------------------
    
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
    
    
    # --- Tableau des clusters -----------------------------------------------
    
    output$cluster_table <- renderTable({
        obj <- model_obj()
        req(obj)
        
        df       <- dataset()
        vars     <- input$active_vars
        clusters <- obj$get_clusters()
        
        if (is.null(vars) || length(vars) == 0L) {
            return(data.frame())
        }
        if (is.null(clusters) || length(clusters) == 0L) {
            return(data.frame())
        }
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
    
    output$cluster_sizes <- renderTable({
        obj <- model_obj()
        req(obj)
        
        clusters <- obj$get_clusters()
        if (is.null(clusters) || length(clusters) == 0L) {
            return(data.frame())
        }
        
        as.data.frame(table(cluster = clusters))
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
    
    
    # --- Graphiques ----------------------------------------------------------
    
    output$plot <- renderPlot({
        obj <- model_obj()
        req(obj)
        
        if (input$plot_type == "inertia") {
            # Courbe du coude : on recalcule le chemin d'inertie pour K = 2..p
            df <- dataset()
            req(input$active_vars)
            X <- df[, input$active_vars, drop = FALSE]
            p <- ncol(X)
            
            if (p >= 2L) {
                obj$compute_inertia_path(K_seq = 2:p, X = X)
                obj$plot(type = "inertia")
            } else {
                plot.new()
                title("Au moins 2 variables actives sont nécessaires pour tracer la courbe d'inertie.")
            }
            
        } else if (input$plot_type == "clusters") {
            obj$plot(type = "clusters")
            
        } else if (input$plot_type == "membership") {
            obj$plot(type = "membership")
            
        } else if (input$plot_type == "profiles") {
            obj$plot(type = "profiles")
        }
    })
    
    
    # --- Diagnostics ---------------------------------------------------------
    
    output$diag_infos <- renderText({
        obj <- model_obj()
        req(obj)
        
        conv  <- tryCatch(obj$get_convergence(), error = function(e) NA)
        inert <- tryCatch(obj$get_inertia(),      error = function(e) NA)
        
        paste(
            "Méthode    :", obj$get_method(), "\n",
            "Convergence:", conv,             "\n",
            "Inertie    :", inert
        )
    })
    
    
    # --- Export clusters (CSV) -----------------------------------------------
    
    output$download_clusters <- downloadHandler(
        filename = function() {
            paste0("mmrClustVar_clusters_", Sys.Date(), ".csv")
        },
        content = function(file) {
            obj <- model_obj()
            validate(need(!is.null(obj), "Aucun modèle appris : lancez d'abord le clustering."))
            
            vars     <- input$active_vars
            clusters <- obj$get_clusters()
            
            if (is.null(vars) || is.null(clusters) || length(vars) != length(clusters)) {
                stop("Clusters indisponibles ou incohérents : vérifiez que le modèle a bien été appris.")
            }
            
            df <- data.frame(
                variable = vars,
                cluster  = clusters,
                stringsAsFactors = FALSE
            )
            
            utils::write.csv(df, file, row.names = FALSE, fileEncoding = "UTF-8")
        }
    )
    
    
    # --- Export rapport complet (texte) --------------------------------------
    
    output$download_report <- downloadHandler(
        filename = function() {
            paste0("mmrClustVar_report_", Sys.Date(), ".txt")
        },
        content = function(file) {
            obj <- model_obj()
            validate(need(!is.null(obj), "Aucun modèle appris : lancez d'abord le clustering."))
            
            df       <- dataset()
            vars     <- input$active_vars
            clusters <- obj$get_clusters()
            
            # 1) Résumé (print + summary)
            summary_txt <- c(
                "===== Résumé du modèle =====",
                "",
                capture.output(obj$print()),
                "",
                capture.output(obj$summary()),
                ""
            )
            
            # 2) Clusters
            clusters_txt <- "===== Affectations des variables ====="
            clust_df <- NULL
            if (!is.null(vars) && !is.null(clusters) && length(vars) == length(clusters)) {
                clust_df <- data.frame(
                    variable = vars,
                    cluster  = clusters,
                    stringsAsFactors = FALSE
                )
            }
            
            # 3) Variables supplémentaires (si présentes)
            suppl_df <- NULL
            if (!is.null(input$suppl_vars) && length(input$suppl_vars) > 0L) {
                X_suppl <- df[, input$suppl_vars, drop = FALSE]
                suppl_df <- tryCatch(
                    obj$predict(X_suppl),
                    error = function(e) NULL
                )
            }
            
            # 4) Diagnostics simples
            conv  <- tryCatch(obj$get_convergence(), error = function(e) NA)
            inert <- tryCatch(obj$get_inertia(),      error = function(e) NA)
            
            diag_txt <- c(
                "===== Diagnostics =====",
                paste("Méthode    :", obj$get_method()),
                paste("Convergence:", conv),
                paste("Inertie    :", inert),
                ""
            )
            
            # 5) Écriture dans le fichier texte
            con <- file(file, open = "wt", encoding = "UTF-8")
            on.exit(close(con), add = TRUE)
            
            writeLines("===== Rapport mmrClustVar =====", con)
            writeLines("", con)
            
            writeLines(summary_txt, con)
            
            writeLines(clusters_txt, con)
            if (!is.null(clust_df)) {
                utils::write.table(clust_df, con, sep = "\t",
                                   row.names = FALSE, quote = FALSE)
                writeLines("", con)
            } else {
                writeLines("(Aucune affectation disponible)", con)
                writeLines("", con)
            }
            
            if (!is.null(suppl_df)) {
                writeLines("===== Variables supplémentaires =====", con)
                utils::write.table(suppl_df, con, sep = "\t",
                                   row.names = FALSE, quote = FALSE)
                writeLines("", con)
            }
            
            writeLines(diag_txt, con)
            
            writeLines("===== Graphiques =====", con)
            writeLines(
                "Les graphiques (inertie, répartition des clusters, adhésion, profils) doivent être sauvegardés séparément depuis l'onglet 'Graphiques'.",
                con
            )
        }
    )
    
    
    # --- Export complet ZIP --------------------------------------------------
    
    output$download_bundle <- downloadHandler(
        filename = function() {
            paste0("mmrClustVar_results_", Sys.Date(), ".zip")
        },
        content = function(file) {
            obj <- model_obj()
            validate(need(!is.null(obj), "Aucun modèle appris : lancez d'abord le clustering."))
            
            df   <- dataset()
            vars <- input$active_vars
            req(vars)
            X    <- df[, vars, drop = FALSE]
            p    <- ncol(X)
            
            # dossier temporaire
            bundle_dir <- file.path(tempdir(), paste0("mmrClustVar_bundle_", Sys.getpid()))
            if (!dir.exists(bundle_dir)) dir.create(bundle_dir, recursive = TRUE)
            
            # 1) Résumé texte
            summary_path <- file.path(bundle_dir, "summary.txt")
            summary_conn <- file(summary_path, open = "wt", encoding = "UTF-8")
            writeLines(capture.output(obj$print()),   summary_conn)
            writeLines("", summary_conn)
            writeLines(capture.output(obj$summary()), summary_conn)
            close(summary_conn)
            
            # 2) Clusters
            clusters <- obj$get_clusters()
            if (!is.null(clusters) && length(clusters) == length(vars)) {
                clusters_df <- data.frame(
                    variable = vars,
                    cluster  = clusters,
                    stringsAsFactors = FALSE
                )
                utils::write.csv(clusters_df,
                                 file = file.path(bundle_dir, "clusters.csv"),
                                 row.names = FALSE, fileEncoding = "UTF-8")
            }
            
            # 3) Prototypes / centres / medoids
            centers <- obj$get_centers()
            saveRDS(centers, file = file.path(bundle_dir, "centers.rds"))
            
            # 4) Courbe d'inertie (2..p)
            if (p >= 2L) {
                inertia_df <- obj$compute_inertia_path(K_seq = 2:p, X = X)
                utils::write.csv(inertia_df,
                                 file = file.path(bundle_dir, "inertia_path.csv"),
                                 row.names = FALSE, fileEncoding = "UTF-8")
            }
            
            # 5) Variables supplémentaires (si présentes)
            if (!is.null(input$suppl_vars) && length(input$suppl_vars) > 0L) {
                X_suppl <- df[, input$suppl_vars, drop = FALSE]
                suppl_df <- tryCatch(obj$predict(X_suppl), error = function(e) NULL)
                if (!is.null(suppl_df)) {
                    utils::write.csv(suppl_df,
                                     file = file.path(bundle_dir, "suppl_predict.csv"),
                                     row.names = FALSE, fileEncoding = "UTF-8")
                }
            }
            
            # 6) README interne
            readme_path <- file.path(bundle_dir, "README.txt")
            readme_conn <- file(readme_path, open = "wt", encoding = "UTF-8")
            writeLines(c(
                "mmrClustVar — Export complet des résultats",
                "",
                "- summary.txt          : résumé textuel du modèle (print + summary)",
                "- clusters.csv         : affectation des variables actives aux clusters",
                "- centers.rds          : prototypes/centres/medoids des clusters (objet R)",
                "- inertia_path.csv     : inertie intra-cluster pour K = 2..p (méthode du coude)",
                "- suppl_predict.csv    : rattachement des variables supplémentaires (si présent)",
                "",
                "Vous pouvez charger centers.rds dans R via :",
                "  centers <- readRDS('centers.rds')"
            ), readme_conn)
            close(readme_conn)
            
            # 7) Création du ZIP
            old_wd <- getwd()
            setwd(bundle_dir)
            on.exit(setwd(old_wd), add = TRUE)
            
            files_to_zip <- list.files(bundle_dir, all.files = FALSE)
            utils::zip(zipfile = file, files = files_to_zip)
        }
    )
}

shinyApp(ui, server)