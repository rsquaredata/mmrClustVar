library(shiny)
library(mmrClustVar)

ui <- fluidPage(
    
    titlePanel("mmrClustVar — Variable clustering"),
    
    sidebarLayout(
        sidebarPanel(
            
            h4("1. Data"),
            
            selectInput(
                "data_source", "Data source",
                choices  = c("Built-in dataset", "User file (CSV/XLSX)"),
                selected = "Built-in dataset"
            ),
            
            conditionalPanel(
                "input.data_source == 'Built-in dataset'",
                selectInput(
                    "builtin_dataset", "Built-in dataset:",
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
                "input.data_source == 'User file (CSV/XLSX)'",
                fileInput("file", "Select a file"),
                checkboxInput("header", "Header? (CSV only)", TRUE),
                selectInput("sep", "Separator (CSV only)",
                            choices = c("," = ",", ";" = ";", "\t" = "\t")),
                selectInput("encoding", "Encoding (CSV only)",
                            choices = c("UTF-8", "latin1", "CP1252"))
            ),
            
            hr(),
            
            h4("2. Variable selection"),
            uiOutput("var_select"),
            
            hr(),
            
            h4("3. Model parameters"),
            
            selectInput(
                "method", "Method",
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
                "K", "Number of clusters K", value = 3, min = 2, step = 1
            ),
            
            checkboxInput(
                "scale", "Standardize numeric variables?", TRUE
            ),
            
            conditionalPanel(
                "input.method == 'kprototypes' || input.method == 'kmedoids'",
                numericInput(
                    "lambda",
                    "λ (weight for categorical part)",
                    value = 1, min = 0.1, step = 0.1
                )
            ),
            
            hr(),
            
            actionButton("run_fit",     "Run clustering", class = "btn-primary"),
            actionButton("run_predict", "Attach supplementary variables"),
            width = 3
        ),
        
        mainPanel(
            tabsetPanel(
                
                tabPanel(
                    "Summary",
                    verbatimTextOutput("model_print"),
                    hr(),
                    verbatimTextOutput("model_summary")
                ),
                
                tabPanel(
                    "Clusters",
                    tableOutput("cluster_table"),
                    tableOutput("cluster_sizes"),
                    downloadButton("download_clusters", "Export clusters (CSV)")
                ),
                
                tabPanel(
                    "Plots",
                    selectInput(
                        "plot_type", "Plot type:",
                        choices = c(
                            "Inertia (elbow method)" = "inertia",
                            "Clusters"               = "clusters",
                            "Membership"             = "membership",
                            "Profiles / heatmaps"    = "profiles"
                        )
                    ),
                    plotOutput("plot")
                ),
                
                tabPanel(
                    "Supplementary variables",
                    tableOutput("predict_table")
                ),
                
                tabPanel(
                    "Diagnostics & export",
                    verbatimTextOutput("diag_infos"),
                    br(),
                    downloadButton("download_report",  "Export text report"),
                    br(), br(),
                    downloadButton("download_bundle",  "Export all results (.zip)")
                )
            )
        )
    )
)

server <- function(input, output, session) {
    
    # --- Data loading --------------------------------------------------------
    
    dataset <- reactive({
        
        if (input$data_source == "Built-in dataset") {
            
            env <- new.env()
            data(
                list    = input$builtin_dataset,
                package = "mmrClustVar",
                envir   = env
            )
            
            get(input$builtin_dataset, envir = env)
            
        } else {
            
            req(input$file)
            ext <- tools::file_ext(input$file$name)
            
            if (tolower(ext) == "csv") {
                utils::read.csv(
                    input$file$datapath,
                    header       = isTRUE(input$header),
                    sep          = input$sep,
                    fileEncoding = input$encoding,
                    stringsAsFactors = FALSE,
                    check.names  = TRUE
                )
            } else if (tolower(ext) %in% c("xls", "xlsx")) {
                if (!requireNamespace("readxl", quietly = TRUE)) {
                    stop("The 'readxl' package is required to read Excel files.")
                }
                as.data.frame(readxl::read_excel(input$file$datapath))
            } else {
                stop("Unsupported file extension. Please use CSV or XLS/XLSX.")
            }
        }
    })
    
    
    # --- Variable selection --------------------------------------------------
    
    output$var_select <- renderUI({
        df <- dataset()
        choices <- names(df)
        
        tagList(
            selectInput(
                "active_vars", "Active variables",
                choices = choices, multiple = TRUE
            ),
            selectInput(
                "suppl_vars", "Supplementary variables",
                choices = choices, multiple = TRUE
            )
        )
    })
    
    
    # --- Clustering ----------------------------------------------------------
    
    model_obj <- reactiveVal(NULL)
    
    observeEvent(input$run_fit, {
        
        df <- dataset()
        req(input$active_vars)
        
        X <- df[, input$active_vars, drop = FALSE]
        p_actives <- ncol(X)
        
        # 1) Check K <= nb of active variables
        if (input$K > p_actives) {
            showNotification(
                "K cannot exceed the number of selected active variables.",
                type = "error"
            )
            return()
        }
        
        # 2) Types of active variables
        is_num  <- vapply(X, is.numeric, logical(1L))
        has_num <- any(is_num)
        
        # 3) Handle scale: ignored if there are no numeric variables
        scale_arg <- input$scale
        if (!has_num && isTRUE(scale_arg)) {
            scale_arg <- FALSE
            showNotification(
                "Standardization ignored: all active variables are categorical.",
                type = "message"
            )
        }
        
        # 4) Create facade object
        lambda_arg <- if (!is.null(input$lambda)) input$lambda else 1
        
        obj <- mmrClustVar$new(
            method = input$method,
            K      = input$K,
            scale  = scale_arg,
            lambda = lambda_arg
        )
        
        # 5) Fit model
        res <- tryCatch({
            obj$fit(X)
            obj
        }, error = function(e) {
            showNotification(
                paste("Error during clustering:", conditionMessage(e)),
                type = "error"
            )
            NULL
        })
        
        if (!is.null(res)) {
            model_obj(res)
            showNotification("Clustering completed.", type = "message")
        }
    })
    
    
    # --- Summary display -----------------------------------------------------
    
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
    
    
    # --- Cluster tables ------------------------------------------------------
    
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
                message = "Mismatch between the number of active variables and the number of cluster assignments."
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
    
    
    # --- Attach supplementary variables --------------------------------------
    
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
                paste("Error during prediction:", conditionMessage(e)),
                type = "error"
            )
            return(NULL)
        })
        
        if (!is.null(pred)) {
            output$predict_table <- renderTable(pred)
            showNotification("Supplementary variables attached.", type = "message")
        }
    })
    
    
    # --- Plots ---------------------------------------------------------------
    
    output$plot <- renderPlot({
        obj <- model_obj()
        req(obj)
        
        if (input$plot_type == "inertia") {
            # Elbow curve: recompute inertia path for K = 2..p
            df <- dataset()
            req(input$active_vars)
            X <- df[, input$active_vars, drop = FALSE]
            p <- ncol(X)
            
            if (p >= 2L) {
                obj$compute_inertia_path(K_seq = 2:p, X = X)
                obj$plot(type = "inertia")
            } else {
                plot.new()
                title("At least 2 active variables are required to draw the inertia curve.")
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
            "Method    :", obj$get_method(), "\n",
            "Converged :", conv,            "\n",
            "Inertia   :", inert
        )
    })
    
    
    # --- Export clusters (CSV) -----------------------------------------------
    
    output$download_clusters <- downloadHandler(
        filename = function() {
            paste0("mmrClustVar_clusters_", Sys.Date(), ".csv")
        },
        content = function(file) {
            obj <- model_obj()
            validate(need(!is.null(obj), "No fitted model: please run clustering first."))
            
            vars     <- input$active_vars
            clusters <- obj$get_clusters()
            
            if (is.null(vars) || is.null(clusters) || length(vars) != length(clusters)) {
                stop("Clusters unavailable or inconsistent: please check that the model has been fitted.")
            }
            
            df <- data.frame(
                variable = vars,
                cluster  = clusters,
                stringsAsFactors = FALSE
            )
            
            utils::write.csv(df, file, row.names = FALSE, fileEncoding = "UTF-8")
        }
    )
    
    
    # --- Export full text report ---------------------------------------------
    
    output$download_report <- downloadHandler(
        filename = function() {
            paste0("mmrClustVar_report_", Sys.Date(), ".txt")
        },
        content = function(file) {
            obj <- model_obj()
            validate(need(!is.null(obj), "No fitted model: please run clustering first."))
            
            df       <- dataset()
            vars     <- input$active_vars
            clusters <- obj$get_clusters()
            
            # 1) Summary (print + summary)
            summary_txt <- c(
                "===== Model summary =====",
                "",
                capture.output(obj$print()),
                "",
                capture.output(obj$summary()),
                ""
            )
            
            # 2) Clusters
            clusters_txt <- "===== Variable-to-cluster assignments ====="
            clust_df <- NULL
            if (!is.null(vars) && !is.null(clusters) && length(vars) == length(clusters)) {
                clust_df <- data.frame(
                    variable = vars,
                    cluster  = clusters,
                    stringsAsFactors = FALSE
                )
            }
            
            # 3) Supplementary variables (if any)
            suppl_df <- NULL
            if (!is.null(input$suppl_vars) && length(input$suppl_vars) > 0L) {
                X_suppl <- df[, input$suppl_vars, drop = FALSE]
                suppl_df <- tryCatch(
                    obj$predict(X_suppl),
                    error = function(e) NULL
                )
            }
            
            # 4) Simple diagnostics
            conv  <- tryCatch(obj$get_convergence(), error = function(e) NA)
            inert <- tryCatch(obj$get_inertia(),      error = function(e) NA)
            
            diag_txt <- c(
                "===== Diagnostics =====",
                paste("Method   :", obj$get_method()),
                paste("Converged:", conv),
                paste("Inertia  :", inert),
                ""
            )
            
            # 5) Write to text file
            con <- file(file, open = "wt", encoding = "UTF-8")
            on.exit(close(con), add = TRUE)
            
            writeLines("===== mmrClustVar report =====", con)
            writeLines("", con)
            
            writeLines(summary_txt, con)
            
            writeLines(clusters_txt, con)
            if (!is.null(clust_df)) {
                utils::write.table(
                    clust_df, con, sep = "\t",
                    row.names = FALSE, quote = FALSE
                )
                writeLines("", con)
            } else {
                writeLines("(No cluster assignment available)", con)
                writeLines("", con)
            }
            
            if (!is.null(suppl_df)) {
                writeLines("===== Supplementary variables =====", con)
                utils::write.table(
                    suppl_df, con, sep = "\t",
                    row.names = FALSE, quote = FALSE
                )
                writeLines("", con)
            }
            
            writeLines(diag_txt, con)
            
            writeLines("===== Plots =====", con)
            writeLines(
                "Plots (inertia, cluster distribution, membership, profiles) must be saved separately from the 'Plots' tab.",
                con
            )
        }
    )
    
    
    # --- Full ZIP export -----------------------------------------------------
    
    output$download_bundle <- downloadHandler(
        filename = function() {
            paste0("mmrClustVar_results_", Sys.Date(), ".zip")
        },
        content = function(file) {
            obj <- model_obj()
            validate(need(!is.null(obj), "No fitted model: please run clustering first."))
            
            df   <- dataset()
            vars <- input$active_vars
            req(vars)
            X    <- df[, vars, drop = FALSE]
            p    <- ncol(X)
            
            # Temporary directory
            bundle_dir <- file.path(tempdir(), paste0("mmrClustVar_bundle_", Sys.getpid()))
            if (!dir.exists(bundle_dir)) dir.create(bundle_dir, recursive = TRUE)
            
            # 1) Text summary
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
                utils::write.csv(
                    clusters_df,
                    file = file.path(bundle_dir, "clusters.csv"),
                    row.names = FALSE, fileEncoding = "UTF-8"
                )
            }
            
            # 3) Prototypes / centers / medoids
            centers <- obj$get_centers()
            saveRDS(centers, file = file.path(bundle_dir, "centers.rds"))
            
            # 4) Inertia path (2..p)
            if (p >= 2L) {
                inertia_df <- obj$compute_inertia_path(K_seq = 2:p, X = X)
                utils::write.csv(
                    inertia_df,
                    file = file.path(bundle_dir, "inertia_path.csv"),
                    row.names = FALSE, fileEncoding = "UTF-8"
                )
            }
            
            # 5) Supplementary variables (if any)
            if (!is.null(input$suppl_vars) && length(input$suppl_vars) > 0L) {
                X_suppl <- df[, input$suppl_vars, drop = FALSE]
                suppl_df <- tryCatch(obj$predict(X_suppl), error = function(e) NULL)
                if (!is.null(suppl_df)) {
                    utils::write.csv(
                        suppl_df,
                        file = file.path(bundle_dir, "suppl_predict.csv"),
                        row.names = FALSE, fileEncoding = "UTF-8"
                    )
                }
            }
            
            # 6) Internal README
            readme_path <- file.path(bundle_dir, "README.txt")
            readme_conn <- file(readme_path, open = "wt", encoding = "UTF-8")
            writeLines(c(
                "mmrClustVar — Full export of results",
                "",
                "- summary.txt      : text summary of the model (print + summary)",
                "- clusters.csv     : active variables and their cluster assignments",
                "- centers.rds      : prototypes/centers/medoids of clusters (R object)",
                "- inertia_path.csv : within-cluster inertia for K = 2..p (elbow method)",
                "- suppl_predict.csv: cluster attachment for supplementary variables (if any)",
                "",
                "You can load centers.rds in R via:",
                "  centers <- readRDS('centers.rds')"
            ), readme_conn)
            close(readme_conn)
            
            # 7) Create ZIP
            old_wd <- getwd()
            setwd(bundle_dir)
            on.exit(setwd(old_wd), add = TRUE)
            
            files_to_zip <- list.files(bundle_dir, all.files = FALSE)
            utils::zip(zipfile = file, files = files_to_zip)
        }
    )
}

shinyApp(ui, server)