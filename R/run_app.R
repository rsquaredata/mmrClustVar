#' Launch the mmrClustVar Shiny application
#'
#' This function launches the Shiny app included in the package,
#' located under \code{inst/shiny/mmrClustVar_app}.
#'
#' @return
#' This function is called for its side-effect (launching the app).
#' @export
#'
#' @examples
#' \dontrun{
#' run_mmrClustVar_app()
#' }
run_mmrClustVar_app <- function() {
    
    app_dir <- system.file("shiny/mmrClustVar_app", package = "mmrClustVar")
    
    message("Chemin Shiny utilisé : ", app_dir)
    
    if (app_dir == "" || !file.exists(file.path(app_dir, "app.R"))) {
        stop("Impossible de trouver l'application Shiny mmrClustVar_app dans inst/shiny/mmrClustVar_app.")
    }
    
    shiny::runApp(app_dir, display.mode = "normal")
}

