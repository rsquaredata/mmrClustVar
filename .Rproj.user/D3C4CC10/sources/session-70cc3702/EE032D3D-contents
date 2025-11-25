#' Launch the demo Shiny application for mmrClustVar
#'
#' @description
#' Launches a Shiny application demonstrating the use of the \code{mmrClustVar} R6 class for clustering variables.
#'
#' @return
#' Nothing. Opens a Shiny app.
#'
#' @examples
#' \dontrun{
#' run_mmrClustVar_app()
#' }
#'
#' @export
run_mmrClustVar_app <- function() {
  app_dir <- system.file("shiny/mmrClustVar_app", package = "mmrClustVar")
  if (app_dir == "") {
    stop("Cannot find app directory. Reinstall the package.", call. = FALSE)
  }
  shiny::runApp(app_dir, display.mode = "normal")
}
