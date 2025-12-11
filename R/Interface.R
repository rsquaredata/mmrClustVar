#' Interface class for variable clustering (Interface)
#'
#' High-level R6 interface to multiple variable clustering algorithms
#' (k-means, k-modes, k-prototypes, k-medoids) implemented in dedicated
#' internal engines inheriting from \code{ClusterBase}.
#'
#' Users can interact with this interface:
#' \enumerate{
#'   \item choose a method (or \code{"auto"}),
#'   \item call \code{$fit(X)} on a data.frame of active variables,
#'   \item inspect \code{$print()}, \code{$summary()}, \code{$plot()}, etc.
#' }
#'
#' @docType class
#' @name Interface
#' @export
Interface <- R6::R6Class(
  "Interface",

  public = list(

    #' @description
    #' Create a new Interface object.
    initialize = function(method = c("auto", "kmeans", "kmodes", "kprototypes", "kmedoids"),
                          K      = 3,
                          scale  = TRUE,
                          lambda = 1) {

      method <- match.arg(method)

      if (!is.numeric(K) || length(K) != 1L || K < 2)
        stop("[Interface] K must be a numeric >= 2.")
      if (!is.logical(scale) || length(scale) != 1L)
        stop("[Interface] 'scale' must be TRUE or FALSE.")
      if (!is.numeric(lambda) || length(lambda) != 1L || lambda <= 0)
        stop("[Interface] 'lambda' must be a numeric > 0.")

      private$FRequestedMethod <- method
      private$FK               <- as.integer(K)
      private$FScale           <- scale
      private$FLambda          <- lambda

      private$FEngine          <- NULL
      private$FEffectiveMethod <- NA_character_
      private$FInertiaPath     <- NULL
    },

    #' @description Fit the model
    fit = function(X) {

      if (!is.data.frame(X))
        X <- as.data.frame(X, stringsAsFactors = FALSE)

      if (ncol(X) < 2)
        stop("[Interface] At least 2 active variables are required.")

      method_eff <- private$decide_effective_method(X)
      engine     <- private$create_engine(method_eff)

      engine$fit(X)

      private$FEngine          <- engine
      private$FEffectiveMethod <- method_eff
      private$FInertiaPath     <- NULL

      invisible(self)
    },

    #' @description Predict supplementary variables
    predict = function(X_new) {
      if (is.null(private$FEngine))
        stop("[Interface] No fitted model: call fit() first.")
      private$FEngine$predict(X_new)
    },

    #' @description Print summary
    print = function(...) {
      cat("Classe 'Interface'\n")
      cat("  Requested method :", private$FRequestedMethod, "\n")
      cat("  Effective method :", private$FEffectiveMethod, "\n")
      cat("  K (clusters)     :", private$FK, "\n")
      cat("  scale            :", private$FScale, "\n")
      cat("  lambda           :", private$FLambda, "\n\n")

      if (!is.null(private$FEngine)) {
        cat("--- Engine summary ---\n")
        private$FEngine$print(...)
      } else {
        cat("(no fitted engine yet)\n")
      }

      invisible(self)
    },

    #' @description Detailed summary
    summary = function(...) {
      cat("=== Interface summary ===\n")
      cat("Requested method :", private$FRequestedMethod, "\n")
      cat("Effective method :", private$FEffectiveMethod, "\n")
      cat("K (clusters)     :", private$FK, "\n")
      cat("scale            :", private$FScale, "\n")
      cat("lambda           :", private$FLambda, "\n\n")

      if (is.null(private$FEngine)) {
        cat("(no fitted engine yet)\n")
        return(invisible(NULL))
      }

      cat("=== Engine-level summary ===\n\n")
      private$FEngine$summary(...)
    },

    #' @description Plot interface
    plot = function(type = c("inertia", "clusters", "membership", "profiles"), ...) {

      type <- match.arg(type)

      # -------------------------------------------------------
      # INERTIA CASE: elbow curve OR inertia object
      # -------------------------------------------------------
      if (type == "inertia") {

        # 1) Inertia path exists → elbow plot
        if (!is.null(private$FInertiaPath)) {

          df <- private$FInertiaPath
          if (!all(c("K", "inertia") %in% names(df)))
            stop("[Interface] Inertia path has an invalid structure.")

          graphics::plot(
            df$K,
            df$inertia,
            type = "b",
            xlab = "Number of clusters K",
            ylab = "Within-cluster inertia",
            main = "Inertia vs K (elbow)",
            ...
          )

          return(invisible(NULL))
        }

        # 2) No path — use the inertia object of the fitted engine
        if (!is.null(private$FEngine)) {
          obj <- private$FEngine$get_inertia_object()
          if (!is.null(obj)) {
            obj$plot(...)
            return(invisible(NULL))
          }
        }

        stop("[Interface] No inertia or inertia path available.")
      }

      # -------------------------------------------------------
      # OTHER PLOT TYPES
      # -------------------------------------------------------
      if (is.null(private$FEngine))
        stop("[Interface] No fitted model: call fit() first.")

      private$FEngine$plot(type = type, ...)
    },

    #' @description Compute inertia path
    compute_inertia_path = function(K_seq, X = NULL) {

      if (missing(K_seq) || length(K_seq) == 0L)
        stop("[Interface] K_seq must be a non-empty vector.")

      K_seq <- sort(unique(as.integer(K_seq)))
      K_seq <- K_seq[K_seq >= 2]

      if (length(K_seq) == 0)
        stop("[Interface] K_seq must contain integers >= 2.")

      if (is.null(X)) {
        if (is.null(private$FEngine))
          stop("[Interface] No fitted engine and no X provided.")
        X <- private$FEngine$get_X_descr()$X_active
      }

      if (!is.data.frame(X))
        X <- as.data.frame(X, stringsAsFactors = FALSE)

      p <- ncol(X)
      K_seq <- K_seq[K_seq <= p]
      if (length(K_seq) == 0)
        stop("[Interface] All K in K_seq exceed number of variables.")

      method_eff <- private$decide_effective_method(X)

      inertias <- numeric(length(K_seq))

      for (i in seq_along(K_seq)) {
        Ki <- K_seq[i]
        engine_i <- private$create_engine(method_eff, K_override = Ki)

        inertias[i] <- tryCatch({
          engine_i$fit(X)
          engine_i$get_inertia()
        }, error = function(e) NA_real_)
      }

      df <- data.frame(K = K_seq, inertia = inertias)
      private$FInertiaPath <- df

      invisible(df)
    },

    #' @description Interpret clusters
    interpret_clusters = function(style = c("compact", "detailed")) {
      if (is.null(private$FEngine))
        stop("[Interface] No fitted model: call fit() first.")
      private$FEngine$interpret_clusters(style = style)
    },

    #' @description Get clusters
    get_clusters = function() {
      if (is.null(private$FEngine)) return(NULL)
      private$FEngine$get_clusters()
    },

    #' @description Get inertia
    get_inertia = function() {
      if (is.null(private$FEngine)) return(NA_real_)
      private$FEngine$get_inertia()
    },

    #' @description Get inertia object
    get_inertia_object = function() {
      if (is.null(private$FEngine)) return(NULL)
      private$FEngine$get_inertia_object()
    },

    #' @description Get method
    get_method = function() private$FEffectiveMethod,

    #' @description Get centers
    get_centers = function() {
      if (is.null(private$FEngine)) return(NULL)
      private$FEngine$get_centers()
    },

    #' @description Get convergence flag
    get_convergence = function() {
      if (is.null(private$FEngine)) return(NA)
      private$FEngine$get_convergence()
    }
  ),

  private = list(

    # Requested params
    FRequestedMethod = NA_character_,
    FK               = NA_integer_,
    FScale           = TRUE,
    FLambda          = 1,

    # Underlying R6 engine
    FEngine          = NULL,

    # Effective method after auto-detection
    FEffectiveMethod = NA_character_,

    # Elbow curve inertia path
    FInertiaPath     = NULL,

    # ------------------------------------------------------------------
    # Decide effective method
    # ------------------------------------------------------------------
    decide_effective_method = function(X) {

      method_req <- private$FRequestedMethod

      is_num <- vapply(X, is.numeric, logical(1L))
      is_cat <- vapply(
        X,
        function(col) is.factor(col) || is.character(col),
        logical(1L)
      )

      has_num <- any(is_num)
      has_cat <- any(is_cat)

      if (method_req != "auto") {
        if (method_req == "kmeans" && !all(is_num))
          stop("[Interface] k-means requires all numeric variables.")
        if (method_req == "kmodes" && !all(is_cat))
          stop("[Interface] k-modes requires all categorical variables.")
        return(method_req)
      }

      if (has_num && !has_cat) return("kmeans")
      if (!has_num && has_cat) return("kmodes")
      return("kprototypes")
    },

    # ------------------------------------------------------------------
    # Create clustering engine
    # ------------------------------------------------------------------
    create_engine = function(method_effective, K_override = NULL) {

      K_val <- if (is.null(K_override)) private$FK else as.integer(K_override)

      if (method_effective == "kmeans") {
        return(Kmeans$new(
          K      = K_val,
          scale  = private$FScale
        ))
      }

      if (method_effective == "kmodes") {
        return(Kmodes$new(
          K      = K_val
        ))
      }

      if (method_effective == "kprototypes") {
        return(Kprototypes$new(
          K      = K_val,
          scale  = private$FScale,
          lambda = private$FLambda
        ))
      }

      if (method_effective == "kmedoids") {
        return(Kmedoids$new(
          K      = K_val,
          lambda = private$FLambda
        ))
      }

      stop("[Interface] Unsupported method: ", method_effective)
    }
  )
)
