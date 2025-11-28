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
    #'
    #' @param method Clustering method: one of \code{"auto"}, \code{"kmeans"},
    #'   \code{"kmodes"}, \code{"kprototypes"}, \code{"kmedoids"}.
    #' @param K Number of clusters (K ≥ 2).
    #' @param scale Logical, standardize numeric variables before modelling?
    #' @param lambda Positive numeric weight for the categorical part (used
    #'   by k-prototypes and k-medoids).
    initialize = function(method = c("auto", "kmeans", "kmodes", "kprototypes", "kmedoids"),
                          K      = 3,
                          scale  = TRUE,
                          lambda = 1) {

      method <- match.arg(method)

      if (!is.numeric(K) || length(K) != 1L || K < 2) {
        stop("[Interface] K must be a numeric >= 2.")
      }
      if (!is.logical(scale) || length(scale) != 1L) {
        stop("[Interface] 'scale' must be TRUE or FALSE.")
      }
      if (!is.numeric(lambda) || length(lambda) != 1L || lambda <= 0) {
        stop("[Interface] 'lambda' must be a numeric > 0.")
      }

      private$FRequestedMethod <- method
      private$FK               <- as.integer(K)
      private$FScale           <- scale
      private$FLambda          <- lambda

      private$FEngine          <- NULL
      private$FEffectiveMethod <- NA_character_
      private$FInertiaPath     <- NULL
    },

    #' @description
    #' Fit the variable clustering model on active variables \code{X}.
    #'
    #' @param X data.frame (or coercible) of active variables to be clustered.
    fit = function(X) {
      if (!is.data.frame(X)) {
        X <- as.data.frame(X, stringsAsFactors = FALSE)
      }

      if (ncol(X) < 2L) {
        stop("[Interface] At least 2 active variables are required.")
      }

      method_eff <- private$decide_effective_method(X)
      engine     <- private$create_engine(method_eff)

      # Fit underlying engine (inherits ClusterBase)
      engine$fit(X)

      private$FEngine          <- engine
      private$FEffectiveMethod <- method_eff
      private$FInertiaPath     <- NULL  # reset any previous inertia path

      invisible(self)
    },

    #' @description
    #' Attach supplementary variables to existing clusters.
    #'
    #' @param X_new data.frame of supplementary variables.
    #' @return A data.frame with cluster assignment and membership indicators.
    predict = function(X_new) {
      if (is.null(private$FEngine)) {
        stop("[Interface] No fitted model: call fit() first.")
      }
      private$FEngine$predict(X_new)
    },

    #' @description
    #' Print a concise summary of the fitted model.
    #' @param ... Regular print params
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

    #' @description
    #' Detailed summary: facade info + engine-level summary.
    #' @param ... Regular summary params
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

    #' @description
    #' High-level plotting interface.
    #'
    #' @param type Plot type:
    #'   \itemize{
    #'     \item \code{"inertia"}        : inertia vs K (requires
    #'           \code{compute_inertia_path()} to have been called),
    #'     \item \code{"clusters"}       : cluster sizes / distribution,
    #'     \item \code{"membership"}     : membership indicators per variable,
    #'     \item \code{"profiles"}       : profiles / heatmaps on individuals.
    #'   }
    #' @param ... Regular plot params
    plot = function(type = c("inertia", "clusters", "membership", "profiles"), ...) {
      type <- match.arg(type)

      if (type == "inertia") {
        if (is.null(private$FInertiaPath)) {
          stop("[Interface] No inertia path available. Call compute_inertia_path() first.")
        }

        df <- private$FInertiaPath
        if (!all(c("K", "inertia") %in% names(df))) {
          stop("[Interface] Inertia path has an invalid structure.")
        }

        graphics::plot(
          df$K, df$inertia,
          type = "b",
          xlab = "Number of clusters K",
          ylab = "Within-cluster inertia",
          main = "Inertia vs K (elbow method)",
          ...
        )
        return(invisible(NULL))
      }

      if (is.null(private$FEngine)) {
        stop("[Interface] No fitted model: call fit() first.")
      }

      private$FEngine$plot(type = type, ...)
    },

    #' @description
    #' Compute an inertia path for a sequence of K values, keeping the
    #' same method / scale / lambda as the current facade.
    #'
    #' @param K_seq Integer vector of K values (e.g. 2:10).
    #' @param X Optional data.frame of active variables. If \code{NULL},
    #'   reuse the active X stored in the fitted engine.
    #' @return Invisibly, a data.frame with columns \code{K} and \code{inertia}.
    compute_inertia_path = function(K_seq, X = NULL) {

      if (missing(K_seq) || length(K_seq) == 0L) {
        stop("[Interface] K_seq must be a non-empty vector of integers.")
      }
      K_seq <- unique(as.integer(K_seq))
      K_seq <- sort(K_seq)
      K_seq <- K_seq[K_seq >= 2L]

      if (length(K_seq) == 0L) {
        stop("[Interface] K_seq must contain values >= 2.")
      }

      # Active data
      if (is.null(X)) {
        if (is.null(private$FEngine)) {
          stop("[Interface] No fitted model and no X provided.")
        }
        X <- private$FEngine$get_X_descr()$X_active
      }

      if (!is.data.frame(X)) {
        X <- as.data.frame(X, stringsAsFactors = FALSE)
      }

      p <- ncol(X)
      K_seq <- K_seq[K_seq <= p]
      if (length(K_seq) == 0L) {
        stop("[Interface] All K in K_seq exceed the number of variables.")
      }

      method_eff <- private$decide_effective_method(X)

      inertias <- numeric(length(K_seq))

      for (i in seq_along(K_seq)) {
        Ki <- K_seq[i]
        engine_i <- private$create_engine(method_eff, K_override = Ki)

        # we ignore potential errors here? better to catch them and set NA
        inertias[i] <- tryCatch({
          engine_i$fit(X)
          engine_i$get_inertia()
        }, error = function(e) {
          NA_real_
        })
      }

      df <- data.frame(
        K       = K_seq,
        inertia = inertias
      )

      private$FInertiaPath <- df
      invisible(df)
    },

    #' @description
    #' Textual interpretation of the clustering solution.
    #'
    #' Delegates to the underlying engine's \code{interpret_clusters()} method.
    #'
    #' @param style Either \code{"compact"} or \code{"detailed"}.
    interpret_clusters = function(style = c("compact", "detailed")) {
      if (is.null(private$FEngine)) {
        stop("[Interface] No fitted model: call fit() first.")
      }
      private$FEngine$interpret_clusters(style = style)
    },

    # ---- Simple getters used by Shiny and external code -------------------

    #' @description Get cluster assignments (one per active variable).
    get_clusters = function() {
      if (is.null(private$FEngine)) return(NULL)
      private$FEngine$get_clusters()
    },

    #' @description Get within-cluster inertia of the fitted model.
    get_inertia = function() {
      if (is.null(private$FEngine)) return(NA_real_)
      private$FEngine$get_inertia()
    },

    #' @description Get the effective method actually used by the engine.
    get_method = function() {
      private$FEffectiveMethod
    },

    #' @description Get cluster prototypes / centers / medoids.
    #'
    #' For k-means and k-modes, this is an internal representation.
    #' For k-medoids, this is a data.frame of medoid variables.
    get_centers = function() {
      if (is.null(private$FEngine)) return(NULL)
      private$FEngine$get_centers()
    },

    #' @description Get the convergence flag of the fitted engine.
    get_convergence = function() {
      if (is.null(private$FEngine)) return(NA)
      private$FEngine$get_convergence()
    }
  ),

  private = list(

    # requested parameters
    FRequestedMethod = NA_character_,
    FK               = NA_integer_,
    FScale           = TRUE,
    FLambda          = 1,

    # underlying fitted engine (ClusterBase subclass)
    FEngine          = NULL,

    # effective method used (after "auto" decision)
    FEffectiveMethod = NA_character_,

    # inertia path data.frame (K, inertia)
    FInertiaPath     = NULL,

    # ------------------------------------------------------------------
    # Decide effective method from requested method + data types
    # ------------------------------------------------------------------
    decide_effective_method = function(X) {
      method_req <- private$FRequestedMethod

      # Determine variable types
      is_num <- vapply(X, is.numeric, logical(1L))
      is_cat <- vapply(
        X,
        function(col) is.factor(col) || is.character(col),
        logical(1L)
      )

      has_num <- any(is_num)
      has_cat <- any(is_cat)

      if (method_req != "auto") {
        # basic consistency checks
        if (method_req == "kmeans" && !all(is_num)) {
          stop("[Interface] k-means requires all active variables to be numeric.")
        }
        if (method_req == "kmodes" && !all(is_cat)) {
          stop("[Interface] k-modes requires all active variables to be categorical.")
        }
        # k-prototypes / k-medoids can handle mixed or homogeneous types
        return(method_req)
      }

      # method = "auto"
      if (has_num && !has_cat) {
        return("kmeans")
      }
      if (!has_num && has_cat) {
        return("kmodes")
      }

      # mixed case: choose k-prototypes by default
      return("kprototypes")
    },

    # ------------------------------------------------------------------
    # Create an engine from a method name
    # ------------------------------------------------------------------
    create_engine = function(method_effective, K_override = NULL) {
      K_val <- if (is.null(K_override)) private$FK else as.integer(K_override)

      if (method_effective == "kmeans") {
        return(Kmeans$new(
          K      = K_val,
          scale  = private$FScale,
          lambda = private$FLambda
        ))
      }

      if (method_effective == "kmodes") {
        return(Kmodes$new(
          K      = K_val,
          scale  = private$FScale,
          lambda = private$FLambda
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
          scale  = private$FScale,
          lambda = private$FLambda
        ))
      }

      stop("[Interface] Unsupported effective method: ", method_effective)
    }
  )
)
