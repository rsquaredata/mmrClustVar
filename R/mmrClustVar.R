#' Variable Clustering (R6 Facade)
#'
#' R6 "facade" class that encapsulates the different implementations of
#' variable clustering (k-means, k-modes, k-prototypes, k-medoids),
#' and provides a unified interface for the end user.
#'
#' The object is created using \code{mmrClustVar$new(...)} and exposes
#' the main methods: \code{$fit()}, \code{$predict()},
#' \code{$print()}, \code{$summary()}, and \code{$plot()}.
#'
#' @docType class
#' @name mmrClustVar
#' @export
NULL

mmrClustVar <- R6::R6Class(
    "mmrClustVar",
    
    public = list(
        
        #' @description
        #' Constructor of the facade class.
        #'
        #' @param method Variable clustering method:
        #'   "kmeans", "kmodes", "kprototypes", "kmedoids", or "auto".
        #' @param K Number of clusters (integer >= 2).
        #' @param scale If TRUE, standardizes the numeric variables.
        #' @param lambda Weighting parameter for categorical distances
        #'   (used in k-prototypes and mixed k-medoids).
        initialize = function(
        method = c("kmeans", "kmodes", "kprototypes", "kmedoids", "auto"),
        K,
        scale  = TRUE,
        lambda = 1,
        ...
        ) {
            method            <- match.arg(method)
            private$method    <- method
            private$K         <- K
            private$scale     <- scale
            private$lambda    <- lambda
            private$engine    <- NULL
            private$X_last    <- NULL
            private$inertia_grid <- NULL
        },
        
        #' @description
        #' Fits the model on a set of active variables.
        #'
        #' @param X A data.frame where columns are variables.
        fit = function(X) {
            if (!is.data.frame(X)) {
                X <- as.data.frame(X)
            }
            
            # Build engine if not already created (or if method = "auto")
            if (is.null(private$engine) || identical(private$method, "auto")) {
                private$engine <- private$build_engine(X)
            }
            
            private$engine$fit(X)
            
            # Store X for inertia path computation
            private$X_last       <- X
            private$inertia_grid <- NULL
            
            invisible(self)
        },
        
        #' @description
        #' Assigns new variables to the existing clusters.
        #'
        #' @param X_new A data.frame of new variables.
        predict = function(X_new) {
            if (is.null(private$engine)) {
                stop("[mmrClustVar] predict() called before fit().")
            }
            private$engine$predict(X_new)
        },
        
        #' @description
        #' Short printed summary of the object + internal engine.
        print = function(...) {
            cat("Class 'mmrClustVar'\n")
            cat("  Requested method :", private$method, "\n")
            cat("  K                :", private$K, "\n")
            cat("  scale            :", private$scale, "\n")
            cat("  lambda           :", private$lambda, "\n\n")
            
            if (!is.null(private$engine)) {
                cat("Internal engine summary:\n")
                private$engine$print()
            } else {
                cat("(No model fitted yet — call $fit(X)).\n")
            }
            invisible(self)
        },
        
        #' @description
        #' Detailed summary: cluster sizes, inertia, and indicators.
        summary = function(...) {
            if (is.null(private$engine)) {
                cat("[mmrClustVar] summary(): no fitted model.\n")
                return(invisible(NULL))
            }
            private$engine$summary(...)
        },
        
        #' @description
        #' Main plots (inertia, clusters, membership, profiles).
        #'
        #' @param type One of:
        #'   "inertia", "clusters", "membership", "profiles".
        plot = function(type = c("inertia", "clusters", "membership", "profiles"),
                        ...) {
            
            type <- match.arg(type)
            
            if (is.null(private$engine)) {
                stop("[mmrClustVar] plot(): no fitted model.")
            }
            
            # If inertia grid already computed
            if (identical(type, "inertia") && !is.null(private$inertia_grid)) {
                df <- private$inertia_grid
                graphics::plot(
                    df$K, df$inertia,
                    type = "b",
                    xlab = "K",
                    ylab = "Within-cluster inertia",
                    main = sprintf("Inertia curve (%s)", private$method)
                )
                return(invisible(NULL))
            }
            
            # Otherwise, delegate to engine
            private$engine$plot(type = type, ...)
        },
        
        #' @description
        #' Computes the intra-cluster inertia for multiple values of K.
        #'
        #' @param K_seq Vector of integers.
        #' @param X Optional data.frame. If NULL, uses the X from fit().
        #' @return data.frame(K, inertia)
        compute_inertia_path = function(K_seq, X = NULL) {
            
            # Retrieve X if not provided
            if (is.null(X)) {
                X <- private$X_last
            }
            if (is.null(X)) {
                stop("[mmrClustVar] compute_inertia_path(): no X available. ",
                     "Provide X or call $fit(X) first.")
            }
            if (!is.data.frame(X)) {
                X <- as.data.frame(X)
            }
            
            p <- ncol(X)
            
            # Keep only valid K (between 2 and p)
            K_seq <- sort(unique(as.integer(K_seq)))
            K_seq <- K_seq[!is.na(K_seq) & K_seq >= 2L & K_seq <= p]
            
            if (length(K_seq) == 0L) {
                stop("[mmrClustVar] compute_inertia_path(): K_seq must contain ",
                     "at least one integer between 2 and ", p, ".")
            }
            
            old_K <- private$K
            
            res <- data.frame(
                K       = K_seq,
                inertia = NA_real_
            )
            
            # Loop to compute inertia for each K
            for (i in seq_along(K_seq)) {
                K_i <- K_seq[i]
                private$K <- K_i
                engine_i  <- private$build_engine(X)
                engine_i$fit(X)
                res$inertia[i] <- engine_i$get_inertia()
            }
            
            private$K            <- old_K
            private$inertia_grid <- res
            
            res
        },
        
        #' @description
        #' Returns the variable cluster assignments.
        get_clusters = function() {
            if (is.null(private$engine)) return(NULL)
            private$engine$get_clusters()
        },
        
        #' @description
        #' Returns the cluster prototypes (centers / modes / medoids).
        get_centers = function() {
            if (is.null(private$engine)) return(NULL)
            private$engine$get_centers()
        },
        
        #' @description
        #' Returns the intra-cluster inertia.
        get_inertia = function() {
            if (is.null(private$engine)) return(NA_real_)
            private$engine$get_inertia()
        },
        
        #' @description
        #' Returns the convergence status.
        get_convergence = function() {
            if (is.null(private$engine)) return(NA)
            private$engine$get_convergence()
        },
        
        #' @description
        #' Returns the requested method.
        get_method = function() private$method,
        
        #' @description
        #' Returns the requested number of clusters K.
        get_K = function() private$K,
        
        #' @description
        #' Human-readable interpretation of variable clusters.
        interpret_clusters = function(...) {
            if (is.null(private$engine)) {
                cat("[mmrClustVar] interpret_clusters(): no fitted model.\n")
                return(invisible(NULL))
            }
            
            # check that the engine exposes this method
            if (!("interpret_clusters" %in% names(private$engine))) {
                cat("[mmrClustVar] interpret_clusters(): not available for this engine.\n")
                return(invisible(NULL))
            }
            
            private$engine$interpret_clusters(...)
        }
    ),
    
    private = list(
        
        method       = NULL,
        K            = NULL,
        scale        = NULL,
        lambda       = NULL,
        
        engine       = NULL,
        X_last       = NULL,
        inertia_grid = NULL,
        
        #' Internal helper: build the appropriate engine
        #' depending on the selected method.
        build_engine = function(X) {
            
            method <- private$method
            K      <- private$K
            scale  <- private$scale
            lambda <- private$lambda
            
            # Automatic method selection based on variable types
            if (identical(method, "auto")) {
                
                is_num <- vapply(X, is.double, logical(1L))
                
                is_cat <- vapply(
                    X,
                    function(col) is.factor(col) || is.character(col) || is.integer(col),
                    logical(1L)
                )
                
                if (all(is_num)) {
                    method_effective <- "kmeans"
                } else if (all(is_cat)) {
                    method_effective <- "kmodes"
                } else {
                    method_effective <- "kprototypes"
                }
                
                cat("[mmrClustVar] method = 'auto' → selected:", method_effective, "\n")
                
            } else {
                method_effective <- method
            }
            
            # Create appropriate engine
            engine <- switch(
                method_effective,
                "kmeans"      = mmrClustVarKMeans$new(K = K, scale = scale, lambda = lambda),
                "kmodes"      = mmrClustVarKModes$new(K = K, scale = scale, lambda = lambda),
                "kprototypes" = mmrClustVarKPrototypes$new(K = K, scale = scale, lambda = lambda),
                "kmedoids"    = mmrClustVarKMedoids$new(K = K, scale = scale, lambda = lambda),
                stop("[mmrClustVar] Unsupported method:", method_effective)
            )
            
            engine
        }
    )
)