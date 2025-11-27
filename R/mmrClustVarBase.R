#' Base Class for Variable Clustering
#'
#' Abstract R6 class that factors out the common behavior of the different
#' variable clustering algorithms (k-means, k-modes, k-prototypes, k-medoids).
#'
#' This class is not intended to be instantiated directly by end users.
#' It serves as the parent class for all specialized algorithm classes.
#'
#' @docType class
#' @name mmrClustVarBase
#' @keywords internal
#' @noRd
#'
#' @section Methods:
#' \describe{
#'   \item{\code{$initialize(method_name, K, scale = TRUE, lambda = 1, ...)}}{
#'     Initializes the general model parameters.
#'   }
#'   \item{\code{$fit(X)}}{
#'     Handles the overall training workflow (data preparation,
#'     call to the specific algorithm, storage of results).
#'   }
#'   \item{\code{$predict(X_new)}}{
#'     Assigns additional variables to the learned clusters.
#'   }
#'   \item{\code{$print()}}{
#'     Provides a short printed summary of the object.
#'   }
#'   \item{\code{$summary()}}{
#'     Detailed summary of the results.
#'   }
#'   \item{\code{$plot(type)}}{
#'     Visualizations associated with the model.
#'   }
#'   \item{\code{$get_clusters()}}{
#'     Returns the cluster assignment of variables.
#'   }
#'   \item{\code{$get_centers()}}{
#'     Returns the cluster prototypes.
#'   }
#'   \item{\code{$get_method()}}{
#'     Returns the name of the clustering method used.
#'   }
#'   \item{\code{$get_K()}}{
#'     Returns the number of clusters.
#'   }
#'   \item{\code{$get_convergence()}}{
#'     Returns the convergence indicator of the algorithm.
#'   }
#' }
NULL

mmrClustVarBase <- R6::R6Class(
    "mmrClustVarBase",
    
    public = list(
        
        initialize = function(K,
                              scale = TRUE,
                              lambda = 1,
                              method_name = "base") {
            
            if (missing(K) || !is.numeric(K) || length(K) != 1L || K < 2) {
                stop("[mmrClustVarBase] K must be an integer >= 2")
            }
            
            private$FNbGroupes   <- as.integer(K)
            private$FScale       <- isTRUE(scale)
            private$FLambda      <- as.numeric(lambda)
            private$FMethod      <- method_name
            private$FConvergence <- FALSE
            private$FInertia     <- NA_real_
            private$FAlgorithme  <- method_name
        },
        
        fit = function(X) {
            
            # 1) Basic checks and internal structure definition
            X <- private$check_and_prepare_X(X, update_structure = TRUE)
            
            # 2) Filter out rows with missing or non-finite values
            #    (important for prcomp() in k-means / numeric part of k-prototypes)
            num_idx <- private$FNumCols
            ok <- stats::complete.cases(X)
            
            if (length(num_idx) > 0L) {
                X_num <- as.matrix(X[, num_idx, drop = FALSE])
                ok_num <- apply(X_num, 1L, function(row) {
                    all(is.finite(row) | is.na(row))
                })
                ok <- ok & ok_num
            }
            
            if (!all(ok)) {
                X <- X[ok, , drop = FALSE]
            }
            
            # 3) Optional standardization of numeric variables
            if (private$FScale) {
                X <- private$scale_active_variables(X)
            }
            
            # 4) Run the algorithm (implemented in child classes)
            res <- private$run_clustering(X)
            
            # 5) Store all shared results
            private$FX_active    <- X
            private$FClusters    <- res$clusters
            private$FCenters     <- res$centers
            private$FInertia     <- res$inertia
            private$FConvergence <- isTRUE(res$converged)
            
            invisible(self)
        },
        
        predict = function(X_new) {
            
            if (is.null(private$FX_active)) {
                stop("[mmrClustVarBase] fit() must be called before predict()")
            }
            
            X_new <- private$check_and_prepare_X(X_new, update_structure = FALSE)
            
            res <- lapply(seq_along(X_new), function(j) {
                var_name <- colnames(X_new)[j]
                private$predict_one_variable(X_new[[j]], var_name)
            })
            
            out <- do.call(rbind, res)
            rownames(out) <- NULL
            private$FX_new <- X_new
            
            out
        },
        
        print = function(...) {
            cat("Class", class(self)[1L], "\n")
            cat("  Method           :", private$FMethod, "\n")
            cat("  K (n groups)     :", private$FNbGroupes, "\n")
            cat("  # active vars    :", length(private$FClusters), "\n")
            cat("  Convergence      :", private$FConvergence, "\n")
            cat("  Within inertia   :", private$FInertia, "\n")
            invisible(self)
        },
        
        summary = function(...) {
            if (is.null(private$FClusters)) {
                cat("Model not fitted yet (fit() not called).\n")
                return(invisible(NULL))
            }
            
            cat("=== Global Summary ===\n")
            cat("Method      :", private$FMethod, "\n")
            cat("K           :", private$FNbGroupes, "\n")
            cat("Convergence :", private$FConvergence, "\n")
            cat("Inertia     :", private$FInertia, "\n\n")
            
            clusters   <- private$FClusters
            tab_sizes  <- table(clusters)
            cat("=== Cluster Sizes ===\n")
            print(tab_sizes)
            cat("\n")
            
            # Hook for child classes: specific membership metrics
            private$summary_membership()
            
            invisible(NULL)
        },
        
        plot = function(type = c("inertia", "clusters", "membership", "profiles"), ...) {
            
            type <- match.arg(type)
            
            if (is.null(private$FClusters)) {
                stop("[mmrClustVarBase] fit() must be called before plot().")
            }
            
            if (type == "clusters") {
                
                barplot(
                    table(private$FClusters),
                    main = sprintf("Variable distribution (%s)", private$FMethod),
                    xlab = "Cluster",
                    ylab = "Number of variables"
                )
                
            } else if (type == "inertia") {
                
                barplot(
                    private$FInertia,
                    main = sprintf("Within-cluster inertia (%s)", private$FMethod),
                    ylab = "Inertia"
                )
                
            } else if (type == "membership") {
                
                private$plot_membership()
                
            } else if (type == "profiles") {
                
                private$plot_profiles()
            }
            
            invisible(NULL)
        },
        
        get_clusters    = function() private$FClusters,
        get_centers     = function() private$FCenters,
        get_method      = function() private$FMethod,
        get_K           = function() private$FNbGroupes,
        get_inertia     = function() private$FInertia,
        get_convergence = function() private$FConvergence
    ),
    
    private = list(
        
        FMethod      = NULL,
        FNbGroupes   = NULL,
        FScale       = NULL,
        FLambda      = NULL,
        FX_active    = NULL,
        FX_new       = NULL,
        FClusters    = NULL,
        FCenters     = NULL,
        FInertia     = NULL,
        FConvergence = NULL,
        FAlgorithme  = NULL,
        FNumCols     = NULL,
        FCatCols     = NULL,
        
        check_and_prepare_X = function(X, update_structure = TRUE) {
            
            if (!is.data.frame(X)) {
                X <- as.data.frame(X)
            }
            
            if (ncol(X) == 0L) {
                stop("[mmrClustVarBase] X must contain at least one column.")
            }
            
            num_idx <- which(vapply(X, is.numeric, logical(1L)))
            cat_idx <- which(vapply(
                X,
                function(col) is.factor(col) || is.character(col),
                logical(1L)
            ))
            
            if (update_structure) {
                private$FNumCols <- num_idx
                private$FCatCols <- cat_idx
            }
            
            # Convert character → factor
            for (j in cat_idx) {
                if (!is.factor(X[[j]])) {
                    X[[j]] <- factor(X[[j]])
                }
            }
            
            X
        },
        
        scale_active_variables = function(X) {
            num_idx <- private$FNumCols
            if (length(num_idx) == 0L) return(X)
            
            X[num_idx] <- lapply(X[num_idx], function(col) {
                m <- mean(col, na.rm = TRUE)
                s <- stats::sd(col, na.rm = TRUE)
                if (is.na(s) || s == 0) return(rep(0, length(col)))
                (col - m) / s
            })
            
            X
        },
        
        run_clustering = function(X) {
            stop("[mmrClustVarBase] run_clustering() must be implemented in a child class.")
        },
        
        predict_one_variable = function(x_new, var_name) {
            stop("[mmrClustVarBase] predict_one_variable() must be implemented in a child class.")
        },
        
        summary_membership = function() {
            cat("(No membership indicators defined for this class.)\n")
        },
        
        plot_membership = function() {
            stop("[mmrClustVarBase] plot(type = 'membership') not implemented for this class.")
        },
        
        plot_profiles = function() {
            warning("[mmrClustVarBase] plot(type = 'profiles') not implemented; no plot produced.")
        }
    )
)