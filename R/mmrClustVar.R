#' mmrClustVar: R6 Class for Clustering of Variables
#'
#' @description
#' An R6 class implementing clustering of variables using:
#' \itemize{
#'   \item \strong{k-means} for numeric variables (correlation-based distance),
#'   \item \strong{k-modes} for categorical variables (simple matching dissimilarity),
#'   \item \strong{k-prototypes} for mixed data (weighted combination of both),
#'   \item \strong{k-medoids} as a robust alternative based on a dissimilarity
#'         matrix between variables.
#' }
#'
#' Variables are treated as objects: each variable is represented by its
#' vector of values across individuals. The algorithm partitions the \eqn{p}
#' variables into \eqn{K} clusters and, for each cluster, computes a latent
#' component, a categorical prototype or a medoid profile used to measure
#' similarity.
#'
#' @section Usage:
#' \preformatted{
#' library(mmrClustVar)
#'
#' obj <- mmrClustVar$new(method = "kmeans", K = 3, scale = TRUE)
#' obj$fit(iris[, 1:4])
#' obj$summary()
#' obj$plot(type = "membership")
#' }
#'
#' @section Public methods:
#' \describe{
#'   \item{\code{$new(method, K, scale, lambda)}}{Create a new clustering object.}
#'   \item{\code{$fit(X)}}{Fit the model on active variables (columns = variables,
#'         rows = individuals).}
#'   \item{\code{$predict(X_new)}}{Attach supplementary variables to existing
#'         clusters (same individuals as in \code{X}).}
#'   \item{\code{$summary()}}{Display a summary with cluster sizes and membership
#'         scores.}
#'   \item{\code{$plot(type, Ks)}}{Plot diagnostics: cluster sizes, membership
#'         scores, or inertia curve as a function of \code{K}.}
#'   \item{\code{$get_clusters()}}{Return the vector of cluster assignments
#'         for the active variables.}
#'   \item{\code{$get_inertia()}}{Return the total within-cluster inertia.}
#'   \item{\code{$get_method()}}{Return the fitted method ("kmeans", "kmodes",
#'         "kprototypes" or "kmedoids").}
#'   \item{\code{$get_centers()}}{Return the list of latent components,
#'         prototypes or medoid profiles.}
#'   \item{\code{$get_convergence()}}{Return a logical indicating whether the
#'         algorithm converged.}
#'   \item{\code{$get_X_new()}}{Return the last set of supplementary variables
#'         passed to \code{$predict()}.}
#'   \item{\code{$get_algorithm()}}{Return a textual description of the
#'         exact algorithm used.}
#' }
#'
#' @param method Character. One of:
#'   \itemize{
#'     \item \code{"kmeans"} – numeric variables only,
#'     \item \code{"kmodes"} – categorical variables only,
#'     \item \code{"kprototypes"} – mixed data,
#'     \item \code{"kmedoids"} – robust clustering on a dissimilarity between variables,
#'     \item \code{"auto"} – automatic selection based on variable types.
#'   }
#'
#' @param K Integer. Number of clusters (must be >= 2).
#' @param scale Logical. If \code{TRUE}, numeric variables are standardized
#'   before clustering.
#' @param lambda Positive numeric. Weight of the categorical dissimilarity in
#'   k-prototypes (ignored for k-means and k-modes).
#'
#' @return
#' An R6Class generator object. Use \code{$new()} to create a clustering object,
#' then \code{$fit()}, \code{$summary()}, \code{$plot()} and the various getters
#' to interpret the results.
#'
#' @examples
#' \dontrun{
#' ## 1) k-means on numeric variables -----------------------------------------
#'
#' # Use only numeric variables from iris
#' X_num <- iris[, 1:4]
#'
#' obj_km <- mmrClustVar$new(method = "kmeans", K = 2, scale = TRUE)
#' obj_km$fit(X_num)
#'
#' # Summary of clusters and membership scores
#' obj_km$summary()
#'
#' # Plot membership scores by variable
#' obj_km$plot(type = "membership")
#'
#' # Attach a supplementary numeric variable
#' res_pred_num <- obj_km$predict(iris["Sepal.Length"])
#'
#'
#' ## 2) k-modes on categorical variables ------------------------------------
#'
#' set.seed(123)
#' df_cat <- data.frame(
#'   Color = factor(sample(c("red", "blue", "green"), 50, replace = TRUE)),
#'   Shape = factor(sample(c("circle", "square"), 50, replace = TRUE)),
#'   Size  = factor(sample(c("S", "M", "L"), 50, replace = TRUE))
#' )
#'
#' obj_kmodes <- mmrClustVar$new(method = "kmodes", K = 2)
#' obj_kmodes$fit(df_cat)
#'
#' # Cluster assignments for the active variables
#' obj_kmodes$get_clusters()
#'
#' # Summary and membership
#' obj_kmodes$summary()
#'
#'
#' ## 3) k-prototypes on mixed data ------------------------------------------
#'
#' set.seed(123)
#' df_mixed <- data.frame(
#'   Income  = rnorm(80, mean = 2000, sd = 500),
#'   Age     = rnorm(80, mean = 35,   sd = 10),
#'   City    = factor(sample(c("Paris", "Lyon", "Marseille"), 80, replace = TRUE)),
#'   Segment = factor(sample(c("A", "B"), 80, replace = TRUE))
#' )
#'
#' obj_kprot <- mmrClustVar$new(method = "kprototypes", K = 3, lambda = 1)
#' obj_kprot$fit(df_mixed)
#'
#' # Inspect fitted method and inertia
#' obj_kprot$get_method()
#' obj_kprot$get_inertia()
#'
#' # Attach supplementary mixed variables
#' df_suppl <- data.frame(
#'   Bonus = rnorm(80),
#'   Label = factor(sample(c("yes", "no"), 80, replace = TRUE))
#' )
#' res_pred_mixed <- obj_kprot$predict(df_suppl)
#'
#'
#' ## 4) Inertia curve to help choose K --------------------------------------
#'
#' obj_elbow <- mmrClustVar$new(method = "kmeans", K = 2)
#' obj_elbow$fit(X_num)
#'
#' # Compute inertia for several K values and plot the curve
#' obj_elbow$plot(type = "inertia", Ks = 2:6)
#' }
#'
#' @export

mmrClustVar <- R6::R6Class(
    "mmrClustVar",
    
    private = list(
        # --- Attributs internes ---
        FMethod      = NULL,  # "kmeans", "kmodes", "kprototypes", "kmedoids", "auto"
        FNbGroupes   = NULL,  # K
        FScale       = NULL,  # TRUE / FALSE (standardisation des variables quantitatives actives)
        FLambda      = NULL,  # pondération partie catégorielle (k-prototypes)
        
        FX_active    = NULL,  # data.frame des variables actives (pré-traitées)
        FClusters    = NULL,  # vecteur d'affectation des variables actives
        FCenters     = NULL,  # centres / modes / prototypes / medoids selon la méthode
        FInertia     = NULL,  # inertie intra-cluster totale
        FConvergence = NULL,  # booléen
        FX_new       = NULL,  # variables supplémentaires fournies à predict()
        FAlgorithme  = NULL,  # nom textuel de la variante exacte de l'algo
        
        # Indices de structure (pour les données d'apprentissage)
        FNumCols     = NULL,  # indices des variables quantitatives dans FX_active
        FCatCols     = NULL,  # indices des variables qualitatives dans FX_active
        
        # --- Helpers internes ---
        
        check_and_prepare_X = function(X, update_structure = TRUE) {
            if (!is.data.frame(X)) {
                stop("X doit être un data.frame.")
            }
            if (nrow(X)  < 2L) {
                stop("X doit contenir au moins 2 individus.")
            }
            
            # Si on met à jour la structure (cas fit()) -> au moins 2 variables actives
            if (update_structure && ncol(X) < 2L) {
                stop("X doit contenir au moins 2 variables actives.")
            }
            
            # Si on ne met pas à jour la structure (cas predict()) -> au moins 1 variable supplémentaire
            if (!update_structure && ncol(X) < 1L) {
                stop("X_new doit contenir au moins 1 variable supplémentaire.")
            }
            
            # Identification des colonnes numériques / catégorielles
            num_idx <- which(vapply(X, is.numeric, logical(1L)))
            cat_idx <- which(vapply(X, function(col) is.factor(col) || is.character(col), logical(1L)))
            
            # Conversion des qualitatives en facteur
            if (length(cat_idx) > 0L) {
                for (j in cat_idx) {
                    if (!is.factor(X[[j]])) {
                        X[[j]] <- factor(X[[j]])
                    }
                }
            }
            
            if (update_structure) {
                private$FNumCols <- num_idx
                private$FCatCols <- cat_idx
                
                if (length(num_idx) == 0L && private$FMethod == "kmeans") {
                    stop("Aucune variable quantitative trouvée pour la méthode 'kmeans'.")
                }
                if (length(cat_idx) == 0L && private$FMethod == "kmodes") {
                    stop("Aucune variable qualitative trouvée pour la méthode 'kmodes'.")
                }
            }
            
            return(X)
        },
        
        scale_active_variables = function(X) {
            if (!isTRUE(private$FScale)) return(X)
            num_idx <- private$FNumCols
            if (length(num_idx) == 0L) return(X)
            
            for (j in num_idx) {
                col <- X[[j]]
                mu  <- mean(col, na.rm = TRUE)
                sdv <- stats::sd(col, na.rm = TRUE)
                if (is.na(sdv) || sdv == 0) {
                    X[[j]] <- col
                } else {
                    X[[j]] <- (col - mu) / sdv
                }
            }
            return(X)
        },
        
        compute_num_profile = function(X, vars_k, num_idx) {
            num_in_cluster <- intersect(vars_k, num_idx)
            if (length(num_in_cluster) == 0L) return(NULL)
            
            Xk_num <- as.matrix(X[, num_in_cluster, drop = FALSE])
            rowMeans(Xk_num, na.rm = TRUE)  # profil numérique : moyenne par individu
        },
        
        compute_cat_profile = function(X, vars_k, cat_idx) {
            cat_in_cluster <- intersect(vars_k, cat_idx)
            if (length(cat_in_cluster) == 0L) return(NULL)
            
            Xk_cat <- X[, cat_in_cluster, drop = FALSE]
            apply(Xk_cat, 1L, function(row_i) {
                row_i <- as.character(row_i)
                tab <- table(row_i, useNA = "no")
                if (length(tab) == 0L) NA_character_ else names(tab)[which.max(tab)]
            })
        },
        
        # --- Cœurs d'algorithmes ---
        
        run_kmeans = function(X) {
            K <- private$FNbGroupes
            n <- nrow(X)
            p <- ncol(X)
            
            is_num <- vapply(X, is.numeric, logical(1L))
            if (!all(is_num)) {
                stop("run_kmeans() : X doit contenir uniquement des variables quantitatives.")
            }
            
            X_mat <- as.matrix(X)
            
            clusters <- rep(seq_len(K), length.out = p)
            clusters <- sample(clusters)
            
            max_iter    <- 100L
            tol         <- 1e-6
            inertia_old <- Inf
            converged   <- FALSE
            Z_list      <- vector("list", K)
            
            for (iter in seq_len(max_iter)) {
                for (k in seq_len(K)) {
                    vars_k <- which(clusters == k)
                    
                    if (length(vars_k) == 0L) {
                        cand <- sample(seq_len(p), 1L)
                        clusters[cand] <- k
                        vars_k <- cand
                    }
                    
                    Xk <- X_mat[, vars_k, drop = FALSE]
                    
                    if (ncol(Xk) == 1L) {
                        Z_list[[k]] <- as.numeric(Xk[, 1L])
                    } else {
                        pc <- stats::prcomp(Xk, center = FALSE, scale. = FALSE)
                        Z_list[[k]] <- as.numeric(pc$x[, 1L])
                    }
                }
                
                dist_mat <- matrix(NA_real_, nrow = p, ncol = K)
                for (j in seq_len(p)) {
                    xj <- X_mat[, j]
                    for (k in seq_len(K)) {
                        zk <- Z_list[[k]]
                        r  <- suppressWarnings(stats::cor(xj, zk, use = "pairwise.complete.obs"))
                        if (is.na(r)) r <- 0
                        d2 <- 1 - r^2
                        dist_mat[j, k] <- d2
                    }
                }
                
                new_clusters <- apply(dist_mat, 1L, which.min)
                
                min_d2   <- dist_mat[cbind(seq_len(p), new_clusters)]
                inertia  <- sum(min_d2, na.rm = TRUE)
                
                if (all(new_clusters == clusters) || abs(inertia_old - inertia) < tol) {
                    converged   <- TRUE
                    clusters    <- new_clusters
                    inertia_old <- inertia
                    break
                }
                
                clusters    <- new_clusters
                inertia_old <- inertia
            }
            
            res <- list(
                clusters  = clusters,
                centers   = Z_list,
                inertia   = inertia_old,
                converged = converged
            )
            
            return(res)
        },
        
        run_kmodes = function(X) {
            K <- private$FNbGroupes
            n <- nrow(X)
            p <- ncol(X)
            
            is_cat <- vapply(X, is.factor, logical(1L))
            if (!all(is_cat)) {
                stop("run_kmodes() : X doit contenir uniquement des variables qualitatives.")
            }
            
            X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)
            
            clusters <- rep(seq_len(K), length.out = p)
            clusters <- sample(clusters)
            
            max_iter    <- 100L
            tol         <- 1e-6
            inertia_old <- Inf
            converged   <- FALSE
            Z_list      <- vector("list", K)
            
            compute_prototype <- function(Xk_char) {
                apply(Xk_char, 1L, function(row_i) {
                    tab <- table(row_i, useNA = "no")
                    if (length(tab) == 0L) {
                        return(NA_character_)
                    } else {
                        names(tab)[which.max(tab)]
                    }
                })
            }
            
            for (iter in seq_len(max_iter)) {
                for (k in seq_len(K)) {
                    vars_k <- which(clusters == k)
                    
                    if (length(vars_k) == 0L) {
                        cand <- sample(seq_len(p), 1L)
                        clusters[cand] <- k
                        vars_k <- cand
                    }
                    
                    Xk_char <- X_char[, vars_k, drop = FALSE]
                    Z_list[[k]] <- compute_prototype(Xk_char)
                }
                
                dist_mat <- matrix(NA_real_, nrow = p, ncol = K)
                for (j in seq_len(p)) {
                    xj <- X_char[[j]]
                    for (k in seq_len(K)) {
                        zk <- Z_list[[k]]
                        if (length(zk) != n) {
                            stop("run_kmodes() : prototype de longueur incompatible.")
                        }
                        mismatch <- xj != zk
                        d <- mean(mismatch, na.rm = TRUE)
                        if (is.na(d)) d <- 1
                        dist_mat[j, k] <- d
                    }
                }
                
                new_clusters <- apply(dist_mat, 1L, which.min)
                
                min_d <- dist_mat[cbind(seq_len(p), new_clusters)]
                inertia <- sum(min_d, na.rm = TRUE)
                
                if (all(new_clusters == clusters) || abs(inertia_old - inertia) < tol) {
                    converged   <- TRUE
                    clusters    <- new_clusters
                    inertia_old <- inertia
                    break
                }
                
                clusters    <- new_clusters
                inertia_old <- inertia
            }
            
            res <- list(
                clusters  = clusters,
                centers   = Z_list,
                inertia   = inertia_old,
                converged = converged
            )
            
            return(res)
        },
        
        run_kprototypes = function(X) {
            n <- nrow(X)
            p <- ncol(X)
            
            num_idx <- private$FNumCols
            cat_idx <- private$FCatCols
            
            K       <- private$FNbGroupes
            lambda  <- private$FLambda
            
            max_iter <- 100L
            tol      <- 1e-6
            
            if (K < 1L || K > p) {
                stop("K must be between 1 and the number of active variables")
            }
            
            dist_num <- function(j, proto_num) {
                if (is.null(proto_num) || !(j %in% num_idx)) return(Inf)
                xj <- as.numeric(X[[j]])
                r  <- suppressWarnings(cor(xj, proto_num,
                                           use = "pairwise.complete.obs"))
                if (is.na(r)) 1 else 1 - r^2
            }
            
            dist_cat <- function(j, proto_cat) {
                if (is.null(proto_cat) || !(j %in% cat_idx)) return(Inf)
                xj <- as.character(X[[j]])
                mismatch <- xj != proto_cat
                d <- mean(mismatch, na.rm = TRUE)
                if (is.na(d)) d <- 1
                lambda * d
            }
            
            init_idx <- sample.int(p, size = K, replace = FALSE)
            clusters <- rep(NA_integer_, p)
            clusters[init_idx] <- seq_len(K)
            
            prototypes <- vector("list", K)
            for (k in seq_len(K)) {
                vars_k <- which(clusters == k)
                prototypes[[k]] <- list(
                    num = private$compute_num_profile(X, vars_k, num_idx),
                    cat = private$compute_cat_profile(X, vars_k, cat_idx)
                )
            }
            
            for (j in which(is.na(clusters))) {
                d_k <- numeric(K)
                if (j %in% num_idx) {
                    for (k in seq_len(K)) d_k[k] <- dist_num(j, prototypes[[k]]$num)
                } else {
                    for (k in seq_len(K)) d_k[k] <- dist_cat(j, prototypes[[k]]$cat)
                }
                clusters[j] <- which.min(d_k)
            }
            
            inertia   <- Inf
            converged <- FALSE
            iter      <- 0L
            
            while (iter < max_iter && !converged) {
                iter <- iter + 1L
                
                prototypes <- vector("list", K)
                for (k in seq_len(K)) {
                    vars_k <- which(clusters == k)
                    prototypes[[k]] <- list(
                        num = private$compute_num_profile(X, vars_k, num_idx),
                        cat = private$compute_cat_profile(X, vars_k, cat_idx)
                    )
                }
                
                new_clusters <- clusters
                for (j in seq_len(p)) {
                    d_k <- numeric(K)
                    if (j %in% num_idx) {
                        for (k in seq_len(K)) d_k[k] <- dist_num(j, prototypes[[k]]$num)
                    } else {
                        for (k in seq_len(K)) d_k[k] <- dist_cat(j, prototypes[[k]]$cat)
                    }
                    new_clusters[j] <- which.min(d_k)
                }
                
                new_inertia <- 0
                for (j in seq_len(p)) {
                    k <- new_clusters[j]
                    if (j %in% num_idx) {
                        new_inertia <- new_inertia + dist_num(j, prototypes[[k]]$num)
                    } else {
                        new_inertia <- new_inertia + dist_cat(j, prototypes[[k]]$cat)
                    }
                }
                
                if (all(new_clusters == clusters) ||
                    abs(inertia - new_inertia) < tol) {
                    converged <- TRUE
                    clusters  <- new_clusters
                    inertia   <- new_inertia
                    break
                }
                
                clusters <- new_clusters
                inertia  <- new_inertia
            }
            
            res <- list(
                clusters  = clusters,
                centers   = prototypes,
                inertia   = inertia,
                converged = converged
            )
            
            return(res)
        },
        
        run_kmedoids = function(X) {
            # k-medoids sur une matrice de dissimilarités entre variables
            if (!requireNamespace("cluster", quietly = TRUE)) {
                stop("Le package 'cluster' est requis pour la méthode 'kmedoids'.")
            }
            
            K <- private$FNbGroupes
            p <- ncol(X)
            
            is_num <- vapply(X, is.numeric, logical(1L))
            is_cat <- vapply(X, function(col) is.factor(col) || is.character(col), logical(1L))
            
            if (all(is_num)) {
                X_mat <- as.matrix(X)
                cor_mat <- suppressWarnings(stats::cor(X_mat, use = "pairwise.complete.obs"))
                cor_mat[is.na(cor_mat)] <- 0
                D <- 1 - cor_mat^2
            } else if (all(is_cat)) {
                X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)
                D <- matrix(0, nrow = p, ncol = p)
                for (j in seq_len(p)) {
                    for (l in seq_len(p)) {
                        if (j == l) {
                            D[j, l] <- 0
                        } else {
                            mismatch <- X_char[[j]] != X_char[[l]]
                            d <- mean(mismatch, na.rm = TRUE)
                            if (is.na(d)) d <- 1
                            D[j, l] <- d
                        }
                    }
                }
            } else {
                stop("k-medoids n'est actuellement implémenté que pour des ensembles de variables toutes numériques ou toutes qualitatives.")
            }
            
            dimnames(D) <- list(colnames(X), colnames(X))
            D_dist <- stats::as.dist(D)
            
            pam_res <- cluster::pam(D_dist,
                                    k   = K,
                                    diss = TRUE)
            
            clusters <- pam_res$clustering          # longueur p
            medoid_indices <- pam_res$medoids       # indices 1..K
            
            centers <- vector("list", K)
            for (k in seq_len(K)) {
                j_med <- medoid_indices[k]
                centers[[k]] <- X[[j_med]]
            }
            
            inertia <- 0
            for (j in seq_len(p)) {
                k  <- clusters[j]
                jm <- medoid_indices[k]
                inertia <- inertia + D[j, jm]
            }
            
            private$FAlgorithme <- "kmedoids.pam"
            
            res <- list(
                clusters  = clusters,
                centers   = centers,
                inertia   = inertia,
                converged = TRUE
            )
            
            return(res)
        },
        
        # Rattachement d'une seule variable supplémentaire
        predict_one_variable = function(x_new, var_name) {
            method <- private$FMethod
            n <- nrow(private$FX_active)
            K <- private$FNbGroupes
            
            if (length(x_new) != n) {
                stop("predict_one_variable() : la variable '", var_name,
                     "' n'a pas le même nombre d'individus que les données actives.")
            }
            
            # --- K-MEANS ---
            if (method == "kmeans") {
                if (!is.numeric(x_new)) {
                    stop("predict_one_variable() : pour la méthode 'kmeans', la variable '",
                         var_name, "' doit être numérique.")
                }
                
                centers <- private$FCenters
                if (is.null(centers)) {
                    stop("predict_one_variable() : les prototypes ne sont pas disponibles. Appelez fit() d'abord.")
                }
                
                distances <- numeric(K)
                adhesions <- numeric(K)
                
                for (k in seq_len(K)) {
                    zk <- centers[[k]]
                    if (is.null(zk)) {
                        distances[k] <- NA_real_
                        adhesions[k] <- NA_real_
                    } else {
                        r <- suppressWarnings(stats::cor(x_new, zk, use = "pairwise.complete.obs"))
                        if (is.na(r)) r <- 0
                        d2 <- 1 - r^2
                        distances[k] <- d2
                        adhesions[k] <- r^2
                    }
                }
                
                best_k   <- which.min(distances)
                best_d2  <- distances[best_k]
                best_r2  <- adhesions[best_k]
                
                res <- data.frame(
                    variable      = var_name,
                    cluster       = best_k,
                    adhesion      = best_r2,
                    distance_best = best_d2,
                    stringsAsFactors = FALSE
                )
                
                for (k in seq_len(K)) {
                    col_name <- paste0("d", k)
                    res[[col_name]] <- distances[k]
                }
                
                return(res)
            }
            
            # --- K-MODES ---
            if (method == "kmodes") {
                if (!(is.factor(x_new) || is.character(x_new))) {
                    stop("predict_one_variable() : pour la méthode 'kmodes', la variable '",
                         var_name, "' doit être qualitative (factor ou character).")
                }
                x_char <- as.character(x_new)
                
                centers <- private$FCenters
                if (is.null(centers)) {
                    stop("predict_one_variable() : les prototypes ne sont pas disponibles. Appelez fit() d'abord.")
                }
                
                distances <- numeric(K)
                adhesions <- numeric(K)
                
                for (k in seq_len(K)) {
                    zk <- centers[[k]]
                    if (is.null(zk) || length(zk) != n) {
                        distances[k] <- NA_real_
                        adhesions[k] <- NA_real_
                    } else {
                        mismatch <- x_char != zk
                        d <- mean(mismatch, na.rm = TRUE)
                        if (is.na(d)) d <- 1
                        distances[k] <- d
                        adhesions[k] <- 1 - d
                    }
                }
                
                best_k  <- which.min(distances)
                best_d  <- distances[best_k]
                best_ad <- adhesions[best_k]
                
                res <- data.frame(
                    variable      = var_name,
                    cluster       = best_k,
                    adhesion      = best_ad,
                    distance_best = best_d,
                    stringsAsFactors = FALSE
                )
                
                for (k in seq_len(K)) {
                    col_name <- paste0("d", k)
                    res[[col_name]] <- distances[k]
                }
                
                return(res)
            }
            
            # --- K-PROTOTYPES ---
            if (method == "kprototypes") {
                centers    <- private$FCenters
                X_active   <- private$FX_active
                clusters   <- private$FClusters
                num_idx    <- private$FNumCols
                cat_idx    <- private$FCatCols
                lambda     <- private$FLambda
                
                if (is.null(centers) || is.null(clusters)) {
                    stop("predict_one_variable() : les prototypes ne sont pas disponibles. Appelez fit() d'abord.")
                }
                
                if (is.numeric(x_new)) {
                    distances <- numeric(K)
                    adhesions <- numeric(K)
                    
                    for (k in seq_len(K)) {
                        zk_num <- centers[[k]]$num
                        if (is.null(zk_num)) {
                            distances[k] <- NA_real_
                            adhesions[k] <- NA_real_
                        } else {
                            r <- suppressWarnings(stats::cor(x_new, zk_num, use = "pairwise.complete.obs"))
                            if (is.na(r)) r <- 0
                            d2 <- 1 - r^2
                            distances[k] <- d2
                            adhesions[k] <- r^2
                        }
                    }
                    
                    best_k  <- which.min(distances)
                    best_d  <- distances[best_k]
                    best_ad <- adhesions[best_k]
                    
                    res <- data.frame(
                        variable      = var_name,
                        cluster       = best_k,
                        adhesion      = best_ad,
                        distance_best = best_d,
                        stringsAsFactors = FALSE
                    )
                    for (k in seq_len(K)) {
                        col_name <- paste0("d", k)
                        res[[col_name]] <- distances[k]
                    }
                    return(res)
                }
                
                if (!(is.factor(x_new) || is.character(x_new))) {
                    stop("predict_one_variable() : pour 'kprototypes', la variable '",
                         var_name, "' doit être soit numérique, soit qualitative (factor/character).")
                }
                x_char <- as.character(x_new)
                
                distances <- numeric(K)
                adhesions <- numeric(K)
                
                for (k in seq_len(K)) {
                    vars_k <- which(clusters == k)
                    cat_in_cluster <- intersect(vars_k, cat_idx)
                    
                    if (length(cat_in_cluster) == 0L) {
                        distances[k] <- lambda * 1
                        adhesions[k] <- 0
                    } else {
                        Xk_cat <- X_active[, cat_in_cluster, drop = FALSE]
                        zk_cat <- apply(Xk_cat, 1L, function(row_i) {
                            row_i <- as.character(row_i)
                            tab <- table(row_i, useNA = "no")
                            if (length(tab) == 0L) NA_character_ else names(tab)[which.max(tab)]
                        })
                        
                        mismatch <- x_char != zk_cat
                        d_raw <- mean(mismatch, na.rm = TRUE)
                        if (is.na(d_raw)) d_raw <- 1
                        distances[k] <- lambda * d_raw
                        adhesions[k] <- 1 - d_raw
                    }
                }
                
                best_k  <- which.min(distances)
                best_d  <- distances[best_k]
                best_ad <- adhesions[best_k]
                
                res <- data.frame(
                    variable      = var_name,
                    cluster       = best_k,
                    adhesion      = best_ad,
                    distance_best = best_d,
                    stringsAsFactors = FALSE
                )
                for (k in seq_len(K)) {
                    col_name <- paste0("d", k)
                    res[[col_name]] <- distances[k]
                }
                
                return(res)
            }
            
            # --- K-MEDOIDS ---
            if (method == "kmedoids") {
                X_active <- private$FX_active
                centers  <- private$FCenters
                if (is.null(centers)) {
                    stop("predict_one_variable() : les medoids ne sont pas disponibles. Appelez fit() d'abord.")
                }
                
                is_num <- is.numeric(X_active[[1]])
                is_cat <- is.factor(X_active[[1]]) || is.character(X_active[[1]])
                
                distances <- numeric(K)
                adhesions <- numeric(K)
                
                if (is_num && !is_cat) {
                    if (!is.numeric(x_new)) {
                        stop("predict_one_variable() : pour 'kmedoids' numérique, '", var_name,
                             "' doit être numérique.")
                    }
                    for (k in seq_len(K)) {
                        zk <- centers[[k]]
                        r <- suppressWarnings(stats::cor(x_new, zk, use = "pairwise.complete.obs"))
                        if (is.na(r)) r <- 0
                        d <- 1 - r^2
                        distances[k] <- d
                        adhesions[k] <- r^2
                    }
                } else if (is_cat && !is_num) {
                    if (!(is.factor(x_new) || is.character(x_new))) {
                        stop("predict_one_variable() : pour 'kmedoids' qualitatif, '",
                             var_name, "' doit être factor/character.")
                    }
                    x_char <- as.character(x_new)
                    for (k in seq_len(K)) {
                        zk <- as.character(centers[[k]])
                        mismatch <- x_char != zk
                        d <- mean(mismatch, na.rm = TRUE)
                        if (is.na(d)) d <- 1
                        distances[k] <- d
                        adhesions[k] <- 1 - d
                    }
                } else {
                    stop("k-medoids n'est implémenté que pour des ensembles de variables toutes numériques ou toutes qualitatives.")
                }
                
                best_k  <- which.min(distances)
                best_d  <- distances[best_k]
                best_ad <- adhesions[best_k]
                
                res <- data.frame(
                    variable      = var_name,
                    cluster       = best_k,
                    adhesion      = best_ad,
                    distance_best = best_d,
                    stringsAsFactors = FALSE
                )
                for (k in seq_len(K)) {
                    col_name <- paste0("d", k)
                    res[[col_name]] <- distances[k]
                }
                
                return(res)
            }
            
            stop("predict_one_variable() : méthode '", method,
                 "' non implémentée pour la prédiction d'une variable.")
        }
        
    ),
    
    public = list(
        # --- Constructeur ---
        initialize = function(method = c("kmeans", "kmodes", "kprototypes", "kmedoids", "auto"),
                              K      = 2L,
                              scale  = TRUE,
                              lambda = 1) {
            method <- match.arg(method)
            if (K < 2L) stop("K doit être supérieur ou égal à 2.")
            if (lambda <= 0) stop("lambda doit être strictement positif.")
            
            private$FMethod    <- method
            private$FNbGroupes <- as.integer(K)
            private$FScale     <- isTRUE(scale)
            private$FLambda    <- lambda
            
            private$FX_active    <- NULL
            private$FClusters    <- NULL
            private$FCenters     <- NULL
            private$FInertia     <- NULL
            private$FConvergence <- FALSE
            private$FNumCols     <- NULL
            private$FCatCols     <- NULL
            private$FX_new       <- NULL
            private$FAlgorithme  <- NULL
        },
        
        # --- Apprentissage ---
        fit = function(X) {
            X <- private$check_and_prepare_X(X, update_structure = TRUE)
            X <- private$scale_active_variables(X)
            private$FX_active <- X
            
            method <- private$FMethod
            if (method == "auto") {
                has_num <- length(private$FNumCols) > 0L
                has_cat <- length(private$FCatCols) > 0L
                if (has_num && !has_cat) {
                    method <- "kmeans"
                } else if (!has_num && has_cat) {
                    method <- "kmodes"
                } else if (has_num && has_cat) {
                    method <- "kprototypes"
                } else {
                    stop("Aucune variable utilisable pour le clustering.")
                }
                private$FMethod <- method
                
                if (method == "kmeans") {
                    private$FAlgorithme <- "kmeans.correlation"
                } else if (method == "kmodes") {
                    private$FAlgorithme <- "kmodes.simple-matching"
                } else if (method == "kprototypes") {
                    private$FAlgorithme <- "kprototypes.corr+lambda"
                }
            }
            
            res <- switch(
                method,
                "kmeans"      = private$run_kmeans(X),
                "kmodes"      = private$run_kmodes(X),
                "kprototypes" = private$run_kprototypes(X),
                "kmedoids"    = private$run_kmedoids(X),
                stop("Méthode non supportée : ", method)
            )
            
            private$FClusters    <- res$clusters
            private$FCenters     <- res$centers
            private$FInertia     <- res$inertia
            private$FConvergence <- isTRUE(res$converged)
            
            invisible(self)
        },
        
        # --- Rattachement des variables supplémentaires ---
        predict = function(X_new) {
            if (is.null(private$FX_active)) {
                stop("Le modèle n'a pas encore été appris. Appelez fit() d'abord.")
            }
            if (!is.data.frame(X_new)) {
                stop("X_new doit être un data.frame.")
            }
            if (nrow(X_new) != nrow(private$FX_active)) {
                stop("X_new doit avoir le même nombre d'individus (lignes) que les données actives.")
            }
            
            X_new <- private$check_and_prepare_X(X_new, update_structure = FALSE)
            
            res_list <- vector("list", length = ncol(X_new))
            for (j in seq_len(ncol(X_new))) {
                res_list[[j]] <- private$predict_one_variable(
                    x_new    = X_new[[j]],
                    var_name = colnames(X_new)[j]
                )
            }
            
            private$FX_new <- X_new
            
            res <- do.call(rbind, res_list)
            rownames(res) <- NULL
            return(res)
        },
        
        # --- Print succinct ---
        print = function(...) {
            cat("Classe 'mmrClustVar'\n")
            cat("  Méthode        :", private$FMethod, "\n")
            cat("  K              :", private$FNbGroupes, "\n")
            if (!is.null(private$FX_active)) {
                cat("  Nb variables   :", ncol(private$FX_active), "\n")
                cat("  Nb individus   :", nrow(private$FX_active), "\n")
            }
            if (!is.null(private$FInertia)) {
                cat("  Inertie intra  :", format(private$FInertia, digits = 4), "\n")
            }
            cat("  Convergence    :", private$FConvergence, "\n")
            invisible(self)
        },
        
        # --- Résumé détaillé ---
        summary = function(...) {
            if (is.null(private$FClusters) || is.null(private$FX_active)) {
                cat("Aucun modèle appris. Appelez fit() d'abord.\n")
                return(invisible(NULL))
            }
            
            method   <- private$FMethod
            K        <- private$FNbGroupes
            clusters <- private$FClusters
            X        <- private$FX_active
            centers  <- private$FCenters
            
            cat("Résumé du modèle 'mmrClustVar'\n")
            cat("  Méthode        :", method, "\n")
            cat("  K              :", K, "\n")
            cat("  Nb variables   :", length(clusters), "\n")
            cat("  Inertie intra  :", format(private$FInertia, digits = 4), "\n")
            cat("  Convergence    :", private$FConvergence, "\n\n")
            
            p <- length(clusters)
            membership <- rep(NA_real_, p)
            
            if (method == "kmeans") {
                X_mat <- as.matrix(X)
                for (j in seq_len(p)) {
                    kj <- clusters[j]
                    zk <- centers[[kj]]
                    xj <- X_mat[, j]
                    r  <- suppressWarnings(stats::cor(xj, zk, use = "pairwise.complete.obs"))
                    if (is.na(r)) r <- 0
                    membership[j] <- r^2
                }
                membership_label <- "r^2 (correlation avec la composante latente du cluster)"
            } else if (method == "kmodes") {
                X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)
                for (j in seq_len(p)) {
                    kj <- clusters[j]
                    xj <- X_char[[j]]
                    zk <- centers[[kj]]
                    mismatch <- xj != zk
                    d <- mean(mismatch, na.rm = TRUE)
                    if (is.na(d)) d <- 1
                    membership[j] <- 1 - d
                }
                membership_label <- "1 - dissimilarité (simple matching) avec le mode du cluster"
            } else if (method == "kprototypes") {
                num_idx <- private$FNumCols
                cat_idx <- private$FCatCols
                X_used  <- X
                
                for (j in seq_len(p)) {
                    kj      <- clusters[j]
                    proto_k <- centers[[kj]]
                    
                    if (j %in% num_idx) {
                        xj     <- X_used[[j]]
                        zk_num <- proto_k$num
                        if (is.null(zk_num)) {
                            membership[j] <- NA_real_
                        } else {
                            r <- suppressWarnings(
                                stats::cor(xj, zk_num, use = "pairwise.complete.obs")
                            )
                            if (is.na(r)) r <- 0
                            membership[j] <- r^2
                        }
                    } else if (j %in% cat_idx) {
                        xj_char <- as.character(X_used[[j]])
                        zk_cat  <- proto_k$cat
                        
                        if (is.null(zk_cat) || length(zk_cat) != length(xj_char)) {
                            membership[j] <- 0
                        } else {
                            mismatch <- xj_char != zk_cat
                            d_raw <- mean(mismatch, na.rm = TRUE)
                            if (is.na(d_raw)) d_raw <- 1
                            membership[j] <- 1 - d_raw
                        }
                    } else {
                        membership[j] <- NA_real_
                    }
                }
                
                membership_label <- "Score d'adhésion (r^2 pour les num., 1 - dissimilarité pour les cat.)"
            } else if (method == "kmedoids") {
                X_active <- X
                is_num <- is.numeric(X_active[[1]])
                is_cat <- is.factor(X_active[[1]]) || is.character(X_active[[1]])
                
                if (is_num && !is_cat) {
                    X_mat <- as.matrix(X_active)
                    for (j in seq_len(p)) {
                        kj <- clusters[j]
                        zk <- centers[[kj]]
                        xj <- X_mat[, j]
                        r  <- suppressWarnings(stats::cor(xj, zk, use = "pairwise.complete.obs"))
                        if (is.na(r)) r <- 0
                        membership[j] <- r^2
                    }
                    membership_label <- "r^2 (corrélation avec le medoid du cluster)"
                } else if (is_cat && !is_num) {
                    X_char <- as.data.frame(lapply(X_active, as.character), stringsAsFactors = FALSE)
                    for (j in seq_len(p)) {
                        kj <- clusters[j]
                        xj <- X_char[[j]]
                        zk <- as.character(centers[[kj]])
                        mismatch <- xj != zk
                        d <- mean(mismatch, na.rm = TRUE)
                        if (is.na(d)) d <- 1
                        membership[j] <- 1 - d
                    }
                    membership_label <- "1 - dissimilarité (simple matching) avec le medoid du cluster"
                } else {
                    membership_label <- "Adhésion non définie (données mixtes non supportées par k-medoids)."
                }
            } else {
                cat("Méthode non reconnue pour summary().\n")
                return(invisible(NULL))
            }
            
            tab_size <- table(clusters)
            cl_ids   <- as.integer(names(tab_size))
            mean_mem <- tapply(membership, clusters, mean, na.rm = TRUE)
            
            cat("Résumé par cluster :\n")
            cl_df <- data.frame(
                cluster          = cl_ids,
                taille           = as.integer(tab_size[as.character(cl_ids)]),
                adhesion_moyenne = as.numeric(mean_mem[as.character(cl_ids)])
            )
            print(cl_df, row.names = FALSE)
            cat("\n")
            
            var_df <- data.frame(
                variable  = colnames(X),
                cluster   = clusters,
                adhesion  = membership,
                stringsAsFactors = FALSE
            )
            
            cat("Indicateur d'adhésion (", membership_label, ")\n", sep = "")
            ord <- order(var_df$cluster, -var_df$adhesion)
            print(utils::head(var_df[ord, ], n = min(10L, nrow(var_df))), row.names = FALSE)
            cat("\n(Le data.frame complet est renvoyé invisiblement.)\n")
            
            invisible(var_df)
        },
        
        # --- Graphiques ---
        plot = function(type = c("inertia", "clusters", "membership", "profiles"),
                        Ks = NULL, ...) {
            type <- match.arg(type)
            
            if (is.null(private$FX_active) || is.null(private$FClusters)) {
                stop("Aucun modèle appris. Appelez fit() d'abord.")
            }
            
            X        <- private$FX_active
            method   <- private$FMethod
            clusters <- private$FClusters
            centers  <- private$FCenters
            K        <- private$FNbGroupes
            
            if (type == "clusters") {
                tab <- table(clusters)
                graphics::barplot(
                    tab,
                    xlab = "Cluster",
                    ylab = "Nombre de variables",
                    main = "Taille des clusters de variables"
                )
                return(invisible(NULL))
            }
            
            if (type == "membership") {
                p <- length(clusters)
                membership <- rep(NA_real_, p)
                ylab <- ""
                
                if (method == "kmeans") {
                    X_mat <- as.matrix(X)
                    for (j in seq_len(p)) {
                        kj <- clusters[j]
                        zk <- centers[[kj]]
                        xj <- X_mat[, j]
                        r  <- suppressWarnings(stats::cor(xj, zk, use = "pairwise.complete.obs"))
                        if (is.na(r)) r <- 0
                        membership[j] <- r^2
                    }
                    ylab <- "r^2 (corrélation avec la composante latente)"
                } else if (method == "kmodes") {
                    X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)
                    for (j in seq_len(p)) {
                        kj <- clusters[j]
                        xj <- X_char[[j]]
                        zk <- centers[[kj]]
                        mismatch <- xj != zk
                        d <- mean(mismatch, na.rm = TRUE)
                        if (is.na(d)) d <- 1
                        membership[j] <- 1 - d
                    }
                    ylab <- "1 - dissimilarité (simple matching)"
                } else if (method == "kprototypes") {
                    num_idx <- private$FNumCols
                    cat_idx <- private$FCatCols
                    X_used  <- X
                    
                    for (j in seq_len(p)) {
                        kj      <- clusters[j]
                        proto_k <- centers[[kj]]
                        
                        if (j %in% num_idx) {
                            xj     <- X_used[[j]]
                            zk_num <- proto_k$num
                            r <- suppressWarnings(
                                stats::cor(xj, zk_num, use = "pairwise.complete.obs")
                            )
                            if (is.na(r)) r <- 0
                            membership[j] <- r^2
                        } else if (j %in% cat_idx) {
                            xj_char <- as.character(X_used[[j]])
                            zk_cat  <- proto_k$cat
                            
                            if (is.null(zk_cat) || length(zk_cat) != length(xj_char)) {
                                membership[j] <- 0
                            } else {
                                mismatch <- xj_char != zk_cat
                                d_raw <- mean(mismatch, na.rm = TRUE)
                                if (is.na(d_raw)) d_raw <- 1
                                membership[j] <- 1 - d_raw
                            }
                        } else {
                            membership[j] <- NA_real_
                        }
                    }
                    ylab <- "Score d'adhésion (r^2 num., 1 - dissimilarité cat.)"
                } else if (method == "kmedoids") {
                    X_active <- X
                    is_num <- is.numeric(X_active[[1]])
                    is_cat <- is.factor(X_active[[1]]) || is.character(X_active[[1]])
                    
                    if (is_num && !is_cat) {
                        X_mat <- as.matrix(X_active)
                        for (j in seq_len(p)) {
                            kj <- clusters[j]
                            zk <- centers[[kj]]
                            xj <- X_mat[, j]
                            r  <- suppressWarnings(stats::cor(xj, zk, use = "pairwise.complete.obs"))
                            if (is.na(r)) r <- 0
                            membership[j] <- r^2
                        }
                        ylab <- "r^2 (corrélation avec le medoid)"
                    } else if (is_cat && !is_num) {
                        X_char <- as.data.frame(lapply(X_active, as.character), stringsAsFactors = FALSE)
                        for (j in seq_len(p)) {
                            kj <- clusters[j]
                            xj <- X_char[[j]]
                            zk <- as.character(centers[[kj]])
                            mismatch <- xj != zk
                            d <- mean(mismatch, na.rm = TRUE)
                            if (is.na(d)) d <- 1
                            membership[j] <- 1 - d
                        }
                        ylab <- "1 - dissimilarité (simple matching) avec le medoid"
                    } else {
                        stop("k-medoids n'est implémenté que pour des ensembles de variables toutes numériques ou toutes qualitatives.")
                    }
                } else {
                    stop("Méthode non reconnue pour plot(type = 'membership').")
                }
                
                ord <- order(clusters, -membership)
                graphics::barplot(
                    membership[ord],
                    names.arg = colnames(X)[ord],
                    las = 2,
                    cex.names = 0.6,
                    xlab = "Variables (triées par cluster et adhésion)",
                    ylab = ylab,
                    main = "Degré d'adhésion des variables à leur cluster",
                    ...
                )
                return(invisible(NULL))
            }
            
            if (type == "inertia") {
                p <- ncol(X)
                if (is.null(Ks)) {
                    Ks <- seq_len(min(10L, p))
                }
                Ks <- Ks[Ks >= 2]
                if (length(Ks) == 0L) {
                    stop("Impossible de tracer la courbe d'inertie : K doit être >= 2.")
                }
                
                old_K    <- private$FNbGroupes
                inertias <- numeric(length(Ks))
                
                for (i in seq_along(Ks)) {
                    private$FNbGroupes <- as.integer(Ks[i])
                    res_i <- switch(
                        method,
                        "kmeans"      = private$run_kmeans(X),
                        "kmodes"      = private$run_kmodes(X),
                        "kprototypes" = private$run_kprototypes(X),
                        "kmedoids"    = private$run_kmedoids(X),
                        stop("Méthode non supportée pour le calcul de l'inertie.")
                    )
                    inertias[i] <- res_i$inertia
                }
                
                private$FNbGroupes <- old_K
                
                graphics::plot(
                    Ks, inertias, type = "b",
                    xlab = "K (nombre de clusters)",
                    ylab = "Inertie intra-cluster",
                    main = "Courbe de l'inertie en fonction de K",
                    ...
                )
                return(invisible(NULL))
            }
            
            if (type == "profiles") {
                n <- nrow(X)
                
                if (method == "kmeans") {
                    Z_mat <- matrix(NA_real_, nrow = n, ncol = K)
                    for (k in seq_len(K)) {
                        zk <- centers[[k]]
                        if (!is.null(zk) && length(zk) == n) {
                            Z_mat[, k] <- zk
                        }
                    }
                    graphics::matplot(
                        seq_len(n), Z_mat, type = "l", lty = 1,
                        xlab = "Individus",
                        ylab = "Score de composante latente",
                        main = "Profils moyens par cluster (composantes latentes)",
                        ...
                    )
                    graphics::legend(
                        "topright",
                        legend = paste("Cluster", seq_len(K)),
                        lty = 1,
                        col = seq_len(K),
                        bty = "n",
                        cex = 0.8
                    )
                    return(invisible(NULL))
                }
                
                if (method == "kmodes") {
                    old_par <- graphics::par(no.readonly = TRUE)
                    on.exit(graphics::par(old_par))
                    
                    n_panel <- min(K, 3L)
                    graphics::par(mfrow = c(1, n_panel))
                    
                    for (k in seq_len(K)) {
                        zk <- centers[[k]]
                        if (is.null(zk)) next
                        
                        tab <- sort(table(zk), decreasing = TRUE)
                        tab <- head(tab, 6L)
                        
                        graphics::barplot(
                            tab,
                            main = paste("Cluster", k),
                            ylab = "Effectif du prototype",
                            las  = 2,
                            ...
                        )
                    }
                    return(invisible(NULL))
                }
                
                if (method == "kprototypes") {
                    Z_mat <- matrix(NA_real_, nrow = n, ncol = K)
                    for (k in seq_len(K)) {
                        proto_k <- centers[[k]]
                        zk_num  <- proto_k$num
                        if (!is.null(zk_num) && length(zk_num) == n) {
                            Z_mat[, k] <- zk_num
                        }
                    }
                    graphics::matplot(
                        seq_len(n), Z_mat, type = "l", lty = 1,
                        xlab = "Individus",
                        ylab = "Profil numérique prototype",
                        main = "Profils numériques des prototypes (k-prototypes)",
                        ...
                    )
                    graphics::legend(
                        "topright",
                        legend = paste("Cluster", seq_len(K)),
                        lty = 1,
                        col = seq_len(K),
                        bty = "n",
                        cex = 0.8
                    )
                    return(invisible(NULL))
                }
                
                if (method == "kmedoids") {
                    Z_mat <- matrix(NA_real_, nrow = n, ncol = K)
                    for (k in seq_len(K)) {
                        zk <- centers[[k]]
                        if (!is.null(zk) && length(zk) == n && is.numeric(zk)) {
                            Z_mat[, k] <- zk
                        }
                    }
                    graphics::matplot(
                        seq_len(n), Z_mat, type = "l", lty = 1,
                        xlab = "Individus",
                        ylab = "Profil du medoid (numérique)",
                        main = "Profils des medoids par cluster (k-medoids)",
                        ...
                    )
                    graphics::legend(
                        "topright",
                        legend = paste("Cluster", seq_len(K)),
                        lty = 1,
                        col = seq_len(K),
                        bty = "n",
                        cex = 0.8
                    )
                    return(invisible(NULL))
                }
                
                stop("Méthode non reconnue pour plot(type = 'profiles').")
            }
        },
        
        # --- Getters publics ---
        get_clusters = function() {
            private$FClusters
        },
        
        get_inertia = function() {
            private$FInertia
        },
        
        get_method = function() {
            private$FMethod
        },
        
        get_centers = function() {
            private$FCenters
        },
        
        get_convergence = function() {
            private$FConvergence
        },
        
        get_X_new = function() {
            private$FX_new
        },
        
        get_algorithm = function() {
            private$FAlgorithme
        }
    )
)
