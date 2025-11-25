#' mmrClustVar: R6 Class for Clustering of Variables
#'
#' @description
#' An R6 class implementing clustering of variables using:
#' \itemize{
#'   \item \strong{k-means} for numeric variables (correlation-based distance),
#'   \item \strong{k-modes} for categorical variables (simple matching dissimilarity),
#'   \item \strong{k-prototypes} for mixed data (weighted combination of both).
#' }
#'
#' Variables are treated as objects: each variable is represented by its
#' vector of values across individuals. The algorithm partitions the \eqn{p}
#' variables into \eqn{K} clusters and, for each cluster, computes a latent
#' component or a categorical prototype used to measure similarity.
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
#'         or "kprototypes").}
#'   \item{\code{$get_centers()}}{Return the list of latent components or
#'         prototypes.}
#'   \item{\code{$get_convergence()}}{Return a logical indicating whether the
#'         algorithm converged.}
#'   \item{\code{$get_X_descr()}}{Return the last set of supplementary variables
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
#' df_descr <- data.frame(
#'   Bonus = rnorm(80),
#'   Label = factor(sample(c("yes", "no"), 80, replace = TRUE))
#' )
#' res_pred_mixed <- obj_kprot$predict(df_descr)
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
        FMethod      = NULL,  # "k-means", "k-modes", "k-prototypes", "auto"
        FNbGroupes   = NULL,  # K
        FScale       = NULL,  # TRUE / FALSE (standardisation des variables quantitatives actives)
        FLambda      = NULL,  # pondération partie catégorielle (k-prototypes)

        FX_active    = NULL,  # data.frame des variables actives (pré-traitées)
        FClusters    = NULL,  # vecteur d'affectation des variables actives
        FCenters     = NULL,  # centres / modes / prototypes selon la méthode
        FInertia     = NULL,  # inertie intra-cluster totale
        FConvergence = NULL,  # booléen
        FX_descr     = NULL,  # variables descriptives fournies à predict()
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

            # Si on ne met pas à jour la structre (cas predict()) -> au moins 1 variable descriptive
            if (!update_structure && ncol(X) < 1L) {
                stop("X_new doit contenir au moins 1 variable descriptive.")
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
                    # Variable constante (on la laisse telle quelle mais on pourra l'exclure plus tard)
                    X[[j]] <- col
                } else {
                    X[[j]] <- (col - mu) / sdv
                }
            }
            return(X)
        },

        compute_num_profile = function(X, vars_k, num_idx) {
            # vars_k : indices de colonnes dans le cluster k
            num_in_cluster <- intersect(vars_k, num_idx)
            if (length(num_in_cluster) == 0L) return(NULL)

            Xk_num <- as.matrix(X[, num_in_cluster, drop = FALSE])
            rowMeans(Xk_num, na.rm = TRUE)  # profil numérique : moyenne par individu
        },

        compute_cat_profile = function(X, vars_k, cat_idx) {
            # vars_k : indices de colonnes dans le cluster k
            cat_in_cluster <- intersect(vars_k, cat_idx)
            if (length(cat_in_cluster) == 0L) return(NULL)

            Xk_cat <- X[, cat_in_cluster, drop = FALSE]
            # profil catégoriel : mode par individu
            apply(Xk_cat, 1L, function(row_i) {
                row_i <- as.character(row_i)
                tab <- table(row_i, useNA = "no")
                if (length(tab) == 0L) NA_character_ else names(tab)[which.max(tab)]
            })
        },

        # --- Cœurs d'algorithmes ---

        run_kmeans = function(X) {
            # X : data.frame ou matrice n x p, uniquement des variables quantitatives
            K <- private$FNbGroupes
            n <- nrow(X)
            p <- ncol(X)

            # Sécurité : vérifier que tout est numérique
            is_num <- vapply(X, is.numeric, logical(1L))
            if (!all(is_num)) {
                stop("run_kmeans() : X doit contenir uniquement des variables quantitatives.")
            }

            X_mat <- as.matrix(X)

            # Initialisation des clusters (partition des variables)
            # On essaie de répartir à peu près équitablement
            clusters <- rep(seq_len(K), length.out = p)
            clusters <- sample(clusters)  # randomisation

            max_iter    <- 100L
            tol         <- 1e-6
            inertia_old <- Inf
            converged   <- FALSE
            Z_list      <- vector("list", K)  # composantes latentes

            for (iter in seq_len(max_iter)) {
                # 1) Calcul des composantes latentes Z_k pour chaque cluster
                for (k in seq_len(K)) {
                    vars_k <- which(clusters == k)

                    if (length(vars_k) == 0L) {
                        # Cluster vide : on choisit une variable au hasard et on la force dans ce cluster
                        cand <- sample(seq_len(p), 1L)
                        clusters[cand] <- k
                        vars_k <- cand
                    }

                    Xk <- X_mat[, vars_k, drop = FALSE]

                    if (ncol(Xk) == 1L) {
                        # Une seule variable : la composante latente est le profil de cette variable
                        Z_list[[k]] <- as.numeric(Xk[, 1L])
                    } else {
                        # Plusieurs variables : on prend la 1ère composante principale
                        # (scores sur les individus)
                        pc <- stats::prcomp(Xk, center = FALSE, scale. = FALSE)
                        Z_list[[k]] <- as.numeric(pc$x[, 1L])
                    }
                }

                # 2) Calcul des distances d^2(j,k) = 1 - r^2(X_j, Z_k)
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

                # 3) Nouvelle affectation : pour chaque variable, le cluster de distance minimale
                new_clusters <- apply(dist_mat, 1L, which.min)

                # Inertie intra-cluster = somme des distances au meilleur cluster
                min_d2   <- dist_mat[cbind(seq_len(p), new_clusters)]
                inertia  <- sum(min_d2, na.rm = TRUE)

                # 4) Critère de convergence
                if (all(new_clusters == clusters) || abs(inertia_old - inertia) < tol) {
                    converged   <- TRUE
                    clusters    <- new_clusters
                    inertia_old <- inertia
                    break
                }

                clusters    <- new_clusters
                inertia_old <- inertia
            }

            # Résultat : on retourne la partition et les composantes latentes
            res <- list(
                clusters = clusters,       # vecteur de longueur p
                centers  = Z_list,         # liste de K vecteurs (longueur n)
                inertia  = inertia_old,    # inertie intra totale
                converged = converged
            )

            return(res)
        },

        run_kmodes = function(X) {
            # X : data.frame n x p, uniquement des variables qualitatives
            K <- private$FNbGroupes
            n <- nrow(X)
            p <- ncol(X)

            # Sécurité : toutes les colonnes doivent être des facteurs
            is_cat <- vapply(X, is.factor, logical(1L))
            if (!all(is_cat)) {
                stop("run_kmodes() : X doit contenir uniquement des variables qualitatives.")
            }

            # On travaille en caractère pour simplifier les comparaisons
            X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)

            # Initialisation des clusters (partition des variables)
            clusters <- rep(seq_len(K), length.out = p)
            clusters <- sample(clusters)  # randomisation

            max_iter    <- 100L
            tol         <- 1e-6
            inertia_old <- Inf
            converged   <- FALSE
            Z_list      <- vector("list", K)  # prototypes catégoriels (vecteurs de longueur n)

            # Fonction locale pour calculer le prototype d'un cluster
            compute_prototype <- function(Xk_char) {
                # Xk_char : data.frame n x p_k (caractères)
                # Pour chaque individu i, on prend la modalité la plus fréquente sur les variables du cluster
                apply(Xk_char, 1L, function(row_i) {
                    tab <- table(row_i, useNA = "no")
                    if (length(tab) == 0L) {
                        return(NA_character_)
                    } else {
                        # mode (on prend la première modalité en cas d'égalité)
                        names(tab)[which.max(tab)]
                    }
                })
            }

            for (iter in seq_len(max_iter)) {
                # 1) Calcul des prototypes Z_k pour chaque cluster
                for (k in seq_len(K)) {
                    vars_k <- which(clusters == k)

                    if (length(vars_k) == 0L) {
                        # Cluster vide : on choisit une variable au hasard et on la force dans ce cluster
                        cand <- sample(seq_len(p), 1L)
                        clusters[cand] <- k
                        vars_k <- cand
                    }

                    Xk_char <- X_char[, vars_k, drop = FALSE]
                    Z_list[[k]] <- compute_prototype(Xk_char)  # vecteur de longueur n
                }

                # 2) Calcul des dissimilarités d(j,k) pour chaque variable / cluster
                dist_mat <- matrix(NA_real_, nrow = p, ncol = K)
                for (j in seq_len(p)) {
                    xj <- X_char[[j]]  # vecteur de longueur n
                    for (k in seq_len(K)) {
                        zk <- Z_list[[k]]
                        if (length(zk) != n) {
                            stop("run_kmodes() : prototype de longueur incompatible.")
                        }
                        mismatch <- xj != zk
                        d <- mean(mismatch, na.rm = TRUE)
                        if (is.na(d)) d <- 1  # si tout est NA, on prend la dissimilarité maximale
                        dist_mat[j, k] <- d
                    }
                }

                # 3) Nouvelle affectation : chaque variable va au cluster de dissimilarité minimale
                new_clusters <- apply(dist_mat, 1L, which.min)

                # Inertie intra-cluster = somme des dissimilarités au meilleur cluster
                min_d <- dist_mat[cbind(seq_len(p), new_clusters)]
                inertia <- sum(min_d, na.rm = TRUE)

                # 4) Critère de convergence
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
                clusters  = clusters,   # vecteur de longueur p : cluster de chaque variable
                centers   = Z_list,     # liste de K prototypes (vecteurs de longueur n)
                inertia   = inertia_old,
                converged = converged
            )

            return(res)
        },

        run_kprototypes = function(X) {
            # X : data.frame des variables actives
            n <- nrow(X)
            p <- ncol(X)

            num_idx <- private$FNumCols   # indices des variables quantitatives
            cat_idx <- private$FCatCols   # indices des variables qualitatives

            K       <- private$FNbGroupes
            lambda  <- private$FLambda

            max_iter <- 100L
            tol      <- 1e-6

            if (K < 1L || K > p) {
                stop("K must be between 1 and the number of active variables")
            }

            # --- helpers distances ---

            dist_num <- function(j, proto_num) {
                # j : indice global de la variable
                if (is.null(proto_num) || !(j %in% num_idx)) return(Inf)
                xj <- as.numeric(X[[j]])
                r  <- suppressWarnings(cor(xj, proto_num,
                                           use = "pairwise.complete.obs"))
                if (is.na(r)) 1 else 1 - r^2
            }

            dist_cat <- function(j, proto_cat) {
                # j : indice global de la variable
                if (is.null(proto_cat) || !(j %in% cat_idx)) return(Inf)
                xj <- as.character(X[[j]])
                mismatch <- xj != proto_cat
                d <- mean(mismatch, na.rm = TRUE)
                if (is.na(d)) d <- 1
                lambda * d
            }

            # --- initialisation des clusters ---

            # on choisit K variables différentes comme noyaux initiaux
            init_idx <- sample.int(p, size = K, replace = FALSE)
            clusters <- rep(NA_integer_, p)
            clusters[init_idx] <- seq_len(K)

            # prototypes initiaux (à partir des noyaux uniquement)
            prototypes <- vector("list", K)
            for (k in seq_len(K)) {
                vars_k <- which(clusters == k)
                prototypes[[k]] <- list(
                    num = private$compute_num_profile(X, vars_k, num_idx),
                    cat = private$compute_cat_profile(X, vars_k, cat_idx)
                )
            }

            # affectation initiale des variables non noyaux
            for (j in which(is.na(clusters))) {
                d_k <- numeric(K)
                if (j %in% num_idx) {
                    for (k in seq_len(K)) d_k[k] <- dist_num(j, prototypes[[k]]$num)
                } else {
                    for (k in seq_len(K)) d_k[k] <- dist_cat(j, prototypes[[k]]$cat)
                }
                clusters[j] <- which.min(d_k)
            }

            # --- boucle de réallocation ---

            inertia <- Inf
            converged <- FALSE
            iter <- 0L

            while (iter < max_iter && !converged) {
                iter <- iter + 1L

                # 1) recalcul des prototypes pour chaque cluster
                prototypes <- vector("list", K)
                for (k in seq_len(K)) {
                    vars_k <- which(clusters == k)
                    prototypes[[k]] <- list(
                        num = private$compute_num_profile(X, vars_k, num_idx),
                        cat = private$compute_cat_profile(X, vars_k, cat_idx)
                    )
                }

                # 2) réaffectation de toutes les variables
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

                # 3) calcul de l'inertie intra-cluster
                new_inertia <- 0
                for (j in seq_len(p)) {
                    k <- new_clusters[j]
                    if (j %in% num_idx) {
                        new_inertia <- new_inertia + dist_num(j, prototypes[[k]]$num)
                    } else {
                        new_inertia <- new_inertia + dist_cat(j, prototypes[[k]]$cat)
                    }
                }

                # 4) test de convergence
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
                clusters   = clusters,
                centers = prototypes,
                inertia    = inertia,
                converged  = converged
            )

            return(res)
        },

        # Rattachement d'une seule variable descriptive
        predict_one_variable = function(x_new, var_name) {
            # x_new : vecteur (numérique ou factor/character)
            # var_name : nom de la variable descriptive

            method <- private$FMethod
            n <- nrow(private$FX_active)
            K <- private$FNbGroupes

            if (length(x_new) != n) {
                stop("predict_one_variable() : la variable '", var_name,
                     "' n'a pas le même nombre d'individus que les données actives.")
            }

            # --- MÉTHODE K-MEANS : variable numérique, distance 1 - r^2 ---
            if (method == "kmeans") {
                if (!is.numeric(x_new)) {
                    stop("predict_one_variable() : pour la méthode 'kmeans', la variable '",
                         var_name, "' doit être numérique.")
                }

                centers <- private$FCenters  # liste de K composantes latentes Z_k
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

            # --- MÉTHODE K-MODES : variable qualitative, simple matching ---
            if (method == "kmodes") {
                # On accepte factor ou character
                if (!(is.factor(x_new) || is.character(x_new))) {
                    stop("predict_one_variable() : pour la méthode 'kmodes', la variable '",
                         var_name, "' doit être qualitative (factor ou character).")
                }
                x_char <- as.character(x_new)

                centers <- private$FCenters  # liste de K prototypes catégoriels (vecteurs de longueur n)
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
                        d <- mean(mismatch, na.rm = TRUE) # simple matching dissimilarity
                        if (is.na(d)) d <- 1
                        distances[k] <- d
                        adhesions[k] <- 1 - d             # proportion de matches
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

            # --- MÉTHODE K-PROTOTYPES : variable num OU quali, distance mixte ---
            if (method == "kprototypes") {
                centers    <- private$FCenters   # liste de K prototypes mixtes (num + cat)
                X_active   <- private$FX_active
                clusters   <- private$FClusters
                num_idx    <- private$FNumCols
                cat_idx    <- private$FCatCols
                lambda     <- private$FLambda

                if (is.null(centers) || is.null(clusters)) {
                    stop("predict_one_variable() : les prototypes ne sont pas disponibles. Appelez fit() d'abord.")
                }

                # --- Cas numérique : même logique que k-means, avec prototypes$num ---
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

                # --- Cas qualitatif : simple matching sur profil catégoriel par cluster ---
                if (!(is.factor(x_new) || is.character(x_new))) {
                    stop("predict_one_variable() : pour 'kprototypes', la variable '",
                         var_name, "' doit être soit numérique, soit qualitative (factor/character).")
                }
                x_char <- as.character(x_new)

                distances <- numeric(K)
                adhesions <- numeric(K)

                for (k in seq_len(K)) {
                    vars_k <- which(clusters == k)               # variables actives dans le cluster k
                    cat_in_cluster <- intersect(vars_k, cat_idx) # indices globaux des quali dans ce cluster

                    if (length(cat_in_cluster) == 0L) {
                        # Aucun support qualitatif dans ce cluster → dissimilarité max
                        distances[k] <- lambda * 1
                        adhesions[k] <- 0
                    } else {
                        Xk_cat <- X_active[, cat_in_cluster, drop = FALSE]
                        # Profil catégoriel par individu : modalité la plus fréquente dans le cluster
                        zk_cat <- apply(Xk_cat, 1L, function(row_i) {
                            row_i <- as.character(row_i)
                            tab <- table(row_i, useNA = "no")
                            if (length(tab) == 0L) NA_character_ else names(tab)[which.max(tab)]
                        })

                        mismatch <- x_char != zk_cat
                        d_raw <- mean(mismatch, na.rm = TRUE)  # simple matching dans [0,1]
                        if (is.na(d_raw)) d_raw <- 1
                        distances[k] <- lambda * d_raw         # même échelle que dans l'algo
                        adhesions[k] <- 1 - d_raw              # proportion de matches, non pondérée
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

            # --- Méthode non gérée ---
            stop("predict_one_variable() : méthode '", method,
                 "' non implémentée pour la prédiction d'une variable.")
        }

    ),

    public = list(
        # --- Constructeur ---
        initialize = function(method = c("kmeans", "kmodes", "kprototypes", "auto"),
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
                stop("Méthode non supportée : ", method)
            )

            private$FClusters    <- res$clusters
            private$FCenters     <- res$centers
            private$FInertia     <- res$inertia
            private$FConvergence <- isTRUE(res$converged)

            invisible(self)
        },


        # --- Rattachement des variables descriptives ---
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

            # Préparation des données descriptives : on convertit les qualitatives en facteurs
            X_new <- private$check_and_prepare_X(X_new, update_structure = FALSE)

            res_list <- vector("list", length = ncol(X_new))
            for (j in seq_len(ncol(X_new))) {
                res_list[[j]] <- private$predict_one_variable(
                    x_new    = X_new[[j]],
                    var_name = colnames(X_new)[j]
                )
            }

            private$FX_descr <- X_new # mémorisation des variables descriptives

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

            # --- 1) Degré d'adhésion par variable ---
            p <- length(clusters)
            membership <- rep(NA_real_, p)

            if (method == "kmeans") {
                # r^2(X_j, Z_k(j))
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
                # 1 - dissimilarité simple matching
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
                # r^2 pour les variables numériques,
                # 1 - dissimilarité simple matching pour les variables qualitatives
                num_idx <- private$FNumCols
                cat_idx <- private$FCatCols
                X_used  <- X

                for (j in seq_len(p)) {
                    kj      <- clusters[j]
                    proto_k <- centers[[kj]]

                    if (j %in% num_idx) {
                        # Partie numérique : r^2(X_j, P_k^num)
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
                        # Partie catégorielle : 1 - d_raw, avec d_raw dans [0,1]
                        xj_char <- as.character(X_used[[j]])
                        zk_cat  <- proto_k$cat   # vecteur de longueur n

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
            } else {
                cat("Méthode non reconnue pour summary().\n")
                return(invisible(NULL))
            }

            # --- 2) Résumé par cluster ---
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

            # --- 3) Tableau variables / cluster / adhésion ---
            var_df <- data.frame(
                variable  = colnames(X),
                cluster   = clusters,
                adhesion  = membership,
                stringsAsFactors = FALSE
            )

            cat("Indicateur d'adhésion (", membership_label, ")\n", sep = "")
            # On affiche seulement les 10 premières lignes triées par cluster puis par adhésion décroissante
            ord <- order(var_df$cluster, -var_df$adhesion)
            print(utils::head(var_df[ord, ], n = min(10L, nrow(var_df))), row.names = FALSE)
            cat("\n(Le data.frame complet est renvoyé invisiblement.)\n")

            invisible(var_df)
        },

        # --- Graphiques ---
        plot = function(type = c("inertia", "clusters", "membership"), Ks = NULL, ...) {
            type <- match.arg(type)

            if (is.null(private$FX_active) || is.null(private$FClusters)) {
                stop("Aucun modèle appris. Appelez fit() d'abord.")
            }

            X        <- private$FX_active
            method   <- private$FMethod
            clusters <- private$FClusters
            centers  <- private$FCenters

            if (type == "clusters") {
                # Barplot des tailles de clusters
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
                # On réutilise la logique de summary() pour calculer l'adhésion
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
                    lambda  <- private$FLambda
                    X_used  <- X

                    for (j in seq_len(p)) {
                        kj <- clusters[j]
                        proto_k <- centers[[kj]]

                        if (j %in% num_idx) {
                            xj <- X_used[[j]]
                            zk_num <- proto_k$num
                            r <- suppressWarnings(stats::cor(xj, zk_num, use = "pairwise.complete.obs"))
                            if (is.na(r)) r <- 0
                            membership[j] <- r^2
                        } else if (j %in% cat_idx) {
                            xj_char <- as.character(X_used[[j]])
                            zk_cat <- proto_k$cat  # vecteur longueur n
                            ...
                            mismatch <- xj_char != zk_cat
                            d_raw <- mean(mismatch, na.rm = TRUE)
                            ...
                            membership[j] <- 1 - d_raw
                        } else {
                            membership[j] <- NA_real_
                        }
                    }
                    ylab <- "Score d'adhésion (r^2 num., 1 - dissimilarité cat.)"
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
                # Courbe inertie(K) en relançant l'algo pour plusieurs valeurs de K
                p <- ncol(X)
                if (is.null(Ks)) {
                    Ks <- seq_len(min(10L, p))
                }
                Ks <- Ks[Ks >= 2]
                if (length(Ks) == 0L) {
                    stop("Impossible de tracer la courbe d'inertie : K doit être >= 2.")
                }

                old_K <- private$FNbGroupes
                inertias <- numeric(length(Ks))

                for (i in seq_along(Ks)) {
                    private$FNbGroupes <- as.integer(Ks[i])
                    res_i <- switch(
                        method,
                        "kmeans"      = private$run_kmeans(X),
                        "kmodes"      = private$run_kmodes(X),
                        "kprototypes" = private$run_kprototypes(X),
                        stop("Méthode non supportée pour le calcul de l'inertie.")
                    )
                    inertias[i] <- res_i$inertia
                }

                private$FNbGroupes <- old_K  # on restaure

                graphics::plot(
                    Ks, inertias, type = "b",
                    xlab = "K (nombre de clusters)",
                    ylab = "Inertie intra-cluster",
                    main = "Courbe de l'inertie en fonction de K",
                    ...
                )
                return(invisible(NULL))
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

        get_X_descr = function() {
            private$FX_descr
        },

        get_algorithm = function() {
            private$FAlgorithme
        }
    )
)
