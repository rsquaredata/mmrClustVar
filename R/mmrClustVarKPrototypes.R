mmrClustVarKPrototypes <- R6::R6Class(
    "mmrClustVarKPrototypes",
    inherit = mmrClustVarBase,
    
    public = list(
        initialize = function(K, scale = TRUE, lambda = 1, ...) {
            super$initialize(
                K           = K,
                scale       = scale,
                lambda      = lambda,
                method_name = "kprototypes"
            )
        }
    ),
    
    private = list(
        
        # ==========================
        # 1. ALGORITHME K-PROTOTYPES
        # ==========================
        run_clustering = function(X) {
            # X : data.frame mixte (num + quali)
            n <- nrow(X)
            p <- ncol(X)
            K <- private$FNbGroupes
            lambda <- private$FLambda
            
            if (K > p) {
                stop("[mmrClustVarKPrototypes] K ne peut pas dépasser le nombre de variables.")
            }
            
            # --- Détection des types ---
            is_num <- vapply(X, is.numeric, logical(1L))
            is_cat <- (!is_num) & vapply(
                X,
                function(col) is.factor(col) || is.character(col),
                logical(1L)
            )
            
            if (!any(is_num) && !any(is_cat)) {
                stop("[mmrClustVarKPrototypes] Aucune variable utilisable (numérique ou qualitative).")
            }
            
            # on autorise des jeux purement num ou purement quali ici
            X_num <- NULL
            X_cat <- NULL
            idx_num <- which(is_num)
            idx_cat <- which(is_cat)
            
            if (length(idx_num) > 0L) {
                X_num <- as.matrix(X[, idx_num, drop = FALSE])
            }
            if (length(idx_cat) > 0L) {
                X_cat <- as.data.frame(
                    lapply(X[, idx_cat, drop = FALSE], as.character),
                    stringsAsFactors = FALSE
                )
                X_cat <- as.matrix(X_cat)
            }
            
            max_iter <- 50L
            converged <- FALSE
            
            # =====================
            # Initialisation simple
            # =====================
            set.seed(123)
            clusters <- sample(seq_len(K), size = p, replace = TRUE)
            
            # éviter les clusters vides
            for (k in seq_len(K)) {
                if (!any(clusters == k)) {
                    j_free <- which.max(tabulate(clusters))
                    clusters[j_free] <- k
                }
            }
            
            # prototypes : liste de K éléments, chacun étant une liste(num = ..., cat = ...)
            centers <- vector("list", K)
            
            # --- helpers prototypes ---
            compute_num_profile <- function(cols_k_global) {
                if (is.null(X_num) || length(cols_k_global) == 0L) return(NULL)
                local_idx <- match(cols_k_global, idx_num)
                local_idx <- local_idx[!is.na(local_idx)]
                if (length(local_idx) == 0L) return(NULL)
                prof <- rowMeans(X_num[, local_idx, drop = FALSE], na.rm = TRUE)
                # si tout NA sur une ligne, rowMeans renvoie NA → on garde, la corr gèrera.
                prof
            }
            
            compute_cat_profile <- function(cols_k_global) {
                if (is.null(X_cat) || length(cols_k_global) == 0L) return(NULL)
                local_idx <- match(cols_k_global, idx_cat)
                local_idx <- local_idx[!is.na(local_idx)]
                if (length(local_idx) == 0L) return(NULL)
                
                mode_vec <- character(n)
                for (i in seq_len(n)) {
                    vals <- X_cat[i, local_idx]
                    vals_no_na <- vals[!is.na(vals)]
                    if (length(vals_no_na) == 0L) {
                        mode_vec[i] <- NA_character_
                    } else {
                        tab <- table(vals_no_na)
                        mode_vec[i] <- names(tab)[which.max(tab)]
                    }
                }
                mode_vec
            }
            
            # --- helpers distance / adhésion ---
            dist_num <- function(x, z) {
                r <- suppressWarnings(
                    stats::cor(x, z, use = "pairwise.complete.obs")
                )
                if (is.na(r)) r <- 0
                d <- 1 - r^2
                list(distance = d, adhesion = r^2)
            }
            
            dist_cat <- function(x_char, z_cat, lambda) {
                mismatch <- (x_char != z_cat)
                d_raw <- mean(mismatch, na.rm = TRUE)
                if (is.na(d_raw)) d_raw <- 1
                list(
                    distance = lambda * d_raw, # utilisé dans l'algo
                    adhesion = 1 - d_raw       # indicateur d'adhésion affiché
                )
            }
            
            # --- Boucle principale ---
            for (iter in seq_len(max_iter)) {
                
                # 1) recalcul des prototypes pour chaque cluster
                for (k in seq_len(K)) {
                    cols_k <- which(clusters == k)
                    if (length(cols_k) == 0L) {
                        # cluster vide → récupère une variable d'un gros cluster
                        j_free <- which.max(tabulate(clusters))
                        clusters[j_free] <- k
                        cols_k <- which(clusters == k)
                    }
                    
                    cols_num_k <- intersect(cols_k, idx_num)
                    cols_cat_k <- intersect(cols_k, idx_cat)
                    
                    centers[[k]] <- list(
                        num = compute_num_profile(cols_num_k),
                        cat = compute_cat_profile(cols_cat_k)
                    )
                }
                
                # 2) réaffectation des variables
                new_clusters <- clusters
                
                for (j in seq_len(p)) {
                    if (is_num[j]) {
                        # variable numérique → distance_num à chaque prototype.num
                        xj <- X[, j]
                        d_all <- numeric(K)
                        for (k in seq_len(K)) {
                            zk <- centers[[k]]$num
                            if (is.null(zk)) {
                                d_all[k] <- 1 # prototype numérique absent → distance max
                            } else {
                                da <- dist_num(xj, zk)
                                d_all[k] <- da$distance
                            }
                        }
                        new_clusters[j] <- which.min(d_all)
                        
                    } else if (is_cat[j]) {
                        # variable qualitative → lambda * simple matching à prototype.cat
                        xj <- as.character(X[, j])
                        d_all <- numeric(K)
                        for (k in seq_len(K)) {
                            zk <- centers[[k]]$cat
                            if (is.null(zk)) {
                                d_all[k] <- lambda # cluster sans partie cat → distance max lambda
                            } else {
                                da <- dist_cat(xj, zk, lambda)
                                d_all[k] <- da$distance
                            }
                        }
                        new_clusters[j] <- which.min(d_all)
                        
                    } else {
                        stop("[mmrClustVarKPrototypes] Type de variable non géré (ni num, ni factor/character).")
                    }
                }
                
                if (all(new_clusters == clusters)) {
                    converged <- TRUE
                    clusters  <- new_clusters
                    break
                }
                
                clusters <- new_clusters
            }
            
            # 3) inertie intra-cluster finale
            inertia <- 0
            for (j in seq_len(p)) {
                k <- clusters[j]
                if (is_num[j]) {
                    xj <- X[, j]
                    zk <- centers[[k]]$num
                    if (is.null(zk)) {
                        inertia <- inertia + 1
                    } else {
                        da <- dist_num(xj, zk)
                        inertia <- inertia + da$distance
                    }
                } else if (is_cat[j]) {
                    xj <- as.character(X[, j])
                    zk <- centers[[k]]$cat
                    if (is.null(zk)) {
                        inertia <- inertia + lambda
                    } else {
                        da <- dist_cat(xj, zk, lambda)
                        inertia <- inertia + da$distance
                    }
                }
            }
            
            list(
                clusters  = clusters,
                centers   = centers,  # liste de K, chacun list(num=..., cat=...)
                inertia   = inertia,
                converged = converged
            )
        },
        
        # =====================
        # 2. PREDICT UNE VAR
        # =====================
        predict_one_variable = function(x_new, var_name) {
            centers <- private$FCenters
            if (is.null(centers)) {
                stop("[mmrClustVarKPrototypes] Aucun prototype disponible (fit() non appelé ?)")
            }
            
            lambda <- private$FLambda
            is_num_var <- is.numeric(x_new)
            
            # helpers cohérents avec run_clustering
            dist_num <- function(x, z) {
                r <- suppressWarnings(
                    stats::cor(x, z, use = "pairwise.complete.obs")
                )
                if (is.na(r)) r <- 0
                d <- 1 - r^2
                list(distance = d, adhesion = r^2)
            }
            
            dist_cat <- function(x_char, z_cat, lambda) {
                mismatch <- (x_char != z_cat)
                d_raw <- mean(mismatch, na.rm = TRUE)
                if (is.na(d_raw)) d_raw <- 1
                list(
                    distance = lambda * d_raw,
                    adhesion = 1 - d_raw
                )
            }
            
            K <- length(centers)
            d_all        <- numeric(K)
            adhesion_all <- numeric(K)
            
            if (is_num_var) {
                for (k in seq_len(K)) {
                    zk <- centers[[k]]$num
                    if (is.null(zk)) {
                        d_all[k]        <- 1
                        adhesion_all[k] <- 0
                    } else {
                        da <- dist_num(x_new, zk)
                        d_all[k]        <- da$distance
                        adhesion_all[k] <- da$adhesion
                    }
                }
            } else {
                x_char <- as.character(x_new)
                for (k in seq_len(K)) {
                    zk <- centers[[k]]$cat
                    if (is.null(zk)) {
                        d_all[k]        <- lambda
                        adhesion_all[k] <- 0
                    } else {
                        da <- dist_cat(x_char, zk, lambda)
                        d_all[k]        <- da$distance
                        adhesion_all[k] <- da$adhesion
                    }
                }
            }
            
            k_best <- which.min(d_all)
            
            data.frame(
                variable = var_name,
                cluster  = as.integer(k_best),
                distance = d_all[k_best],
                adhesion = adhesion_all[k_best],
                stringsAsFactors = FALSE
            )
        },
        
        # ===========================
        # 3. SUMMARY : adhésions
        # ===========================
        summary_membership = function() {
            X <- private$FX_active
            if (is.null(X)) {
                cat("(k-prototypes) Aucune variable active stockée.\n")
                return(invisible(NULL))
            }
            
            n <- nrow(X)
            p <- ncol(X)
            clusters <- private$FClusters
            centers  <- private$FCenters
            lambda   <- private$FLambda
            
            is_num <- vapply(X, is.numeric, logical(1L))
            is_cat <- (!is_num) & vapply(
                X,
                function(col) is.factor(col) || is.character(col),
                logical(1L)
            )
            
            dist_num <- function(x, z) {
                r <- suppressWarnings(
                    stats::cor(x, z, use = "pairwise.complete.obs")
                )
                if (is.na(r)) r <- 0
                d <- 1 - r^2
                list(distance = d, adhesion = r^2)
            }
            
            dist_cat <- function(x_char, z_cat, lambda) {
                mismatch <- (x_char != z_cat)
                d_raw <- mean(mismatch, na.rm = TRUE)
                if (is.na(d_raw)) d_raw <- 1
                list(
                    distance = lambda * d_raw,
                    adhesion = 1 - d_raw
                )
            }
            
            distance   <- numeric(p)
            adhesion   <- numeric(p)
            type_var   <- character(p)
            
            for (j in seq_len(p)) {
                k <- clusters[j]
                if (is_num[j]) {
                    xj <- X[, j]
                    zk <- centers[[k]]$num
                    if (is.null(zk)) {
                        distance[j] <- 1
                        adhesion[j] <- 0
                    } else {
                        da <- dist_num(xj, zk)
                        distance[j] <- da$distance
                        adhesion[j] <- da$adhesion
                    }
                    type_var[j] <- "numérique"
                } else if (is_cat[j]) {
                    xj <- as.character(X[, j])
                    zk <- centers[[k]]$cat
                    if (is.null(zk)) {
                        distance[j] <- lambda
                        adhesion[j] <- 0
                    } else {
                        da <- dist_cat(xj, zk, lambda)
                        distance[j] <- da$distance
                        adhesion[j] <- da$adhesion
                    }
                    type_var[j] <- "qualitative"
                } else {
                    distance[j] <- NA_real_
                    adhesion[j] <- NA_real_
                    type_var[j] <- "autre"
                }
            }
            
            df <- data.frame(
                variable = colnames(X),
                type     = type_var,
                cluster  = as.integer(clusters),
                distance = distance,
                adhesion = adhesion,
                stringsAsFactors = FALSE
            )
            
            cat("=== Indicateurs d'adhésion (k-prototypes) ===\n")
            cat("- distance :\n")
            cat("    * numériques : 1 - r^2(X_j, prototype_num)\n")
            cat("    * qualitatives : lambda × proportion de mismatches\n")
            cat("- adhesion :\n")
            cat("    * numériques : r^2\n")
            cat("    * qualitatives : 1 - proportion de mismatches (non pondéré par lambda)\n\n")
            print(df)
            
            invisible(df)
        },
        
        # ===========================
        # 4. PLOT : adhésion
        # ===========================
        plot_membership = function() {
            X <- private$FX_active
            if (is.null(X)) {
                stop("[mmrClustVarKPrototypes] plot(type='membership') : aucun X actif.")
            }
            
            n <- nrow(X)
            p <- ncol(X)
            clusters <- private$FClusters
            centers  <- private$FCenters
            lambda   <- private$FLambda
            
            is_num <- vapply(X, is.numeric, logical(1L))
            is_cat <- (!is_num) & vapply(
                X,
                function(col) is.factor(col) || is.character(col),
                logical(1L)
            )
            
            dist_num <- function(x, z) {
                r <- suppressWarnings(
                    stats::cor(x, z, use = "pairwise.complete.obs")
                )
                if (is.na(r)) r <- 0
                1 - r^2
            }
            
            dist_cat_raw <- function(x_char, z_cat) {
                mismatch <- (x_char != z_cat)
                d_raw <- mean(mismatch, na.rm = TRUE)
                if (is.na(d_raw)) d_raw <- 1
                d_raw
            }
            
            membership <- numeric(p)
            
            for (j in seq_len(p)) {
                k <- clusters[j]
                if (is_num[j]) {
                    xj <- X[, j]
                    zk <- centers[[k]]$num
                    if (is.null(zk)) {
                        membership[j] <- 0
                    } else {
                        d <- dist_num(xj, zk)
                        membership[j] <- 1 - d  # = r^2
                    }
                } else if (is_cat[j]) {
                    xj <- as.character(X[, j])
                    zk <- centers[[k]]$cat
                    if (is.null(zk)) {
                        membership[j] <- 0
                    } else {
                        d_raw <- dist_cat_raw(xj, zk)
                        membership[j] <- 1 - d_raw
                    }
                } else {
                    membership[j] <- NA_real_
                }
            }
            
            o <- order(membership, decreasing = TRUE)
            barplot(
                membership[o],
                names.arg = colnames(X)[o],
                las = 2,
                main = "Adhésion des variables aux clusters (k-prototypes)",
                ylab = "adhésion (r^2 ou 1 - mismatches)",
                cex.names = 0.7
            )
        }
    )
)
