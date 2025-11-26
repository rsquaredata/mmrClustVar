mmrClustVarKMeans <- R6::R6Class(
    "mmrClustVarKMeans",
    inherit = mmrClustVarBase,
    
    public = list(
        initialize = function(K, scale = TRUE, lambda = 1, ...) {
            # lambda est ignoré pour k-means, mais on garde la signature homogène
            super$initialize(
                K           = K,
                scale       = scale,
                lambda      = lambda,
                method_name = "kmeans"
            )
        }
    ),
    
    private = list(
        
        # ==========================
        # 1. ALGORITHME K-MEANS VAR
        # ==========================
        run_clustering = function(X) {
            # X : data.frame numérique, déjà passé par check_and_prepare_X()
            # et éventuellement standardisé
            
            # Vérifier qu'on n'a QUE du numérique
            if (!all(vapply(X, is.numeric, logical(1L)))) {
                stop("[mmrClustVarKMeans] Toutes les variables doivent être numériques.")
            }
            
            n <- nrow(X)
            p <- ncol(X)
            K <- private$FNbGroupes
            
            if (K > p) {
                stop("[mmrClustVarKMeans] K ne peut pas dépasser le nombre de variables.")
            }
            
            # matrice n x p (individus x variables)
            X_mat <- as.matrix(X)
            
            # paramètres algo
            max_iter <- 50L
            tol      <- 1e-6
            
            # =====================
            # Initialisation simple
            # =====================
            set.seed(123)
            noyaux <- sample(seq_len(p), size = K, replace = FALSE)
            clusters <- rep(NA_integer_, p)
            clusters[noyaux] <- seq_len(K)
            
            # affectation initiale des autres colonnes au hasard
            idx_other <- setdiff(seq_len(p), noyaux)
            if (length(idx_other) > 0L) {
                clusters[idx_other] <- sample(seq_len(K), size = length(idx_other), replace = TRUE)
            }
            
            # Assurer aucun cluster vide
            for (k in seq_len(K)) {
                if (!any(clusters == k)) {
                    j_free <- which.max(tabulate(clusters))
                    clusters[j_free] <- k
                }
            }
            
            # liste des composantes latentes Z_k (vecteurs de longueur n)
            centers <- vector("list", K)
            
            # fonction interne : calcule Z_k = 1ère composante principale sur les variables du cluster k
            compute_Zk <- function(cols_k) {
                if (length(cols_k) == 1L) {
                    # cas trivial : une seule variable → profil = variable centrée
                    zk <- X_mat[, cols_k]
                    zk <- zk - mean(zk, na.rm = TRUE)
                    return(as.numeric(zk))
                } else {
                    # plusieurs variables : ACP sur les colonnes du cluster
                    Xk <- X_mat[, cols_k, drop = FALSE]
                    
                    # on ne garde que les lignes complètes pour l'ACP
                    idx_ok <- stats::complete.cases(Xk)
                    
                    # si pas assez de lignes complètes -> 1ère variable centrée
                    if (sum(idx_ok) < 2L) {
                        zk <- Xk[, 1]
                        zk <- zk - mean(zk, na.rm = TRUE)
                        return(as.numeric(zk))
                    }
                    
                    # ACP sur les lignes complètes uniquement
                    pc <- stats::prcomp(
                        Xk[idx_ok, , drop = FALSE],
                        center  = FALSE,
                        scale.  = FALSE
                    )
                    
                    z_short <- pc$x[, 1]              # scores sur les lignes complètes
                    zk <- rep(NA_real_, nrow(Xk))     # vecteur de longueur n
                    zk[idx_ok] <- z_short             # on remplit uniquement là où on a l'ACP
                    
                    return(as.numeric(zk))
                }
            }
            
            # fonction interne : r^2(X_j, Z_k)
            r2_corr <- function(x, z) {
                r <- suppressWarnings(stats::cor(x, z, use = "pairwise.complete.obs"))
                if (is.na(r)) r <- 0
                r^2
            }
            
            # fonction pour calculer l'objectif W = somme r^2(X_j, Z_{cluster(j)})
            compute_W <- function(clusters, centers) {
                W <- 0
                for (j in seq_len(p)) {
                    k  <- clusters[j]
                    zk <- centers[[k]]
                    xj <- X_mat[, j]
                    W <- W + r2_corr(xj, zk)
                }
                W
            }
            
            # =============================
            # Boucle de réallocation
            # =============================
            old_W    <- -Inf
            converged <- FALSE
            
            for (iter in seq_len(max_iter)) {
                # 1) recalcul des composantes latentes Z_k
                for (k in seq_len(K)) {
                    cols_k <- which(clusters == k)
                    if (length(cols_k) == 0L) {
                        cols_k <- sample(seq_len(p), size = 1L)
                        clusters[cols_k] <- k
                    }
                    centers[[k]] <- compute_Zk(cols_k)
                }
                
                # 2) réaffectation des variables
                new_clusters <- clusters
                
                for (j in seq_len(p)) {
                    xj <- X_mat[, j]
                    r2_all <- vapply(centers,
                                     function(zk) r2_corr(xj, zk),
                                     numeric(1L))
                    k_best <- which.max(r2_all)
                    new_clusters[j] <- k_best
                }
                
                # 3) calcul de W et critère d'arrêt
                W <- compute_W(new_clusters, centers)
                
                if (all(new_clusters == clusters)) {
                    converged <- TRUE
                    clusters  <- new_clusters
                    break
                }
                
                if (abs(W - old_W) < tol) {
                    converged <- TRUE
                    clusters  <- new_clusters
                    break
                }
                
                clusters <- new_clusters
                old_W    <- W
            }
            
            # Inertie intra définie comme somme (1 - r^2) pour la partition finale
            inertia <- 0
            for (j in seq_len(p)) {
                k  <- clusters[j]
                zk <- centers[[k]]
                xj <- X_mat[, j]
                inertia <- inertia + (1 - r2_corr(xj, zk))
            }
            
            list(
                clusters  = clusters,
                centers   = centers,   # liste de Z_k (vecteur de longueur n)
                inertia   = inertia,
                converged = converged
            )
        },
        
        # =====================
        # 2. PREDICT UNE VAR
        # =====================
        predict_one_variable = function(x_new, var_name) {
            # x_new : vecteur numérique de longueur n (même nb d'individus)
            if (!is.numeric(x_new)) {
                stop("[mmrClustVarKMeans] predict() pour k-means requiert une variable numérique.")
            }
            
            centers <- private$FCenters
            if (is.null(centers)) {
                stop("[mmrClustVarKMeans] Aucun prototype disponible (fit() non appelé ?)")
            }
            
            r2_corr <- function(x, z) {
                r <- suppressWarnings(stats::cor(x, z, use = "pairwise.complete.obs"))
                if (is.na(r)) r <- 0
                r^2
            }
            
            K <- length(centers)
            r2_all <- vapply(
                centers,
                function(zk) r2_corr(x_new, zk),
                numeric(1L)
            )
            
            k_best   <- which.max(r2_all)
            adhesion <- r2_all[k_best]
            
            data.frame(
                variable = var_name,
                cluster  = as.integer(k_best),
                adhesion = adhesion,
                stringsAsFactors = FALSE
            )
        },
        
        # ===========================
        # 3. SUMMARY : adhésions
        # ===========================
        summary_membership = function() {
            X <- private$FX_active
            if (is.null(X)) {
                cat("(k-means) Aucune variable active stockée.\n")
                return(invisible(NULL))
            }
            
            X_mat    <- as.matrix(X)
            p        <- ncol(X_mat)
            clusters <- private$FClusters
            centers  <- private$FCenters
            
            r2_corr <- function(x, z) {
                r <- suppressWarnings(stats::cor(x, z, use = "pairwise.complete.obs"))
                if (is.na(r)) r <- 0
                r^2
            }
            
            membership <- numeric(p)
            for (j in seq_len(p)) {
                k  <- clusters[j]
                zk <- centers[[k]]
                xj <- X_mat[, j]
                membership[j] <- r2_corr(xj, zk)
            }
            
            df <- data.frame(
                variable = colnames(X_mat),
                cluster  = as.integer(clusters),
                r2       = membership,
                stringsAsFactors = FALSE
            )
            
            cat("=== Indicateurs d'adhésion (k-means) ===\n")
            cat("r^2 = corr^2(variable, composante latente du cluster)\n\n")
            
            # Part d'inertie expliquée = moyenne des r^2 (entre 0 et 1)
            part_expliquee <- mean(df$r2, na.rm = TRUE)
            
            # Statistiques par cluster (moyenne / min / max des r^2)
            stats_list <- lapply(split(df, df$cluster), function(dsub) {
                c(
                    cluster = dsub$cluster[1],
                    r2_mean = mean(dsub$r2, na.rm = TRUE),
                    r2_min  = min(dsub$r2, na.rm = TRUE),
                    r2_max  = max(dsub$r2, na.rm = TRUE)
                )
            })
            stats_par_cluster <- as.data.frame(do.call(rbind, stats_list))
            stats_par_cluster$cluster <- as.integer(stats_par_cluster$cluster)
            
            cat(sprintf("Part d'inertie expliquée (moyenne des r^2) : %.3f\n\n",
                        part_expliquee))
            
            cat("--- Statistiques par cluster ---\n")
            print(stats_par_cluster)
            
            cat("\n--- Détail par variable ---\n")
            print(df)
            
            invisible(df)
        },
        
        # ===========================
        # 4. PLOT : adhésion
        # ===========================
        plot_membership = function() {
            X <- private$FX_active
            if (is.null(X)) {
                stop("[mmrClustVarKMeans] plot(type = 'membership') : aucun X actif.")
            }
            
            X_mat    <- as.matrix(X)
            p        <- ncol(X_mat)
            clusters <- private$FClusters
            centers  <- private$FCenters
            
            r2_corr <- function(x, z) {
                r <- suppressWarnings(stats::cor(x, z, use = "pairwise.complete.obs"))
                if (is.na(r)) r <- 0
                r^2
            }
            
            membership <- numeric(p)
            for (j in seq_len(p)) {
                k  <- clusters[j]
                zk <- centers[[k]]
                xj <- X_mat[, j]
                membership[j] <- r2_corr(xj, zk)
            }
            
            o <- order(membership, decreasing = TRUE)
            barplot(
                membership[o],
                names.arg = colnames(X_mat)[o],
                las = 2,
                main = "Adhésion des variables aux clusters (k-means)",
                ylab = "r^2",
                cex.names = 0.7
            )
        },
        
        # ===========================
        # 5. PLOT : profils moyens
        # ===========================
        plot_profiles = function() {
            X <- private$FX_active
            if (is.null(X)) {
                warning("[mmrClustVarKMeans] plot(type = 'profiles') : aucun X actif.")
                return(invisible(NULL))
            }
            
            X_mat    <- as.matrix(X)
            n        <- nrow(X_mat)
            clusters <- private$FClusters
            K        <- private$FNbGroupes
            
            if (is.null(clusters)) {
                warning("[mmrClustVarKMeans] Pas de clusters disponibles.")
                return(invisible(NULL))
            }
            
            # profil moyen par individu et par cluster :
            # moyenne des variables du cluster
            prof_mat <- matrix(NA_real_, nrow = n, ncol = K)
            colnames(prof_mat) <- paste0("Cluster ", seq_len(K))
            rownames(prof_mat) <- seq_len(n)
            
            for (k in seq_len(K)) {
                vars_k <- which(clusters == k)
                if (length(vars_k) == 0L) next
                prof_mat[, k] <- rowMeans(X_mat[, vars_k, drop = FALSE], na.rm = TRUE)
            }
            
            op <- par(no.readonly = TRUE)
            on.exit(par(op))
            
            image(
                x = seq_len(K),
                y = seq_len(n),
                z = t(prof_mat),
                xlab = "Clusters",
                ylab = "Individus",
                main = "Profils moyens par cluster (k-means)",
                axes = FALSE
            )
            axis(1, at = seq_len(K), labels = colnames(prof_mat))
            axis(2, at = pretty(seq_len(n)))
            box()
            
            invisible(NULL)
        }
    )
)
