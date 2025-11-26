mmrClustVarKModes <- R6::R6Class(
    "mmrClustVarKModes",
    inherit = mmrClustVarBase,
    
    public = list(
        initialize = function(K, scale = TRUE, lambda = 1, ...) {
            # scale et lambda sont ignorés pour k-modes, mais on garde la même API
            super$initialize(
                K           = K,
                scale       = scale,
                lambda      = lambda,
                method_name = "kmodes"
            )
        }
    ),
    
    private = list(
        
        # ==========================
        # 1. ALGORITHME K-MODES VAR
        # ==========================
        run_clustering = function(X) {
            # X : data.frame catégoriel (factor/character),
            # déjà passé par check_and_prepare_X()
            
            # --- 1) Vérifier que tout est qualitatif ---
            is_cat <- vapply(
                X,
                function(col) is.factor(col) || is.character(col),
                logical(1L)
            )
            if (!all(is_cat)) {
                stop("[mmrClustVarKModes] Toutes les variables doivent être qualitatives (factor ou character).")
            }
            
            # On travaille en matrice de caractères pour simplifier
            X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)
            X_mat  <- as.matrix(X_char)
            
            n <- nrow(X_mat)   # individus
            p <- ncol(X_mat)   # variables
            K <- private$FNbGroupes
            
            if (K > p) {
                stop("[mmrClustVarKModes] K ne peut pas dépasser le nombre de variables.")
            }
            
            max_iter <- 50L
            converged <- FALSE
            
            # =====================
            # Initialisation simple
            # =====================
            set.seed(123)
            noyaux <- sample(seq_len(p), size = K, replace = FALSE)
            clusters <- rep(NA_integer_, p)
            clusters[noyaux] <- seq_len(K)
            
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
            
            # prototypes : liste de K vecteurs de longueur n (un "mode" par individu)
            centers <- vector("list", K)
            
            # --- helper : calcul du mode par individu pour un cluster ---
            compute_mode_profile <- function(cols_k) {
                # cols_k : indices des colonnes dans le cluster
                if (length(cols_k) == 1L) {
                    # avec une seule variable, le "mode" est juste cette variable
                    return(X_mat[, cols_k])
                } else {
                    # pour chaque individu i, on prend la modalité la plus fréquente parmi les variables du cluster
                    mode_vec <- character(n)
                    for (i in seq_len(n)) {
                        vals <- X_mat[i, cols_k]
                        # table des modalités (on ignore les NA)
                        vals_no_na <- vals[!is.na(vals)]
                        if (length(vals_no_na) == 0L) {
                            mode_vec[i] <- NA_character_
                        } else {
                            tab <- table(vals_no_na)
                            # en cas d'égalité, on prend la première dans l'ordre de table()
                            mode_vec[i] <- names(tab)[which.max(tab)]
                        }
                    }
                    mode_vec
                }
            }
            
            # --- helper : dissimilarité simple matching ---
            # d = proportion de positions où x != m
            simple_matching <- function(x, m) {
                mismatch <- (x != m)
                # on ignore les NA pour le calcul, et si tout est NA, on renvoie 1 (max dissimilarité)
                d <- mean(mismatch, na.rm = TRUE)
                if (is.na(d)) d <- 1
                d
            }
            
            # --- Boucle principale ---
            for (iter in seq_len(max_iter)) {
                
                # 1) recalcul des "modes" pour chaque cluster
                for (k in seq_len(K)) {
                    cols_k <- which(clusters == k)
                    if (length(cols_k) == 0L) {
                        # cluster vide : on lui donne une variable au hasard
                        cols_k <- sample(seq_len(p), size = 1L)
                        clusters[cols_k] <- k
                    }
                    centers[[k]] <- compute_mode_profile(cols_k)
                }
                
                # 2) réaffectation des variables
                new_clusters <- clusters
                
                for (j in seq_len(p)) {
                    xj <- X_mat[, j]
                    d_all <- vapply(
                        centers,
                        function(mk) simple_matching(xj, mk),
                        numeric(1L)
                    )
                    # on cherche la dissimilarité minimale
                    k_best <- which.min(d_all)
                    new_clusters[j] <- k_best
                }
                
                # 3) critère de convergence : plus aucun changement d'affectation
                if (all(new_clusters == clusters)) {
                    converged <- TRUE
                    clusters  <- new_clusters
                    break
                }
                
                clusters <- new_clusters
            }
            
            # Calcul de l'inertie intra-cluster :
            # somme des dissimilarités simple matching pour la partition finale
            inertia <- 0
            for (j in seq_len(p)) {
                k  <- clusters[j]
                mk <- centers[[k]]
                xj <- X_mat[, j]
                inertia <- inertia + simple_matching(xj, mk)
            }
            
            list(
                clusters  = clusters,
                centers   = centers,   # liste de vecteurs de caractères (modes par individu)
                inertia   = inertia,
                converged = converged
            )
        },
        
        # =====================
        # 2. PREDICT UNE VAR
        # =====================
        predict_one_variable = function(x_new, var_name) {
            # x_new : vecteur factor/character de longueur n
            # var_name : nom de la variable
            
            if (!(is.factor(x_new) || is.character(x_new))) {
                stop("[mmrClustVarKModes] predict() pour k-modes requiert une variable qualitative.")
            }
            
            centers <- private$FCenters
            if (is.null(centers)) {
                stop("[mmrClustVarKModes] Aucun prototype disponible (fit() non appelé ?)")
            }
            
            x_char <- as.character(x_new)
            
            simple_matching <- function(x, m) {
                mismatch <- (x != m)
                d <- mean(mismatch, na.rm = TRUE)
                if (is.na(d)) d <- 1
                d
            }
            
            K <- length(centers)
            d_all <- vapply(
                centers,
                function(mk) simple_matching(x_char, mk),
                numeric(1L)
            )
            
            k_best   <- which.min(d_all)
            d_raw    <- d_all[k_best]
            adhesion <- 1 - d_raw  # proportion de matches
            
            data.frame(
                variable = var_name,
                cluster  = as.integer(k_best),
                distance = d_raw,
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
                cat("(k-modes) Aucune variable active stockée.\n")
                return(invisible(NULL))
            }
            
            X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)
            X_mat  <- as.matrix(X_char)
            
            p        <- ncol(X_mat)
            clusters <- private$FClusters
            centers  <- private$FCenters
            
            simple_matching <- function(x, m) {
                mismatch <- (x != m)
                d <- mean(mismatch, na.rm = TRUE)
                if (is.na(d)) d <- 1
                d
            }
            
            membership <- numeric(p)
            distance   <- numeric(p)
            
            for (j in seq_len(p)) {
                k  <- clusters[j]
                mk <- centers[[k]]
                xj <- X_mat[, j]
                d_raw <- simple_matching(xj, mk)
                distance[j]   <- d_raw
                membership[j] <- 1 - d_raw
            }
            
            df <- data.frame(
                variable = colnames(X_mat),
                cluster  = as.integer(clusters),
                distance = distance,
                adhesion = membership,
                stringsAsFactors = FALSE
            )
            
            cat("=== Indicateurs d'adhésion (k-modes) ===\n")
            cat("distance = proportion de mismatches avec le mode du cluster\n")
            cat("adhesion = 1 - distance (proportion de matches)\n\n")
            
            # Indicateurs globaux
            dist_globale <- mean(df$distance, na.rm = TRUE)
            adh_globale  <- mean(df$adhesion, na.rm = TRUE)
            
            # Stats par cluster
            stats_list <- lapply(split(df, df$cluster), function(dsub) {
                c(
                    cluster   = dsub$cluster[1],
                    dist_mean = mean(dsub$distance, na.rm = TRUE),
                    dist_min  = min(dsub$distance, na.rm = TRUE),
                    dist_max  = max(dsub$distance, na.rm = TRUE),
                    adh_mean  = mean(dsub$adhesion, na.rm = TRUE),
                    adh_min   = min(dsub$adhesion, na.rm = TRUE),
                    adh_max   = max(dsub$adhesion, na.rm = TRUE)
                )
            })
            stats_par_cluster <- as.data.frame(do.call(rbind, stats_list))
            stats_par_cluster$cluster <- as.integer(stats_par_cluster$cluster)
            
            cat(sprintf("Distance moyenne globale : %.3f\n", dist_globale))
            cat(sprintf("Adhésion moyenne globale : %.3f\n\n", adh_globale))
            
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
                stop("[mmrClustVarKModes] plot(type = 'membership') : aucun X actif.")
            }
            
            X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)
            X_mat  <- as.matrix(X_char)
            
            p        <- ncol(X_mat)
            clusters <- private$FClusters
            centers  <- private$FCenters
            
            simple_matching <- function(x, m) {
                mismatch <- (x != m)
                d <- mean(mismatch, na.rm = TRUE)
                if (is.na(d)) d <- 1
                d
            }
            
            membership <- numeric(p)
            for (j in seq_len(p)) {
                k  <- clusters[j]
                mk <- centers[[k]]
                xj <- X_mat[, j]
                d_raw       <- simple_matching(xj, mk)
                membership[j] <- 1 - d_raw
            }
            
            o <- order(membership, decreasing = TRUE)
            barplot(
                membership[o],
                names.arg = colnames(X_mat)[o],
                las = 2,
                main = "Adhésion des variables aux clusters (k-modes)",
                ylab = "adhésion (1 - proportion de mismatches)",
                cex.names = 0.7
            )
        },
        
        # ===========================
        # 5. PLOT : profils moyens
        # ===========================
        plot_profiles = function() {
            X <- private$FX_active
            if (is.null(X)) {
                warning("[mmrClustVarKModes] plot(type = 'profiles') : aucun X actif.")
                return(invisible(NULL))
            }
            
            X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)
            X_mat  <- as.matrix(X_char)
            
            n        <- nrow(X_mat)
            p        <- ncol(X_mat)
            clusters <- private$FClusters
            centers  <- private$FCenters
            K        <- private$FNbGroupes
            
            if (is.null(clusters) || is.null(centers)) {
                warning("[mmrClustVarKModes] Pas de clusters / prototypes disponibles.")
                return(invisible(NULL))
            }
            
            # prof_mat[i, k] = proportion de variables du cluster k
            # pour lesquelles x_{i,j} == mode_{i}^{(k)}
            prof_mat <- matrix(NA_real_, nrow = n, ncol = K)
            colnames(prof_mat) <- paste0("Cluster ", seq_len(K))
            rownames(prof_mat) <- seq_len(n)
            
            for (k in seq_len(K)) {
                vars_k <- which(clusters == k)
                if (length(vars_k) == 0L) next
                
                mk <- centers[[k]]          # vecteur de longueur n (mode par individu)
                
                for (i in seq_len(n)) {
                    vals <- X_mat[i, vars_k]
                    # on compare aux modes pour cet individu
                    matches <- (vals == mk[i])
                    v <- mean(matches, na.rm = TRUE)
                    if (is.na(v)) v <- NA_real_
                    prof_mat[i, k] <- v
                }
            }
            
            op <- par(no.readonly = TRUE)
            on.exit(par(op))
            
            image(
                x = seq_len(K),
                y = seq_len(n),
                z = t(prof_mat),
                xlab = "Clusters",
                ylab = "Individus",
                main = "Profils moyens par cluster (k-modes)\n(proportion de matches au mode)",
                axes = FALSE
            )
            axis(1, at = seq_len(K), labels = colnames(prof_mat))
            axis(2, at = pretty(seq_len(n)))
            box()
            
            invisible(NULL)
        }
    )
)
