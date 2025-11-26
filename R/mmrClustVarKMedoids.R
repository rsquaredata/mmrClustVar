mmrClustVarKMedoids <- R6::R6Class(
    "mmrClustVarKMedoids",
    inherit = mmrClustVarBase,
    
    public = list(
        initialize = function(K, scale = TRUE, lambda = 1, ...) {
            # lambda est conservé dans l'API mais peu utilisé ici
            super$initialize(
                K           = K,
                scale       = scale,
                lambda      = lambda,
                method_name = "kmedoids"
            )
        }
    ),
    
    private = list(
        
        # ==========================
        # 1. ALGORITHME K-MEDOIDS
        # ==========================
        run_clustering = function(X) {
            n <- nrow(X)
            p <- ncol(X)
            K <- private$FNbGroupes
            
            if (K > p) {
                stop("[mmrClustVarKMedoids] K ne peut pas dépasser le nombre de variables.")
            }
            
            # --- types ---
            is_num <- vapply(X, is.numeric, logical(1L))
            is_cat <- (!is_num) & vapply(
                X,
                function(col) is.factor(col) || is.character(col),
                logical(1L)
            )
            
            if (!any(is_num) && !any(is_cat)) {
                stop("[mmrClustVarKMedoids] Aucune variable utilisable (numérique ou qualitative).")
            }
            
            # version character pour les quali
            X_cat <- NULL
            if (any(is_cat)) {
                X_cat <- as.data.frame(
                    lapply(X[, is_cat, drop = FALSE], as.character),
                    stringsAsFactors = FALSE
                )
                X_cat <- as.matrix(X_cat)
            }
            X_num <- NULL
            if (any(is_num)) {
                X_num <- as.matrix(X[, is_num, drop = FALSE])
            }
            
            idx_num <- which(is_num)
            idx_cat <- which(is_cat)
            
            # --- helper distance entre deux variables X_j, X_l ---
            dist_var <- function(j, l) {
                # Determiner les types
                if (is_num[j] && is_num[l]) {
                    # 1 - r^2
                    xj <- X_num[, match(j, idx_num)]
                    xl <- X_num[, match(l, idx_num)]
                    r <- suppressWarnings(
                        stats::cor(xj, xl, use = "pairwise.complete.obs")
                    )
                    if (is.na(r)) r <- 0
                    return(1 - r^2)
                } else if (is_cat[j] && is_cat[l]) {
                    # simple matching
                    xj <- X_cat[, match(j, idx_cat)]
                    xl <- X_cat[, match(l, idx_cat)]
                    mismatch <- (xj != xl)
                    d <- mean(mismatch, na.rm = TRUE)
                    if (is.na(d)) d <- 1
                    return(d)
                } else {
                    # types différents : distance max
                    return(1)
                }
            }
            
            # --- matrice de dissimilarités D ---
            D <- matrix(0, nrow = p, ncol = p)
            for (j in seq_len(p)) {
                for (l in seq_len(j - 1L)) {
                    djl <- dist_var(j, l)
                    D[j, l] <- djl
                    D[l, j] <- djl
                }
            }
            
            max_iter <- 50L
            converged <- FALSE
            
            # =====================
            # Initialisation medoids
            # =====================
            set.seed(123)
            medoids <- sample(seq_len(p), size = K, replace = FALSE)
            
            # affectations initiales
            clusters <- apply(D[, medoids, drop = FALSE], 1, which.min)
            
            # éviter les clusters vides
            for (k in seq_len(K)) {
                if (!any(clusters == k)) {
                    j_free <- which.max(tabulate(clusters))
                    clusters[j_free] <- k
                    medoids[k] <- j_free
                }
            }
            
            # --- Boucle principale ---
            for (iter in seq_len(max_iter)) {
                
                # 1) mise à jour des medoids
                new_medoids <- medoids
                for (k in seq_len(K)) {
                    Gk <- which(clusters == k)
                    if (length(Gk) == 0L) {
                        # cluster vide → on récupère une variable d'un gros cluster
                        j_free <- which.max(tabulate(clusters))
                        clusters[j_free] <- k
                        Gk <- which(clusters == k)
                    }
                    
                    # coût pour chaque candidat j ∈ Gk : somme des distances D[j, l] pour l ∈ Gk
                    cost <- numeric(length(Gk))
                    for (idx in seq_along(Gk)) {
                        j_candidate <- Gk[idx]
                        cost[idx] <- sum(D[j_candidate, Gk])
                    }
                    new_medoids[k] <- Gk[which.min(cost)]
                }
                
                # 2) réaffectation des variables aux nouveaux medoids
                clusters_new <- apply(D[, new_medoids, drop = FALSE], 1, which.min)
                
                # 3) convergence ?
                if (all(new_medoids == medoids) && all(clusters_new == clusters)) {
                    converged <- TRUE
                    medoids   <- new_medoids
                    clusters  <- clusters_new
                    break
                }
                
                medoids  <- new_medoids
                clusters <- clusters_new
            }
            
            # inertie intra-cluster : somme des distances au medoid du cluster
            inertia <- 0
            for (j in seq_len(p)) {
                k <- clusters[j]
                m <- medoids[k]
                inertia <- inertia + D[j, m]
            }
            
            # centres = liste de K medoids (position + nom)
            centers <- vector("list", K)
            for (k in seq_len(K)) {
                m <- medoids[k]
                centers[[k]] <- list(
                    medoid_index = m,
                    medoid_name  = colnames(X)[m]
                )
            }
            
            list(
                clusters  = clusters,
                centers   = centers,
                inertia   = inertia,
                converged = converged
            )
        },
        
        # =====================
        # 2. PREDICT UNE VAR
        # =====================
        predict_one_variable = function(x_new, var_name) {
            X <- private$FX_active
            centers <- private$FCenters
            if (is.null(X) || is.null(centers)) {
                stop("[mmrClustVarKMedoids] fit() doit être appelé avant predict().")
            }
            
            n <- nrow(X)
            p <- ncol(X)
            
            is_num <- vapply(X, is.numeric, logical(1L))
            is_cat <- (!is_num) & vapply(
                X,
                function(col) is.factor(col) || is.character(col),
                logical(1L)
            )
            
            # helpers distance cohérents avec D
            dist_pair <- function(x1, x2, type1, type2) {
                if (type1 == "num" && type2 == "num") {
                    r <- suppressWarnings(
                        stats::cor(x1, x2, use = "pairwise.complete.obs")
                    )
                    if (is.na(r)) r <- 0
                    return(1 - r^2)
                } else if (type1 == "cat" && type2 == "cat") {
                    x1c <- as.character(x1)
                    x2c <- as.character(x2)
                    mismatch <- (x1c != x2c)
                    d <- mean(mismatch, na.rm = TRUE)
                    if (is.na(d)) d <- 1
                    return(d)
                } else {
                    return(1)
                }
            }
            
            type_new <- if (is.numeric(x_new)) "num" else "cat"
            
            # distances à chaque medoid
            K <- length(centers)
            d_all <- numeric(K)
            
            for (k in seq_len(K)) {
                m <- centers[[k]]$medoid_index
                xm <- X[, m]
                type_m <- if (is_num[m]) "num" else if (is_cat[m]) "cat" else "autre"
                d_all[k] <- dist_pair(x_new, xm, type_new, type_m)
            }
            
            k_best   <- which.min(d_all)
            d_best   <- d_all[k_best]
            adhesion <- 1 - d_best
            
            data.frame(
                variable = var_name,
                cluster  = as.integer(k_best),
                distance = d_best,
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
                cat("(k-medoids) Aucune variable active stockée.\n")
                return(invisible(NULL))
            }
            
            n <- nrow(X)
            p <- ncol(X)
            clusters <- private$FClusters
            centers  <- private$FCenters
            
            is_num <- vapply(X, is.numeric, logical(1L))
            is_cat <- (!is_num) & vapply(
                X,
                function(col) is.factor(col) || is.character(col),
                logical(1L)
            )
            
            # helper distance pairwise
            dist_pair <- function(x1, x2, type1, type2) {
                if (type1 == "num" && type2 == "num") {
                    r <- suppressWarnings(
                        stats::cor(x1, x2, use = "pairwise.complete.obs")
                    )
                    if (is.na(r)) r <- 0
                    return(1 - r^2)
                } else if (type1 == "cat" && type2 == "cat") {
                    x1c <- as.character(x1)
                    x2c <- as.character(x2)
                    mismatch <- (x1c != x2c)
                    d <- mean(mismatch, na.rm = TRUE)
                    if (is.na(d)) d <- 1
                    return(d)
                } else {
                    return(1)
                }
            }
            
            distance <- numeric(p)
            adhesion <- numeric(p)
            type_var <- character(p)
            
            # retrouver pour chaque cluster son medoid_index
            medoid_index <- vapply(centers, function(c) c$medoid_index, integer(1L))
            
            for (j in seq_len(p)) {
                k <- clusters[j]
                m <- medoid_index[k]
                
                xj <- X[, j]
                xm <- X[, m]
                
                type_j <- if (is_num[j]) "num" else if (is_cat[j]) "cat" else "autre"
                type_m <- if (is_num[m]) "num" else if (is_cat[m]) "cat" else "autre"
                
                d <- dist_pair(xj, xm, type_j, type_m)
                distance[j] <- d
                adhesion[j] <- 1 - d
                type_var[j] <- if (is_num[j]) "numérique" else if (is_cat[j]) "qualitative" else "autre"
            }
            
            df <- data.frame(
                variable = colnames(X),
                type     = type_var,
                cluster  = as.integer(clusters),
                distance = distance,
                adhesion = adhesion,
                stringsAsFactors = FALSE
            )
            
            cat("=== Indicateurs d'adhésion (k-medoids) ===\n")
            cat("distance = dissimilarité à la variable medoid du cluster (dans [0,1])\n")
            cat("adhesion = 1 - distance\n\n")
            
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
                stop("[mmrClustVarKMedoids] plot(type='membership') : aucun X actif.")
            }
            
            n <- nrow(X)
            p <- ncol(X)
            clusters <- private$FClusters
            centers  <- private$FCenters
            
            is_num <- vapply(X, is.numeric, logical(1L))
            is_cat <- (!is_num) & vapply(
                X,
                function(col) is.factor(col) || is.character(col),
                logical(1L)
            )
            
            dist_pair <- function(x1, x2, type1, type2) {
                if (type1 == "num" && type2 == "num") {
                    r <- suppressWarnings(
                        stats::cor(x1, x2, use = "pairwise.complete.obs")
                    )
                    if (is.na(r)) r <- 0
                    return(1 - r^2)
                } else if (type1 == "cat" && type2 == "cat") {
                    x1c <- as.character(x1)
                    x2c <- as.character(x2)
                    mismatch <- (x1c != x2c)
                    d <- mean(mismatch, na.rm = TRUE)
                    if (is.na(d)) d <- 1
                    return(d)
                } else {
                    return(1)
                }
            }
            
            medoid_index <- vapply(centers, function(c) c$medoid_index, integer(1L))
            
            membership <- numeric(p)
            for (j in seq_len(p)) {
                k <- clusters[j]
                m <- medoid_index[k]
                
                xj <- X[, j]
                xm <- X[, m]
                
                type_j <- if (is_num[j]) "num" else if (is_cat[j]) "cat" else "autre"
                type_m <- if (is_num[m]) "num" else if (is_cat[m]) "cat" else "autre"
                
                d <- dist_pair(xj, xm, type_j, type_m)
                membership[j] <- 1 - d
            }
            
            o <- order(membership, decreasing = TRUE)
            barplot(
                membership[o],
                names.arg = colnames(X)[o],
                las = 2,
                main = "Adhésion des variables aux clusters (k-medoids)",
                ylab = "adhésion (1 - distance au medoid)",
                cex.names = 0.7
            )
        },
        
        # ===========================
        # 5. PLOT : profils
        # ===========================
        plot_profiles = function() {
            X <- private$FX_active
            if (is.null(X)) {
                stop("[mmrClustVarKMedoids] plot(type = 'profiles') : aucun X actif.")
            }
            
            clusters <- private$FClusters
            if (is.null(clusters)) {
                stop("[mmrClustVarKMedoids] plot(type = 'profiles') : aucune partition disponible.")
            }
            
            K <- private$FNbGroupes
            num_idx <- private$FNumCols
            cat_idx <- private$FCatCols
            
            has_num <- length(num_idx) > 0L
            has_cat <- length(cat_idx) > 0L
            
            op <- par(no.readonly = TRUE)
            on.exit(par(op))
            if (has_num && has_cat) {
                par(mfrow = c(1, 2))
            } else {
                par(mfrow = c(1, 1))
            }
            
            # --- PARTIE NUMERIQUE : profils moyens par cluster ---
            if (has_num) {
                X_num <- as.matrix(X[, num_idx, drop = FALSE])
                
                # Standardisation colonne par colonne
                X_num_std <- scale(X_num)
                
                n <- nrow(X_num_std)
                prof_num <- matrix(NA_real_, nrow = n, ncol = K)
                colnames(prof_num) <- paste0("Cluster ", seq_len(K))
                rownames(prof_num) <- seq_len(n)
                
                # pour chaque cluster : moyenne des variables numériques du cluster
                for (k in seq_len(K)) {
                    vars_k <- which(clusters == k & seq_len(ncol(X)) %in% num_idx)
                    if (length(vars_k) == 0L) next
                    
                    # indices locaux dans X_num_std
                    cols_k_num <- match(vars_k, num_idx)
                    prof_num[, k] <- rowMeans(
                        X_num_std[, cols_k_num, drop = FALSE],
                        na.rm = TRUE
                    )
                }
                
                keep <- which(colSums(!is.na(prof_num)) > 0L)
                if (length(keep) > 0L) {
                    image(
                        x = seq_len(n),
                        y = seq_along(keep),
                        z = t(prof_num[, keep, drop = FALSE]),
                        xlab = "Individus",
                        ylab = "Clusters (partie numérique)",
                        main = "Profils moyens (k-medoids, variables numériques)",
                        axes = FALSE
                    )
                    axis(1, at = pretty(seq_len(n)))
                    axis(2,
                         at = seq_along(keep),
                         labels = colnames(prof_num)[keep])
                    box()
                } else {
                    plot.new()
                    title("Aucun profil numérique (k-medoids)")
                }
            }
            
            # --- PARTIE CATEGORIELLE : proportion de matches avec la modalité majoritaire ---
            if (has_cat) {
                X_cat <- as.data.frame(
                    lapply(X[, cat_idx, drop = FALSE], as.character),
                    stringsAsFactors = FALSE
                )
                X_cat_mat <- as.matrix(X_cat)
                
                n <- nrow(X_cat_mat)
                prof_cat <- matrix(NA_real_, nrow = n, ncol = K)
                colnames(prof_cat) <- paste0("Cluster ", seq_len(K))
                rownames(prof_cat) <- seq_len(n)
                
                # pour chaque cluster, pour chaque individu :
                # proportion de variables du cluster prenant la modalité majoritaire (par variable)
                for (k in seq_len(K)) {
                    vars_k <- which(clusters == k & seq_len(ncol(X)) %in% cat_idx)
                    if (length(vars_k) == 0L) next
                    
                    cols_k_cat <- match(vars_k, cat_idx)
                    
                    mat_k <- X_cat_mat[, cols_k_cat, drop = FALSE]
                    # mode par variable
                    modes <- apply(mat_k, 2, function(col) {
                        col_no_na <- col[!is.na(col)]
                        if (length(col_no_na) == 0L) return(NA_character_)
                        tab <- table(col_no_na)
                        names(tab)[which.max(tab)]
                    })
                    
                    # matrice binaire : 1 si modalité == mode de la variable, 0 sinon
                    bin_mat <- mapply(
                        function(col, m) as.numeric(col == m),
                        as.data.frame(mat_k, stringsAsFactors = FALSE),
                        modes,
                        SIMPLIFY = FALSE
                    )
                    bin_mat <- as.matrix(as.data.frame(bin_mat, stringsAsFactors = FALSE))
                    
                    # profil : proportion moyenne de matches par individu
                    prof_cat[, k] <- rowMeans(bin_mat, na.rm = TRUE)
                }
                
                keep <- which(colSums(!is.na(prof_cat)) > 0L)
                if (length(keep) > 0L) {
                    image(
                        x = seq_len(n),
                        y = seq_along(keep),
                        z = t(prof_cat[, keep, drop = FALSE]),
                        xlab = "Individus",
                        ylab = "Clusters (partie catégorielle)",
                        main = "Profils moyens (k-medoids, variables qualitatives)\n(proportion de matches avec la modalité majoritaire)",
                        axes = FALSE
                    )
                    axis(1, at = pretty(seq_len(n)))
                    axis(2,
                         at = seq_along(keep),
                         labels = colnames(prof_cat)[keep])
                    box()
                } else {
                    plot.new()
                    title("Aucun profil catégoriel (k-medoids)")
                }
            }
            
            invisible(NULL)
        }
    )
)
