mmrClustVarKMedoids <- R6::R6Class(
    "mmrClustVarKMedoids",
    inherit = mmrClustVarBase,
    
    public = list(
        initialize = function(K, scale = TRUE, lambda = 1, ...) {
            # lambda is kept in the API for consistency, but barely used here
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
        # 1. VARIABLE K-MEDOIDS ALGORITHM
        # ==========================
        run_clustering = function(X) {
            n <- nrow(X)
            p <- ncol(X)
            K <- private$FNbGroupes
            
            if (K > p) {
                stop("[mmrClustVarKMedoids] K cannot exceed the number of variables.")
            }
            
            # --- types ---
            is_num <- vapply(X, is.numeric, logical(1L))
            is_cat <- (!is_num) & vapply(
                X,
                function(col) is.factor(col) || is.character(col),
                logical(1L)
            )
            
            if (!any(is_num) && !any(is_cat)) {
                stop("[mmrClustVarKMedoids] No usable variables (numeric or categorical).")
            }
            
            # character version for categorical variables
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
            
            # --- helper: distance between two variables X_j, X_l ---
            dist_var <- function(j, l) {
                # Decide types
                if (is_num[j] && is_num[l]) {
                    # numeric vs numeric: 1 - r^2
                    xj <- X_num[, match(j, idx_num)]
                    xl <- X_num[, match(l, idx_num)]
                    r <- suppressWarnings(
                        stats::cor(xj, xl, use = "pairwise.complete.obs")
                    )
                    if (is.na(r)) r <- 0
                    return(1 - r^2)
                } else if (is_cat[j] && is_cat[l]) {
                    # categorical vs categorical: simple matching
                    xj <- X_cat[, match(j, idx_cat)]
                    xl <- X_cat[, match(l, idx_cat)]
                    mismatch <- (xj != xl)
                    d <- mean(mismatch, na.rm = TRUE)
                    if (is.na(d)) d <- 1
                    return(d)
                } else {
                    # mixed types: maximal distance
                    return(1)
                }
            }
            
            # --- dissimilarity matrix D ---
            D <- matrix(0, nrow = p, ncol = p)
            for (j in seq_len(p)) {
                for (l in seq_len(j - 1L)) {
                    djl <- dist_var(j, l)
                    D[j, l] <- djl
                    D[l, j] <- djl
                }
            }
            
            max_iter  <- 50L
            converged <- FALSE
            
            # =====================
            # Medoid initialization
            # =====================
            set.seed(123)
            medoids <- sample(seq_len(p), size = K, replace = FALSE)
            
            # initial assignments
            clusters <- apply(D[, medoids, drop = FALSE], 1, which.min)
            
            # avoid empty clusters
            for (k in seq_len(K)) {
                if (!any(clusters == k)) {
                    j_free <- which.max(tabulate(clusters))
                    clusters[j_free] <- k
                    medoids[k] <- j_free
                }
            }
            
            # --- main loop ---
            for (iter in seq_len(max_iter)) {
                
                # 1) update medoids
                new_medoids <- medoids
                for (k in seq_len(K)) {
                    Gk <- which(clusters == k)
                    if (length(Gk) == 0L) {
                        # empty cluster → take a variable from the largest cluster
                        j_free <- which.max(tabulate(clusters))
                        clusters[j_free] <- k
                        Gk <- which(clusters == k)
                    }
                    
                    # cost for each candidate j ∈ Gk: sum of distances D[j, l] for l ∈ Gk
                    cost <- numeric(length(Gk))
                    for (idx in seq_along(Gk)) {
                        j_candidate <- Gk[idx]
                        cost[idx] <- sum(D[j_candidate, Gk])
                    }
                    new_medoids[k] <- Gk[which.min(cost)]
                }
                
                # 2) reassign variables to updated medoids
                clusters_new <- apply(D[, new_medoids, drop = FALSE], 1, which.min)
                
                # 3) convergence?
                if (all(new_medoids == medoids) && all(clusters_new == clusters)) {
                    converged <- TRUE
                    medoids   <- new_medoids
                    clusters  <- clusters_new
                    break
                }
                
                medoids  <- new_medoids
                clusters <- clusters_new
            }
            
            # within-cluster inertia: sum of distances to the medoid of each cluster
            inertia <- 0
            for (j in seq_len(p)) {
                k <- clusters[j]
                m <- medoids[k]
                inertia <- inertia + D[j, m]
            }
            
            # centers = list of K medoids (index + name)
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
        # 2. PREDICT ONE VARIABLE
        # =====================
        predict_one_variable = function(x_new, var_name) {
            X <- private$FX_active
            centers <- private$FCenters
            if (is.null(X) || is.null(centers)) {
                stop("[mmrClustVarKMedoids] fit() must be called before predict().")
            }
            
            n <- nrow(X)
            p <- ncol(X)
            
            is_num <- vapply(X, is.numeric, logical(1L))
            is_cat <- (!is_num) & vapply(
                X,
                function(col) is.factor(col) || is.character(col),
                logical(1L)
            )
            
            # pairwise distance helper consistent with D
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
            
            # distances to each medoid
            K <- length(centers)
            d_all <- numeric(K)
            
            for (k in seq_len(K)) {
                m  <- centers[[k]]$medoid_index
                xm <- X[, m]
                type_m <- if (is_num[m]) "num" else if (is_cat[m]) "cat" else "other"
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
        # 3. SUMMARY: membership indicators
        # ===========================
        summary_membership = function() {
            X <- private$FX_active
            if (is.null(X)) {
                cat("(k-medoids) No active variables stored.\n")
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
            
            # pairwise distance helper
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
            
            # retrieve medoid_index for each cluster
            medoid_index <- vapply(centers, function(c) c$medoid_index, integer(1L))
            
            for (j in seq_len(p)) {
                k <- clusters[j]
                m <- medoid_index[k]
                
                xj <- X[, j]
                xm <- X[, m]
                
                type_j <- if (is_num[j]) "num" else if (is_cat[j]) "cat" else "other"
                type_m <- if (is_num[m]) "num" else if (is_cat[m]) "cat" else "other"
                
                d <- dist_pair(xj, xm, type_j, type_m)
                distance[j] <- d
                adhesion[j] <- 1 - d
                type_var[j] <- if (is_num[j]) "numeric" else if (is_cat[j]) "categorical" else "other"
            }
            
            df <- data.frame(
                variable = colnames(X),
                type     = type_var,
                cluster  = as.integer(clusters),
                distance = distance,
                adhesion = adhesion,
                stringsAsFactors = FALSE
            )
            
            cat("=== Membership indicators (k-medoids) ===\n")
            cat("distance = dissimilarity to the cluster medoid variable (in [0, 1])\n")
            cat("adhesion = 1 - distance\n\n")
            
            # global indicators
            dist_global <- mean(df$distance, na.rm = TRUE)
            adh_global  <- mean(df$adhesion, na.rm = TRUE)
            
            # cluster-level stats
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
            stats_by_cluster <- as.data.frame(do.call(rbind, stats_list))
            stats_by_cluster$cluster <- as.integer(stats_by_cluster$cluster)
            
            cat(sprintf("Global mean distance : %.3f\n", dist_global))
            cat(sprintf("Global mean adhesion : %.3f\n\n", adh_global))
            
            cat("--- Cluster-level statistics ---\n")
            print(stats_by_cluster)
            
            cat("\n--- Variable-level details ---\n")
            print(df)
            
            invisible(df)
        },
        
        # ===========================
        # 4. PLOT: membership barplot
        # ===========================
        plot_membership = function() {
            X <- private$FX_active
            if (is.null(X)) {
                stop("[mmrClustVarKMedoids] plot(type = 'membership'): no active X available.")
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
                
                type_j <- if (is_num[j]) "num" else if (is_cat[j]) "cat" else "other"
                type_m <- if (is_num[m]) "num" else if (is_cat[m]) "cat" else "other"
                
                d <- dist_pair(xj, xm, type_j, type_m)
                membership[j] <- 1 - d
            }
            
            o <- order(membership, decreasing = TRUE)
            barplot(
                membership[o],
                names.arg = colnames(X)[o],
                las = 2,
                main = "Variable–Cluster Membership (k-medoids)",
                ylab = "adhesion (1 - distance to medoid)",
                cex.names = 0.7
            )
        },
        
        # ===========================
        # 5. PLOT: profiles
        # ===========================
        plot_profiles = function() {
            X <- private$FX_active
            if (is.null(X)) {
                stop("[mmrClustVarKMedoids] plot(type = 'profiles'): no active X available.")
            }
            
            clusters <- private$FClusters
            if (is.null(clusters)) {
                stop("[mmrClustVarKMedoids] plot(type = 'profiles'): no partition available.")
            }
            
            K      <- private$FNbGroupes
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
            
            # --- NUMERIC PART: mean profiles per cluster ---
            if (has_num) {
                X_num <- as.matrix(X[, num_idx, drop = FALSE])
                
                # column-wise standardisation
                X_num_std <- scale(X_num)
                
                n <- nrow(X_num_std)
                prof_num <- matrix(NA_real_, nrow = n, ncol = K)
                colnames(prof_num) <- paste0("Cluster ", seq_len(K))
                rownames(prof_num) <- seq_len(n)
                
                # for each cluster: mean of numeric variables belonging to that cluster
                for (k in seq_len(K)) {
                    vars_k <- which(clusters == k & seq_len(ncol(X)) %in% num_idx)
                    if (length(vars_k) == 0L) next
                    
                    # local indices in X_num_std
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
                        xlab = "Individuals",
                        ylab = "Clusters (numeric part)",
                        main = "Mean profiles (k-medoids, numeric variables)",
                        axes = FALSE
                    )
                    axis(1, at = pretty(seq_len(n)))
                    axis(2,
                         at = seq_along(keep),
                         labels = colnames(prof_num)[keep])
                    box()
                } else {
                    plot.new()
                    title("No numeric profiles (k-medoids)")
                }
            }
            
            # --- CATEGORICAL PART: proportion of matches with majority modality ---
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
                
                # for each cluster, for each individual:
                # proportion of cluster variables taking the majority modality (per variable)
                for (k in seq_len(K)) {
                    vars_k <- which(clusters == k & seq_len(ncol(X)) %in% cat_idx)
                    if (length(vars_k) == 0L) next
                    
                    cols_k_cat <- match(vars_k, cat_idx)
                    
                    mat_k <- X_cat_mat[, cols_k_cat, drop = FALSE]
                    # mode for each variable
                    modes <- apply(mat_k, 2, function(col) {
                        col_no_na <- col[!is.na(col)]
                        if (length(col_no_na) == 0L) return(NA_character_)
                        tab <- table(col_no_na)
                        names(tab)[which.max(tab)]
                    })
                    
                    # binary matrix: 1 if modality == mode of the variable, 0 otherwise
                    bin_mat <- mapply(
                        function(col, m) as.numeric(col == m),
                        as.data.frame(mat_k, stringsAsFactors = FALSE),
                        modes,
                        SIMPLIFY = FALSE
                    )
                    bin_mat <- as.matrix(as.data.frame(bin_mat, stringsAsFactors = FALSE))
                    
                    # profile: mean match proportion per individual
                    prof_cat[, k] <- rowMeans(bin_mat, na.rm = TRUE)
                }
                
                keep <- which(colSums(!is.na(prof_cat)) > 0L)
                if (length(keep) > 0L) {
                    image(
                        x = seq_len(n),
                        y = seq_along(keep),
                        z = t(prof_cat[, keep, drop = FALSE]),
                        xlab = "Individuals",
                        ylab = "Clusters (categorical part)",
                        main = "Mean profiles (k-medoids, categorical variables)\n(match proportion with majority modality)",
                        axes = FALSE
                    )
                    axis(1, at = pretty(seq_len(n)))
                    axis(2,
                         at = seq_along(keep),
                         labels = colnames(prof_cat)[keep])
                    box()
                } else {
                    plot.new()
                    title("No categorical profiles (k-medoids)")
                }
            }
            
            invisible(NULL)
        }
    )
)