mmrClustVarKModes <- R6::R6Class(
    "mmrClustVarKModes",
    inherit = mmrClustVarBase,
    
    public = list(
        initialize = function(K, scale = TRUE, lambda = 1, ...) {
            # scale and lambda are ignored for k-modes, but kept for API consistency
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
        # 1. VARIABLE K-MODES ALGORITHM
        # ==========================
        run_clustering = function(X) {
            # X: categorical data.frame (factor/character), already validated
            
            # --- 1) Check all variables are categorical ---
            is_cat <- vapply(
                X,
                function(col) is.factor(col) || is.character(col),
                logical(1L)
            )
            if (!all(is_cat)) {
                stop("[mmrClustVarKModes] All variables must be categorical.")
            }
            
            # Work with a character matrix for simplicity
            X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)
            X_mat  <- as.matrix(X_char)
            
            n <- nrow(X_mat)
            p <- ncol(X_mat)
            K <- private$FNbGroupes
            
            if (K > p) {
                stop("[mmrClustVarKModes] K cannot exceed the number of variables.")
            }
            
            max_iter  <- 50L
            converged <- FALSE
            
            # =====================
            # Simple initialization
            # =====================
            set.seed(123)
            seeds <- sample(seq_len(p), size = K, replace = FALSE)
            clusters <- rep(NA_integer_, p)
            clusters[seeds] <- seq_len(K)
            
            others <- setdiff(seq_len(p), seeds)
            if (length(others) > 0L) {
                clusters[others] <- sample(seq_len(K), size = length(others), replace = TRUE)
            }
            
            # Ensure no empty cluster
            for (k in seq_len(K)) {
                if (!any(clusters == k)) {
                    j_free <- which.max(tabulate(clusters))
                    clusters[j_free] <- k
                }
            }
            
            # Each center is a vector of length n (mode profile per individual)
            centers <- vector("list", K)
            
            # --- helper: compute mode profile for a cluster ---
            compute_mode_profile <- function(cols_k) {
                if (length(cols_k) == 1L) {
                    return(X_mat[, cols_k])
                } else {
                    mode_vec <- character(n)
                    for (i in seq_len(n)) {
                        vals <- X_mat[i, cols_k]
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
            }
            
            # Simple matching dissimilarity: proportion of mismatches
            simple_matching <- function(x, m) {
                mismatch <- (x != m)
                d <- mean(mismatch, na.rm = TRUE)
                if (is.na(d)) d <- 1
                d
            }
            
            # ===========================
            # Main loop
            # ===========================
            for (iter in seq_len(max_iter)) {
                
                # 1) recompute cluster modes
                for (k in seq_len(K)) {
                    cols_k <- which(clusters == k)
                    if (length(cols_k) == 0L) {
                        cols_k <- sample(seq_len(p), size = 1L)
                        clusters[cols_k] <- k
                    }
                    centers[[k]] <- compute_mode_profile(cols_k)
                }
                
                # 2) reassign variables
                new_clusters <- clusters
                
                for (j in seq_len(p)) {
                    xj <- X_mat[, j]
                    d_all <- vapply(
                        centers,
                        function(mk) simple_matching(xj, mk),
                        numeric(1L)
                    )
                    k_best <- which.min(d_all)
                    new_clusters[j] <- k_best
                }
                
                # 3) convergence criterion
                if (all(new_clusters == clusters)) {
                    converged <- TRUE
                    clusters  <- new_clusters
                    break
                }
                
                clusters <- new_clusters
            }
            
            # Final inertia = sum of simple-matching dissimilarities
            inertia <- 0
            for (j in seq_len(p)) {
                k  <- clusters[j]
                mk <- centers[[k]]
                xj <- X_mat[, j]
                inertia <- inertia + simple_matching(xj, mk)
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
            
            if (!(is.factor(x_new) || is.character(x_new))) {
                stop("[mmrClustVarKModes] predict() requires a categorical variable.")
            }
            
            centers <- private$FCenters
            if (is.null(centers)) {
                stop("[mmrClustVarKModes] No prototypes available (did you run fit()?).")
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
            adhesion <- 1 - d_raw   # proportion of matches
            
            data.frame(
                variable = var_name,
                cluster  = as.integer(k_best),
                distance = d_raw,
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
                cat("(k-modes) No active variables stored.\n")
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
            
            cat("=== Membership Indicators (k-modes) ===\n")
            cat("distance  = proportion of mismatches with the cluster mode\n")
            cat("adhesion  = 1 - distance (proportion of matches)\n\n")
            
            # Global metrics
            dist_global <- mean(df$distance, na.rm = TRUE)
            adh_global  <- mean(df$adhesion, na.rm = TRUE)
            
            # Cluster-level statistics
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
            
            cat(sprintf("Global average distance : %.3f\n", dist_global))
            cat(sprintf("Global average adhesion : %.3f\n\n", adh_global))
            
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
                stop("[mmrClustVarKModes] plot(type='membership'): no active X available.")
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
                d_raw <- simple_matching(xj, mk)
                membership[j] <- 1 - d_raw
            }
            
            o <- order(membership, decreasing = TRUE)
            barplot(
                membership[o],
                names.arg = colnames(X_mat)[o],
                las = 2,
                main = "Variable–Cluster Membership (k-modes)",
                ylab = "adhesion (1 − mismatch proportion)",
                cex.names = 0.7
            )
        },
        
        # ===========================
        # 5. PLOT: mean profiles heatmap
        # ===========================
        plot_profiles = function() {
            X <- private$FX_active
            if (is.null(X)) {
                warning("[mmrClustVarKModes] plot(type='profiles'): no active X available.")
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
                warning("[mmrClustVarKModes] No clusters or prototypes available.")
                return(invisible(NULL))
            }
            
            # prof_mat[i, k] = proportion of matches with mode profile for cluster k
            prof_mat <- matrix(NA_real_, nrow = n, ncol = K)
            colnames(prof_mat) <- paste0("Cluster ", seq_len(K))
            rownames(prof_mat) <- seq_len(n)
            
            for (k in seq_len(K)) {
                vars_k <- which(clusters == k)
                if (length(vars_k) == 0L) next
                
                mk <- centers[[k]]
                
                for (i in seq_len(n)) {
                    vals <- X_mat[i, vars_k]
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
                ylab = "Individuals",
                main = "Mean Profiles by Cluster (k-modes)\n(proportion of matches to cluster mode)",
                axes = FALSE
            )
            axis(1, at = seq_len(K), labels = colnames(prof_mat))
            axis(2, at = pretty(seq_len(n)))
            box()
            
            invisible(NULL)
        }
    )
)