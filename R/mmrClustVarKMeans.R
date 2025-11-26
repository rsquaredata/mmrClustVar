mmrClustVarKMeans <- R6::R6Class(
    "mmrClustVarKMeans",
    inherit = mmrClustVarBase,
    
    public = list(
        initialize = function(K, scale = TRUE, lambda = 1, ...) {
            # lambda is ignored for k-means, but kept for signature consistency
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
        # 1. VARIABLE K-MEANS ALGORITHM
        # ==========================
        run_clustering = function(X) {
            # X: numeric data.frame already validated by check_and_prepare_X()
            # and possibly standardized
            
            # Ensure all variables are numeric
            if (!all(vapply(X, is.numeric, logical(1L)))) {
                stop("[mmrClustVarKMeans] All variables must be numeric.")
            }
            
            n <- nrow(X)
            p <- ncol(X)
            K <- private$FNbGroupes
            
            if (K > p) {
                stop("[mmrClustVarKMeans] K cannot exceed the number of variables.")
            }
            
            # Convert to n × p matrix
            X_mat <- as.matrix(X)
            
            # Algorithm parameters
            max_iter <- 50L
            tol      <- 1e-6
            
            # =====================
            # Simple initialization
            # =====================
            set.seed(123)
            nuclei <- sample(seq_len(p), size = K, replace = FALSE)
            clusters <- rep(NA_integer_, p)
            clusters[nuclei] <- seq_len(K)
            
            # Random assignment for the remaining variables
            idx_other <- setdiff(seq_len(p), nuclei)
            if (length(idx_other) > 0L) {
                clusters[idx_other] <- sample(seq_len(K), size = length(idx_other), replace = TRUE)
            }
            
            # Ensure no empty cluster
            for (k in seq_len(K)) {
                if (!any(clusters == k)) {
                    j_free <- which.max(tabulate(clusters))
                    clusters[j_free] <- k
                }
            }
            
            # List of latent components Z_k (each a vector of length n)
            centers <- vector("list", K)
            
            # Internal function: compute Z_k = 1st principal component
            compute_Zk <- function(cols_k) {
                
                if (length(cols_k) == 1L) {
                    # Trivial case: single variable → centered variable
                    zk <- X_mat[, cols_k]
                    zk <- zk - mean(zk, na.rm = TRUE)
                    return(as.numeric(zk))
                    
                } else {
                    # Multiple variables → PCA
                    Xk <- X_mat[, cols_k, drop = FALSE]
                    
                    # Keep only complete rows
                    idx_ok <- stats::complete.cases(Xk)
                    
                    # If too few complete rows → fallback to centered first variable
                    if (sum(idx_ok) < 2L) {
                        zk <- Xk[, 1]
                        zk <- zk - mean(zk, na.rm = TRUE)
                        return(as.numeric(zk))
                    }
                    
                    # PCA on complete rows
                    pc <- stats::prcomp(
                        Xk[idx_ok, , drop = FALSE],
                        center  = FALSE,
                        scale.  = FALSE
                    )
                    
                    z_short <- pc$x[, 1]          # scores on complete rows
                    zk <- rep(NA_real_, nrow(Xk)) # full-length vector
                    zk[idx_ok] <- z_short
                    
                    return(as.numeric(zk))
                }
            }
            
            # r^2(X_j, Z_k)
            r2_corr <- function(x, z) {
                r <- suppressWarnings(stats::cor(x, z, use = "pairwise.complete.obs"))
                if (is.na(r)) r <- 0
                r^2
            }
            
            # Objective function: W = sum_j r^2(X_j, Z_cluster(j))
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
            # Reallocation loop
            # =============================
            old_W    <- -Inf
            converged <- FALSE
            
            for (iter in seq_len(max_iter)) {
                
                # 1) Recompute latent components Z_k
                for (k in seq_len(K)) {
                    cols_k <- which(clusters == k)
                    
                    if (length(cols_k) == 0L) {
                        cols_k <- sample(seq_len(p), size = 1L)
                        clusters[cols_k] <- k
                    }
                    
                    centers[[k]] <- compute_Zk(cols_k)
                }
                
                # 2) Reassign variables
                new_clusters <- clusters
                
                for (j in seq_len(p)) {
                    xj <- X_mat[, j]
                    r2_all <- vapply(
                        centers,
                        function(zk) r2_corr(xj, zk),
                        numeric(1L)
                    )
                    k_best <- which.max(r2_all)
                    new_clusters[j] <- k_best
                }
                
                # 3) Compute W and stopping criteria
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
            
            # Inertia = sum_j (1 - r^2)
            inertia <- 0
            for (j in seq_len(p)) {
                k  <- clusters[j]
                zk <- centers[[k]]
                xj <- X_mat[, j]
                inertia <- inertia + (1 - r2_corr(xj, zk))
            }
            
            list(
                clusters  = clusters,
                centers   = centers,   # list of Z_k vectors
                inertia   = inertia,
                converged = converged
            )
        },
        
        # =====================
        # 2. PREDICT ONE VARIABLE
        # =====================
        predict_one_variable = function(x_new, var_name) {
            
            if (!is.numeric(x_new)) {
                stop("[mmrClustVarKMeans] predict() for k-means requires a numeric variable.")
            }
            
            centers <- private$FCenters
            if (is.null(centers)) {
                stop("[mmrClustVarKMeans] No prototypes available (was fit() called?).")
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
        # 3. SUMMARY: membership indicators
        # ===========================
        summary_membership = function() {
            X <- private$FX_active
            if (is.null(X)) {
                cat("(k-means) No active variables stored.\n")
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
            
            cat("=== Membership Indicators (k-means) ===\n")
            cat("r^2 = squared correlation between variable and cluster latent component\n\n")
            
            # Explained inertia (mean r^2)
            explained <- mean(df$r2, na.rm = TRUE)
            
            # Stats by cluster
            stats_list <- lapply(split(df, df$cluster), function(dsub) {
                c(
                    cluster = dsub$cluster[1],
                    r2_mean = mean(dsub$r2, na.rm = TRUE),
                    r2_min  = min(dsub$r2, na.rm = TRUE),
                    r2_max  = max(dsub$r2, na.rm = TRUE)
                )
            })
            stats_by_cluster <- as.data.frame(do.call(rbind, stats_list))
            stats_by_cluster$cluster <- as.integer(stats_by_cluster$cluster)
            
            cat(sprintf("Explained inertia (mean r^2): %.3f\n\n", explained))
            
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
                stop("[mmrClustVarKMeans] plot(type = 'membership') : no active X available.")
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
                main = "Variable–Cluster Membership (k-means)",
                ylab = "r^2",
                cex.names = 0.7
            )
        },
        
        # ===========================
        # 5. PLOT: mean profiles heatmap
        # ===========================
        plot_profiles = function() {
            X <- private$FX_active
            if (is.null(X)) {
                warning("[mmrClustVarKMeans] plot(type = 'profiles') : no active X available.")
                return(invisible(NULL))
            }
            
            X_mat    <- as.matrix(X)
            n        <- nrow(X_mat)
            clusters <- private$FClusters
            K        <- private$FNbGroupes
            
            if (is.null(clusters)) {
                warning("[mmrClustVarKMeans] No clusters available.")
                return(invisible(NULL))
            }
            
            # Average profile per individual per cluster
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
                ylab = "Individuals",
                main = "Mean Profiles by Cluster (k-means)",
                axes = FALSE
            )
            axis(1, at = seq_len(K), labels = colnames(prof_mat))
            axis(2, at = pretty(seq_len(n)))
            box()
            
            invisible(NULL)
        }
    )
)