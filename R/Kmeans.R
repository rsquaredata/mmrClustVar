#' K-means for variable clustering (numeric variables only)
#'
#' R6 class implementing a k-means-like algorithm for clustering
#' \emph{variables} (columns) when all active variables are numeric.
#'
#' @docType class
#' @name Kmeans
#' @export
Kmeans <- R6::R6Class(
  "Kmeans",
  inherit = .ClusterBase,

  public = list(
    
    #' @description
    #' Create a new Kmeans instance
    #' @param K Number of cluster
    #' @param scale A boolean defining whether to scale the data or not
    #' @param random_state The random seed to use
    initialize = function(K, scale = TRUE, random_state = NULL) {
      super$initialize(
        K           = K,
        scale       = scale,
        method_name = "kmeans",
        random_state = random_state
      )
    },

    #' @description
    #' Textual interpretation of the k-means variable clustering solution.
    #'
    #' For each cluster, this method reports:
    #'   - the number of variables (all numeric),
    #'   - a few most representative variables (highest r^2 to the local PC),
    #'   - simple membership statistics based on r^2.
    #'
    #' The "distance" is defined as 1 - r^2(x_j, PC_k),
    #' and the "adhesion" as r^2(x_j, PC_k).
    #'
    #' @param style Character, either "compact" or "detailed".
    interpret_clusters = function(style = c("compact", "detailed")) {
      style <- match.arg(style)

      X <- private$FX_active
      clusters <- private$FClusters

      if (is.null(X) || is.null(clusters)) {
        stop("[Kmeans] interpret_clusters(): no fitted model or missing state.")
      }

      p <- ncol(X)
      K <- private$FNbGroupes

      # type check
      is_num <- vapply(X, is.numeric, logical(1L))
      if (!all(is_num)) {
        warning("[Kmeans] Some non-numeric variables detected; ",
                "only numeric variables are used for interpretation.")
      }

      distance <- rep(NA_real_, p)
      adhesion <- rep(NA_real_, p)
      type_var <- ifelse(is_num, "numeric", "other")

      # work on numeric part only
      X_num  <- as.matrix(X[, is_num, drop = FALSE])
      idx_num <- which(is_num)

      # for each cluster, compute a local PC1, then r^2 to that PC
      for (k in seq_len(K)) {
        vars_k <- which(clusters == k & is_num)
        if (length(vars_k) == 0L) next

        # local numeric matrix for this cluster
        cols_k_num <- match(vars_k, idx_num)
        Xk <- X_num[, cols_k_num, drop = FALSE]

        # standardise by column
        Xk_std <- scale(Xk)

        # local PCA
        pc <- stats::prcomp(Xk_std, center = FALSE, scale. = FALSE)
        score1 <- pc$x[, 1]

        # compute r^2(x_j, PC1_k) for variables in this cluster
        for (j in vars_k) {
          xj <- scale(X[, j])
          r2 <- private$r2_corr(as.numeric(xj), score1)
          distance[j] <- 1 - r2
          adhesion[j] <- r2
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

      n_num <- sum(type_var == "numeric")

      cat("=== Global overview (k-means, variable clustering) ===\n")
      cat("Number of clusters :", K, "\n")
      cat("Number of variables:", p,
          sprintf("(%d numeric, %d other)\n", n_num, p - n_num))
      cat("\n")

      # per-cluster interpretation
      for (k in seq_len(K)) {
        df_k <- df[df$cluster == k & df$type == "numeric", , drop = FALSE]

        cat(sprintf("--- Cluster %d ---\n", k))

        if (nrow(df_k) == 0L) {
          cat("Cluster contains no numeric variables.\n\n")
          next
        }

        size_k <- nrow(df_k)
        cat(sprintf("Size: %d numeric variables.\n", size_k))

        adh_mean_k <- mean(df_k$adhesion, na.rm = TRUE)
        adh_min_k  <- min(df_k$adhesion, na.rm = TRUE)
        adh_max_k  <- max(df_k$adhesion, na.rm = TRUE)

        # top variables by adhesion
        o_k   <- order(df_k$adhesion, decreasing = TRUE)
        top_k <- df_k[o_k, , drop = FALSE]
        if (nrow(top_k) > 3L) {
          top_k <- top_k[seq_len(3L), , drop = FALSE]
        }

        if (style == "compact") {
          cat("Most representative variables (top by r^2 to the local PC):\n")
          for (i in seq_len(nrow(top_k))) {
            cat(sprintf("  - %s (adhesion r^2 = %.3f)\n",
                        top_k$variable[i],
                        top_k$adhesion[i]))
          }
          cat("Membership statistics (based on r^2):\n")
          cat(sprintf("  - Mean adhesion: %.3f\n", adh_mean_k))
          cat(sprintf("  - Min / Max adhesion: %.3f / %.3f\n",
                      adh_min_k, adh_max_k))
          cat("Interpretation:\n")
          cat("  This cluster groups numeric variables that share a strong\n")
          cat("  common latent component (local principal component).\n\n")

        } else {  # detailed
          cat("Most representative variables (top 3 by r^2 to the local PC):\n")
          for (i in seq_len(nrow(top_k))) {
            cat(sprintf("  - %s (adhesion r^2 = %.3f)\n",
                        top_k$variable[i],
                        top_k$adhesion[i]))
          }
          cat("Membership statistics (based on r^2):\n")
          cat(sprintf("  - Mean adhesion: %.3f\n", adh_mean_k))
          cat(sprintf("  - Min / Max adhesion: %.3f / %.3f\n",
                      adh_min_k, adh_max_k))
          cat("Interpretation:\n")
          cat("  This cluster is characterised by a local principal component\n")
          cat("  that summarises the behaviour of its variables; those with\n")
          cat("  high r^2 are the most typical members of the group.\n\n")
        }
      }

      invisible(df)
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
        stop("[Kmeans] All variables must be numeric.")
      }

      n <- nrow(X)
      p <- ncol(X)
      K <- private$FNbGroupes

      if (K > p) {
        stop("[Kmeans] K cannot exceed the number of variables.")
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

      # Objective function: W = sum_j r^2(X_j, Z_cluster(j))
      compute_W <- function(clusters, centers) {
        W <- 0
        for (j in seq_len(p)) {
          k  <- clusters[j]
          zk <- centers[[k]]
          xj <- X_mat[, j]
          W <- W + private$r2_corr(xj, zk)
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
            function(zk) private$r2_corr(xj, zk),
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
        inertia <- inertia + (1 - private$r2_corr(xj, zk))
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
        stop("[Kmeans] predict() for k-means requires a numeric variable.")
      }

      centers <- private$FCenters
      if (is.null(centers)) {
        stop("[Kmeans] No prototypes available (was fit() called?).")
      }

      K <- length(centers)
      r2_all <- vapply(
        centers,
        function(zk) private$r2_corr(x_new, zk),
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
    summary_membership_impl = function() {
      X <- private$FX_active
      if (is.null(X)) {
        cat("(k-means) No active variables stored.\n")
        return(invisible(NULL))
      }

      X_mat    <- as.matrix(X)
      p        <- ncol(X_mat)
      clusters <- private$FClusters
      centers  <- private$FCenters

      if (is.null(clusters) || is.null(centers)) {
        cat("(k-means) No clusters or centers available.\n")
        return(invisible(NULL))
      }

      r2_vec   <- numeric(p)
      dist_vec <- numeric(p)

      for (j in seq_len(p)) {
        k  <- clusters[j]
        zk <- centers[[k]]
        xj <- X_mat[, j]
        r2 <- private$r2_corr(xj, zk)
        r2_vec[j]   <- r2
        dist_vec[j] <- 1 - r2
      }

      df <- data.frame(
        variable = colnames(X_mat),
        cluster  = as.integer(clusters),
        r2       = r2_vec,
        distance = dist_vec,
        stringsAsFactors = FALSE
      )

      cat("=== Membership indicators (k-means) ===\n")
      cat("r^2      = squared correlation between variable and cluster latent component\n")
      cat("distance = 1 - r^2 (within-cluster dissimilarity)\n\n")

      ## ---- Global indicators ----
      dist_global <- mean(df$distance, na.rm = TRUE)
      r2_global   <- mean(df$r2,       na.rm = TRUE)

      explained_inertia_global <- 1 - dist_global  # == r2_global

      cat(sprintf("Global mean distance              : %.3f\n", dist_global))
      cat(sprintf("Global mean r^2 (adhesion)       : %.3f\n", r2_global))
      cat(sprintf("Explained inertia (global, ratio) : %.3f\n", explained_inertia_global))
      cat(sprintf("Explained inertia (global, %% )    : %.1f%%\n",
                  100 * explained_inertia_global))
      cat("\n")

      ## ---- Cluster-level statistics ----
      stats_list <- lapply(split(df, df$cluster), function(dsub) {
        dist_mean <- mean(dsub$distance, na.rm = TRUE)
        r2_mean   <- mean(dsub$r2,       na.rm = TRUE)

        c(
          cluster              = dsub$cluster[1],
          dist_mean            = dist_mean,
          dist_min             = min(dsub$distance, na.rm = TRUE),
          dist_max             = max(dsub$distance, na.rm = TRUE),
          r2_mean              = r2_mean,
          r2_min               = min(dsub$r2,       na.rm = TRUE),
          r2_max               = max(dsub$r2,       na.rm = TRUE),
          explained_inertia_cl = 1 - dist_mean      # proportion of explained inertia for this cluster
        )
      })

      stats_by_cluster <- as.data.frame(do.call(rbind, stats_list))
      stats_by_cluster$cluster <- as.integer(stats_by_cluster$cluster)

      cat("--- Cluster-level statistics ---\n")
      print(stats_by_cluster)

      cat("\n--- Variable-level details ---\n")
      print(df)

      invisible(df)
    },

    # ===========================
    # 4. PLOT: membership barplot
    # ===========================
    plot_membership_impl = function() {
      X <- private$FX_active
      if (is.null(X)) {
        stop("[Kmeans] plot(type = 'membership') : no active X available.")
      }

      X_mat    <- as.matrix(X)
      p        <- ncol(X_mat)
      clusters <- private$FClusters
      centers  <- private$FCenters

      membership <- numeric(p)
      for (j in seq_len(p)) {
        k  <- clusters[j]
        zk <- centers[[k]]
        xj <- X_mat[, j]
        membership[j] <- private$r2_corr(xj, zk)
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
    plot_profiles_impl = function() {
      X <- private$FX_active
      if (is.null(X)) {
        warning("[Kmeans] plot(type = 'profiles') : no active X available.")
        return(invisible(NULL))
      }

      X_mat    <- as.matrix(X)
      n        <- nrow(X_mat)
      clusters <- private$FClusters
      K        <- private$FNbGroupes

      if (is.null(clusters)) {
        warning("[Kmeans] No clusters available.")
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
