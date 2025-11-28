#' K-prototypes for variable clustering (mixed variables)
#'
#' R6 class implementing a k-prototypes-like algorithm for clustering
#' \emph{variables} (columns) when active variables are mixed (numeric +
#' categorical).
#'
#' @docType class
#' @name Kprototypes
#' @export
Kprototypes <- R6::R6Class(
  "Kprototypes",
  inherit = .ClusterBase,

  public = list(

    #' @description
    #' Create a new Kprototypes instance
    #' @param K Number of cluster
    #' @param scale A boolean defining whether to scale the data or not
    #' @param lambda Weight for qualitative data
    #' @param random_state The random seed to use
    initialize = function(K, scale = TRUE, lambda = 1, random_state = NULL) {
      # lambda > 0: weight for the categorical part
      if (!is.numeric(lambda) || length(lambda) != 1L || lambda <= 0) {
        stop("[Kprototypes] lambda must be a numeric > 0.")
      }

      super$initialize(
        K           = K,
        scale       = scale,
        lambda      = lambda,
        method_name = "kprototypes",
        random_state = random_state
      )
    },

    #' @description
    #' Textual interpretation of the k-prototypes variable clustering solution.
    #'
    #' For each cluster, this method reports:
    #'   - the number of numeric and categorical variables,
    #'   - a few most representative variables (highest adhesion to the
    #'     numeric and/or categorical prototype),
    #'   - simple membership statistics.
    #'
    #' For numeric variables, adhesion is based on r^2 to a local PC (as in k-means).
    #' For categorical variables, adhesion is based on agreement with a cluster-level
    #' modal profile (as in k-modes).
    #'
    #' @param style Character, either "compact" or "detailed".
    interpret_clusters = function(style = c("compact", "detailed")) {
      style <- match.arg(style)

      X        <- private$FX_active
      clusters <- private$FClusters

      if (is.null(X) || is.null(clusters)) {
        stop("[Kprototypes] interpret_clusters(): no fitted model or missing state.")
      }

      p <- ncol(X)
      K <- private$FNbGroupes

      is_num <- vapply(X, is.numeric, logical(1L))
      is_cat <- (!is_num) & vapply(
        X,
        function(col) is.factor(col) || is.character(col),
        logical(1L)
      )

      type_var <- ifelse(is_num, "numeric",
                         ifelse(is_cat, "categorical", "other"))

      distance <- rep(NA_real_, p)
      adhesion <- rep(NA_real_, p)

      ## --- Numeric part (k-means-like, r^2 to local PC) ---
      if (any(is_num)) {
        X_num   <- as.matrix(X[, is_num, drop = FALSE])
        idx_num <- which(is_num)

        for (k in seq_len(K)) {
          vars_k <- which(clusters == k & is_num)
          if (length(vars_k) == 0L) next

          cols_k_num <- match(vars_k, idx_num)
          Xk <- X_num[, cols_k_num, drop = FALSE]

          Xk_std <- scale(Xk)
          pc     <- stats::prcomp(Xk_std, center = FALSE, scale. = FALSE)
          score1 <- pc$x[, 1]

          for (j in vars_k) {
            xj <- scale(X[, j])
            r2 <- private$r2_corr(as.numeric(xj), as.numeric(score1))
            adhesion[j] <- r2
            distance[j] <- 1 - r2
          }
        }
      }

      ## --- Categorical part (k-modes-like, agreement with modal profile) ---
      if (any(is_cat)) {
        X_cat <- as.data.frame(
          lapply(X[, is_cat, drop = FALSE], as.character),
          stringsAsFactors = FALSE
        )
        X_cat_mat     <- as.matrix(X_cat)
        idx_cat_global <- which(is_cat)
        n <- nrow(X_cat_mat)

        for (k in seq_len(K)) {
          vars_k <- which(clusters == k & is_cat)
          if (length(vars_k) == 0L) next

          cols_k_cat <- match(vars_k, idx_cat_global)
          mat_k      <- X_cat_mat[, cols_k_cat, drop = FALSE]

          # cluster-level modal profile by individual
          modal_profile <- apply(mat_k, 1, function(row) {
            row_no_na <- row[!is.na(row)]
            if (length(row_no_na) == 0L) return(NA_character_)
            tab <- table(row_no_na)
            names(tab)[which.max(tab)]
          })

          for (j in vars_k) {
            cpos <- match(j, idx_cat_global)
            if (is.na(cpos)) next
            col_j <- X_cat_mat[, cpos]
            agree <- (col_j == modal_profile)
            prop_agree <- mean(agree, na.rm = TRUE)
            if (is.na(prop_agree)) prop_agree <- 0

            adhesion[j] <- prop_agree
            distance[j] <- 1 - prop_agree
          }
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
      n_cat <- sum(type_var == "categorical")

      cat("=== Global overview (k-prototypes, variable clustering) ===\n")
      cat("Number of clusters :", K, "\n")
      cat("Number of variables:", p,
          sprintf("(%d numeric, %d categorical, %d other)\n",
                  n_num, n_cat, p - n_num - n_cat))
      cat("\n")

      for (k in seq_len(K)) {
        df_k <- df[df$cluster == k & df$type != "other", , drop = FALSE]

        cat(sprintf("--- Cluster %d ---\n", k))

        if (nrow(df_k) == 0L) {
          cat("Cluster contains no numeric or categorical variables.\n\n")
          next
        }

        size_k  <- nrow(df_k)
        n_num_k <- sum(df_k$type == "numeric")
        n_cat_k <- sum(df_k$type == "categorical")

        cat(sprintf("Size: %d variables (%d numeric, %d categorical).\n",
                    size_k, n_num_k, n_cat_k))

        adh_mean_k <- mean(df_k$adhesion, na.rm = TRUE)
        adh_min_k  <- min(df_k$adhesion, na.rm = TRUE)
        adh_max_k  <- max(df_k$adhesion, na.rm = TRUE)

        o_k   <- order(df_k$adhesion, decreasing = TRUE)
        top_k <- df_k[o_k, , drop = FALSE]
        if (nrow(top_k) > 3L) {
          top_k <- top_k[seq_len(3L), , drop = FALSE]
        }

        if (style == "compact") {
          cat("Most representative variables (top by adhesion to the numeric/categorical prototype):\n")
          for (i in seq_len(nrow(top_k))) {
            cat(sprintf("  - %s (%s, adhesion = %.3f)\n",
                        top_k$variable[i],
                        top_k$type[i],
                        top_k$adhesion[i]))
          }
          cat("Membership statistics:\n")
          cat(sprintf("  - Mean adhesion: %.3f\n", adh_mean_k))
          cat(sprintf("  - Min / Max adhesion: %.3f / %.3f\n",
                      adh_min_k, adh_max_k))
          cat("Interpretation:\n")
          cat("  This cluster groups variables that share both a common\n")
          cat("  numeric pattern (through a local PC) and, when relevant,\n")
          cat("  a common categorical pattern (through a modal profile).\n\n")

        } else {  # detailed
          cat("Most representative variables (top 3 by adhesion to the prototype):\n")
          for (i in seq_len(nrow(top_k))) {
            cat(sprintf("  - %s (%s, adhesion = %.3f)\n",
                        top_k$variable[i],
                        top_k$type[i],
                        top_k$adhesion[i]))
          }
          cat("Membership statistics:\n")
          cat(sprintf("  - Mean adhesion: %.3f\n", adh_mean_k))
          cat(sprintf("  - Min / Max adhesion: %.3f / %.3f\n",
                      adh_min_k, adh_max_k))
          cat("Interpretation:\n")
          cat("  This cluster is driven by a mixed prototype combining\n")
          cat("  a numeric component and a categorical profile; variables\n")
          cat("  with high adhesion follow one or both of these aspects closely.\n\n")
        }
      }

      invisible(df)
    }
  ),

  private = list(

    # ================================
    # 1. VARIABLE K-PROTOTYPES ALGORITHM
    # ================================
    # Prototypes representation:
    # centers[[k]] is a list with:
    #   $num : numeric vector of length n (latent numeric profile, k-means-like)
    #   $cat : character vector of length n (categorical profile, mode per individual)
    #
    run_clustering = function(X) {
      # X: mixed data.frame (numeric + factor/character),
      # already checked by check_and_prepare_X() and possibly scaled.

      n <- nrow(X)
      p <- ncol(X)
      K <- private$FNbGroupes

      if (K > p) {
        stop("[Kprototypes] K cannot exceed the number of variables.")
      }

      num_idx <- private$FNumCols
      cat_idx <- private$FCatCols

      # numeric matrix (possibly empty)
      X_num <- if (length(num_idx) > 0L) {
        as.matrix(X[, num_idx, drop = FALSE])
      } else {
        NULL
      }

      # categorical matrix (possibly empty)
      X_cat <- if (length(cat_idx) > 0L) {
        as.matrix(as.data.frame(
          lapply(X[, cat_idx, drop = FALSE], as.character),
          stringsAsFactors = FALSE
        ))
      } else {
        NULL
      }

      max_iter  <- 50L
      converged <- FALSE
      lambda    <- private$FLambda

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

      # prototypes: list of K lists (num, cat)
      centers <- vector("list", K)

      # --- internal helpers ---

      # Numeric part: same idea as variable k-means
      compute_Zk <- function(cols_k_num) {
        if (is.null(X_num) || length(cols_k_num) == 0L) return(NULL)

        if (length(cols_k_num) == 1L) {
          zk <- X_num[, cols_k_num]
          zk <- zk - mean(zk, na.rm = TRUE)
          return(as.numeric(zk))
        } else {
          Xk <- X_num[, cols_k_num, drop = FALSE]
          pc <- stats::prcomp(Xk, center = FALSE, scale. = FALSE)
          zk <- pc$x[, 1]
          return(as.numeric(zk))
        }
      }

      # Categorical part: same idea as variable k-modes
      compute_mode_profile <- function(cols_k_cat) {
        if (is.null(X_cat) || length(cols_k_cat) == 0L) return(NULL)

        if (length(cols_k_cat) == 1L) {
          return(X_cat[, cols_k_cat])
        } else {
          n_loc <- nrow(X_cat)
          mode_vec <- character(n_loc)
          for (i in seq_len(n_loc)) {
            vals       <- X_cat[i, cols_k_cat]
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

      # distance d(variable j, cluster k)
      # - if var j is numeric: d = 1 - r^2(x_j, Z_k)
      # - if var j is categorical: d = lambda * simple_matching(x_j, mode_k)
      compute_distance <- function(j, k) {
        is_num_j <- j %in% num_idx
        is_cat_j <- j %in% cat_idx
        c_k      <- centers[[k]]

        if (is_num_j && !is.null(c_k$num) && !is.null(X_num)) {
          col_num <- which(num_idx == j)
          xj      <- X_num[, col_num]
          r2      <- private$r2_corr(xj, c_k$num)
          d       <- 1 - r2
          if (is.na(d)) d <- 1
          return(d)
        }

        if (is_cat_j && !is.null(c_k$cat) && !is.null(X_cat)) {
          col_cat <- which(cat_idx == j)
          xj      <- X_cat[, col_cat]
          d_raw   <- private$simple_matching(xj, c_k$cat)
          d       <- lambda * d_raw
          return(d)
        }

        # fallback: no meaningful prototype → very large distance
        return(1e6)
      }

      # Inertia = sum of final distances
      compute_inertia <- function(clusters_vec) {
        tot <- 0
        for (j in seq_len(p)) {
          k <- clusters_vec[j]
          tot <- tot + compute_distance(j, k)
        }
        tot
      }

      old_inertia <- Inf

      # =============================
      # Main loop
      # =============================
      for (iter in seq_len(max_iter)) {

        # 1) recompute prototypes
        for (k in seq_len(K)) {
          vars_k <- which(clusters == k)
          if (length(vars_k) == 0L) {
            # empty cluster → force a random variable
            vars_k <- sample(seq_len(p), size = 1L)
            clusters[vars_k] <- k
          }

          vars_k_num <- intersect(vars_k, num_idx)
          vars_k_cat <- intersect(vars_k, cat_idx)

          zk <- if (!is.null(X_num) && length(vars_k_num) > 0L) {
            cols_k_num <- match(vars_k_num, num_idx)
            compute_Zk(cols_k_num)
          } else {
            NULL
          }

          mk <- if (!is.null(X_cat) && length(vars_k_cat) > 0L) {
            cols_k_cat <- match(vars_k_cat, cat_idx)
            compute_mode_profile(cols_k_cat)
          } else {
            NULL
          }

          centers[[k]] <- list(
            num = zk,
            cat = mk
          )
        }

        # 2) reassign variables
        new_clusters <- clusters

        for (j in seq_len(p)) {
          d_all <- vapply(
            X   = seq_len(K),
            FUN = function(k) compute_distance(j, k),
            FUN.VALUE = numeric(1L)
          )
          k_best <- which.min(d_all)
          new_clusters[j] <- k_best
        }

        inertia <- compute_inertia(new_clusters)

        # 3) stopping criteria
        if (all(new_clusters == clusters)) {
          converged <- TRUE
          clusters  <- new_clusters
          break
        }

        if (abs(old_inertia - inertia) < 1e-6) {
          converged <- TRUE
          clusters  <- new_clusters
          break
        }

        clusters    <- new_clusters
        old_inertia <- inertia
      }

      inertia <- compute_inertia(clusters)

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
      num_idx <- private$FNumCols
      cat_idx <- private$FCatCols
      centers <- private$FCenters
      lambda  <- private$FLambda

      if (is.null(centers)) {
        stop("[Kprototypes] No prototypes available (did you run fit()?).")
      }

      is_num <- is.numeric(x_new)
      is_cat <- is.factor(x_new) || is.character(x_new)

      if (!is_num && !is_cat) {
        stop("[Kprototypes] New variable must be numeric or categorical.")
      }

      K          <- length(centers)
      d_all      <- numeric(K)
      type_metric <- NA_character_

      if (is_num) {
        x <- as.numeric(x_new)
        for (k in seq_len(K)) {
          zk <- centers[[k]]$num
          if (is.null(zk)) {
            d_all[k] <- 1e6
          } else {
            r2      <- private$r2_corr(x, zk)
            d_all[k] <- 1 - r2
          }
        }

        k_best   <- which.min(d_all)
        d_raw    <- d_all[k_best]
        adhesion <- 1 - d_raw     # ≈ r^2
        type_metric <- "r2_num"

      } else {
        x <- as.character(x_new)
        for (k in seq_len(K)) {
          mk <- centers[[k]]$cat
          if (is.null(mk)) {
            d_all[k] <- 1e6
          } else {
            d_all[k] <- lambda * private$simple_matching(x, mk)
          }
        }

        k_best   <- which.min(d_all)
        d_raw    <- d_all[k_best]
        # d_raw = lambda * (1 - match_prop)
        match_prop <- 1 - (d_raw / lambda)
        adhesion   <- match_prop
        type_metric <- "match_cat"
      }

      data.frame(
        variable    = var_name,
        type        = if (is_num) "numeric" else "categorical",
        cluster     = as.integer(k_best),
        distance    = d_raw,
        adhesion    = adhesion,
        metric_type = type_metric,
        stringsAsFactors = FALSE
      )
    },

    # ===========================
    # 3. SUMMARY: membership indicators
    # ===========================
    summary_membership_impl = function() {
      X <- private$FX_active
      if (is.null(X)) {
        cat("(k-prototypes) No active variables stored.\n")
        return(invisible(NULL))
      }

      num_idx <- private$FNumCols
      cat_idx <- private$FCatCols
      centers <- private$FCenters
      lambda  <- private$FLambda

      X_num <- if (length(num_idx) > 0L) as.matrix(X[, num_idx, drop = FALSE]) else NULL
      X_cat <- if (length(cat_idx) > 0L) {
        as.matrix(as.data.frame(
          lapply(X[, cat_idx, drop = FALSE], as.character),
          stringsAsFactors = FALSE
        ))
      } else NULL

      p        <- ncol(X)
      clusters <- private$FClusters

      var_names   <- colnames(X)
      type_vec    <- character(p)
      dist_vec    <- numeric(p)
      adh_vec     <- numeric(p)
      metric_type <- character(p)

      for (j in seq_len(p)) {
        k   <- clusters[j]
        c_k <- centers[[k]]

        if (j %in% num_idx && !is.null(X_num)) {
          type_vec[j] <- "numeric"
          col_num <- which(num_idx == j)
          xj      <- X_num[, col_num]
          r2      <- private$r2_corr(xj, c_k$num)
          dist_vec[j]    <- 1 - r2
          adh_vec[j]     <- r2
          metric_type[j] <- "r2_num"

        } else if (j %in% cat_idx && !is.null(X_cat)) {
          type_vec[j] <- "categorical"
          col_cat <- which(cat_idx == j)
          xj      <- X_cat[, col_cat]
          d_raw   <- private$simple_matching(xj, c_k$cat)
          dist_vec[j]    <- lambda * d_raw
          adh_vec[j]     <- 1 - d_raw
          metric_type[j] <- "match_cat"

        } else {
          type_vec[j]    <- "unknown"
          dist_vec[j]    <- NA_real_
          adh_vec[j]     <- NA_real_
          metric_type[j] <- NA_character_
        }
      }

      df <- data.frame(
        variable    = var_names,
        type        = type_vec,
        cluster     = as.integer(clusters),
        distance    = dist_vec,
        adhesion    = adh_vec,
        metric_type = metric_type,
        stringsAsFactors = FALSE
      )

      cat("=== Membership indicators (k-prototypes) ===\n")
      cat("numeric     : adhesion = r^2(variable, numeric cluster profile)\n")
      cat("categorical : adhesion = match proportion with the cluster categorical profile\n\n")

      # Global indicators
      dist_global <- mean(df$distance, na.rm = TRUE)
      adh_global  <- mean(df$adhesion, na.rm = TRUE)
      explained_inertia_global <- 1 - dist_global   # same as adh_global

      cat(sprintf("Global mean distance              : %.3f\n", dist_global))
      cat(sprintf("Global mean adhesion              : %.3f\n", adh_global))
      cat(sprintf("Explained inertia (global, ratio) : %.3f\n", explained_inertia_global))
      cat(sprintf("Explained inertia (global, %% )    : %.1f%%\n",
                  100 * explained_inertia_global))
      cat("\n")

      # Separate numeric / categorical subsets (optional but informative)
      df_num <- df[df$metric_type == "r2_num", ]
      df_cat <- df[df$metric_type == "match_cat", ]

      adh_num <- if (nrow(df_num) > 0) mean(df_num$adhesion, na.rm = TRUE) else NA_real_
      adh_cat <- if (nrow(df_cat) > 0) mean(df_cat$adhesion, na.rm = TRUE) else NA_real_

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

      cat(sprintf("Global mean distance       : %.3f\n", dist_global))
      cat(sprintf("Global mean adhesion       : %.3f\n", adh_global))
      if (!is.na(adh_num)) {
        cat(sprintf("Mean adhesion (numeric)    : %.3f\n", adh_num))
      }
      if (!is.na(adh_cat)) {
        cat(sprintf("Mean adhesion (categorical): %.3f\n", adh_cat))
      }
      cat("\n")

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
        stop("[Kprototypes] plot(type = 'membership'): no active X available.")
      }

      num_idx <- private$FNumCols
      cat_idx <- private$FCatCols
      centers <- private$FCenters

      X_num <- if (length(num_idx) > 0L) as.matrix(X[, num_idx, drop = FALSE]) else NULL
      X_cat <- if (length(cat_idx) > 0L) {
        as.matrix(as.data.frame(
          lapply(X[, cat_idx, drop = FALSE], as.character),
          stringsAsFactors = FALSE
        ))
      } else NULL

      p        <- ncol(X)
      clusters <- private$FClusters

      adh_vec   <- numeric(p)
      var_names <- colnames(X)

      for (j in seq_len(p)) {
        k   <- clusters[j]
        c_k <- centers[[k]]

        if (j %in% num_idx && !is.null(X_num)) {
          col_num <- which(num_idx == j)
          xj      <- X_num[, col_num]
          adh_vec[j] <- private$r2_corr(xj, c_k$num)
        } else if (j %in% cat_idx && !is.null(X_cat)) {
          col_cat <- which(cat_idx == j)
          xj      <- X_cat[, col_cat]
          d_raw   <- private$simple_matching(xj, c_k$cat)
          adh_vec[j] <- 1 - d_raw
        } else {
          adh_vec[j] <- NA_real_
        }
      }

      o <- order(adh_vec, decreasing = TRUE, na.last = NA)
      barplot(
        adh_vec[o],
        names.arg = var_names[o],
        las = 2,
        main = "Variable–Cluster Membership (k-prototypes)",
        ylab = "adhesion",
        cex.names = 0.7
      )
    },

    # ===========================
    # 5. PLOT: profiles (numeric + categorical)
    # ===========================
    plot_profiles_impl = function() {
      X <- private$FX_active
      if (is.null(X)) {
        warning("[Kprototypes] plot(type = 'profiles'): no active X available.")
        return(invisible(NULL))
      }

      num_idx <- private$FNumCols
      cat_idx <- private$FCatCols
      centers <- private$FCenters

      n <- nrow(X)
      K <- private$FNbGroupes

      # numeric part: heatmap of Z_k (like k-means)
      has_num <- length(num_idx) > 0L &&
        any(vapply(centers, function(c_k) !is.null(c_k$num), logical(1L)))

      # categorical part: mean match proportion per individual (like k-modes)
      has_cat <- length(cat_idx) > 0L &&
        any(vapply(centers, function(c_k) !is.null(c_k$cat), logical(1L)))

      if (!has_num && !has_cat) {
        warning("[Kprototypes] No usable profiles for 'profiles' plot.")
        return(invisible(NULL))
      }

      op <- par(no.readonly = TRUE)
      on.exit(par(op))

      if (has_num && has_cat) {
        par(mfrow = c(1, 2))
      }

      # --- NUMERIC PROFILES ---
      if (has_num) {
        Z_mat <- matrix(NA_real_, nrow = n, ncol = K)
        colnames(Z_mat) <- paste0("Cluster ", seq_len(K))
        rownames(Z_mat) <- seq_len(n)

        for (k in seq_len(K)) {
          zk <- centers[[k]]$num
          if (!is.null(zk)) {
            Z_mat[, k] <- zk
          }
        }

        image(
          x = seq_len(K),
          y = seq_len(n),
          z = t(Z_mat),
          xlab = "Clusters",
          ylab = "Individuals",
          main = "Numeric profiles (k-prototypes)",
          axes = FALSE
        )
        axis(1, at = seq_len(K), labels = colnames(Z_mat))
        axis(2, at = pretty(seq_len(n)))
        box()
      }

      # --- CATEGORICAL PROFILES ---
      if (has_cat) {
        X_cat <- as.matrix(as.data.frame(
          lapply(X[, cat_idx, drop = FALSE], as.character),
          stringsAsFactors = FALSE
        ))
        clusters <- private$FClusters

        prof_cat <- matrix(NA_real_, nrow = n, ncol = K)
        colnames(prof_cat) <- paste0("Cluster ", seq_len(K))
        rownames(prof_cat) <- seq_len(n)

        for (k in seq_len(K)) {
          vars_k     <- which(clusters == k)
          vars_k_cat <- intersect(vars_k, cat_idx)
          if (length(vars_k_cat) == 0L) next

          mk <- centers[[k]]$cat
          if (is.null(mk)) next

          # match proportion per individual
          for (i in seq_len(n)) {
            cols_k_cat <- match(vars_k_cat, cat_idx)
            vals       <- X_cat[i, cols_k_cat]
            matches    <- (vals == mk[i])
            v <- mean(matches, na.rm = TRUE)
            if (is.na(v)) v <- NA_real_
            prof_cat[i, k] <- v
          }
        }

        image(
          x = seq_len(K),
          y = seq_len(n),
          z = t(prof_cat),
          xlab = "Clusters",
          ylab = "Individuals",
          main = "Categorical profiles (k-prototypes)\n(match proportion with cluster mode)",
          axes = FALSE
        )
        axis(1, at = seq_len(K), labels = colnames(prof_cat))
        axis(2, at = pretty(seq_len(n)))
        box()
      }

      invisible(NULL)
    }
  )
)
