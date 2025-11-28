#' K-modes for variable clustering (categorical variables only)
#'
#' R6 class implementing a k-modes-like algorithm for clustering
#' \emph{variables} (columns) when all active variables are categorical.
#'
#' @docType class
#' @name Kmodes
#' @export
Kmodes <- R6::R6Class(
  "Kmodes",
  inherit = .ClusterBase,

  public = list(
    
    #' @description
    #' Create a new Kmodes instance
    #' @param K Number of cluster
    #' @param scale A boolean defining whether to scale the data or not
    #' @param random_state The random seed to use
    initialize = function(K, random_state = NULL) {
      super$initialize(
        K           = K,
        method_name = "kmodes",
        random_state = random_state
      )
    },

    #' @description
    #' Textual interpretation of the k-modes variable clustering solution.
    #'
    #' For each cluster, this method reports:
    #'   - the number of variables (categorical),
    #'   - a few most representative variables (highest agreement with the
    #'     cluster modal profile),
    #'   - simple membership statistics based on this agreement.
    #'
    #' The "distance" is defined as 1 - proportion of matches with the
    #' cluster-level modal profile, and "adhesion" as that proportion.
    #'
    #' @param style Character, either "compact" or "detailed".
    interpret_clusters = function(style = c("compact", "detailed")) {
      style <- match.arg(style)

      X <- private$FX_active
      clusters <- private$FClusters

      if (is.null(X) || is.null(clusters)) {
        stop("[Kmodes] interpret_clusters(): no fitted model or missing state.")
      }

      p <- ncol(X)
      K <- private$FNbGroupes

      is_cat <- vapply(
        X,
        function(col) is.factor(col) || is.character(col),
        logical(1L)
      )
      if (!all(is_cat)) {
        warning("[Kmodes] Some non-categorical variables detected; ",
                "only categorical variables are used for interpretation.")
      }

      type_var <- ifelse(is_cat, "categorical", "other")

      X_cat <- as.data.frame(
        lapply(X[, is_cat, drop = FALSE], as.character),
        stringsAsFactors = FALSE
      )
      X_cat_mat     <- as.matrix(X_cat)
      idx_cat_global <- which(is_cat)

      distance <- rep(NA_real_, p)
      adhesion <- rep(NA_real_, p)

      n <- nrow(X_cat_mat)

      # for each cluster, build a cluster-level modal profile by individual,
      # then compute agreement of each variable with that profile
      for (k in seq_len(K)) {
        vars_k <- which(clusters == k & is_cat)
        if (length(vars_k) == 0L) next

        cols_k_cat <- match(vars_k, idx_cat_global)
        mat_k <- X_cat_mat[, cols_k_cat, drop = FALSE]

        # cluster-level modal profile by individual (majority category)
        modal_profile <- apply(mat_k, 1, function(row) {
          row_no_na <- row[!is.na(row)]
          if (length(row_no_na) == 0L) return(NA_character_)
          tab <- table(row_no_na)
          names(tab)[which.max(tab)]
        })

        # agreement for each variable of the cluster
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

      df <- data.frame(
        variable = colnames(X),
        type     = type_var,
        cluster  = as.integer(clusters),
        distance = distance,
        adhesion = adhesion,
        stringsAsFactors = FALSE
      )

      n_cat <- sum(type_var == "categorical")

      cat("=== Global overview (k-modes, variable clustering) ===\n")
      cat("Number of clusters :", K, "\n")
      cat("Number of variables:", p,
          sprintf("(%d categorical, %d other)\n", n_cat, p - n_cat))
      cat("\n")

      for (k in seq_len(K)) {
        df_k <- df[df$cluster == k & df$type == "categorical", , drop = FALSE]

        cat(sprintf("--- Cluster %d ---\n", k))

        if (nrow(df_k) == 0L) {
          cat("Cluster contains no categorical variables.\n\n")
          next
        }

        size_k <- nrow(df_k)
        cat(sprintf("Size: %d categorical variables.\n", size_k))

        adh_mean_k <- mean(df_k$adhesion, na.rm = TRUE)
        adh_min_k  <- min(df_k$adhesion, na.rm = TRUE)
        adh_max_k  <- max(df_k$adhesion, na.rm = TRUE)

        o_k   <- order(df_k$adhesion, decreasing = TRUE)
        top_k <- df_k[o_k, , drop = FALSE]
        if (nrow(top_k) > 3L) {
          top_k <- top_k[seq_len(3L), , drop = FALSE]
        }

        if (style == "compact") {
          cat("Most representative variables (top by agreement with the cluster modal profile):\n")
          for (i in seq_len(nrow(top_k))) {
            cat(sprintf("  - %s (adhesion = %.3f)\n",
                        top_k$variable[i],
                        top_k$adhesion[i]))
          }
          cat("Membership statistics:\n")
          cat(sprintf("  - Mean adhesion: %.3f\n", adh_mean_k))
          cat(sprintf("  - Min / Max adhesion: %.3f / %.3f\n",
                      adh_min_k, adh_max_k))
          cat("Interpretation:\n")
          cat("  This cluster groups categorical variables that tend to\n")
          cat("  share the same dominant categories across individuals.\n\n")

        } else {  # detailed
          cat("Most representative variables (top 3 by agreement with the modal profile):\n")
          for (i in seq_len(nrow(top_k))) {
            cat(sprintf("  - %s (adhesion = %.3f)\n",
                        top_k$variable[i],
                        top_k$adhesion[i]))
          }
          cat("Membership statistics:\n")
          cat(sprintf("  - Mean adhesion: %.3f\n", adh_mean_k))
          cat(sprintf("  - Min / Max adhesion: %.3f / %.3f\n",
                      adh_min_k, adh_max_k))
          cat("Interpretation:\n")
          cat("  This cluster is characterised by a shared modal pattern\n")
          cat("  across individuals; variables with high adhesion follow\n")
          cat("  this pattern closely.\n\n")
        }
      }

      invisible(df)
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
        stop("[Kmodes] All variables must be categorical.")
      }

      # Work with a character matrix for simplicity
      X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)
      X_mat  <- as.matrix(X_char)

      n <- nrow(X_mat)
      p <- ncol(X_mat)
      K <- private$FNbGroupes

      if (K > p) {
        stop("[Kmodes] K cannot exceed the number of variables.")
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
            function(mk) private$simple_matching(xj, mk),
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
        inertia <- inertia + private$simple_matching(xj, mk)
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
        stop("[Kmodes] predict() requires a categorical variable.")
      }

      centers <- private$FCenters
      if (is.null(centers)) {
        stop("[Kmodes] No prototypes available (did you run fit()?).")
      }

      x_char <- as.character(x_new)

      K <- length(centers)
      d_all <- vapply(
        centers,
        function(mk) private$simple_matching(x_char, mk),
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
    summary_membership_impl = function() {
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

      if (is.null(clusters) || is.null(centers)) {
        cat("(k-modes) No clusters or modes available.\n")
        return(invisible(NULL))
      }

      distance <- numeric(p)
      adhesion <- numeric(p)

      for (j in seq_len(p)) {
        k  <- clusters[j]
        mk <- centers[[k]]  # mode profile for cluster k (length n)
        xj <- X_mat[, j]

        d_raw <- private$simple_matching(xj, mk)  # proportion of mismatches
        distance[j] <- d_raw
        adhesion[j] <- 1 - d_raw                  # proportion of matches
      }

      df <- data.frame(
        variable = colnames(X_mat),
        cluster  = as.integer(clusters),
        distance = distance,
        adhesion = adhesion,
        stringsAsFactors = FALSE
      )

      cat("=== Membership indicators (k-modes) ===\n")
      cat("distance = proportion of mismatches with the cluster mode profile\n")
      cat("adhesion = 1 - distance (proportion of matches)\n\n")

      ## ---- Global indicators ----
      dist_global <- mean(df$distance, na.rm = TRUE)
      adh_global  <- mean(df$adhesion, na.rm = TRUE)

      explained_inertia_global <- 1 - dist_global  # same as adh_global

      cat(sprintf("Global mean distance              : %.3f\n", dist_global))
      cat(sprintf("Global mean adhesion              : %.3f\n", adh_global))
      cat(sprintf("Explained inertia (global, ratio) : %.3f\n", explained_inertia_global))
      cat(sprintf("Explained inertia (global, %% )    : %.1f%%\n",
                  100 * explained_inertia_global))
      cat("\n")

      ## ---- Cluster-level statistics ----
      stats_list <- lapply(split(df, df$cluster), function(dsub) {
        dist_mean <- mean(dsub$distance, na.rm = TRUE)
        adh_mean  <- mean(dsub$adhesion, na.rm = TRUE)

        c(
          cluster              = dsub$cluster[1],
          dist_mean            = dist_mean,
          dist_min             = min(dsub$distance, na.rm = TRUE),
          dist_max             = max(dsub$distance, na.rm = TRUE),
          adh_mean             = adh_mean,
          adh_min              = min(dsub$adhesion, na.rm = TRUE),
          adh_max              = max(dsub$adhesion, na.rm = TRUE),
          explained_inertia_cl = 1 - dist_mean
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
        stop("[Kmodes] plot(type='membership'): no active X available.")
      }

      X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)
      X_mat  <- as.matrix(X_char)

      p        <- ncol(X_mat)
      clusters <- private$FClusters
      centers  <- private$FCenters

      membership <- numeric(p)

      for (j in seq_len(p)) {
        k  <- clusters[j]
        mk <- centers[[k]]
        xj <- X_mat[, j]
        d_raw <- private$simple_matching(xj, mk)
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
    plot_profiles_impl = function() {
      X <- private$FX_active
      if (is.null(X)) {
        warning("[Kmodes] plot(type='profiles'): no active X available.")
        return(invisible(NULL))
      }

      X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)
      X_mat  <- as.matrix(X_char)

      n        <- nrow(X_mat)
      clusters <- private$FClusters
      centers  <- private$FCenters
      K        <- private$FNbGroupes

      if (is.null(clusters) || is.null(centers)) {
        warning("[Kmodes] No clusters or prototypes available.")
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
        if (is.null(mk)) next

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
