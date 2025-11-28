#' K-medoids for variable clustering (mixed variables)
#'
#' R6 class implementing a k-medoids-like algorithm for clustering
#' \emph{variables} (columns) based on a dissimilarity matrix:
#'   - numeric vs numeric: 1 - r^2(x_j, x_l)
#'   - categorical vs categorical: simple matching dissimilarity
#'   - mixed (numeric vs categorical): fixed dissimilarity = 1
#'
#' @docType class
#' @name Kmedoids
#' @export
Kmedoids <- R6::R6Class(
  "Kmedoids",
  inherit = .ClusterBase,

  public = list(
    
    #' @description
    #' Create a new Kmedoids instance
    #' @param K Number of cluster
    #' @param lambda Weight for qualitative data
    #' @param random_state The random seed to use
    initialize = function(K, lambda = 1, random_state = NULL) {
      # lambda kept for signature consistency (can be used later to weight cat part)
      if (!is.numeric(lambda) || length(lambda) != 1L || lambda <= 0) {
        stop("[Kmedoids] lambda must be a numeric > 0.")
      }

      super$initialize(
        K           = K,
        scale       = scale,
        lambda      = lambda,
        method_name = "kmedoids",
        random_state = random_state
      )
    },

    #' @description
    #' Textual interpretation of the k-medoids variable clustering solution.
    #'
    #' For each cluster, this method reports:
    #'   - the number of variables,
    #'   - a few most central variables (closest to the medoid),
    #'   - simple membership statistics based on distance to the medoid.
    #'
    #' @param style Character, either "compact" or "detailed".
    interpret_clusters = function(style = c("compact", "detailed")) {
      style <- match.arg(style)

      X        <- private$FX_active
      clusters <- private$FClusters
      centers  <- private$FCenters

      if (is.null(X) || is.null(clusters) || is.null(centers)) {
        stop("[Kmedoids] interpret_clusters(): no fitted model or missing state.")
      }

      p <- ncol(X)
      K <- private$FNbGroupes

      # helper: distance entre deux variables (j, l)
      compute_dist_vars <- private$make_var_var_distance_fun(X)

      var_names <- colnames(X)

      cat("=== Global overview (k-medoids, variable clustering) ===\n")
      cat("Number of clusters :", K, "\n")
      cat("Number of variables:", p, "\n\n")

      for (k in seq_len(K)) {
        vars_k <- which(clusters == k)
        cat(sprintf("--- Cluster %d ---\n", k))

        if (length(vars_k) == 0L) {
          cat("Cluster contains no variables.\n\n")
          next
        }

        medoid_j <- centers[k]
        d_k <- vapply(
          vars_k,
          function(j) compute_dist_vars(j, medoid_j),
          numeric(1L)
        )

        # stats de distance à la médoïde
        d_mean <- mean(d_k, na.rm = TRUE)
        d_min  <- min(d_k, na.rm = TRUE)
        d_max  <- max(d_k, na.rm = TRUE)

        cat(sprintf("Size: %d variables.\n", length(vars_k)))
        cat(sprintf("Distances to medoid (cluster %d):\n", k))
        cat(sprintf("  - Mean: %.3f\n", d_mean))
        cat(sprintf("  - Min : %.3f\n", d_min))
        cat(sprintf("  - Max : %.3f\n", d_max))

        # top variables les plus proches de la médoïde
        o_local <- order(d_k, decreasing = FALSE)
        top_idx <- vars_k[o_local]
        if (length(top_idx) > 3L) {
          top_idx <- top_idx[seq_len(3L)]
        }

        cat("Most central variables (closest to the medoid):\n")
        for (j in top_idx) {
          d_j <- compute_dist_vars(j, medoid_j)
          cat(sprintf("  - %s (distance to medoid = %.3f)\n",
                      var_names[j], d_j))
        }

        if (style == "detailed") {
          cat("Interpretation:\n")
          cat("  This cluster groups variables that are mutually similar according\n")
          cat("  to the mixed dissimilarity (1 - r^2 for numeric pairs, simple\n")
          cat("  matching for categorical pairs, 1 for mixed pairs). The medoid is\n")
          cat("  the most central variable of the cluster.\n")
        }
        cat("\n")
      }

      invisible(NULL)
    }
  ),

  private = list(

    # ================================
    # 1. VARIABLE K-MEDOIDS ALGORITHM
    # ================================
    # centers: vector of medoid indices (1..p)
    #
    run_clustering = function(X) {
      # X: mixed data.frame (numeric + factor/character),
      # already checked by check_and_prepare_X() and possibly scaled.

      if (!requireNamespace("cluster", quietly = TRUE)) {
        stop("[Kmedoids] The 'cluster' package is required for k-medoids.")
      }

      n <- nrow(X)
      p <- ncol(X)
      K <- private$FNbGroupes

      if (K > p) {
        stop("[Kmedoids] K cannot exceed the number of variables.")
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

      # --- helper: distance entre deux variables j et l (indices 1..p) ----
      compute_dist_vars <- function(j, l) {
        if (j == l) return(0)

        j_is_num <- j %in% num_idx
        l_is_num <- l %in% num_idx
        j_is_cat <- j %in% cat_idx
        l_is_cat <- l %in% cat_idx

        # numeric vs numeric : 1 - r^2
        if (j_is_num && l_is_num && !is.null(X_num)) {
          col_j <- which(num_idx == j)
          col_l <- which(num_idx == l)
          xj    <- X_num[, col_j]
          xl    <- X_num[, col_l]
          r2    <- private$r2_corr(xj, xl)
          d     <- 1 - r2
          if (is.na(d)) d <- 1
          return(d)
        }

        # categorical vs categorical : simple matching
        if (j_is_cat && l_is_cat && !is.null(X_cat)) {
          col_j <- which(cat_idx == j)
          col_l <- which(cat_idx == l)
          xj    <- X_cat[, col_j]
          xl    <- X_cat[, col_l]
          d     <- private$simple_matching(xj, xl)
          if (is.na(d)) d <- 1
          return(d)
        }

        # Mixed types (numeric vs categorical) → dissimilarité maximale 1
        return(1)
      }

      # On stocke la fonction dans l'objet pour réutilisation (summary, interpret)
      private$make_var_var_distance_fun <- function(X_internal) {
        # capture environment: num_idx, cat_idx, X_num, X_cat, private$r2_corr, private$simple_matching
        function(j, l) compute_dist_vars(j, l)
      }

      # --- Construction de la matrice de dissimilarité entre variables ----
      D <- matrix(0, nrow = p, ncol = p)
      colnames(D) <- rownames(D) <- colnames(X)

      for (j in seq_len(p)) {
        for (l in seq_len(j)) {
          d_jl      <- compute_dist_vars(j, l)
          D[j, l]   <- d_jl
          D[l, j]   <- d_jl
        }
      }

      # --- k-medoids (PAM) sur la matrice de dissimilarité ------------------
      pam_fit <- cluster::pam(D, k = K, diss = TRUE)

      clusters    <- as.integer(pam_fit$clustering)   # longueur = p
      medoid_pos  <- pam_fit$medoids                  # positions (1..p)
      medoid_idx  <- as.integer(medoid_pos)           # indices de variables

      # pam$objective[1] = sum des distances aux médoïdes
      inertia <- as.numeric(pam_fit$objective[1])
      
      names(clusters) <- names(X)

      list(
        clusters  = clusters,
        centers   = medoid_idx,  # indices des médoïdes
        inertia   = inertia,
        converged = TRUE
      )
    },

    # petit slot pour exposer la fonction distance (définie dans run_clustering)
    make_var_var_distance_fun = NULL,

    # =====================
    # 2. PREDICT ONE VARIABLE
    # =====================
    # Ici, on affecte une nouvelle variable en comparant sa dissimilarité
    # aux médoïdes des clusters existants.
    predict_one_variable = function(x_new, var_name) {
      X_ref   <- private$FX_active
      centers <- private$FCenters
      num_idx <- private$FNumCols
      cat_idx <- private$FCatCols

      if (is.null(X_ref) || is.null(centers)) {
        stop("[Kmedoids] No prototypes available (did you run fit()?).")
      }

      is_num <- is.numeric(x_new)
      is_cat <- is.factor(x_new) || is.character(x_new)

      if (!is_num && !is_cat) {
        stop("[Kmedoids] New variable must be numeric or categorical.")
      }

      # Préparer x_new selon son type
      if (is_num) {
        x_new_vec <- as.numeric(x_new)
      } else {
        x_new_vec <- as.character(x_new)
      }

      # Construit un petit helper distance(x_new, var_j)
      compute_new_to_var <- function(j) {
        j_is_num <- j %in% num_idx
        j_is_cat <- j %in% cat_idx

        if (is_num && j_is_num) {
          col_j <- which(num_idx == j)
          xj    <- as.numeric(X_ref[, col_j])
          r2    <- private$r2_corr(x_new_vec, xj)
          d     <- 1 - r2
          if (is.na(d)) d <- 1
          return(d)
        }

        if (is_cat && j_is_cat) {
          col_j <- which(cat_idx == j)
          xj    <- as.character(X_ref[, col_j])
          d     <- private$simple_matching(x_new_vec, xj)
          if (is.na(d)) d <- 1
          return(d)
        }

        # mixte → dissimilarité max
        return(1)
      }

      K    <- length(centers)
      d_all <- numeric(K)

      for (k in seq_len(K)) {
        med_j <- centers[k]
        d_all[k] <- compute_new_to_var(med_j)
      }

      k_best   <- which.min(d_all)
      d_raw    <- d_all[k_best]
      adhesion <- 1 - d_raw

      data.frame(
        variable = var_name,
        type     = if (is_num) "numeric" else "categorical",
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
        cat("(k-medoids) No active variables stored.\n")
        return(invisible(NULL))
      }

      clusters <- private$FClusters
      centers  <- private$FCenters
      p        <- ncol(X)

      if (is.null(clusters) || is.null(centers)) {
        cat("(k-medoids) No clustering result available.\n")
        return(invisible(NULL))
      }

      # récupère la fonction distance(j, l) définie dans run_clustering
      if (is.null(private$make_var_var_distance_fun)) {
        # au cas où, on la recrée à partir de X
        private$make_var_var_distance_fun <- function(X_internal) {
          num_idx <- private$FNumCols
          cat_idx <- private$FCatCols

          X_num <- if (length(num_idx) > 0L) {
            as.matrix(X_internal[, num_idx, drop = FALSE])
          } else NULL

          X_cat <- if (length(cat_idx) > 0L) {
            as.matrix(as.data.frame(
              lapply(X_internal[, cat_idx, drop = FALSE], as.character),
              stringsAsFactors = FALSE
            ))
          } else NULL

          function(j, l) {
            if (j == l) return(0)

            j_is_num <- j %in% num_idx
            l_is_num <- l %in% num_idx
            j_is_cat <- j %in% cat_idx
            l_is_cat <- l %in% cat_idx

            if (j_is_num && l_is_num && !is.null(X_num)) {
              col_j <- which(num_idx == j)
              col_l <- which(num_idx == l)
              xj    <- X_num[, col_j]
              xl    <- X_num[, col_l]
              r2    <- private$r2_corr(xj, xl)
              d     <- 1 - r2
              if (is.na(d)) d <- 1
              return(d)
            }

            if (j_is_cat && l_is_cat && !is.null(X_cat)) {
              col_j <- which(cat_idx == j)
              col_l <- which(cat_idx == l)
              xj    <- X_cat[, col_j]
              xl    <- X_cat[, col_l]
              d     <- private$simple_matching(xj, xl)
              if (is.na(d)) d <- 1
              return(d)
            }

            return(1)
          }
        }
      }

      compute_dist_vars <- private$make_var_var_distance_fun(X)

      var_names <- colnames(X)
      dist_vec  <- numeric(p)
      adh_vec   <- numeric(p)

      for (j in seq_len(p)) {
        k       <- clusters[j]
        medoid  <- centers[k]
        d_j     <- compute_dist_vars(j, medoid)
        dist_vec[j] <- d_j
        adh_vec[j]  <- 1 - d_j
      }

      df <- data.frame(
        variable = var_names,
        cluster  = as.integer(clusters),
        distance = dist_vec,
        adhesion = adh_vec,
        stringsAsFactors = FALSE
      )

      cat("=== Membership indicators (k-medoids) ===\n")
      cat("distance = dissimilarity to cluster medoid\n")
      cat("adhesion = 1 - distance (closeness to medoid)\n\n")

      dist_global <- mean(df$distance, na.rm = TRUE)
      adh_global  <- mean(df$adhesion, na.rm = TRUE)
      explained_inertia_global <- 1 - dist_global

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
      cat(sprintf("Explained inertia (global, ratio) : %.3f\n", explained_inertia_global))
      cat(sprintf("Explained inertia (global, %% )    : %.1f%%\n",
                  100 * explained_inertia_global))
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
        stop("[Kmedoids] plot(type = 'membership'): no active X available.")
      }

      clusters <- private$FClusters
      centers  <- private$FCenters
      p        <- ncol(X)

      if (is.null(clusters) || is.null(centers)) {
        stop("[Kmedoids] No clustering result available.")
      }

      if (is.null(private$make_var_var_distance_fun)) {
        private$summary_membership_impl()  # initialise la distance
      }

      compute_dist_vars <- private$make_var_var_distance_fun(X)

      var_names <- colnames(X)
      adh_vec   <- numeric(p)

      for (j in seq_len(p)) {
        k      <- clusters[j]
        medoid <- centers[k]
        d_j    <- compute_dist_vars(j, medoid)
        adh_vec[j] <- 1 - d_j
      }

      o <- order(adh_vec, decreasing = TRUE, na.last = NA)
      barplot(
        adh_vec[o],
        names.arg = var_names[o],
        las = 2,
        main = "Variable–Cluster Membership (k-medoids)",
        ylab = "adhesion (1 - distance to medoid)",
        cex.names = 0.7
      )
    }

  )
)
