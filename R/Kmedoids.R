#' K-medoids for variable clustering (mixed variables)
#'
#' R6 class implementing mixed-variable k-medoids for clustering columns of X.
#' Dissimilarity rules:
#'   - numeric vs numeric: 1 - r^2
#'   - categorical vs categorical: simple matching
#'   - mixed: 1
#'
#' @export
Kmedoids <- R6::R6Class(
  "Kmedoids",
  inherit = .ClusterBase,

  public = list(

    initialize = function(K, lambda = 1, scale = TRUE, random_state = NULL) {
      # lambda retained for compatibility with prototypes logic
      if (!is.numeric(lambda) || lambda <= 0)
        stop("[Kmedoids] lambda must be > 0.")

      super$initialize(
        K            = K,
        scale        = scale,
        lambda       = lambda,
        method_name  = "kmedoids",
        random_state = random_state
      )
    },

    #' Textual interpretation
    interpret_clusters = function(style = c("compact", "detailed")) {
      style <- match.arg(style)

      X        <- private$FX_active
      clusters <- private$FClusters
      centers  <- private$FCenters

      if (is.null(X) || is.null(clusters) || is.null(centers))
        stop("[Kmedoids] interpret_clusters(): no fitted model.")

      p <- ncol(X)
      K <- private$FNbGroupes

      # ensure distance function exists
      if (is.null(private$F_dist_fun))
        private$F_dist_fun <- private$build_distance_function(X)

      dfun      <- private$F_dist_fun
      var_names <- colnames(X)

      cat("=== Cluster interpretation (k-medoids) ===\n")
      cat("Variables:", p, " | Clusters:", K, "\n\n")

      for (k in seq_len(K)) {

        vars_k <- which(clusters == k)
        cat(sprintf("--- Cluster %d ---\n", k))

        if (length(vars_k) == 0) {
          cat("  (empty)\n\n")
          next
        }

        med <- centers[k]
        dist_k <- vapply(vars_k, function(j) dfun(j, med), numeric(1))

        cat("Size:", length(vars_k), "\n")
        cat("Distances to medoid:\n")
        cat("  mean:", sprintf("%.3f", mean(dist_k)), "\n")
        cat("  min :", sprintf("%.3f", min(dist_k)), "\n")
        cat("  max :", sprintf("%.3f", max(dist_k)), "\n")

        ord <- order(dist_k)
        top <- vars_k[ord][1:min(3, length(vars_k))]

        cat("Most central variables:\n")
        for (j in top) {
          cat(sprintf("  - %s (d=%.3f)\n", var_names[j], dfun(j, med)))
        }

        if (style == "detailed") {
          cat("Interpretation:\n")
          cat("  Variables in this cluster share similar behaviour under the mixed\n")
          cat("  dissimilarity (1 - r^2 for numeric, matching for categorical, 1 for mixed).\n")
          cat("  The medoid is the most central variable.\n")
        }
        cat("\n")
      }

      invisible(NULL)
    }
  ),

  private = list(

    # ===============================================================
    # 1. BUILD DISTANCE FUNCTION BETWEEN VARIABLES ------------------
    # ===============================================================

    F_dist_fun = NULL,

    build_distance_function = function(X) {

      num_idx <- private$FNumCols
      cat_idx <- private$FCatCols

      # numeric matrix
      Xn <- if (length(num_idx))
        as.matrix(X[, num_idx, drop = FALSE]) else NULL

      # categorical as character
      Xc <- if (length(cat_idx))
        as.matrix(as.data.frame(
          lapply(X[, cat_idx, drop = FALSE], as.character),
          stringsAsFactors = FALSE
        )) else NULL

      dfun <- function(j, l) {

        if (j == l) return(0)

        j_num <- j %in% num_idx
        l_num <- l %in% num_idx
        j_cat <- j %in% cat_idx
        l_cat <- l %in% cat_idx

        # numeric vs numeric
        if (j_num && l_num) {
          cj <- which(num_idx == j)
          cl <- which(num_idx == l)
          r2 <- private$r2_corr(Xn[, cj], Xn[, cl])
          return(1 - r2)
        }

        # categorical vs categorical
        if (j_cat && l_cat) {
          cj <- which(cat_idx == j)
          cl <- which(cat_idx == l)
          return(private$simple_matching(Xc[, cj], Xc[, cl]))
        }

        # mixed → max distance
        return(1)
      }

      dfun
    },

    # ===============================================================
    # 2. RUN K-MEDOIDS CLUSTERING ----------------------------------
    # ===============================================================

    run_clustering = function(X) {

      if (!requireNamespace("cluster", quietly = TRUE))
        stop("[Kmedoids] package 'cluster' is required.")

      p <- ncol(X)
      K <- private$FNbGroupes

      if (K > p)
        stop("[Kmedoids] Cannot cluster into more groups than variables.")

      # Build + cache distance function
      private$F_dist_fun <- private$build_distance_function(X)
      dfun <- private$F_dist_fun

      # Build dissimilarity matrix D (p x p)
      D <- matrix(0, p, p)
      colnames(D) <- rownames(D) <- colnames(X)

      for (j in seq_len(p))
        for (l in seq_len(j)) {
          d <- dfun(j, l)
          if (is.na(d)) d <- 1
          D[j, l] <- D[l, j] <- d
        }

      # PAM
      pam <- cluster::pam(D, k = K, diss = TRUE)

      clusters <- as.integer(pam$clustering)
      medoids  <- as.integer(pam$medoids)
      inertia  <- pam$objective[1]

      names(clusters) <- colnames(X)

      # compute per-cluster inertia (REQUIRED for InertiaObject)
      inertia_per_cluster <- tapply(
        seq_len(p),
        clusters,
        function(idx) {
          sum(sapply(idx, function(j) {
            med <- medoids[clusters[j]]
            dfun(j, med)
          }))
        }
      )

      return(list(
        clusters            = clusters,
        centers             = medoids,
        inertia             = inertia,
        inertia_per_cluster = inertia_per_cluster,
        converged           = TRUE
      ))
    },

    # ===============================================================
    # 3. PREDICT ONE VARIABLE ---------------------------------------
    # ===============================================================

    predict_one_variable = function(x_new, var_name) {

      Xref    <- private$FX_active
      centers <- private$FCenters
      num_idx <- private$FNumCols
      cat_idx <- private$FCatCols

      if (is.null(private$F_dist_fun))
        private$F_dist_fun <- private$build_distance_function(Xref)

      dfun <- private$F_dist_fun

      is_num <- is.numeric(x_new)
      is_cat <- is.factor(x_new) || is.character(x_new)

      if (!is_num && !is_cat)
        stop("[Kmedoids] new variable must be numeric or categorical.")

      xvec <- if (is_num) as.numeric(x_new) else as.character(x_new)

      # distance between x_new and each cluster medoid
      dist_to_medoid <- numeric(length(centers))

      for (k in seq_along(centers)) {

        j <- centers[k]  # variable index of medoid

        if (is_num && j %in% num_idx) {
          col_j <- which(num_idx == j)
          xj    <- as.numeric(Xref[, col_j])
          d     <- 1 - private$r2_corr(xvec, xj)

        } else if (is_cat && j %in% cat_idx) {
          col_j <- which(cat_idx == j)
          xj    <- as.character(Xref[, col_j])
          d     <- private$simple_matching(xvec, xj)

        } else {
          d <- 1
        }

        if (is.na(d)) d <- 1
        dist_to_medoid[k] <- d
      }

      k_best   <- which.min(dist_to_medoid)
      d_best   <- dist_to_medoid[k_best]
      adhesion <- 1 - d_best

      data.frame(
        variable = var_name,
        type     = if (is_num) "numeric" else "categorical",
        cluster  = k_best,
        distance = d_best,
        adhesion = adhesion,
        stringsAsFactors = FALSE
      )
    },

    # ===============================================================
    # 4. MEMBERSHIP SUMMARY & PLOT ----------------------------------
    # ===============================================================

    summary_membership_impl = function() {

      X        <- private$FX_active
      clusters <- private$FClusters
      centers  <- private$FCenters
      p        <- ncol(X)

      if (is.null(private$F_dist_fun))
        private$F_dist_fun <- private$build_distance_function(X)

      dfun <- private$F_dist_fun
      vnames <- colnames(X)

      dist_vec <- numeric(p)
      adh_vec  <- numeric(p)

      for (j in seq_len(p)) {
        med   <- centers[clusters[j]]
        d     <- dfun(j, med)
        dist_vec[j] <- d
        adh_vec[j]  <- 1 - d
      }

      df <- data.frame(
        variable = vnames,
        cluster  = clusters,
        distance = dist_vec,
        adhesion = adh_vec,
        stringsAsFactors = FALSE
      )

      cat("=== Membership indicators (k-medoids) ===\n")
      cat("distance = dissimilarity to cluster medoid\n")
      cat("adhesion = 1 - distance\n\n")

      print(df)
      invisible(df)
    },

    plot_membership_impl = function() {

      X <- private$FX_active
      clusters <- private$FClusters
      centers  <- private$FCenters
      p <- ncol(X)

      if (is.null(private$F_dist_fun))
        private$F_dist_fun <- private$build_distance_function(X)

      dfun <- private$F_dist_fun

      adh <- numeric(p)
      for (j in seq_len(p)) {
        med <- centers[clusters[j]]
        adh[j] <- 1 - dfun(j, med)
      }

      ord <- order(adh, decreasing = TRUE)
      barplot(
        adh[ord],
        names.arg = colnames(X)[ord],
        las = 2,
        main = "Variable–Cluster Membership (k-medoids)",
        ylab = "Adhesion (1 - distance)",
        cex.names = 0.7
      )
    }
  )
)
