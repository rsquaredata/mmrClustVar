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
        d_k <- vapply(vars_k, function(j) compute_dist_vars(j, medoid_j), numeric(1))

        cat(sprintf("Size: %d variables.\n", length(vars_k)))
        cat(sprintf("Distances to medoid (cluster %d):\n", k))
        cat(sprintf("  - Mean: %.3f\n", mean(d_k)))
        cat(sprintf("  - Min : %.3f\n", min(d_k)))
        cat(sprintf("  - Max : %.3f\n", max(d_k)))

        # central vars
        o <- order(d_k)
        top_idx <- vars_k[o]
        if (length(top_idx) > 3L) top_idx <- top_idx[1:3]

        cat("Most central variables (closest to the medoid):\n")
        for (j in top_idx) {
          cat(sprintf("  - %s (distance = %.3f)\n",
                      var_names[j],
                      compute_dist_vars(j, medoid_j)))
        }

        if (style == "detailed") {
          cat("Interpretation:\n")
          cat("  Variables grouped by mixed dissimilarity.\n")
          cat("  Medoid = central variable.\n")
        }
        cat("\n")
      }

      invisible(NULL)
    }
  ),

  private = list(

    # =============================================================
    # 1. Core algorithm
    # =============================================================
    run_clustering = function(X) {

      if (!requireNamespace("cluster", quietly = TRUE)) {
        stop("[Kmedoids] The 'cluster' package is required for k-medoids.")
      }

      p <- ncol(X)
      K <- private$FNbGroupes
      num_idx <- private$FNumCols
      cat_idx <- private$FCatCols

      # Pre-processed matrices
      X_num <- if (length(num_idx) > 0) as.matrix(X[, num_idx, drop = FALSE]) else NULL
      X_cat <- if (length(cat_idx) > 0) {
        as.matrix(as.data.frame(
          lapply(X[, cat_idx, drop = FALSE], as.character),
          stringsAsFactors = FALSE
        ))
      } else NULL

      # ------------------------------
      # distance between variables j,l
      # ------------------------------
      compute_dist_vars <- function(j, l) {

        if (j == l) return(0)

        j_is_num <- j %in% num_idx
        l_is_num <- l %in% num_idx
        j_is_cat <- j %in% cat_idx
        l_is_cat <- l %in% cat_idx

        # num–num → 1 - r²
        if (j_is_num && l_is_num && !is.null(X_num)) {
          col_j <- which(num_idx == j)
          col_l <- which(num_idx == l)
          r2 <- private$r2_corr(X_num[, col_j], X_num[, col_l])
          d <- 1 - r2
          if (is.na(d)) d <- 1
          return(d)
        }

        # cat–cat → simple matching
        if (j_is_cat && l_is_cat && !is.null(X_cat)) {
          col_j <- which(cat_idx == j)
          col_l <- which(cat_idx == l)
          d <- private$simple_matching(X_cat[, col_j], X_cat[, col_l])
          if (is.na(d)) d <- 1
          return(d)
        }

        # mixed → max
        return(1)
      }

      private$make_var_var_distance_fun <- function(Xref) {
        function(j, l) compute_dist_vars(j, l)
      }

      # ------------------------------
      # dissimilarity matrix p×p
      # ------------------------------
      D <- matrix(0, p, p)
      colnames(D) <- rownames(D) <- colnames(X)

      for (j in seq_len(p)) {
        for (l in seq_len(j)) {
          d <- compute_dist_vars(j, l)
          D[j, l] <- d
          D[l, j] <- d
        }
      }

      # ------------------------------
      # PAM algorithm
      # ------------------------------
      pam_fit <- cluster::pam(D, k = K, diss = TRUE)

      clusters   <- as.integer(pam_fit$clustering)
      medoid_idx <- as.integer(pam_fit$medoids)
      inertia    <- as.numeric(pam_fit$objective[1])

      # ------------------------------
      # inertia per cluster
      # ------------------------------
      inertia_per_cluster <- tapply(
        seq_len(p),
        clusters,
        function(idx) sum(sapply(idx, function(j) {
          compute_dist_vars(j, medoid_idx[clusters[j]])
        }))
      )

      names(clusters) <- colnames(X)

      return(list(
        clusters  = clusters,
        centers   = medoid_idx,
        inertia   = inertia,
        inertia_per_cluster = inertia_per_cluster,
        converged = TRUE
      ))
    },

    # Required for predict() — unchanged
    make_var_var_distance_fun = NULL,

    predict_one_variable = function(x_new, var_name) {
      # … (identique à ton code original)
      # aucune modification
      stop("predict_one_variable() unchanged here for brevity.")
    },

    summary_membership_impl = function() {
      # … (identique à ton code original)
      # aucune modification
      stop("summary_membership_impl() unchanged here for brevity.")
    },

    plot_membership_impl = function() {
      # … (identique à ton code original)
      # aucune modification
      stop("plot_membership_impl() unchanged here for brevity.")
    }
  )
)
