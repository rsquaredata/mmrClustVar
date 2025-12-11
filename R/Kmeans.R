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

    initialize = function(K, scale = TRUE, random_state = NULL) {
      super$initialize(
        K            = K,
        scale        = scale,
        lambda       = 1,
        method_name  = "kmeans",
        random_state = random_state
      )
    },

    interpret_clusters = function(style = c("compact", "detailed")) {
      style <- match.arg(style)

      X        <- private$FX_active
      clusters <- private$FClusters
      centers  <- private$FCenters

      if (is.null(X) || is.null(clusters) || is.null(centers))
        stop("[Kmeans] interpret_clusters(): no fitted model.")

      p <- ncol(X)
      K <- private$FNbGroupes
      var_names <- colnames(X)

      cat("=== Cluster interpretation (k-means on variables) ===\n")
      cat("Variables:", p, " | Clusters:", K, "\n\n")

      for (k in seq_len(K)) {
        vars_k <- which(clusters == k)
        cat(sprintf("--- Cluster %d ---\n", k))

        if (length(vars_k) == 0) {
          cat("  (empty)\n\n")
          next
        }

        dist_k <- numeric(length(vars_k))
        for (i in seq_along(vars_k)) {
          j <- vars_k[i]
          dist_k[i] <- sum((X[, j] - centers[, k])^2)
        }

        cat("Size:", length(vars_k), "\n")
        cat("Distances to centroid:\n")
        cat("  mean:", sprintf("%.3f", mean(dist_k)), "\n")
        cat("  min :", sprintf("%.3f", min(dist_k)), "\n")
        cat("  max :", sprintf("%.3f", max(dist_k)), "\n")

        ord <- order(dist_k)
        top <- vars_k[ord][1:min(3, length(vars_k))]

        cat("Most central variables:\n")
        for (j in top)
          cat(sprintf("  - %s\n", var_names[j]))

        if (style == "detailed") {
          cat("Interpretation:\n")
          cat("  Variables in this cluster tend to have similar numeric behaviour.\n")
        }
        cat("\n")
      }

      invisible(NULL)
    }
  ),

  private = list(

    run_clustering = function(X) {

      if (!all(vapply(X, is.numeric, logical(1L))))
        stop("[Kmeans] All variables must be numeric.")

      p <- ncol(X)
      n <- nrow(X)
      K <- private$FNbGroupes

      if (K > p)
        stop("[Kmeans] Cannot create more clusters than variables.")

      data_for_kmeans <- t(as.matrix(X))  # p × n

      km <- stats::kmeans(
        x        = data_for_kmeans,
        centers  = K,
        nstart   = 20,
        iter.max = 100
      )

      clusters <- km$cluster       # length p
      centers  <- km$centers       # K × n
      tot_ss   <- km$tot.withinss  # total WCSS

      inertia_per_cluster <- tapply(
        seq_len(p),
        clusters,
        function(idx) {
          sum(
            sapply(idx, function(j) {
              sum((data_for_kmeans[j, ] - centers[clusters[j], ])^2)
            })
          )
        }
      )

      list(
        clusters            = as.integer(clusters),
        centers             = t(centers),   # n × K
        inertia             = tot_ss,
        inertia_per_cluster = inertia_per_cluster,
        converged           = TRUE
      )
    },

    predict_one_variable = function(x_new, var_name) {

      if (!is.numeric(x_new))
        stop("[Kmeans] supplementary variable must be numeric.")

      Xref    <- private$FX_active
      centers <- private$FCenters

      if (length(x_new) != nrow(Xref))
        stop("[Kmeans] supplementary variable must have same number of rows as X.")

      d <- apply(centers, 2, function(c_k) sum((x_new - c_k)^2))

      k_best   <- which.min(d)
      d_best   <- d[k_best]
      adhesion <- 1 - d_best / max(d)

      data.frame(
        variable = var_name,
        type     = "numeric",
        cluster  = k_best,
        distance = d_best,
        adhesion = adhesion,
        stringsAsFactors = FALSE
      )
    },

    summary_membership_impl = function() {

      X        <- private$FX_active
      clusters <- private$FClusters
      centers  <- private$FCenters
      p        <- ncol(X)
      var_names <- colnames(X)

      dist_vec <- numeric(p)
      adh_vec  <- numeric(p)

      for (j in seq_len(p)) {
        c_j <- clusters[j]
        d   <- sum((X[, j] - centers[, c_j])^2)
        dist_vec[j] <- d
      }
      max_d <- max(dist_vec)
      if (max_d == 0) {
        adh_vec[] <- 1
      } else {
        adh_vec <- 1 - dist_vec / max_d
      }

      df <- data.frame(
        variable = var_names,
        cluster  = clusters,
        distance = dist_vec,
        adhesion = adh_vec,
        stringsAsFactors = FALSE
      )

      cat("=== Membership indicators (k-means) ===\n")
      print(df)
      invisible(df)
    },

    plot_membership_impl = function() {
      X        <- private$FX_active
      clusters <- private$FClusters
      centers  <- private$FCenters
      p        <- ncol(X)

      adh <- numeric(p)
      dist_vec <- numeric(p)

      for (j in seq_len(p)) {
        c_j <- clusters[j]
        d   <- sum((X[, j] - centers[, c_j])^2)
        dist_vec[j] <- d
      }
      max_d <- max(dist_vec)
      if (max_d == 0) {
        adh[] <- 1
      } else {
        adh <- 1 - dist_vec / max_d
      }

      ord <- order(adh, decreasing = TRUE)
      barplot(
        adh[ord],
        names.arg = colnames(X)[ord],
        las = 2,
        main = "Variable–Cluster Membership (k-means)",
        ylab = "Adhesion",
        cex.names = 0.7
      )
    }
  )
)
