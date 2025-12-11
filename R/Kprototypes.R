#' K-prototypes for variable clustering (mixed data)
#'
#' R6 class implementing a k-prototypes algorithm for clustering
#' variables mixing numeric and categorical data.
#'
#' @docType class
#' @name Kprototypes
#' @export
Kprototypes <- R6::R6Class(
  "Kprototypes",
  inherit = .ClusterBase,

  public = list(

    initialize = function(K, scale = TRUE, lambda = 1) {

      if (!is.numeric(lambda) || lambda <= 0) {
        stop("[Kprototypes] lambda must be numeric > 0.")
      }

      super$initialize(
        K           = K,
        scale       = scale,
        lambda      = lambda,
        method_name = "kprototypes"
      )
    },

    interpret_clusters = function(style = c("compact", "detailed")) {
      # identique à ta version — non modifiée
      stop("interpret_clusters() left unchanged — paste your version here.")
    }
  ),

  private = list(

    # =============================================================
    # K-prototypes core algorithm
    # =============================================================
    run_clustering = function(X) {

      num_idx <- private$FNumCols
      cat_idx <- private$FCatCols
      lambda  <- private$FLambda
      K       <- private$FNbGroupes

      p <- ncol(X)

      X_num <- if (length(num_idx) > 0) as.matrix(X[, num_idx, drop = FALSE]) else NULL
      X_cat <- if (length(cat_idx) > 0) as.matrix(X[, cat_idx, drop = FALSE]) else NULL

      # -----------------------------------
      # 1. init prototypes
      # -----------------------------------
      centers_num <- if (!is.null(X_num)) {
        X_num[, sample(seq_len(ncol(X_num)), K), drop = FALSE]
      } else NULL

      centers_cat <- if (!is.null(X_cat)) {
        lapply(seq_len(K), function(k) {
          X_cat[, sample(seq_len(ncol(X_cat)), 1)]
        })
      } else NULL

      clusters <- rep(NA_integer_, p)
      converged <- FALSE

      iter <- 0
      max_iter <- 50

      repeat {

        iter <- iter + 1
        old_clusters <- clusters

        # -----------------------------------
        # 2. assign each variable to closest prototype
        # -----------------------------------
        dist_jk <- matrix(0, nrow = p, ncol = K)

        for (j in seq_len(p)) {
          for (k in seq_len(K)) {

            d_num <- if (!is.null(X_num)) {
              sum((X_num[, j] - centers_num[, k])^2)
            } else 0

            d_cat <- if (!is.null(X_cat)) {
              m <- X_cat[, j]
              c <- centers_cat[[k]]
              sum(m != c)
            } else 0

            dist_jk[j, k] <- d_num + lambda * d_cat
          }
        }

        clusters <- max.col(-dist_jk)  # best cluster

        # -----------------------------------
        # 3. recompute prototypes
        # -----------------------------------
        for (k in seq_len(K)) {
          vars_k <- which(clusters == k)
          if (length(vars_k) == 0) next

          if (!is.null(X_num)) {
            centers_num[, k] <- rowMeans(X_num[, vars_k, drop = FALSE])
          }

          if (!is.null(X_cat)) {
            centers_cat[[k]] <- apply(X_cat[, vars_k, drop = FALSE], 1, function(v) {
              names(which.max(table(v)))
            })
          }
        }

        if (identical(clusters, old_clusters) || iter >= max_iter) {
          converged <- TRUE
          break
        }
      }

      # -----------------------------------
      # 4. compute inertia
      # -----------------------------------
      inertia_total <- sum(sapply(seq_len(p), function(j) {
        k <- clusters[j]

        d_num <- if (!is.null(X_num)) {
          sum((X_num[, j] - centers_num[, k])^2)
        } else 0

        d_cat <- if (!is.null(X_cat)) {
          sum(X_cat[, j] != centers_cat[[k]])
        } else 0

        d_num + lambda * d_cat
      }))

      # -----------------------------------
      # inertia per cluster
      # -----------------------------------
      inertia_per_cluster <- tapply(
        seq_len(p),
        clusters,
        function(idx) sum(sapply(idx, function(j) {
          k <- clusters[j]

          d_num <- if (!is.null(X_num)) {
            sum((X_num[, j] - centers_num[, k])^2)
          } else 0

          d_cat <- if (!is.null(X_cat)) {
            sum(X_cat[, j] != centers_cat[[k]])
          } else 0

          d_num + lambda * d_cat
        }))
      )

      return(list(
        clusters  = clusters,
        centers   = list(num = centers_num, cat = centers_cat),
        inertia   = inertia_total,
        inertia_per_cluster = inertia_per_cluster,
        converged = converged
      ))
    }
  )
)
