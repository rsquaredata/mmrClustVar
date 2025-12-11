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

    initialize = function(K, lambda = 1, scale = FALSE, random_state = NULL) {
      super$initialize(
        K            = K,
        scale        = FALSE,   # kmodes ne scale jamais
        lambda       = lambda,
        method_name  = "kmodes",
        random_state = random_state
      )
    },

    interpret_clusters = function(style = c("compact", "detailed")) {
      style <- match.arg(style)

      X        <- private$FX_active
      clusters <- private$FClusters
      centers  <- private$FCenters

      if (is.null(X) || is.null(clusters) || is.null(centers))
        stop("[Kmodes] interpret_clusters(): no fitted model.")

      p      <- ncol(X)
      K      <- private$FNbGroupes
      vnames <- colnames(X)

      cat("=== Cluster interpretation (k-modes) ===\n")
      cat("Variables:", p, " | Clusters:", K, "\n\n")

      Xc <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)

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
          dist_k[i] <- private$simple_matching(Xc[[j]], centers[[k]])
        }

        cat("Size:", length(vars_k), "\n")
        cat("Distances to prototype:\n")
        cat("  mean:", sprintf("%.3f", mean(dist_k)), "\n")
        cat("  min :", sprintf("%.3f", min(dist_k)), "\n")
        cat("  max :", sprintf("%.3f", max(dist_k)), "\n")

        ord <- order(dist_k)
        top <- vars_k[ord][1:min(3, length(vars_k))]

        cat("Most central variables:\n")
        for (j in top) {
          d <- private$simple_matching(Xc[[j]], centers[[k]])
          cat(sprintf("  - %s (d=%.3f)\n", vnames[j], d))
        }

        if (style == "detailed") {
          cat("Interpretation:\n")
          cat("  Variables in this cluster share similar categorical patterns.\n")
        }
        cat("\n")
      }

      invisible(NULL)
    }
  ),

  private = list(

    run_clustering = function(X) {

      if (!all(vapply(X, function(z) is.factor(z) || is.character(z), logical(1L))))
        stop("[Kmodes] All variables must be categorical.")

      p <- ncol(X)
      n <- nrow(X)
      K <- private$FNbGroupes

      if (K > p)
        stop("[Kmodes] Cannot create more clusters than variables.")

      Xc <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)

      # init: K variables aléatoires comme prototypes
      init_idx <- sample(seq_len(p), K)
      prototypes <- vector("list", K)
      for (k in seq_len(K)) {
        prototypes[[k]] <- Xc[[init_idx[k]]]
      }

      clusters <- rep(NA_integer_, p)
      changed  <- TRUE
      iter     <- 0L
      max_iter <- 50L

      while (changed && iter < max_iter) {
        iter <- iter + 1L
        old_clusters <- clusters

        # assignation
        clusters <- sapply(seq_len(p), function(j) {
          d <- sapply(prototypes, function(proto) {
            private$simple_matching(Xc[[j]], proto)
          })
          which.min(d)
        })

        # mise à jour des prototypes : mode ligne par ligne
        for (k in seq_len(K)) {
          vars_k <- which(clusters == k)
          if (length(vars_k) == 0) next

          new_proto <- character(n)
          for (i in seq_len(n)) {
            vals <- Xc[vars_k][i, ]
            tab  <- table(na.omit(unlist(vals)))
            new_proto[i] <- names(which.max(tab))
          }
          prototypes[[k]] <- new_proto
        }

        changed <- !identical(clusters, old_clusters)
      }

      # inertie totale + par cluster
      inertia <- 0
      inertia_per_cluster <- numeric(K)

      for (k in seq_len(K)) {
        vars_k <- which(clusters == k)
        if (length(vars_k) == 0) next

        d_k <- sum(
          sapply(vars_k, function(j) {
            private$simple_matching(Xc[[j]], prototypes[[k]])
          })
        )

        inertia_per_cluster[k] <- d_k
        inertia <- inertia + d_k
      }

      list(
        clusters            = clusters,
        centers             = prototypes,
        inertia             = inertia,
        inertia_per_cluster = inertia_per_cluster,
        converged           = TRUE
      )
    },

    predict_one_variable = function(x_new, var_name) {

      if (!(is.factor(x_new) || is.character(x_new)))
        stop("[Kmodes] supplementary variable must be categorical.")

      x_new   <- as.character(x_new)
      centers <- private$FCenters
      K       <- length(centers)

      d <- sapply(centers, function(proto) {
        private$simple_matching(x_new, proto)
      })

      k_best   <- which.min(d)
      d_best   <- d[k_best]
      adhesion <- 1 - d_best

      data.frame(
        variable = var_name,
        type     = "categorical",
        cluster  = k_best,
        distance = d_best,
        adhesion = adhesion,
        stringsAsFactors = FALSE
      )
    },

    summary_membership_impl = function() {

      Xc       <- as.data.frame(lapply(private$FX_active, as.character), stringsAsFactors = FALSE)
      clusters <- private$FClusters
      centers  <- private$FCenters
      vnames   <- colnames(Xc)
      p        <- ncol(Xc)

      dist_vec <- numeric(p)
      adh_vec  <- numeric(p)

      for (j in seq_len(p)) {
        k <- clusters[j]
        d <- private$simple_matching(Xc[[j]], centers[[k]])
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

      cat("=== Membership indicators (k-modes) ===\n")
      print(df)
      invisible(df)
    },

    plot_membership_impl = function() {

      Xc       <- as.data.frame(lapply(private$FX_active, as.character), stringsAsFactors = FALSE)
      clusters <- private$FClusters
      centers  <- private$FCenters
      p        <- ncol(Xc)

      adh <- numeric(p)
      for (j in seq_len(p)) {
        k  <- clusters[j]
        d  <- private$simple_matching(Xc[[j]], centers[[k]])
        adh[j] <- 1 - d
      }

      ord <- order(adh, decreasing = TRUE)
      barplot(
        adh[ord],
        names.arg = colnames(Xc)[ord],
        las = 2,
        main = "Variable–Cluster Membership (k-modes)",
        ylab = "Adhesion",
        cex.names = 0.7
      )
    }
  )
)
