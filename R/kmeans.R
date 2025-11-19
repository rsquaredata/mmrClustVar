library(R6)
library(ggplot2)


#' @export
#' @title K-means
#' @description Implementation of the K-means variable clustering algorithm.
#' @examples
#' kmeans_model <- Kmeans$new(n_cluster=5)
#' result <- kmeans_model$fit(mtcars)
#' print(result$clusters)
Kmeans <- R6Class("K-means",

    inherit = .Clustering,

    private = list(

      .X = NULL,
      .centers = NULL,
      .clusters = NULL,

      .stop_iter = function(clust_1, clust_2, iter) {
        max_iter <- self$get.max_iter()
        if (iter > max_iter) {
          return(TRUE)
        }
        return(identical(clust_1, clust_2))
      },

      .latent_center = function(X_subset) {
        # compute first latent component of a cluster (first PC score)
        if (ncol(X_subset) == 1L) {
          center <- as.numeric(scale(X_subset[, 1], center=TRUE, scale=TRUE))
        } else {
          pca <- prcomp(X_subset, center=FALSE, scale.=FALSE, rank.=1L)
          center <- as.numeric(scale(pca$x[, 1], center=TRUE, scale=TRUE))
        }
        if (anyNA(center) || sd(center, na.rm=TRUE) == 0) {
          center <- rep(0, nrow(X_subset))
        }
        return(center)
      },

      .compute_centers = function(X, clusters, K) {
        n_obs <- nrow(X)
        centers <- matrix(NA_real_, nrow=n_obs, ncol=K)
        for (k in 1:K) {
          idx <- which(clusters == k)
          if (length(idx) == 0L) {
            next
          }
          centers[, k] <- private$.latent_center(X[, idx, drop=FALSE])
        }
        return(centers)
      },

      .allocate = function(X, n, K, centroids, square=TRUE) {
        scores <- matrix(-Inf, nrow=n, ncol=K)
        for (k in 1:K) {
          centroid <- centroids[, k]
          if (all(is.na(centroid))) next
          corr_vec <- cor(X, centroid, use="pairwise.complete.obs")
          if (square) corr_vec <- corr_vec^2
          scores[, k] <- corr_vec
        }
        return(max.col(scores, ties.method="random"))
      },

      .ensure_non_empty_clusters = function(assignments, K) {
        tab <- tabulate(assignments, nbins=K)
        empty <- which(tab == 0L)
        if (length(empty) == 0L) return(assignments)

        for (k in empty) {
          donor <- which.max(tab)
          donor_idx <- which(assignments == donor)
          move <- sample(donor_idx, 1L)
          assignments[move] <- k
          tab[donor] <- tab[donor] - 1L
          tab[k] <- 1L
        }
        return(assignments)
      },

      .clusterize = function(X, n, K) {
        n_obs <- nrow(X)

        # initial seeds (random variables used)
        seeds <- sample(1:n, K)
        centers <- vapply(seeds,
                          function(j) as.numeric(scale(X[, j], center=TRUE, scale=TRUE)),
                          FUN.VALUE = numeric(n_obs))
        centers <- matrix(centers, nrow=n_obs, ncol=K)

        cluster_assignment <- private$.allocate(X, n, K, centers)
        cluster_assignment <- private$.ensure_non_empty_clusters(cluster_assignment, K)

        previous_cluster_assignment <- rep(NA_integer_, n)
        iter <- 1L

        while (!private$.stop_iter(cluster_assignment, previous_cluster_assignment, iter)) {
          previous_cluster_assignment <- cluster_assignment

          centers <- private$.compute_centers(X, cluster_assignment, K)

          # if some centers are NA (empty cluster), re-init with random variable
          na_centers <- which(apply(centers, 2, function(col) all(is.na(col))))
          if (length(na_centers) > 0L) {
            for (k in na_centers) {
              j <- sample(1:n, 1L)
              centers[, k] <- as.numeric(scale(X[, j], center = TRUE, scale = TRUE))
            }
          }

          cluster_assignment <- private$.allocate(X, n, K, centers)
          cluster_assignment <- private$.ensure_non_empty_clusters(cluster_assignment, K)
          iter <- iter + 1L
        }

        return(list(
          clusters = cluster_assignment,
          centers = centers
        ))
      }
    ),

    public = list(

      #' @description
      #' Create a new K-means object.
      #' @param n_cluster An optional integer specifying the number of clusters. If `NULL`, the model will determine the optimal number of clusters.
      #' @param center A logical indicating whether to center the data (default is `TRUE`).
      #' @param scale A logical indicating whether to scale the data (default is `TRUE`).
      #' @param max_iter An integer specifying the maximum number of iterations for the K-means algorithm (default is 300).
      #' @param random_seed An optional integer to set the random seed for reproducibility.
      initialize = function(n_cluster=NULL, center=TRUE, scale=TRUE, max_iter=300, random_seed=NULL) {
        super$initialize(n_cluster=n_cluster, center=center, scale=scale, max_iter=300, seed=random_seed)
      },

      #' @description
      #' Fits the K-means model to data, automatically finding the best K if not provided
      #' @param X A data.frame or matrix on which to perform clustering
      #' @return A list containing the following components:
      #' - `clusters`: An integer vector indicating the cluster to which each point is allocated.
      #' - `centers`: A matrix with the coordinates of cluster centers.
      #' - `W`: The within-cluster sum of squares for the fitted clusters.
      fit = function(X) {
        if (!is.data.frame(X) && !is.matrix(X)) {
          stop("Wrong type for `X`. Expected data.frame or matrix.")
        }
        if (!is.numeric(as.matrix(X))) {
          stop("Non numeric values found in `X`, expected numeric only.")
        }

        X <- private$.perform_scale(X, set=TRUE)
        n <- ncol(X)
        K <- self$get.n_cluster()

        if (is.null(K)) {

          compute.W_for_k <- function(k) {
            result <- private$.clusterize(X, n, k)
            W <- private$.compute_W(X, clusters=result$clusters, centers=result$centers)
            return(W)
          }

          compute.silhouette_for_K <- function(k) {
            result <- private$.clusterize(X, n, k)
            silhouette <- private$.silhouette(X, result$clusters)
            return(silhouette$mean)
          }

          K_values <- 2:min(ncol(X), 10)
          W_values <- sapply(K_values, compute.W_for_k)

          #silhouette_values <- sapply(K_values, compute.silhouette_for_K)
          elbow_index <- private$.find_elbow(K_values, W_values, plot=TRUE, plot_axis=c("K", "W score"))
          K <- K_values[elbow_index]
          self$set.n_cluster(K)
        }

        result <- private$.clusterize(X, n, K)
        W <- private$.compute_W(X, clusters=result$clusters, centers=result$centers)

        private$.clusters <- result$clusters
        private$.centers <- result$centers
        private$.X <- X

        return(list(
          clusters = private$.clusters,
          centers = result$centers,
          W = W
        ))
      },

      #' @description
      #' Project clusters on factorial plane
      visualize = function() {
        pca <- prcomp(private$.X, center=FALSE, scale.=FALSE)
        loadings <- as.data.frame(pca$rotation[, 1:2])
        colnames(loadings) <- c("Comp.1", "Comp.2")
        loadings$variable <- colnames(private$.X)
        loadings$cluster <- factor(private$.clusters)

        ggplot(loadings, aes(x=Comp.1, y=Comp.2, label=variable, color=cluster)) +
          geom_hline(yintercept=0, linetype="dashed", color="grey70") +
          geom_vline(xintercept=0, linetype="dashed", color="grey70") +
          geom_segment(aes(x=0, y=0, xend=Comp.1, yend=Comp.2), arrow=arrow(length=unit(0.15, "cm")), show.legend=FALSE, alpha=0.5) +
          geom_point(size=2) +
          geom_text(aes(x=1.1*Comp.1, y=1.1*Comp.2, label=variable, color=cluster), size=3) +
          labs(title = "Variable clusters on PCA plane", x="Comp.1", y="Comp.2") +
          theme_minimal()
      },

      silhouette = function() {
        private$.silhouette(private$.X, private$.clusters)
      }
    )

)



