library(R6)


#' @export
#' @title K-means
#' @description Implementation of the K-means variable clustering algorithm.
#' @examples
#' active.idx <- c(1, 2, 3, 4, 5)
#' descriptive.idx <- c(6, 7, 8, 9, 10)
#' model <- Kmeans$new(n_cluster=3)
#' model$fit(mtcars[, actives.idx])
#' kmeans$predict(mtcars[, descriptive.idx])
Kmeans <- R6Class("K-means",

    inherit = .Cluster,

    private = list(

      .stop_iter = function(clust_1, clust_2, iter) {
        max_iter <- self$get.max_iter()
        if (iter > max_iter) {
          warning("Maximum iteration limit reached; you can set the limit with the `max_iter` parameter.")
          return(TRUE)
        }
        return(identical(clust_1, clust_2))
      },

      .latent_center = function(X_subset) {
        # compute first latent component of a cluster (first PC score)
        if (ncol(X_subset) == 1) {
          center <- as.numeric(X_subset[, 1])
        } else {
          pca <- prcomp(X_subset, center=FALSE, scale.=FALSE, rank.=1)
          center <- as.numeric(pca$x[, 1])
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

      .allocate = function(X, n, K, centroids) {
        scores <- matrix(-Inf, nrow=n, ncol=K)
        for (k in 1:K) {
          centroid <- centroids[, k]
          if (all(is.na(centroid))) next
          corr_vec <- cor(X, centroid, use="pairwise.complete.obs")^2
          scores[, k] <- corr_vec
        }
        return(
          unname(max.col(scores, ties.method="random"))
          )
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

        previous_cluster_assignment <- rep(NA_integer_, n)
        iter <- 1

        while (!private$.stop_iter(cluster_assignment, previous_cluster_assignment, iter)) {
          previous_cluster_assignment <- cluster_assignment
          centers <- private$.compute_centers(X, cluster_assignment, K)
          cluster_assignment <- private$.allocate(X, n, K, centers)
          iter <- iter + 1
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
      #' @param max_iter An integer specifying the maximum number of iterations for the K-means algorithm (default is 300).
      #' @param random_seed An optional integer to set the random seed for reproducibility.
      initialize = function(n_cluster=NULL, max_iter=300, random_seed=NULL) {
        super$initialize(n_cluster=n_cluster, max_iter=max_iter, seed=random_seed)
      },

      #' @description
      #' Fits the K-means model to data, automatically finding the best K if not provided
      #' @param X a data.frame or matrix of actives variables to compute clusters
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

        X <- scale(X, center=TRUE, scale=TRUE)
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
          message("Set K to ", K)
        }

        result <- private$.clusterize(X, n, K)
        W <- private$.compute_W(X, clusters=result$clusters, centers=result$centers)

        private$.X <- X
        private$.clusters <- data.frame(
          cluster=result$clusters,
          row.names=colnames(X)
        )
        private$.centers <- result$centers

        # TODO: return clusters as data.frame (variable name as column name, cluster assignment as value)

        return(list(
          clusters = private$.clusters,
          W = W
        ))
      },

      #' @description
      #' Assign descriptive variables to clusters
      #' @param descriptives A data.frame or matrix of descriptive variables to assign to cluster.
      #' @return An integer vector indicating, for each provided variable,
      #'         the index of the cluster to which it is assigned.
      predict = function(descriptives) {

        if (is.null(private$.centers)) {
          stop("You must create clusters with `fit()` before predicting cluster assignments.")
        }

        if (!is.data.frame(descriptives) && !is.matrix(descriptives)) {
          stop("Wrong type for `descriptives`. Expected data.frame or matrix.")
        }
        if (!is.numeric(as.matrix(descriptives))) {
          stop("Non numeric values found in `descriptives`, expected numeric only.")
        }

        n <- ncol(descriptives)
        K <- self$get.n_cluster()

        assignments <- private$.allocate(X=descriptives, n=n, K=K, centroids=private$.centers)

        # Append descriptive values to the .X object
        private$.descriptives <- descriptives
        # Saves predicted clusters
        private$.predicts <- data.frame(
          cluster=assignments,
          row.names=colnames(descriptives)
        )

        return(data.frame(
          cluster=assignments,
          row.names=colnames(descriptives)
        ))
      },

      #' @description
      #' Project clusters on factorial plane
      visualize = function() {
        pca <- prcomp(private$.X, center=FALSE, scale.=FALSE)
        loadings <- as.data.frame(pca$rotation[, 1:2])
        colnames(loadings) <- c("Comp.1", "Comp.2")
        loadings$variable <- colnames(private$.X)
        loadings$cluster <- factor(private$.clusters[,1])

        ggplot2::ggplot(loadings, ggplot2::aes(x=Comp.1, y=Comp.2, label=variable, color=cluster)) +
          ggplot2::geom_hline(yintercept=0, linetype="dashed", color="grey70") +
          ggplot2::geom_vline(xintercept=0, linetype="dashed", color="grey70") +
          ggplot2::geom_segment(ggplot2::aes(x=0, y=0, xend=Comp.1, yend=Comp.2), arrow=ggplot2::arrow(length=ggplot2::unit(0.15, "cm")), show.legend=FALSE, alpha=0.5) +
          ggplot2::geom_point(size=2) +
          ggplot2::geom_text(ggplot2::aes(x=1.1*Comp.1, y=1.1*Comp.2, label=variable, color=cluster), size=3) +
          ggplot2::labs(title = "Variable clusters on PCA plane", x="Comp.1", y="Comp.2") +
          ggplot2::theme_minimal()
      },

      # Getters

      #' @description
      #' Get active variables cluster assignments
      get.clusters = function() return(private$.clusters),

      #' @description
      #' Get predicted descriptive cluster assignments
      get.predicts = function() return(private$.predicts)
    )

)
