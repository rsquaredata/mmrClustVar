library(R6)


#' @export
#' @title K-means
#' @description Implementation of the K-means clustering algorithm.
#' @examples
#' kmeans_model <- Kmeans$new(n_cluster=5)
#' result <- kmeans_model$fit(mtcars)
#' print(result$clusters)
#' kmeans_model$scatterplot(mtcars, 'mpg', 'hp')
Kmeans <- R6Class("K-means",

  inherit = .Clustering,

  private = list(

    .centers = NULL,
    .clusters = NULL,

    .stop_iter = function(clust_1, clust_2, iter) {

      max_iter <- self$get.max_iter()
      if(iter > max_iter) {
        return(TRUE)
      }

      if(identical(clust_1, clust_2)) {
        return(TRUE)
      }

      return(FALSE)
    },

    .allocate = function(X, K, centers) {
      # Compute distance to centroids
      distance_matrix <- as.matrix(dist(rbind(centers, X), method="euclidean"))
      distance_matrix <- distance_matrix[-(1:K), 1:K]
      # Assign to cluster
      cluster_assignment <- apply(distance_matrix, 1, which.min)
      return(cluster_assignment)
    },

    .clusterize = function(X, n, K) {

      # Initialize K random centers
      centers <- X[sample(1:n, K), ]
      rownames(centers) <- NULL
      cluster_assignment <- rep(0, n)
      previous_cluster_assignment <- rep(-1, n)

      # Iterate until convergence
      iter = 0
      while(!private$.stop_iter(cluster_assignment, previous_cluster_assignment, iter=iter)) {
        previous_cluster_assignment <- cluster_assignment
        iter = iter+1
        # Allocation step
        cluster_assignment <- private$.allocate(X, K=K, centers=centers)
        # Recalculate centers
        for (g in 1:K) {
          if (sum(cluster_assignment == g) > 0) {
            centers[g, ] <- colMeans(X[cluster_assignment == g, ])
          }
        }
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
    #' @param latent A logical indicating whether to use latent variables computed on PCA (default is `FALSE`).
    #' @param center A logical indicating whether to center the data (default is `TRUE`).
    #' @param scale A logical indicating whether to scale the data (default is `TRUE`).
    #' @param max_iter An integer specifying the maximum number of iterations for the K-means algorithm (default is 300).
    #' @param random_seed An optional integer to set the random seed for reproducibility.
    initialize = function(n_cluster=NULL, latent=FALSE, center=TRUE, scale=TRUE, max_iter=300, random_seed=NULL) {
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

      # Check if X is a data.frame or matrix
      if(!is.data.frame(X) && !is.matrix(X)){
        print(paste("Wrong type for `X`. Expected data.frame got ", class(X)))
        return(NULL)
      }
      # Check if X is numeric
      if(!is.numeric(as.matrix(X))) {
        print("Non numeric values found in `X`, expected numeric only")
        return(NULL)
      }

      X <- private$.perform_scale(X, set=TRUE) # center and scale X if `center` and `scale` == TRUE
      n <- nrow(X)
      K <- self$get.n_cluster()

      if(is.null(K)) {

        # Find best K number of cluster

        compute.W_for_k <- function(k) {
          result <- private$.clusterize(X, n, k)
          W <- private$.compute_W(X, clusters=result$clusters, centers=result$centers)
          return(W)
        }

        # Compute W score for K = 2 to 10
        K_values <- 2:15
        W_values <- sapply(K_values, compute.W_for_k)
        # Find curve elbow (best result)
        elbow_index <- private$.find_elbow(K_values, W_values, plot=TRUE, plot_axis=c("K", "W score"))
        K <- K_values[elbow_index]
        self$set.n_cluster(K)
      }

      # Compute cluster and W score
      result <- private$.clusterize(X, n, K)
      W <- private$.compute_W(X, clusters=result$clusters, centers=result$centers)
      # Save private attributes
      private$.clusters <- result$clusters
      private$.centers <- result$centers

      # TODO: compute "proportion de variance expliquée par le partitionnement" for each group

      return(list(
        clusters = private$.clusters,
        centers = private$.perform_unscale(result$centers),
        W = W
      ))
    },

    #' @description
    #' Predict the cluster group for each observation of the data
    #' @param X A data.frame or matrix on which to predict cluster groups
    #' @return `clusters`: An integer vector indicating the cluster to which each point is allocated
    predict = function(X) {

      X <- private$.perform_scale(X)
      K <- self$get.n_cluster()
      clusters <- private$.allocate(X, K=K, centers=private$.centers)
      W <- private$.compute_W(X, clusters=clusters, centers=private$.centers)

      private$.clusters <- c(private$.clusters, clusters)

      return(list(
        clusters=clusters,
        W=W
      ))
    },

    #' @description
    #' Draw a scatterplot highlighting the cluster groups on 2 variables of the data
    #' @param X A data.frame or matrix used for computing the clusters
    #' @param var1 A string representing the column name to display on x axis
    #' @param var2 A string representing the column name to display on y axis
    scatterplot = function(X, var1, var2) {

      # convert to a sequence
      cluster_colors <- as.numeric(private$.clusters[rownames(X)])
      # unscale centers
      centers <- private$.perform_unscale(private$.centers)

      plot(X[[var1]], X[[var2]],
         col = cluster_colors,
         pch = 19,
         xlab = var1,
         ylab = var2,
         main = "K-means Clustering")

      # Add cluster centers to the plot
      points(centers[, var1], centers[, var2],
         col = 1:self$get.n_cluster(),
         pch = 4,
         cex = 2,
         lwd = 2)
    }
  )

)




