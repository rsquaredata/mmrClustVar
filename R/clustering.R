library(R6)


.Clustering <- R6Class("Clustering",

    private = list(

      # Attributes
      .n_cluster = NULL,
      .center = NULL,
      .scale = NULL,
      .max_iter = NULL,

      # Util methods
      .perform_scale = function(X) {
        X <- as.data.frame(scale(X, center=private$.center, scale=private$.scale))
        return(X)
      },

      .perform_unscale = function(X, og_means, og_sds) {

        unscaled_X <- X

        if (private$.scale) {
          unscaled_X <- sweep(unscaled_X, 2, og_sds, `*`)
        }
        if (private$.center) {
          unscaled_X <- sweep(unscaled_X, 2, og_means, `+`)
        }

        return(unscaled_X)
      },

      .compute_W = function(X, clusters, centers) {

        W <- 0
        K <- nrow(centers)

        for (g in 1:K) {
          cluster_points <- X[clusters == g, ]
          if (nrow(cluster_points) > 0) {
            # For each point, compute the squared distance to the cluster center
            distances <- sweep(as.matrix(cluster_points), 2, as.matrix(centers[g, ]), `-`)
            W <- W + sum(distances^2)
          }
        }

        return(W)
      },

      .find_elbow = function(x, y, plot=FALSE, plot_axis=c("x", "y")) {

        perpendicular_distance <- function(x, y, slope, intercept) {
          distance <- abs(slope * x - y + intercept) / sqrt(slope^2 + 1)
          return(distance)
        }

        # Compute slope and intercept of the line (y = ax + b)
        slope <- (y[length(y)] - y[1]) / (x[length(x)] - x[1])
        intercept <- y[1] - slope * x[1]
        # Calculate distances
        distances <- perpendicular_distance(x, y, slope, intercept)
        # Find the 'elbow' point (the one with max distance)
        elbow_index <- which.max(distances)

        # Plot curve if `plot` == TRUE
        if(plot) {
          plot(x, y, type = "b", pch = 19, col = "blue",
               main = "Elbow Plot", xlab = "K", ylab = "W score")
          points(x[elbow_index], y[elbow_index], col = "red", pch = 19, cex = 1.5)
          abline(a = intercept, b = slope, lty = 2)
          legend("topright", legend = c("Data points", "Elbow point", "Fitted line"),
                 col = c("blue", "red", "black"), pch = c(19, 19, NA),
                 lty = c(0, 0, 2), cex = 0.8)
        }

        return(elbow_index)
      }
    ),

    public = list(

      # Init
      initialize = function(n_cluster, center, scale, max_iter, seed) {
        private$.n_cluster <- n_cluster
        private$.center <- center
        private$.scale <- scale
        private$.max_iter <- max_iter
        set.seed(seed) # set random seed
      },

      # Getter methods
      get.n_cluster = function() {
        return(private$.n_cluster)
      },
      get.max_iter = function() {
        return(private$.max_iter)
      },

      # Setter methods
      set.n_cluster = function(K) {
        private$.n_cluster <- K
      },

      # Abstract methods
      fit = function() {
        stop("The `fit` method should be implemented in subclasses of Clustering")
      },
      predict = function() {
        stop("The `predict` method should be implemented in subclasses of Clustering")
      }
    )

)
