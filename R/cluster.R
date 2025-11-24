library(R6)


.Cluster <- R6Class("Cluster",

    private = list(

      # Attributes
      .X = NULL,            # active data
      .descriptives = NULL, # descriptive data
      .clusters = NULL,     # cluster assignments
      .predicts = NULL,     # predicted cluster assignments
      .n_cluster = NULL,    # number of clusters
      .max_iter = NULL,     # maximum iteration

      # Getter methods
      get.n_cluster = function() return(private$.n_cluster),
      get.max_iter = function() return(private$.max_iter),

      # Setter methods
      set.n_cluster = function(K) private$.n_cluster <- K,

      # Util methods
      .compute_W = function(X, clusters, centers) {
        if (!is.matrix(X)) {
          X <- as.matrix(X)
        }

        n_vars <- ncol(X)
        if (length(clusters) != n_vars) {
          stop("`clusters` must have the same length as the number of variables in `X`.")
        }

        K <- ncol(centers)
        if (is.null(K) || K == 0L) {
          stop("`centers` must have at least one column.")
        }

        W <- 0
        for (g in 1:K) {
          idx <- which(clusters == g)
          if (length(idx) == 0L) next

          cluster_vars <- X[, idx, drop=FALSE]
          center_vec  <- centers[, g]

          if (all(is.na(center_vec))) next

          distances <- sweep(cluster_vars, 1, center_vec, FUN="-")
          W <- W + sum(distances^2, na.rm=TRUE)
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
          plot(x, y, type="b", pch=19, col="blue", main="Elbow Plot", xlab=plot_axis[1], ylab=plot_axis[2])
          points(x[elbow_index], y[elbow_index], col="red", pch=19, cex=1.5)
          abline(a=intercept, b=slope, lty=2)
          legend("topright", legend=c("Data points", "Elbow point", "Fitted line"), col=c("blue", "red", "black"), pch=c(19, 19, NA), lty=c(0, 0, 2), cex=0.8)
        }

        return(elbow_index)
      },

      .silhouette = function(X, clusters, squared=TRUE, use_abs=TRUE) {
        X <- as.matrix(X)
        n_vars <- ncol(X)
        if (length(clusters) != n_vars) {
          stop("`clusters` must have the same length as the number of variables in `X`.")
        }

        # similarity = |cor| or cor, optionally squared
        corr <- cor(X, use = "pairwise.complete.obs")
        if (use_abs) corr <- abs(corr)
        if (squared) corr <- corr^2

        # convert to dissimilarity
        diss <- 1 - corr
        diag(diss) <- 0

        uniq_clusters <- sort(unique(clusters))
        sil <- numeric(n_vars)

        for (i in seq_len(n_vars)) {
          k_i <- clusters[i]

          # average dissimilarity within the same cluster
          same_idx <- which(clusters == k_i & seq_len(n_vars) != i)
          a_i <- if (length(same_idx) == 0L) 0 else mean(diss[i, same_idx])

          # minimum average dissimilarity to other clusters
          b_i <- Inf
          for (k in uniq_clusters[uniq_clusters != k_i]) {
            other_idx <- which(clusters == k)
            if (length(other_idx) == 0L) next
            b_i <- min(b_i, mean(diss[i, other_idx]))
          }
          if (!is.finite(b_i)) b_i <- 0

          denom <- max(a_i, b_i)
          sil[i] <- if (denom == 0) 0 else (b_i - a_i) / denom
        }

        return(list(
          per_variable = data.frame(
            variable = colnames(X),
            cluster = clusters,
            silhouette = sil,
            row.names = NULL
          ),
          mean = mean(sil)
        ))
      }
    ),

    public = list(

      # Init
      initialize = function(n_cluster, max_iter, seed) {
        private$.n_cluster <- n_cluster
        private$.max_iter <- max_iter
        set.seed(seed) # set random seed
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
