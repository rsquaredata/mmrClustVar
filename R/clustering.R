library(R6)


.Clustering <- R6Class("Clustering",

    private = list(
      # Attributes
      .n_cluster = NULL
    ),

    public = list(

      # Init
      initialize = function(n_cluster = NA) {
        private$.n_cluster <- n_cluster
      },

      # Getter methods
      get_n_cluster = function() {
        return(private$.n_cluster)
      },

      # Abstract methods
      fit = function() {
        stop("The `fit` method should be implemented in subclasses of Clustering")
      }
    )

)
