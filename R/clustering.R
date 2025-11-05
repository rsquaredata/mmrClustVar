library(R6)


.Clustering <- R6Class("Clustering",

    public = list(

      n_cluster = NULL,

      initialize = function(n_cluster = NA) {
        self$n_cluster <- n_cluster
      },

      fit = function() {
        stop()
      }
    )

)
