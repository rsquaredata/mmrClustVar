library(R6)


#' @export
Kmeans <- R6Class("K-means",

  inherit = .Clustering,

  public = list(

    initialize = function(n_cluster = NA) {
      super$initialize()
    },

    fit = function() {
      # test
      return(self$n_cluster)
    }
  )

)
