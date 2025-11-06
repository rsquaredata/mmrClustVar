library(R6)


#' @export
Kmeans <- R6Class("K-means",

  inherit = .Clustering,

  public = list(

    initialize = function(n_cluster = NA) {
      super$initialize(n_cluster = n_cluster)
    },

    fit = function(data) {
      # test
      n_cluster <- self$get_n_cluster()
      return(n_cluster)
    }
  )

)

t <- Kmeans$new(n_cluster=5)
hello <- t$fit()
