library(R6)



.PCA <- R6Class("PCA",

  private = list(
    .X = NA
  ),

  public = list(

    initialize = function(data) {

      df.ok <- is.data.frame(data)
      if(!df.ok) {
        stop(paste("Invalid `data` type. Expected data.frame, got", class(data)))
      }
      nb.ok <- sum(sapply(data,function(x){is.numeric(x)}))
      if(!nb.ok) {
        stop("Non numeric `data` values")
      }

      private$.X <- data
    }

    perform_pca = function(n_comps=NA, center=TRUE) {

      # Center the data
      X <- scale(private$.X, center=center, scale=FALSE)

      # Calculate covariance matrix
      covariance_matrix <- cov(X)

      # Get eigenvalues and eigenvectors
      eigen_decomposition <- eigen(covariance_matrix)
      eigenvalues <- eigen_decomposition$values
      eigenvectors <- eigen_decomposition$vectors
    }

  )

)
