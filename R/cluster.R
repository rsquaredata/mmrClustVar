library(R6)

Clustering <- R6Class("Clustering",
                      
                      private = list(
                        .n_cluster = NULL,
                        .seed = NULL
                      ),
                      
                      public = list(
                        
                        initialize = function(n_cluster = NULL, seed = NULL) {
                          
                          # n_cluster peut être NULL → auto-K dans les classes enfants
                          if (!is.null(n_cluster)) {
                            if (!is.numeric(n_cluster) || n_cluster <= 0) {
                              stop("`n_cluster` must be a strictly positive integer or NULL.")
                            }
                            private$.n_cluster <- as.integer(n_cluster)
                          }
                          
                          # seed optionnel
                          if (!is.null(seed)) {
                            private$.seed <- as.integer(seed)
                            set.seed(private$.seed)
                          }
                        },
                        
                        # Getters
                        get_n_cluster = function() private$.n_cluster,
                        get_seed = function() private$.seed,
                        
                        # Setter
                        set_n_cluster = function(k) {
                          if (!is.numeric(k) || k <= 0) stop("`k` must be > 0")
                          private$.n_cluster <- as.integer(k)
                        },
                        
                        # Abstract methods
                        fit = function(X) stop("fit() must be implemented in child classes."),
                        predict = function(newdata) stop("predict() must be implemented."),
                        
                        print = function(...) {
                          cat("<Clustering object>\n")
                          cat(" - n_cluster:", private$.n_cluster, "\n")
                          if (!is.null(private$.seed)) cat(" - seed     :", private$.seed, "\n")
                          invisible(self)
                        }
                      )
)
