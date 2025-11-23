library(R6)


Kmodes <- R6Class("Kmodes",

   inherit = .Clustering,

   private = list(

     .modes    = NULL,   # centers (variables fictives)
     .clusters = NULL,   # cluster assignment of variables


     # DISTANCE entre deux variables (colonnes)

     .col_dissimilarity = function(x, y) {
       sum(x != y)
     },

     # Init CAO pour variables

     .cao_init = function(X, K) {

       p <- ncol(X)        # nb de variables
       n <- nrow(X)
       modes <- matrix(NA, nrow = K, ncol = n)
       colnames(modes) <- paste0("ind_", 1:n)

       # 1er mode : variable la plus fréquente
       freq_score <- sapply(1:p, function(j) {
         sum(table(X[, j]) / n)
       })

       modes[1, ] <- as.character(X[, which.max(freq_score)])

       # Modes suivants : max distance cumulée aux modes précédents
       for (k in 2:K) {
         dist_to_existing <- sapply(1:p, function(j) {
           sum(sapply(1:(k - 1), function(m) {
             private$.col_dissimilarity(
               as.character(X[, j]),
               as.character(modes[m, ])
             )
           }))
         })

         modes[k, ] <- as.character(X[, which.max(dist_to_existing)])
       }

       return(as.data.frame(modes, stringsAsFactors = FALSE))
     },


     # Mise à jour des modes pour les variables

     .update_modes = function(X, clusters, K) {

       n <- nrow(X)
       new_modes <- matrix(NA, nrow = K, ncol = n)
       new_modes <- as.data.frame(new_modes, stringsAsFactors = FALSE)

       for (k in 1:K) {

         vars_in_cluster <- which(clusters == k)

         if (length(vars_in_cluster) == 0) {
           # garder ancien mode
           new_modes[k, ] <- private$.modes[k, ]
         } else {

           # Pour chaque individu i :
           # mode = modalité la plus fréquente parmi les variables du cluster
           for (i in 1:n) {
             col_vals <- as.character(X[i, vars_in_cluster])
             new_modes[k, i] <- names(which.max(table(col_vals)))
           }
         }
       }

       return(as.data.frame(new_modes, stringsAsFactors = FALSE))
     },

     # Compute cost

     .compute_cost = function(X, modes, clusters) {

       p <- ncol(X)
       n <- nrow(X)

       sum(sapply(1:p, function(j) {
         k <- clusters[j]
         private$.col_dissimilarity(
           as.character(X[, j]),
           as.character(modes[k, ])
         )
       }))
     },

     # Sélection automatique du K optimal

     .auto_select_K = function(X) {

       K_values <- 2:6
       costs <- numeric(length(K_values))

       for (i in seq_along(K_values)) {
         K <- K_values[i]
         modes_K <- private$.cao_init(X, K)

         clusters_K <- sapply(1:ncol(X), function(j) {
           which.min(sapply(1:K, function(m) {
             private$.col_dissimilarity(
               as.character(X[, j]),
               as.character(modes_K[m, ])
             )
           }))
         })

         costs[i] <- private$.compute_cost(X, modes_K, clusters_K)
       }

       return(K_values[which.min(costs)])
     }
   ),

   public = list(

     initialize = function(n_cluster=NULL, center=TRUE, scale=TRUE, max_iter=300, random_seed=NULL) {
       super$initialize(n_cluster=n_cluster, center=center, scale=scale, max_iter=300, seed=random_seed)
     },

     # Fit

     fit = function(X) {

       X <- as.data.frame(X)

       if (!all(sapply(X, is.factor))) {
         stop("Toutes les variables doivent être qualitatives (facteurs).")
       }

       # Choix de K automatique si absent
       K <- self$get.n_cluster()
       if (is.null(K)) {
         K <- private$.auto_select_K(X)
         message("K auto choisi : ", K)
         self$set.n_cluster(K)
       }

       # Init CAO
       private$.modes <- private$.cao_init(X, K)

       # Boucle d'optimisation
       for (iter in 1:self$get.max_iter()) {

         # -Assignation des variables
         clusters <- sapply(1:ncol(X), function(j) {
           which.min(sapply(1:K, function(m) {
             private$.col_dissimilarity(
               as.character(X[, j]),
               as.character(private$.modes[m, ])
             )
           }))
         })

         # Mise à jour des modes
         new_modes <- private$.update_modes(X, clusters, K)

         # Convergence
         if (all(new_modes == private$.modes)) {
           private$.clusters <- clusters
           message("Convergence atteinte en ", iter, " itérations.")
           return(invisible(self))
         }

         private$.modes <- new_modes
       }

       warning("Max_iter atteint sans convergence.")
       private$.clusters <- clusters
       invisible(self)
     },

     # Predict

     predict = function(newvars) {

       newvars <- as.data.frame(newvars)

       sapply(1:ncol(newvars), function(j) {
         which.min(sapply(1:self$get.n_cluster(), function(m) {
           private$.col_dissimilarity(
             as.character(newvars[, j]),
             as.character(private$.modes[m, ])
           )
         }))
       })
     },

     # Getters
     get.modes = function() private$.modes,
     get.clusters = function() private$.clusters
   )
)

