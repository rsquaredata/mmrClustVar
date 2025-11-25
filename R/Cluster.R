



.Cluster <- R6::R6Class(
  "Cluster",
  
  private = list(
    
    # --- Attributs internes ---
    
    FNbGroupes   = NULL,  # K
    FX_active    = NULL,  # data.frame des variables actives (pré-traitées)
    FClusters    = NULL,  # vecteur d'affectation des variables actives
    FCenters     = NULL,  # centres / modes / prototypes selon la méthode
    FInertia     = NULL,  # inertie intra-cluster totale
    FConvergence = NULL,  # booléen
    FX_descr     = NULL,  # variables descriptives fournies à predict()
    FAlgorithme  = NULL,  # nom textuel de la variante exacte de l'algo
    
    # --- Helpers internes ---
    
    check_X = function(X, update_structure = TRUE) {
      if (!is.data.frame(X)) {
        stop("X doit être un data.frame.")
      }
      if (nrow(X)  < 2L) {
        stop("X doit contenir au moins 2 individus.")
      }
      
      # Si on met à jour la structure (cas fit()) -> au moins 2 variables actives
      if (update_structure && ncol(X) < 2L) {
        stop("X doit contenir au moins 2 variables actives.")
      }
      
      # Si on ne met pas à jour la structre (cas predict()) -> au moins 1 variable descriptive
      if (!update_structure && ncol(X) < 1L) {
        stop("X_new doit contenir au moins 1 variable descriptive.")
      }
    },
    
    check_X_new = function() {
      if (!is.data.frame(X_new)) {
        stop("X_new doit être un data.frame.")
      }
      if (nrow(X_new) != nrow(private$FX_active)) {
        stop("X_new doit avoir le même nombre d'individus (lignes) que les données actives.")
      }
    },
    
    scale_active_variables = function(X, num_idx) {
      if (length(num_idx) == 0L) return(X)
      
      for (j in num_idx) {
        col <- X[[j]]
        mu  <- mean(col, na.rm = TRUE)
        sdv <- stats::sd(col, na.rm = TRUE)
        if (is.na(sdv) || sdv == 0) {
          # Variable constante (on la laisse telle quelle mais on pourra l'exclure plus tard)
          X[[j]] <- col
        } else {
          X[[j]] <- (col - mu) / sdv
        }
      }
      return(X)
    }
    
    
  )
)