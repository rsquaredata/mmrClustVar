
#' @export

Kmodes <- R6::R6Class(
  "K-modes",
  
  inherit = .Cluster,
  
  private = list(
    
    # --- Private attributs ---
    
    FCatCols     = NULL,  # indices des variables qualitatives dans FX_active
    
    # --- Helper methods ---
    
    prepare_X = function(X, update_structure = TRUE) {
      # Check X matrix
      private$check_X(X, update_structure)
      # Identification des colonnes numériques / catégorielles
      cat_idx <- which(vapply(X, function(col) is.factor(col) || is.character(col), logical(1L)))
      
      # Conversion des qualitatives en facteur
      if (length(cat_idx) > 0L) {
        for (j in cat_idx) {
          if (!is.factor(X[[j]])) {
            X[[j]] <- factor(X[[j]])
          }
        }
      }
      
      if (update_structure) {
        private$FCatCols <- cat_idx
        if (length(cat_idx) == 0L) {
          stop("No qualitative variables found for the “kmodes” method.")
        }
      }
      
      return(X)
    },
    
    # Rattachement d'une seule variable descriptive
    predict_one_variable = function(x_new, var_name) {
      # x_new : vecteur (numérique ou factor/character)
      # var_name : nom de la variable descriptive
      
      method <- private$FMethod
      n <- nrow(private$FX_active)
      K <- private$FNbGroupes
      
      if (length(x_new) != n) {
        stop("predict_one_variable() : variable '", var_name,
             "' doesn't have the same number of observations as the active data.")
      }
      
      # On accepte factor ou character
      if (!(is.factor(x_new) || is.character(x_new))) {
        stop("predict_one_variable() : for 'kmodes', variable '",
             var_name, "' must be qualitative (factor or character).")
      }
      x_char <- as.character(x_new)
      
      centers <- private$FCenters  # liste de K prototypes catégoriels (vecteurs de longueur n)
      if (is.null(centers)) {
        stop("predict_one_variable(): prototypes are not available. Call fit() first.")
      }
      
      distances <- numeric(K)
      adhesions <- numeric(K)
      
      for (k in seq_len(K)) {
        zk <- centers[[k]]
        if (is.null(zk) || length(zk) != n) {
          distances[k] <- NA_real_
          adhesions[k] <- NA_real_
        } else {
          mismatch <- x_char != zk
          d <- mean(mismatch, na.rm = TRUE) # simple matching dissimilarity
          if (is.na(d)) d <- 1
          distances[k] <- d
          adhesions[k] <- 1 - d             # proportion de matches
        }
      }
      
      best_k  <- which.min(distances)
      best_d  <- distances[best_k]
      best_ad <- adhesions[best_k]
      
      res <- data.frame(
        variable      = var_name,
        cluster       = best_k,
        adhesion      = best_ad,
        distance_best = best_d,
        stringsAsFactors = FALSE
      )
      
      for (k in seq_len(K)) {
        col_name <- paste0("d", k)
        res[[col_name]] <- distances[k]
      }
      
      return(res)
    },
    
    compute_membership = function(X, clusters, centers) {
      p <- length(clusters)
      membership <- rep(NA_real_, p)
      names(membership) <- names(clusters)
      
      # 1 - dissimilarité simple matching
      X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)
      for (j in seq_len(p)) {
        kj <- clusters[j]
        xj <- X_char[[j]]
        zk <- centers[[kj]]
        mismatch <- xj != zk
        d <- mean(mismatch, na.rm = TRUE)
        if (is.na(d)) d <- 1
        membership[j] <- 1 - d
      }
      membership_label <- "1 - dissimilarity (simple matching) with cluster mode"
      
      return(list(
        content = membership,
        label = membership_label
      ))
    },
    
    # --- Main algorithm method ---
    
    clusterize = function(X) {
      # X : data.frame n x p, uniquement des variables qualitatives
      K <- private$FNbGroupes
      n <- nrow(X)
      p <- ncol(X)
      
      # Sécurité : toutes les colonnes doivent être des facteurs
      is_cat <- vapply(X, is.factor, logical(1L))
      if (!all(is_cat)) {
        stop("clusterize() : X must contain only qualitative variables.")
      }
      
      # On travaille en caractère pour simplifier les comparaisons
      X_char <- as.data.frame(lapply(X, as.character), stringsAsFactors = FALSE)
      
      # Initialisation des clusters (partition des variables)
      clusters <- rep(seq_len(K), length.out = p)
      clusters <- sample(clusters)  # randomisation
      
      max_iter    <- 100L
      tol         <- 1e-6
      inertia_old <- Inf
      converged   <- FALSE
      Z_list      <- vector("list", K)  # prototypes catégoriels (vecteurs de longueur n)
      
      # Fonction locale pour calculer le prototype d'un cluster
      compute_prototype <- function(Xk_char) {
        # Xk_char : data.frame n x p_k (caractères)
        # Pour chaque individu i, on prend la modalité la plus fréquente sur les variables du cluster
        apply(Xk_char, 1L, function(row_i) {
          tab <- table(row_i, useNA = "no")
          if (length(tab) == 0L) {
            return(NA_character_)
          } else {
            # mode (on prend la première modalité en cas d'égalité)
            names(tab)[which.max(tab)]
          }
        })
      }
      
      for (iter in seq_len(max_iter)) {
        # 1) Calcul des prototypes Z_k pour chaque cluster
        for (k in seq_len(K)) {
          vars_k <- which(clusters == k)
          
          if (length(vars_k) == 0L) {
            # Cluster vide : on choisit une variable au hasard et on la force dans ce cluster
            cand <- sample(seq_len(p), 1L)
            clusters[cand] <- k
            vars_k <- cand
          }
          
          Xk_char <- X_char[, vars_k, drop = FALSE]
          Z_list[[k]] <- compute_prototype(Xk_char)  # vecteur de longueur n
        }
        
        # 2) Calcul des dissimilarités d(j,k) pour chaque variable / cluster
        dist_mat <- matrix(NA_real_, nrow = p, ncol = K)
        for (j in seq_len(p)) {
          xj <- X_char[[j]]  # vecteur de longueur n
          for (k in seq_len(K)) {
            zk <- Z_list[[k]]
            if (length(zk) != n) {
              stop("clusterize() : prototype of incompatible length.")
            }
            mismatch <- xj != zk
            d <- mean(mismatch, na.rm = TRUE)
            if (is.na(d)) d <- 1  # si tout est NA, on prend la dissimilarité maximale
            dist_mat[j, k] <- d
          }
        }
        
        # 3) Nouvelle affectation : chaque variable va au cluster de dissimilarité minimale
        new_clusters <- apply(dist_mat, 1L, which.min)
        
        # Inertie intra-cluster = somme des dissimilarités au meilleur cluster
        min_d <- dist_mat[cbind(seq_len(p), new_clusters)]
        inertia <- sum(min_d, na.rm = TRUE)
        
        # 4) Critère de convergence
        if (all(new_clusters == clusters) || abs(inertia_old - inertia) < tol) {
          converged   <- TRUE
          clusters    <- new_clusters
          inertia_old <- inertia
          break
        }
        
        clusters    <- new_clusters
        inertia_old <- inertia
      }
      
      names(clusters) <- names(X) # set column names
      
      res <- list(
        clusters  = clusters,   # vecteur de longueur p : cluster de chaque variable
        centers   = Z_list,     # liste de K prototypes (vecteurs de longueur n)
        inertia   = inertia_old,
        converged = converged
      )
      
      return(res)
    }
  ),
  
  public = list(
    
    initialize = function(K = 2L, random_state = NULL) {
      if (K < 2L) stop("K must be greater or egal to 2")
      
      private$FMethod <- "kmodes"
      private$FNbGroupes <- as.integer(K)
      set.seed(random_state)
    },
    
    fit = function(X) {
      X <- private$prepare_X(X, update_structure = TRUE)
      private$FX_active <- X
      
      res <- private$clusterize(X)
      
      private$FClusters    <- res$clusters
      private$FCenters     <- res$centers
      private$FInertia     <- res$inertia
      private$FConvergence <- isTRUE(res$converged)
      
      invisible(res)
    },
    
    predict = function(X_new) {
      if (is.null(private$FX_active)) {
        stop("The model has not yet been trained. Call fit() first.")
      }
      private$check_X_new(X_new)
      
      res_list <- vector("list", length = ncol(X_new))
      for (j in seq_len(ncol(X_new))) {
        res_list[[j]] <- private$predict_one_variable(
          x_new    = X_new[[j]],
          var_name = colnames(X_new)[j]
        )
      }
      
      private$FX_descr <- X_new # mémorisation des variables descriptives
      
      res <- do.call(rbind, res_list)
      rownames(res) <- NULL
      return(res)
    }
  )
)