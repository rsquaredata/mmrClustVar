

# Base class for K-means, K-modes and K-prototypes variable clustering classes




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
    FMethod      = NULL,  # "k-means", "k-modes", "k-prototypes"
    
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
    
    check_X_new = function(X_new) {
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
    },
    
    # --- Abstract interns  ---
    
    compute_membership = function(X, clusters, centers) {
      # Compute cluster membership
      # Arguments:
      #   X:        data matrix
      #   clusters: vecteur d'affectation des variables actives
      #   centers:  centres / modes / prototypes selon la méthode
      # Returns:
      #   list(
      #     content: membership numerics
      #     label:   distance label
      #   )
      stop(paste("Membership method hasn't been implemented on ", private$FMethod, " and is required for summary."))
    }
  ),
  
  public = list(
    
    # --- Interpretability methods ---
    
    print = function(...) {
      cat("Variable clustering model\n")
      cat("  Algorithm      :", private$FMethod, "\n")
      cat("  K              :", private$FNbGroupes, "\n")
      if (!is.null(private$FX_active)) {
        cat("  Nb of variables  :", ncol(private$FX_active), "\n")
        cat("  Nb of obs     :", nrow(private$FX_active), "\n")
      }
      if (!is.null(private$FInertia)) {
        cat("  Intra-inertia  :", format(private$FInertia, digits = 4), "\n")
      }
      cat("  Convergence    :", private$FConvergence, "\n")
      invisible(self)
    },
    
    summary = function(...) {
      X        <- private$FX_active
      K        <- private$FNbGroupes
      clusters <- private$FClusters
      centers  <- private$FCenters
      
      cat("Model summary\n")
      cat("  Algorithm        :", method, "\n")
      cat("  K                :", K, "\n")
      cat("  Nb of variables  :", length(clusters), "\n")
      cat("  Intra-inertia    :", format(private$FInertia, digits = 4), "\n")
      cat("  Convergence      :", private$FConvergence, "\n\n")
      
      # --- 1) Degré d'adhésion par variable ---
      membership <- private$compute_membership(X, clusters, centers)
      
      # --- 2) Résumé par cluster ---
      tab_size <- table(clusters)
      cl_ids   <- as.integer(names(tab_size))
      mean_mem <- tapply(membership$content, clusters, mean, na.rm = TRUE)
      
      cat("Per cluster summary:\n")
      cl_df <- data.frame(
        cluster          = cl_ids,
        taille           = as.integer(tab_size[as.character(cl_ids)]),
        adhesion_moyenne = as.numeric(mean_mem[as.character(cl_ids)])
      )
      print(cl_df, row.names = FALSE)
      cat("\n")
      
      # --- 3) Tableau variables / cluster / adhésion ---
      var_df <- data.frame(
        variable  = colnames(X),
        cluster   = clusters,
        adhesion  = membership$content,
        stringsAsFactors = FALSE
      )
      
      cat("Adherence indicator (", membership$label, ")\n", sep = "")
      # On affiche seulement les 10 premières lignes triées par cluster puis par adhésion décroissante
      ord <- order(var_df$cluster, -var_df$adhesion)
      print(utils::head(var_df[ord, ], n = min(10L, nrow(var_df))), row.names = FALSE)
      cat("\n(The full data.frame is returned as an invisible object)\n")
      
      invisible(var_df)
    },
    
    plot = function(type = "clusters", Ks = NULL, ...) {
      .types <- c("inertia", "clusters", "membership")
      if (!(type %in% .types)) {
        stop("")
      }
      
      if (is.null(private$FX_active) || is.null(private$FClusters)) {
        stop("Aucun modèle appris. Appelez fit() d'abord.")
      }
      
      X        <- private$FX_active
      method   <- private$FMethod
      clusters <- private$FClusters
      centers  <- private$FCenters
      
      if (type == "clusters") {
        # Barplot des tailles de clusters
        tab <- table(clusters)
        graphics::barplot(
          tab,
          xlab = "Cluster",
          ylab = "Nombre de variables",
          main = "Taille des clusters de variables",
          col = names(tab)
        )
        return(invisible(NULL))
      }
      
      if (type == "membership") {
        # On réutilise la logique de summary() pour calculer l'adhésion
        membership <- private$compute_membership(X, clusters, centers)
        
        ord <- order(clusters)
        graphics::barplot(
          membership$content[ord],
          names.arg = colnames(X)[ord],
          las = 2,
          cex.names = 0.6,
          xlab = "Variables",
          ylab = membership$label,
          main = "Degré d'adhésion des variables à leur cluster",
          col = clusters,
          ...
        )
        return(ord)
        return(invisible(NULL))
      }
      
      if (type == "inertia") {
        # Courbe inertie(K) en relançant l'algo pour plusieurs valeurs de K
        p <- ncol(X)
        if (is.null(Ks)) {
          Ks <- seq_len(min(10L, p))
        }
        Ks <- Ks[Ks >= 2]
        if (length(Ks) == 0L) {
          stop("Impossible de tracer la courbe d'inertie : K doit être >= 2.")
        }
        
        old_K <- private$FNbGroupes
        inertias <- numeric(length(Ks))
        
        for (i in seq_along(Ks)) {
          private$FNbGroupes <- as.integer(Ks[i])
          res_i <- self$fit(X)
          inertias[i] <- res_i$inertia
        }
        
        private$FNbGroupes <- old_K  # on restaure
        
        graphics::plot(
          Ks, inertias, type = "b",
          xlab = "K (nombre de clusters)",
          ylab = "Inertie intra-cluster",
          main = "Courbe de l'inertie en fonction de K",
          ...
        )
        return(invisible(NULL))
      }
    },
    
    # --- Getter methods ---
    
    get_clusters = function() {
      private$FClusters
    },
    
    get_inertia = function() {
      private$FInertia
    },
    
    get_method = function() {
      private$FMethod
    },
    
    get_centers = function() {
      private$FCenters
    },
    
    get_convergence = function() {
      private$FConvergence
    },
    
    get_X_descr = function() {
      private$FX_descr
    },
    
    # --- Abstract methods ---
    
    fit = function(X) {
      stop(paste("$fit() method hasn't been implemented on ", private$FMethod))
    },
    predict = function(X_new) {
      stop(paste("$predict() method hasn't been implemented on ", private$FMethod))
    }
  )
)