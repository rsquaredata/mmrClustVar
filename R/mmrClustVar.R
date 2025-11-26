#' Clustering de variables (façade R6)
#'
#' Classe R6 "façade" qui encapsule les différentes implémentations
#' de clustering de variables (k-means, k-modes, k-prototypes, k-medoids)
#' et fournit une interface unique pour l'utilisateur.
#'
#' L'objet se crée avec \code{mmrClustVar$new(...)} puis expose
#' les méthodes principales \code{$fit()}, \code{$predict()},
#' \code{$print()}, \code{$summary()} et \code{$plot()}.
#'
#' @docType class
#' @name mmrClustVar
#' @export
NULL

mmrClustVar <- R6::R6Class(
    "mmrClustVar",
    
    public = list(
        
        #' @description
        #' Constructeur de la façade.
        #'
        #' @param method Méthode de clustering de variables :
        #'   "kmeans", "kmodes", "kprototypes", "kmedoids" ou "auto".
        #' @param K Nombre de clusters (entier >= 2).
        #' @param scale Si TRUE, standardise les variables numériques.
        #' @param lambda Pondération pour les distances catégorielles (k-prototypes, k-medoids mixtes).
        initialize = function(
        method = c("kmeans", "kmodes", "kprototypes", "kmedoids", "auto"),
        K,
        scale  = TRUE,
        lambda = 1,
        ...
        ) {
            method            <- match.arg(method)
            private$method    <- method
            private$K         <- K
            private$scale     <- scale
            private$lambda    <- lambda
            private$engine    <- NULL
            private$X_last    <- NULL
            private$inertia_grid <- NULL
        },
        
        #' @description
        #' Ajuste le modèle sur un jeu de variables actives.
        #'
        #' @param X data.frame de variables (colonnes = variables).
        fit = function(X) {
            if (!is.data.frame(X)) {
                X <- as.data.frame(X)
            }
            
            if (is.null(private$engine) || identical(private$method, "auto")) {
                private$engine <- private$build_engine(X)
            }
            
            private$engine$fit(X)
            
            # mémoriser X pour la courbe d'inertie
            private$X_last       <- X
            private$inertia_grid <- NULL
            
            invisible(self)
        },
        
        #' @description
        #' Affecte des variables supplémentaires aux clusters existants.
        #'
        #' @param X_new data.frame de nouvelles variables (colonnes).
        predict = function(X_new) {
            if (is.null(private$engine)) {
                stop("[mmrClustVar] predict() appelé avant fit().")
            }
            private$engine$predict(X_new)
        },
        
        #' @description
        #' Résumé court de l'objet + moteur interne.
        print = function(...) {
            cat("Classe 'mmrClustVar'\n")
            cat("  Méthode demandée :", private$method, "\n")
            cat("  K                :", private$K, "\n")
            cat("  scale            :", private$scale, "\n")
            cat("  lambda           :", private$lambda, "\n\n")
            
            if (!is.null(private$engine)) {
                cat("Résumé du moteur interne :\n")
                private$engine$print()
            } else {
                cat("(Aucun modèle encore ajusté, appelez $fit(X)).\n")
            }
            invisible(self)
        },
        
        #' @description
        #' Résumé détaillé : tailles des clusters, inertie et indicateurs.
        summary = function(...) {
            if (is.null(private$engine)) {
                cat("[mmrClustVar] summary() : aucun modèle (fit() non appelé).\n")
                return(invisible(NULL))
            }
            private$engine$summary(...)
        },
        
        #' @description
        #' Graphiques principaux (inertie, clusters, adhésions, profils).
        #'
        #' @param type "inertia", "clusters", "membership" ou "profiles".
        plot = function(type = c("inertia", "clusters", "membership", "profiles"),
                        ...) {
            type <- match.arg(type)
            
            if (is.null(private$engine)) {
                stop("[mmrClustVar] plot() : aucun modèle (fit() non appelé).")
            }
            
            # si on a déjà calculé une grille d'inertie et qu'on demande "inertia"
            if (identical(type, "inertia") && !is.null(private$inertia_grid)) {
                df <- private$inertia_grid
                graphics::plot(
                    df$K, df$inertia,
                    type = "b",
                    xlab = "K",
                    ylab = "Inertie intra-cluster",
                    main = sprintf("Courbe d'inertie (%s)", private$method)
                )
                return(invisible(NULL))
            }
            
            # sinon on délègue à la classe fille
            private$engine$plot(type = type, ...)
        },
        
        #' @description
        #' Calcule l'inertie intra pour plusieurs valeurs de K.
        #'
        #' @param K_seq vecteur d'entiers (valeurs de K à tester).
        #' @param X data.frame (optionnel). Si NULL, réutilise le X de fit().
        #' @return data.frame(K, inertia)
        compute_inertia_path = function(K_seq, X = NULL) {
            # récupère X si non fourni
            if (is.null(X)) {
                X <- private$X_last
            }
            if (is.null(X)) {
                stop("[mmrClustVar] compute_inertia_path() : aucun X disponible. ",
                     "Fournis X ou appelle d'abord $fit(X).")
            }
            if (!is.data.frame(X)) {
                X <- as.data.frame(X)
            }
            
            p <- ncol(X)
            
            # K entre 2 et p
            K_seq <- sort(unique(as.integer(K_seq)))
            K_seq <- K_seq[!is.na(K_seq) & K_seq >= 2L & K_seq <= p]
            if (length(K_seq) == 0L) {
                stop("[mmrClustVar] compute_inertia_path() : K_seq doit contenir ",
                     "au moins un entier entre 2 et ", p, ".")
            }
            
            old_K <- private$K
            
            res <- data.frame(
                K       = K_seq,
                inertia = NA_real_
            )
            
            for (i in seq_along(K_seq)) {
                K_i <- K_seq[i]
                private$K <- K_i
                engine_i  <- private$build_engine(X)
                engine_i$fit(X)
                # get_inertia() est défini dans mmrClustVarBase
                res$inertia[i] <- engine_i$get_inertia()
            }
            
            private$K          <- old_K
            private$inertia_grid <- res
            
            res
        },
        
        #' @description
        #' Affectations des variables aux clusters.
        get_clusters = function() {
            if (is.null(private$engine)) return(NULL)
            private$engine$get_clusters()
        },
        
        #' @description
        #' Prototypes (centres/modes/medoids) des clusters.
        get_centers = function() {
            if (is.null(private$engine)) return(NULL)
            private$engine$get_centers()
        },
        
        #' @description
        #' Inertie intra du modèle courant.
        get_inertia = function() {
            if (is.null(private$engine)) return(NA_real_)
            private$engine$get_inertia()
        },
        
        #' @description
        #' Statut de convergence du modèle courant.
        get_convergence = function() {
            if (is.null(private$engine)) return(NA)
            private$engine$get_convergence()
        },
        
        #' @description
        #' Méthode demandée à la création.
        get_method = function() private$method,
        
        #' @description
        #' Nombre de clusters K demandé.
        get_K = function() private$K
    ),
    
    private = list(
        
        method       = NULL,
        K            = NULL,
        scale        = NULL,
        lambda       = NULL,
        
        engine       = NULL,   # instance de mmrClustVar*
        X_last       = NULL,   # dernier X passé à fit()
        inertia_grid = NULL,   # data.frame(K, inertia)
        
        # fabrique le moteur interne en fonction de method / types
        build_engine = function(X) {
            method <- private$method
            K      <- private$K
            scale  <- private$scale
            lambda <- private$lambda
            
            if (identical(method, "auto")) {
                is_num <- vapply(X, is.numeric, logical(1L))
                is_cat <- (!is_num) & vapply(
                    X,
                    function(col) is.factor(col) || is.character(col),
                    logical(1L)
                )
                
                if (all(is_num)) {
                    method_effective <- "kmeans"
                } else if (all(is_cat)) {
                    method_effective <- "kmodes"
                } else {
                    method_effective <- "kprototypes"
                }
                
                cat("[mmrClustVar] method = 'auto' -> méthode choisie :",
                    method_effective, "\n")
            } else {
                method_effective <- method
            }
            
            engine <- switch(
                method_effective,
                "kmeans" = mmrClustVarKMeans$new(
                    K      = K,
                    scale  = scale,
                    lambda = lambda
                ),
                "kmodes" = mmrClustVarKModes$new(
                    K      = K,
                    scale  = scale,
                    lambda = lambda
                ),
                "kprototypes" = mmrClustVarKPrototypes$new(
                    K      = K,
                    scale  = scale,
                    lambda = lambda
                ),
                "kmedoids" = mmrClustVarKMedoids$new(
                    K      = K,
                    scale  = scale,
                    lambda = lambda
                ),
                stop("[mmrClustVar] Méthode non supportée : ", method_effective)
            )
            
            engine
        }
    )
)