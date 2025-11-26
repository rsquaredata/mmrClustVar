#' Classe de base pour le clustering de variables
#'
#' Classe R6 abstraite qui factorise le comportement commun aux
#' différentes variantes de clustering de variables
#' (k-means, k-modes, k-prototypes, k-medoids).
#'
#' Elle n'est pas destinée à être instanciée directement par l'utilisateur
#' final, mais sert de classe mère pour les classes spécialisées.
#'
#' @docType class
#' @name mmrClustVarBase
#' @keywords internal
#' @noRd
#'
#' @section Méthodes :
#' \describe{
#'   \item{\code{$initialize(method_name, K, scale = TRUE, lambda = 1, ...)}}{
#'     Initialise les paramètres généraux du modèle.
#'   }
#'   \item{\code{$fit(X)}}{
#'     Gère le flux global d'apprentissage (préparation des données,
#'     appel à l'algorithme spécifique, stockage des résultats).
#'   }
#'   \item{\code{$predict(X_new)}}{
#'     Rattache des variables supplémentaires aux clusters appris.
#'   }
#'   \item{\code{$print()}}{
#'     Résumé succinct de l'objet.
#'   }
#'   \item{\code{$summary()}}{
#'     Résumé détaillé des résultats.
#'   }
#'   \item{\code{$plot(type)}}{
#'     Visualisations associées au modèle.
#'   }
#'   \item{\code{$get_clusters()}}{
#'     Renvoie l'affectation des variables aux clusters.
#'   }
#'   \item{\code{$get_centers()}}{
#'     Renvoie les prototypes des clusters.
#'   }
#'   \item{\code{$get_method()}}{
#'     Renvoie le nom de la méthode de clustering utilisée.
#'   }
#'   \item{\code{$get_K()}}{
#'     Renvoie le nombre de clusters.
#'   }
#'   \item{\code{$get_convergence()}}{
#'     Renvoie l'indicateur de convergence de l'algorithme.
#'   }
#' }
NULL

mmrClustVarBase <- R6::R6Class(
    "mmrClustVarBase",
    
    public = list(
        
        initialize = function(K,
                              scale = TRUE,
                              lambda = 1,
                              method_name = "base") {
            
            if (missing(K) || !is.numeric(K) || length(K) != 1L || K < 2) {
                stop("[mmrClustVarBase] K doit être un entier >= 2")
            }
            private$FNbGroupes   <- as.integer(K)
            private$FScale       <- isTRUE(scale)
            private$FLambda      <- as.numeric(lambda)
            private$FMethod      <- method_name
            private$FConvergence <- FALSE
            private$FInertia     <- NA_real_
            private$FAlgorithme  <- method_name
        },
        
        fit = function(X) {
            # 1) Vérification / typage / structure
            X <- private$check_and_prepare_X(X, update_structure = TRUE)
            
            # 2) Filtrage des lignes avec valeurs manquantes / non finies
            #    (important pour prcomp() dans k-means / partie numérique de k-prototypes)
            num_idx <- private$FNumCols
            ok <- stats::complete.cases(X)
            
            if (length(num_idx) > 0L) {
                X_num <- as.matrix(X[, num_idx, drop = FALSE])
                ok_num <- apply(X_num, 1L, function(row) {
                    all(is.finite(row) | is.na(row))
                })
                ok <- ok & ok_num
            }
            
            if (!all(ok)) {
                X <- X[ok, , drop = FALSE]
            }
            
            # 3) Standardisation éventuelle des quantitatives
            if (private$FScale) {
                X <- private$scale_active_variables(X)
            }
            
            # 4) Lancement de l'algo spécifique (classe fille)
            res <- private$run_clustering(X)
            
            # 5) Stockage des résultats communs
            private$FX_active    <- X
            private$FClusters    <- res$clusters
            private$FCenters     <- res$centers
            private$FInertia     <- res$inertia
            private$FConvergence <- isTRUE(res$converged)
            
            invisible(self)
        },
        
        predict = function(X_new) {
            if (is.null(private$FX_active)) {
                stop("[mmrClustVarBase] fit() doit être appelé avant predict()")
            }
            X_new <- private$check_and_prepare_X(X_new, update_structure = FALSE)
            
            res <- lapply(seq_along(X_new), function(j) {
                var_name <- colnames(X_new)[j]
                private$predict_one_variable(X_new[[j]], var_name)
            })
            
            out <- do.call(rbind, res)
            rownames(out) <- NULL
            private$FX_new <- X_new
            out
        },
        
        print = function(...) {
            cat("Classe", class(self)[1L], "\n")
            cat("  Méthode          :", private$FMethod, "\n")
            cat("  K (nb groupes)   :", private$FNbGroupes, "\n")
            cat("  Nb variables act :", length(private$FClusters), "\n")
            cat("  Convergence      :", private$FConvergence, "\n")
            cat("  Inertie intra    :", private$FInertia, "\n")
            invisible(self)
        },
        
        summary = function(...) {
            if (is.null(private$FClusters)) {
                cat("Modèle non encore ajusté (fit() non appelé).\n")
                return(invisible(NULL))
            }
            
            cat("=== Résumé global ===\n")
            cat("Méthode      :", private$FMethod, "\n")
            cat("K            :", private$FNbGroupes, "\n")
            cat("Convergence  :", private$FConvergence, "\n")
            cat("Inertie      :", private$FInertia, "\n\n")
            
            clusters   <- private$FClusters
            tab_taille <- table(clusters)
            cat("=== Taille des clusters ===\n")
            print(tab_taille)
            cat("\n")
            
            # hook pour les classes filles : métriques d'adhésion spécifiques
            private$summary_membership()
            
            invisible(NULL)
        },
        
        plot = function(type = c("inertia", "clusters", "membership", "profiles"), ...) {
            type <- match.arg(type)
            
            if (is.null(private$FClusters)) {
                stop("[mmrClustVarBase] fit() doit être appelé avant plot().")
            }
            
            if (type == "clusters") {
                barplot(
                    table(private$FClusters),
                    main = sprintf("Répartition des variables (%s)", private$FMethod),
                    xlab = "Cluster",
                    ylab = "Nombre de variables"
                )
            } else if (type == "inertia") {
                barplot(
                    private$FInertia,
                    main = sprintf("Inertie intra-cluster (%s)", private$FMethod),
                    ylab = "Inertie"
                )
            } else if (type == "membership") {
                private$plot_membership()
            } else if (type == "profiles") {
                private$plot_profiles()
            }
            
            invisible(NULL)
        },
        
        get_clusters    = function() private$FClusters,
        get_centers     = function() private$FCenters,
        get_method      = function() private$FMethod,
        get_K           = function() private$FNbGroupes,
        get_inertia     = function() private$FInertia,
        get_convergence = function() private$FConvergence
    ),
    
    private = list(
        
        FMethod      = NULL,
        FNbGroupes   = NULL,
        FScale       = NULL,
        FLambda      = NULL,
        FX_active    = NULL,
        FX_new       = NULL,
        FClusters    = NULL,
        FCenters     = NULL,
        FInertia     = NULL,
        FConvergence = NULL,
        FAlgorithme  = NULL,
        FNumCols     = NULL,
        FCatCols     = NULL,
        
        check_and_prepare_X = function(X, update_structure = TRUE) {
            if (!is.data.frame(X)) {
                X <- as.data.frame(X)
            }
            if (ncol(X) == 0L) {
                stop("[mmrClustVarBase] X doit contenir au moins une colonne.")
            }
            
            num_idx <- which(vapply(X, is.numeric, logical(1L)))
            cat_idx <- which(vapply(
                X,
                function(col) is.factor(col) || is.character(col),
                logical(1L)
            ))
            
            if (update_structure) {
                private$FNumCols <- num_idx
                private$FCatCols <- cat_idx
            }
            
            # conversion caractères -> facteur
            for (j in cat_idx) {
                if (!is.factor(X[[j]])) {
                    X[[j]] <- factor(X[[j]])
                }
            }
            
            X
        },
        
        scale_active_variables = function(X) {
            num_idx <- private$FNumCols
            if (length(num_idx) == 0L) return(X)
            
            X[num_idx] <- lapply(X[num_idx], function(col) {
                m <- mean(col, na.rm = TRUE)
                s <- stats::sd(col, na.rm = TRUE)
                if (is.na(s) || s == 0) return(rep(0, length(col)))
                (col - m) / s
            })
            
            X
        },
        
        run_clustering = function(X) {
            stop("[mmrClustVarBase] run_clustering() doit être redéfini dans une classe fille.")
        },
        
        predict_one_variable = function(x_new, var_name) {
            stop("[mmrClustVarBase] predict_one_variable() doit être redéfini dans une classe fille.")
        },
        
        summary_membership = function() {
            cat("(Pas d'indicateurs d'adhésion spécifiques définis pour cette classe.)\n")
        },
        
        plot_membership = function() {
            stop("[mmrClustVarBase] plot(type = 'membership') non défini pour cette classe.")
        },
        
        plot_profiles = function() {
            warning("[mmrClustVarBase] plot(type = 'profiles') non défini pour cette classe, aucun graphique produit.")
        }
    )
)
