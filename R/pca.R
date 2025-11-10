library(R6)
library(corrplot)



PCA <- R6Class("PCA",

  public = list(

    eigen = NULL,
    calc = NULL,

    initialize = function(X) {

      df.ok <- is.data.frame(X)
      if(!df.ok) {
        stop(paste("Invalid `data` type. Expected data.frame, got", class(data)))
      }
      nb.ok <- sum(sapply(X,function(x){is.numeric(x)}))
      if(!nb.ok) {
        stop("Non numeric `data` values")
      }

      X <- as.data.frame(X)
      self$calc <- princomp(X, cor=TRUE, scores=TRUE)
      self$eigen <- self$calc$sdev^2
    },

    print = function(X) {
      cat("Variables : ", colnames(X), "\n")
      cat("Valeurs propres : ", self$eigen, "\n")
    },

    summary = function(ncomp=2) {

      if (is.null(ncomp) || ncomp <= 0 || ncomp > length(self$eigen))
      {
        ncomp = min(2, length(self$eigen))
      }

      cat("Valeurs propres : ", self$eigen, "\n")
      cat("Correlations","\n")
      temp <- as.matrix(self$calc$loadings[,1:ncomp])
      correl <- sapply(1:ncomp, function(j){ temp[,j]*sqrt(self$eigen[j]) })
      colnames(correl) <- paste("Comp.", 1:ncomp, sep="")
      print(correl, digits=3)
    },

    screeplot = function(thresholds=FALSE) {
      plot(1:length(self$eigen), self$eigen, main="Eigen values", type="b")
      if(thresholds) {
        bb <- sort(cumsum(1/(length(self$eigen):1)), decreasing=TRUE)
        lines(1:length(self$eigen), bb, type="b", col="red")
      }
    },

    correl_circle = function(X, comp1=1, comp2=2) {

      c1 <- self$calc$loadings[,comp1] * self$calc$sdev[comp1]
      c2 <- self$calc$loadings[,comp2] * self$calc$sdev[comp2]

      plot(c1, c2, xlim=c(-1,+1), ylim=c(-1,+1), type="n", asp=1)
      abline(h=0,v=0)
      text(ifelse(c1>0, c1+0.1, c1-0.1), c2, labels=colnames(X), cex=1, col="blue", font=2)
      symbols(0, 0, circles=1, inches=F, add=TRUE)

      arrows(0, 0, x1=c1, y1=c2, angle=10, col="gray")
    },

    correl_heat = function(X, used_comp=1) {

      if (is.null(used_comp) || used_comp < 1 || used_comp > length(self$eigen)){
        corrplot(cor(X))
      } else
      {
        ord <- order(self$calc$loadings[,used_comp])
        corrplot(cor(X[ord]))
      }
    }

  )

)
