library(R6)


#' @export
TestModule <- R6Class("TestModule",

  public = list(

    value = 3,

    my_func = function() {
      return("Hello world")
    }
  )

)
