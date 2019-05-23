predict.PoissonSpline<- 
  function (object, x, derivative = 0, ...) 
  {
    if (missing(x)) {
      x <- object$x
    }
    c(splint(object$predicted$x, object$predicted$y, 
             x, derivative = derivative, ...))
  }