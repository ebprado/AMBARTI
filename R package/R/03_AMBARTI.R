#' @export
#'
run_AMBARTI <- function(data,
                        ntrees = 50,
                        nburn = 2000,
                        npost = 1000){

  # Some pre-processing
  y = data$y
  x = data$x
  x$g = as.factor(x$g)
  x$e = as.factor(x$e)

  # Run AMBARTI
  fit.ambarti = ambarti(x, y, ntrees = ntrees, nburn = nburn, npost = npost)

  return(fit.ambarti)
}
