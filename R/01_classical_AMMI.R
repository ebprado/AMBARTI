#' @export
#' @importFrom agricolae 'AMMI'
#'
run_classical_AMMI <- function(data){

  x <- data$x
  y <- data$y
  sim_data <- cbind(x, y)

  REP = 1 # no repetitions

  model <- with(sim_data, AMMI(ENV = env,
                               GEN = gen,
                               REP = REP,
                               Y =  y,
                               console = FALSE))
  return(model)
}
