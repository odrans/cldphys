#' Virtual temperature
#' @param T numeric, air temperature (in K)
#' @param qv numeric, specific humidity (in kg kg^-1)
#' @return numeric, virtual temperature (in K)
#' @export
T_virtual <- function(T, qv) {

  Tv <- T * (1 + 0.608 * qv)

  return(Tv)

}
