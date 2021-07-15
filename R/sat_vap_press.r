#' Saturation vapour pressure with respect to water, using the Magnus formula
#' @param T numeric, air temperature (in K)
#' @return numeric, saturation vapour pressure (in Pa)
#' @export
es_w_magnus <- function(T) {

  T <- T - 273.15

  esat <- 6.1094 * exp(17.625 * T / (T + 243.04)) * 1E2

  return(esat)

}


#' Saturation vapour pressure with respect to water, using the empirical formula by Murphy and Koop (2005)
#' @param T numeric, air temperature (in K)
#' @return numeric, saturation vapour pressure (in Pa)
#' @export
es_w_koop <- function(T) {

  esat <- exp(54.842763 - 6763.22 / T - 4.210 * log(T) + 0.000367 * T + tanh(0.0415 * (T - 218.8)) * (53.878 - 1331.22 / T - 9.44523 * log(T) + 0.014025 * T))

  return(esat)

}
