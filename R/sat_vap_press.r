#' Saturation vapour pressure with respect to water, using the Magnus formula
#' @param T numeric, air temperature (in K)
#' @return numeric, saturation vapour pressure (in Pa)
#' @export
es_w_magnus <- function(T) {

  T <- T - 273.15

  esat <- 6.1094 * exp(17.625 * T / (T + 243.04)) * 1E2

  return(esat)

}


#' Saturation vapour pressure with respect to water
#' Parameterization from Murphy and Koop (2005), Eq. 10
#' Considered valid for 123 < T < 332 K
#' @param ta numeric, air temperature (zin K)
#' @return numeric, saturation vapour pressure (in Pa)
#' @export
es_w_koop <- function(ta) {

  esat <- exp(54.842763 - 6763.22 / ta - 4.210 * log(ta) + 0.000367 * ta + tanh(0.0415 * (ta - 218.8)) *
              (53.878 - 1331.22 / ta - 9.44523 * log(ta) + 0.014025 * ta))

  return(esat)

}

#' Saturation vapour pressure with respect to ice
#' Parameterization from Murphy and Koop (2005), Eq. 7
#' Considered valid for T > 110 K
#' @param ta numeric, air temperature (in K)
#' @return numeric, saturation vapour pressure (in Pa)
#' @export
es_i_koop <- function(ta) {

  b0 <- 9.550426
  b1 <- -5723.265
  b2 <- 3.53068
  b3 <- -0.00728332

  esat <- exp(b0 + b1/ta + b2*log(ta) + b3*ta)

  return(esat)

}
