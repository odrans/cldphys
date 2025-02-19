
#' Surface tension of between ice and water vapor
#' Parameterization from Hale and Plummer (1974)
#' @param ta Temperature (K)
#' @return Surface tension (N/m)
#' @export
surf_tens_ice_vap <- function(ta) {
  0.001 * (141 - 0.15 * ta)
}


#' Surface tension of between ice and liquid water
#' Parameterization from Reinhardt and Doye (2013)
#' @param ta Temperature (K)
#' @return Surface tension (N/m)
#' @export
surf_tens_ice_wat <- function(ta) {
  ##23.24E-3 * (ta / 235.8) ^0.35 ## Nemec (2013)
  24 * 1E-3 + 0.18 * 1E-3 * (ta - 240)
}

#' Surface tension of between liquid water and water vapor
#' Parameterization from Pruppacher and Klett (1997); Eq. 5.12
#' @param ta Temperature (K)
#' @return Surface tension (N/m)
#' @export
surf_tens_wat_vap <- function(ta) {

  tc <- ta - 273.15

  a0 <- 75.93
  a1 <- 0.115
  a2 <- 6.618 * 1E-2
  a3 <- 6.511 * 1E-3
  a4 <- 2.993 * 1E-4
  a5 <- 6.283 * 1E-6
  a6 <- 5.285 * 1E-8

  sigma_w_a <- (a0 + a1 * tc + a2 * tc^2 + a3 * tc^3 + a4 * tc^4 + a5 * tc^5 + a6 * tc^6) * 1E-3

  return(sigma_w_a)

}
