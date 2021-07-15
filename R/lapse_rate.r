#' Moist adiabatic lapse rate (from Hemond and Fechner, 2015)
#' @param T numeric, air temperature (in K)
#' @param w numeric, mixing ratio (in kg kg^-1)
#' @return numeric, moist adiabatic lapse rate in K m^-1
#' @export
lapse_rate_moist <- function(T, w) {

  l.const <- cldphys::const()
  g <- l.const$g
  Rd <- l.const$Rd
  Rv <- l.const$Rv
  eps <- l.const$eps
  cp <- l.const$cp

  Lv <- L_vap(T)

  gamma_w <- g * (1 + (Lv * w) / (Rd * T)) / (cp + (Lv^2 * w) / (Rv * T^2))

  return(gamma_w)

}


#' Dry adiabatic lapse rate
#' @return numeric, dry adiabatic lapse rate in K m^-1
#' @export
lapse_rate_dry <- function() {

  l.const <- cldphys::const()
  g <- l.const$g
  cp <- l.const$cp

  gamma_d <- g / cp

  return(gamma_d)

}
