#' Compute the vaporisation latent heat
#' @param T numeric, air temperature in K
#' @return numeric, latent heat of vaporisation in J / kg
#' @export
L_vap <- function(T) {

  if(length(T) !=1 ) stop("Input should be one-dimensional")

  Mw <- cldphys::const()$Mw

  if(T < 273.15) {
    # From Murphy and Koop (2005)
    Lv <- (56579 - 42.212 * T + exp(0.1149 * (281.6 - T))) * 1 / Mw
  } else {
    # From Hess (1959)
    Lv <- 2.5E6 - 2.37E3 * (T - 273.15)
  }

  return(Lv)

}
