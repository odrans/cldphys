#' Compute the vaporisation latent heat
#' @param ta numeric, air temperature (in K)
#' @return numeric, latent heat of vaporisation (in J / kg)
#' @export
L_vap <- function(ta) {

  Mw <- cldphys::const()$Mw

  Lv <- ifelse(
    ta < 273.15,
    (56579 - 42.212 * ta + exp(0.1149 * (281.6 - ta))) * 1 / Mw, # From Murphy and Koop (2005)
    2.5E6 - 2.37E3 * (ta - 273.15) # From Hess (1959)
  )

  return(Lv)
}

#' Compute the sublimation latent heat
#' Parameters are from Murphy and Koop (2005)
#' @param ta numeric, air temperature (in K)
#' @return numeric, latent heat of sublimation (in J / kg)
#' @export
L_sub <- function(ta) {

  Mw <- cldphys::const()$Mw

  Ls <- (46782.5 + 35.8925 * ta - 0.07414 * ta^2 + 541.5 * exp(- (ta / 123.75)^2)) * 1 / Mw

  return(Ls)
}
