#' Cloud droplet number concentration (from Grosvenor et al, 2015).
#' Default values for many constants are taken from Brenguier et al (2000), leading to consistency with Eq. 1 from Quaas et al (2006).
#' @param cod numeric, cloud optical depth
#' @param re numeric, cloud effective radius (in m)
#' @param k numeric, droplet size distribution shape parameter. Default is 0.8.
#' @param fad numeric, adiabaticity factor. Default is 1.
#' @param cw numeric, condensation rate (in kg m^-3 m^-3). Default is 1.9E-6.
#' @param Qext numeric, extinction efficiency factor. Default is 2.
#' @return numeric, cloud droplet number concentration (in # m^-3)
#' @export
cdnc <- function(cod, re, k = 0.8, fad = 1, cw = 1.9E-6, Qext = 2) {

  rho_w <- cldphys::const()$rho_w

  K = sqrt(5) / (2 * pi * k) * sqrt(fad * cw / (Qext * rho_w))

  cdnc <- K * cod^0.5 * re^{-2.5}

  return(cdnc)

}
