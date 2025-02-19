
#' Gibbs free energy for the transition from vapor to water
#' @param r radius of the spherical water embryo (in m)
#' @param ta air temperature (in K)
#' @param Sw Supersaturation ratio with respect to water
#' @return Gibbs free energy (in J)
#' @export
gibbs_vapor_wat <- function(r, ta, Sw) {

  kB <- cldphys::const()$kB
  sigma_wat_vap <- cldphys::surf_tens_wat_vap(ta)

  G <- - (4 * pi * r^3) / (3 * v_wat(ta)) * kB * ta * log(Sw) + 4 * pi * r^2 * sigma_wat_vap

  return(G)
}

#' Critical water embryo radius for the transition from vapor to liquid water
#' @param ta air temperature (in K)
#' @param Sw Supersaturation ratio with respect to water
#' @export
rc_vapor_wat <- function(ta, Sw) {
  kB <- cldphys::const()$kB
  sigma_wat_vap <- cldphys::surf_tens_wat_vap(ta)

  rc <- 2 * v_wat(ta) * sigma_wat_vap / (kB * ta * log(Sw))

  return(rc)
}


#' Volume of water molecules in liquid water
#' @param ta air temperature (in K)
#' @return Volume of water molecules in water (in m^3)
#' @export
v_wat <- function(ta) {

  Na <- cldphys::const()$Na
  Mw <- cldphys::const()$Mw
  rho_wat <- cldphys::density_wat(ta)

  v_wat <- (Mw / Na) / rho_wat

  return(v_wat)
}



#' Number of water molecules forming the liquid water embryo
#' @param r radius of the spherical liquid water embryo (in m)
#' @param ta air temperature (in K)
#' @return Number of water molecules forming the water embryo
#' @export
n_wat <- function(r, ta) {
  nwat <- 4/3 * pi * r^3 / v_wat(ta)

  return(nwat)
}

#' Nuleation rate for water from homogeneous nucleation
#' Paramterization by Miller et al (1983), doi:10.1063/1.445236
#' From laboratory measurements, valid over the temperature range 230-290 K
#' @param S Supersaturation ratio with respect to water
#' @param ta air temperature (in K)
#' @return Nucleation rate (in m^-3 s^-1)
#' @export
nucl_rate_wat_hom_miller <- function(S, ta) {

  A <- 328.124 - 5.58243 * ta + 0.030365 * ta^2 - 5.0319*1E-5 * ta^3
  B <- 999.814 - 4.10087 * ta + 3.01084*1E-3 * ta^2

  ## J, in cm-3 s-1
  J <- S^2 * exp(A - B / log(S)^2 )

  ## Convert to m-3 s-1
  J <- J * 1E6

  return(J)

}

#' Nuleation rate for water from homogeneous nucleation
#' Theory by Pruppacher and Klett (1997), chapter 7
#' @param S Supersaturation ratio with respect to water
#' @param ta air temperature (in K)
#' @return Nucleation rate (in m^-3 s^-1)
#' @export
nucl_rate_wat_hom <- function(S, ta) {

  kB <- cldphys::const()$kB
  Na <- cldphys::const()$Na
  Mw <- cldphys::const()$Mw

  ## Mass of a water molecule
  m_w <- Mw / Na

  ## Surface tension between liquid water and water vapor
  sigma_wv <- cldphys::surf_tens_wat_vap(ta)

  ## Saturation vapor pressure of water
  es <- cldphys::es_w_koop(ta)
  e <- S * es

  ## Critical radius of liquid embryo at metastable equilibrium
  r_c <- cldphys::rc_vapor_wat(ta, S)

  ## Surface of the liquid embryo at critical size
  omega_c <- 4 * pi * r_c^2

  ## Number of water molecules forming the liquid embryo
  n_c <- cldphys::n_wat(r_c, ta)

  ## Gibbs free energy for the transition from vapor to water at critical size
  dG_c <- cldphys::gibbs_vapor_wat(r_c, ta, S)

  ## Zeldovich factor
  Z <- (dG_c / (3 * pi * kB * ta * n_c^2))^0.5

  ## Molecular flux of water molecules to surface of the embryo (m^-2 s^-1)
  w <- e / (2 * pi * m_w * kB * ta)^0.5

  ## Number concentration of embryos at saturation with respect to a flat water surface
  n_c_sat <- es / (kB * ta)

  ## Number concentration of embryos of critical size
  n_c <- n_c_sat * exp(- dG_c / (kB * ta))

  ## Nucleation rate (m^-3 s^-1)
  J <-  n_c * omega_c * w * Z

  return(J)
}
