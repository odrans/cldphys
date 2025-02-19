
#' Gibbs free energy for the transition from vapor to ice
#' @param r radius of the spherical ice embryo (in m)
#' @param ta air temperature (in K)
#' @param Si Supersaturation ratio with respect to ice
#' @return Gibbs free energy (in J)
#' @export
gibbs_vapor_ice <- function(r, ta, Si) {

  kB <- cldphys::const()$kB
  sigma_ice_vap <- cldphys::surf_tens_ice_vap(ta)

  G <- - (4 * pi * r^3) / (3 * v_ice(ta)) * kB * ta * log(Si) + 4 * pi * r^2 * sigma_ice_vap

  return(G)
}


#' Gibbs free energy for the transition from water to ice
#' @param r radius of the spherical ice embryo (in m)
#' @param ta air temperature (in K)
#' @return Gibbs free energy (in J)
#' @export
gibbs_wat_ice <- function(r, ta) {

  kB <- cldphys::const()$kB

  ## Surface tension between ice and water
  sigma_ice_wat <- cldphys::surf_tens_ice_wat(ta)

  ## Saturation vapor pressure of water and ice
  es_w <- cldphys::es_w_koop(ta)
  es_i <- cldphys::es_i_koop(ta)

  ## Gibbs free energy for the transition from water to ice
  G <- - (4 * pi * r^3) / (3 * v_ice(ta)) * kB * ta * log(es_w / es_i) + 4 * pi * r^2 * sigma_ice_wat

  return(G)
}

#' Critical ice embryo radius for the transition from vapor to ice
#' @param ta air temperature (in K)
#' @param Si Supersaturation ratio with respect to ice
#' @export
rc_vapor_ice <- function(ta, Si) {

  kB <- cldphys::const()$kB
  sigma_ice_vap <- cldphys::surf_tens_ice_vap(ta)

  rc <- 2 * v_ice(ta) * sigma_ice_vap / (kB * ta * log(Si))

  return(rc)
}

#' Critical ice embryo radius for the transition from water to ice
#' @param ta air temperature (in K)
#' @export
rc_wat_ice <- function(ta) {

  kB <- cldphys::const()$kB

  ## Surface tension between ice and water
  sigma_ice_wat <- cldphys::surf_tens_ice_wat(ta)

  ## Saturation vapor pressure of water and ice
  es_w <- cldphys::es_w_koop(ta)
  es_i <- cldphys::es_i_koop(ta)
  S <- es_w / es_i

  rc <- 2 * v_ice(ta) * sigma_ice_wat / (kB * ta * log(S))

  return(rc)
}

#' Volume of water molecules in ice
#' @param ta air temperature (in K)
#' @return Volume of water molecules in ice (in m^3)
#' @export
v_ice <- function(ta) {
  Na <- cldphys::const()$Na
  Mw <- cldphys::const()$Mw
  rho_ice <- cldphys::density_ice(ta)

  v_ice <- (Mw / Na) / rho_ice

  return(v_ice)
}

#' Number of water molecules forming the ice embryo
#' @param r radius of the spherical ice embryo (in m)
#' @param ta air temperature (in K)
#' @return Number of water molecules forming the ice embryo
#' @export
n_ice <- function(r, ta) {
  nice <- 4/3 * pi * r^3 / v_ice(ta)

  return(nice)
}


#' Nuleation rate for ice crystals from homogeneous nucleation
#' Theory by Pruppacher and Klett (1997), chapter 7
#' @param Si Supersaturation ratio with respect to ice
#' @param ta air temperature (in K)
#' @return Nucleation rate (in m^-3 s^-1)
#' @export
nucl_rate_ice_hom <- function(Si, ta) {

  kB <- cldphys::const()$kB
  Na <- cldphys::const()$Na
  Mw <- cldphys::const()$Mw

  ## Mass of a water molecule
  m_w <- Mw / Na

  ## Saturation vapor pressure of water
  es <- cldphys::es_i_koop(ta)
  e <- Si * es

  ## Critical radius of liquid embryo at metastable equilibrium
  r_c <- cldphys::rc_vapor_ice(ta, Si)

  ## Surface of the liquid embryo at critical size
  omega_c <- 4 * pi * r_c^2

  ## Number of water molecules forming the liquid embryo
  n_c <- cldphys::n_ice(r_c, ta)

  ## Gibbs free energy for the transition from vapor to water at critical size
  dG_c <- cldphys::gibbs_vapor_ice(r_c, ta, Si)

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


#' Nuleation rate for ice crystals from homogeneous nucleation
#' Theory by Pruppacher and Klett (1997), chapter 7
#' @param Si Supersaturation ratio with respect to ice
#' @param ta air temperature (in K)
#' @param fix_Cpre logical, if TRUE the prefactor is fixed to 10^41 regardless of the temperature
#' @return Nucleation rate (in m^-3 s^-1)
#' @export
nucl_rate_ice_hom_wat <- function(ta, fix_Cpre = FALSE) {

  kB <- cldphys::const()$kB
  Na <- cldphys::const()$Na
  Mw <- cldphys::const()$Mw
  h <- 6.62607015e-34

  ## Saturation vapor pressure of water
  ## es <- cldphys::es_i_koop(ta)
  ## e <- Si * es

  ## Critical radius of liquid embryo at metastable equilibrium
  r_c <- cldphys::rc_wat_ice(ta)

  ## Surface of the liquid embryo at critical size
  omega_c <- 4 * pi * r_c^2

  ## Number of water molecules forming the liquid embryo
  n_c <- cldphys::n_ice(r_c, ta)

  ## Gibbs free energy for the transition from vapor to water at critical size
  dG_c <- cldphys::gibbs_wat_ice(r_c, ta)

  ## Zeldovich factor
  Z <- (dG_c / (3 * pi * kB * ta * n_c^2))^0.5

  ## Molecular Gibbs free energy of activation for diffusion of water moleules across the water-ice boundary
  dg <- cldphys::gibbs_diff_water_ice(ta)

  ## Number density of water molecules in the surrounding water
  N1 <- cldphys::density_wat(ta) * Na / Mw

  ## Number of molecules in contat with the surface of the embryo (here from Fletcher)
  Nc <- 3 * 1E-10 * N1 ## Fletcher

  ## Molecular flux of water molecules to surface of the embryo (m^-2 s^-1)
  w <- Nc * kB * ta / h * exp(-dg / (kB * ta))

  ## Number concentration of embryos of critical size
  N_c <- N1 * exp(- dG_c / (kB * ta))

  ## Prefactor in Ickes et al. (2015), should be around 10^41 regardless of the temperature
  ## Cpre <- omega_c * Nc * Z *  N1 * kB * ta / h

  ## Nucleation rate (m^-3 s^-1)
  J <-  N_c * omega_c * w * Z

  if(fix_Cpre) {
    J <- 10^41 * exp(-dg / (kB * ta)) * exp(-dG_c / (kB * ta))
  }

  return(J)
}


#' Gibbs free energy for the transition from water to ice
#' @param ta air temperature (in K)
#' @return Gibbs free energy (in J)
#' @export
gibbs_diff_water_ice <- function(ta) {

  kB <- cldphys::const()$kB
  E <- 892
  ta_star <- 118

  dg <- kB * ta^2 * E / (ta - ta_star)^2

  return(dg)
}
