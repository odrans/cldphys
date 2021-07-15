#' Adiabatic condensation rate at satutation, following Albrecht et al (1990)
#' @param T numeric, air temperature (in K)
#' @param p numeric, atmospheric pressure (in Pa)
#' @param use_elr boolean, true if the environmental lapse rate (-6.5 K / km) should be used. If false, a moist adiabatic lapse rate will be used.
#' @return numeric, adiabatic condentation rate in kg m-3 m-1 (kg )
#' @export
cw_adiab_albrecht <- function(T, p, use_elr = FALSE) {

  l.const <- cldphys::const()
  eps <- l.const$eps
  Rd <- l.const$Rd
  g <- l.const$g
  elr <- l.const$gamma_env

  Lv <- cldphys::L_vap(T)
  es <- cldphys::es_w_koop(T)
  ws <- cldphys::mixing_ratio(es, p)
  qs <- cldphys::spec_hum(es, p)
  Tv <- cldphys::T_virtual(T, qs)

  H <- Rd * T / g
  rho_air <- p / (Rd * Tv)

  gamma_m <- cldphys::lapse_rate_moist(T, ws)
  if(use_elr) gamma_m <- elr

  mr_rate <- ((eps + ws) * ws * Lv) / (Rd * T^2) * gamma_m - ws * p / ((p - es) * H)  ## rate of change of the saturation mixing ratio

  cw <- mr_rate * rho_air  ## Convertion to a mass condensation rate

  return(cw)

}


#' Adiabatic condensation rate, following Ahmad et al (2013)
#' @param T numeric, air temperature (in K)
#' @param p numeric, atmospheric pressure (in Pa)
#' @param use_elr boolean, true if the environmental lapse rate (-6.5 K / km) should be used. If false, a moist adiabatic lapse rate will be used.
#' @return numeric, adiabatic condentation rate in kg m-3 m-1 (kg )
#' @export
cw_adiab_ahmad <- function(T, p, use_elr = FALSE) {

  l.const <- cldphys::const()
  eps <- l.const$eps
  Rd <- l.const$Rd
  g <- l.const$g
  elr <- l.const$gamma_env
  cp  <- l.const$cp

  Lv <- cldphys::L_vap(T)
  es <- cldphys::es_w_koop(T)
  ws <- cldphys::mixing_ratio(es, p)
  qs <- cldphys::spec_hum(es, p)
  Tv <- cldphys::T_virtual(T, qs)

  rho_air <- p / (Rd * Tv)

  gamma_d <- cldphys::lapse_rate_dry()
  gamma_m <- cldphys::lapse_rate_moist(T, ws)
  if(use_elr) gamma_m <- elr

  mr_rate <- cp / Lv * (gamma_d - gamma_m)

  cw <- mr_rate * rho_air

  return(cw)

}


#' Adiabatic condensation rate, following Bennartz et al (2009)
#' @param T numeric, air temperature (in K)
#' @param p numeric, atmospheric pressure (in Pa)
#' @param use_elr boolean, true if the environmental lapse rate (-6.5 K / km) should be used. If false, a moist adiabatic lapse rate will be used.
#' @return numeric, adiabatic condentation rate in kg m-3 m-1 (kg )
#' @export
cw_adiab_bennartz <- function(T, p, use_elr = FALSE) {

  l.const <- cldphys::const()
  Rv <- l.const$Rv
  elr <- l.const$gamma_env

  Lv <- cldphys::L_vap(T)
  es <- cldphys::es_w_koop(T)
  ws <- cldphys::mixing_ratio(es, p)

  gamma_m <- cldphys::lapse_rate_moist(T, ws)
  if(use_elr) gamma_m <- elr

  cw <- 1 / (Rv * T) * Lv * es / (Rv * T^2) * gamma_m

  return(cw)

}
