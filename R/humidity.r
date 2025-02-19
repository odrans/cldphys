#' Mixing ratio of the mass of water vapour to the mass of dry air
#' @param e numeric, partial pressure of water vapour (in Pa)
#' @param p numeric, atmospheric pressure (in Pa)
#' @return numeric, water vapour mixing ratio (kg kg^-1)
#' @export
mixing_ratio <- function(e, p) {

  w_v <- cldphys::const()$eps * e / (p - e)

  return(w_v)

}


#' Specific humidity
#' @param e numeric, partial pressure of water vapour (in Pa)
#' @param p numeric, atmospheric pressure (in Pa)
#' @return numeric, specitic humidity (kg kg^-1)
#' @export
spec_hum <- function(e, p) {

  eps <- cldphys::const()$eps

  q_v <- eps * e / (p - e + eps * e)

  return(q_v)

}
