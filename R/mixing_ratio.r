#' Water vapour mixing ratio
#' @param e numeric, partial pressure of water vapour (in Pa)
#' @param p numeric, atmospheric pressure (in Pa)
#' @return numeric, water vapour mixing ratio (-)
mixing_ratio <- function(e, p) {

  w <- cldphys::const()$eps * e / (p - e)

  return(w)

}
