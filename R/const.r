#' Load useful constants
#' @return list containing constants
#' @export
const <- function() {

  l <- list(
    rho_w = 1E3, # Water density (kg m^-3)
    cp = 1003.5, # Specific heat of dry air at constant pressure (J kg^-1 K^-1)
    g = 9.81, # Gravitational acceleration constant (m s^-2)
    Rd = 287, # Specific gas constant for dry air (J kg^-1 K^-1)
    Rv = 461.5, # Specific gas constant for water vapor (J kg^-1 K^-1)
    Mw = 18 * 1E-3,  ## Molar mass of water (kg mol^-1)
    gamma_env = 6.49E-3 ## Environmental lapse rate (K m^-1)
  )

  l$eps <- l$Rd / l$Rv # Dry to moist gas constant ratio

  return(l)

}
