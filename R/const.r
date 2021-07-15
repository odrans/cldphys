#' Load useful constants
#' @return list containing constants
#' @export
const <- function() {

  l <- list(
    Rd = 287, # Specific gas constant for dry air (J kg^-1 K^-1)
    Rv = 461.5, # Specific gas constant for water vapor (J kg^-1 K^-1)
    Mw = 18 * 1E-3  ## Molar mass of water (kg mol^-1)
  )

  l$eps <- l$Rd / l$Rv # Dry to moist gas constant ratio

  return(l)

}
