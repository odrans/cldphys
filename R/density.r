
#' Density of ice
#' Parameterization from Pruppacher and Klett (1997)
#' @param t Temperature (in K)
#' @return Density of ice (in kg/m^3)
#' @export
density_ice <- function(ta) {
  tc <- ta - 273.15

  rho_ice <- 916.7 - 0.175*tc - 5E-4*tc^2

  return(rho_ice)
}

#' Density of water
#' Parameterization from Pruppacher and Klett (1997)
#' @param T Temperature (in K)
#' @return Density of water (in kg/m^3)
#' @export
density_wat <- function(ta) {
  # Convert temperature from Kelvin to Celsius
  tc <- ta - 273.15

  # Pre-allocate rho_wat vector with same length as tc
  rho_wat <- numeric(length(tc))

  # Condition to identify values greater than or equal to 0 degrees Celsius
  positive_mask <- tc >= 0

  # For temperatures greater than or equal to 0 degrees Celsius
  if (any(positive_mask)) {
    a0 <- 0.99986
    a1 <- 6.690 * 1E-5
    a2 <- -8.486 * 1E-6
    a3 <- 1.518 * 1E-7
    a4 <- -6.9984 * 1E-9
    a5 <- -3.6449 * 1E-10
    a6 <- -7.497 * 1E-12

    tc_pos <- tc[positive_mask]
    rho_wat[positive_mask] <- (a0 + a1 * tc_pos + a2 * tc_pos^2 + a3 * tc_pos^3 + a4 * tc_pos^4 +
                               a5 * tc_pos^5 + a6 * tc_pos^6) * 1E3
  }

  # For temperatures below 0 degrees Celsius
  if (any(!positive_mask)) {
    A0 <- 999.8396
    A1 <- 18.224944
    A2 <- -7.922210 * 1E-3
    A3 <- -55.44846 * 1E-6
    A4 <- 149.7562 * 1E-9
    A5 <- -393.2952 * 1E-12
    B <- 18.159725 * 1E-3

    tc_neg <- tc[!positive_mask]
    rho_wat[!positive_mask] <- (A0 + A1 * tc_neg + A2 * tc_neg^2 + A3 * tc_neg^3 +
                               A4 * tc_neg^4 + A5 * tc_neg^5) / (1 + B * tc_neg)
  }

  return(rho_wat)
}
