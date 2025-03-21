% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hom_ice_nucleation.r
\name{nucl_rate_ice_hom_wat}
\alias{nucl_rate_ice_hom_wat}
\title{Nuleation rate for ice crystals from homogeneous nucleation
Theory by Pruppacher and Klett (1997), chapter 7}
\usage{
nucl_rate_ice_hom_wat(ta, fix_Cpre = FALSE)
}
\arguments{
\item{ta}{air temperature (in K)}

\item{fix_Cpre}{logical, if TRUE the prefactor is fixed to 10^41 regardless of the temperature}

\item{Si}{Supersaturation ratio with respect to ice}
}
\value{
Nucleation rate (in m^-3 s^-1)
}
\description{
Nuleation rate for ice crystals from homogeneous nucleation
Theory by Pruppacher and Klett (1997), chapter 7
}
