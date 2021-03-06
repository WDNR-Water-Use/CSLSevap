% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radiation_functions.R
\name{R_a}
\alias{R_a}
\title{Radiation, extraterrestrial solar}
\usage{
R_a(dt, datetimes, phi, Lm, Lz, Gsc = 0.082)
}
\arguments{
\item{dt}{a string indicating the timestep of calculation. If not "hourly",
assumes a daily or larger timestep}

\item{datetimes}{timeseries of datetimes to use in calculations \link{POSIXct}}

\item{phi}{latitude of the location (radians). Positive for northern
hemisphere, negative for southern hemisphere.}

\item{Lm}{longitude of the measurement site (degrees west of Greenwich).}

\item{Lz}{longitude of the center of the time zone used for measurements
(degrees west of Greenwich). For example, Lz = 75, 90, 105 and 120°
for the Eastern, Central, Rocky Mountain and Pacific time zones
(United States) and Lz = 0° for Greenwich, 330° for Cairo (Egypt),
and 255° for Bangkok (Thailand).}

\item{Gsc}{solar constant (MJ/m^2/min). Defaults to 0.0820 MJ/m^2/min.}
}
\value{
\item{Ra}{extraterrestrial radiation (MJ/m^2/day) or (MJ/m^2/hour),
depending on length of input time period}
}
\description{
Calculates the extraterrestrial radiation for a given location and timestep.
For daily periods, uses Equation 21 from Allen et al. (1998). For hourly or
shorter periods, uses Equation 28 from Allen et al. (1998).
}
\references{
Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
evapotranspiration: Guidelines for computing crop water requirements. Rome:
FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
}
