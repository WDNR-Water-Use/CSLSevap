% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/radiation_functions.R
\name{declination}
\alias{declination}
\title{Solar declination}
\usage{
declination(J)
}
\arguments{
\item{J}{Julian day, i.e., the number of the day in the year between 1 (1
January) and 365 or 366 (31 December).}
}
\value{
\item{delta}{solar declination (radians)}
}
\description{
Calculates the solar declination (radians) using Equation 24 from Allen et
al. (1998). Similar to Equation 2.1 in Harwell (2012), but with some rounding
differences.
}
\references{
Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
evapotranspiration: Guidelines for computing crop water requirements. Rome:
FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.

Harwell, G.R., 2012, Estimation of evaporation from open water—A
review of selected studies, summary of U.S. Army Corps of Engineers data
collection and methods, and evaluation of two methods for estimation of
evaporation from five reservoirs in Texas: U.S. Geological Survey
Scientific Investigations Report 2012–5202, 96 p.
}
