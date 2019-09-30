# ET_functions.R
# Includes:
# - FAO_ET

#' FAO Penman-Monteith Evapotranspiration
#'
#' Calculated potential evapotranspiration for a reference grass crop using the
#' FAO Penman-Montieth Equation. Calculates reference ET in mm/hr, if given
#' hourly input data, or in mm/day, if given daily or monthly input data.
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param loc a list with location information that includes:
#' \itemize{
#' \item z - elevation above mean sea level (m)
#' \item phi - latitude (radians). Positive for northern hemisphere, negative
#'             for southern hemisphere.
#' \item Lm - longitude of location (degrees west of Greenwich).
#' \item Lz - longitude of location's timezone (degrees west of Greenwich). For
#'            example, Lz = 75, 90, 105 and 120° for the Eastern, Central, Rocky
#'            Mountain and Pacific time zones (United States) and Lz = 0° for
#'            Greenwich, 330° for Cairo (Egypt), and 255° for Bangkok (Thailand)
#' }
#' @param weather a list with weather data that includes:
#' \itemize{
#' \item dt - string indicating the timestep of input weather series. Expects
#'            "hourly", "daily", or "monthly".
#' \item datetimes - datetimes of weather records [POSIXct]. If monthly
#'                   timestep, make sure date is the 15th of each month.
#' \item atmp - If hourly timestep, vector of air temperature (degrees C)
#'            corresponding with datetimes vector. If daily or monthly timestep,
#'            list with two vectors, "min" and "max", with mean daily min and
#'            max air temperature (degrees C) corresponding with datetimes
#'            vector
#' \item RH - If hourly timestep, vector of relative humidity (percent)
#'            corresponding with datetimes vector. If daily or monthly timestep,
#'            list with two vectors, "min" and "max", with mean daily min and
#'            max relative humidity (percent) corresponding with datetimes
#'            vector.
#' \item Rs - vector of incoming solar or shortwave radiation (MJ/m^2/timestep),
#'            corresponding with datetimes vector.
#' \item wind - vector with mean windspead (m/s), corresponding with datetimes
#'              vector.
#' \item wind_elev - atomic number, height at which wind is measured (m)
#' }
#' @param albedo albedo or canopy reflection coefficient (-). Defaults to 0.23
#'               for the hypothetical grass reference crop.
#' @param Gsc solar constant (MJ/m^2/min). Defaults to 0.0820 MJ/m^2/min.
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#'
#' @return \item{ET}{Evapotranspiration (mm/hour or mm/day, depending on dt)}
#'
#' @importFrom NISTunits NISTdegCtOk
#'
#' @export
FAO_ET <- function(loc, weather, albedo = 0.23, Gsc = 0.0820, SBc = 4.903e-9){
  # Mean temperature as deg C and K
  if (weather$dt == "hourly"){
    tmp_C <- weather$atmp
  } else {
    tmp_C <- (weather$atmp$min + weather$atmp$max)/2
  }
  tmp_K <- NISTdegCtOk(tmp_C)

  # Vapour parameters
  Delta <- FAO_slope_es_curve(weather$dt, weather$atmp)
  gamma <- FAO_psychrometric_constant(loc$z)
  vpd   <- FAO_vpd(weather$dt, weather$atmp, weather$RH)

  # Wind at 2m
  if (weather$wind_elev != 2){
    weather$wind <- FAO_u2(weather$wind, weather$wind_elev)
  }

  # Radiation
  Rn <- FAO_Rn(loc, weather, albedo, Gsc, SBc)

  # Soil Heat
  G  <- FAO_G(weather$dt, weather$datetimes, loc$phi, loc$Lm, loc$Lz, Rn,
              weather$atmp)

  ET <- (0.408*Delta*(Rn - G) + gamma*(900/tmp_K)*weather$wind*vpd)/
        (Delta + gamma*(1 + 0.34*weather$wind))

  return(ET)
}
