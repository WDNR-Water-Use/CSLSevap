#lake_radiation_functions.R
# Includes
# - McMahon_Rnwb
# - McMahon_Rolwb
# - McMahon_Ril
# - McJannet_Cf

# ------------------------------------------------------------------------------
#' Net Radiation at Water Temperature
#'
#' Calculates the net radiation for a lake as a function of water temperature
#' based on McMahon et al. (2013) Equation S11.25.
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R.: Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363,
#'   https://doi.org/10.5194/hess-17-1331-2013, 2013.
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
#' @param lake a list with lake data that includes:
#' \itemize{
#' \item A - surface area of the lake (km^2).
#' \item depth_m - depth of the lake (m).
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
#' @param albedo albedo of the lake, defaults to albedo for water, or 0.08.
#' @param Gsc solar constant (MJ/m^2/min). Defaults to 0.0820 MJ/m^2/min.
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#'
#' @return {Rnwb}{net radiation (MJ/m^2/day) at water temperature}
#'
#' @export
McMahon_Rn <- function(loc, lake, weather, albedo = 0.08, Gsc = 0.0820,
                         SBc = 4.903e-9) {

  Rns  <- FAO_Rns(weather$Rs, albedo)
  Ra   <- FAO_Ra(weather$dt, weather$datetimes, loc$phi, loc$Lm, loc$Lz, Gsc)
  Rso  <- FAO_Rso(Ra, loc$z)
  Cf   <- McJannet_Cf(weather$Rs, Rso)
  Ril  <- McMahon_Ril(Cf, weather$atmp, SBc)
  wtmp <- McJannet_wtmp(loc, lake, weather)
  Rol  <- McMahon_Rol(wtmp, SBc)

  Rn  <- Rns + (Ril - Rol)

  return(Rn)
}

# ------------------------------------------------------------------------------
#' Net Radiation at Wet Bulb Temperature
#'
#' Calculates the net radiation for a lake as a function of wet bulb temperature
#' based on McMahon et al. (2013) Equation S11.31.
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R.: Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363,
#'   https://doi.org/10.5194/hess-17-1331-2013, 2013.
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
#' @param albedo albedo of the lake, defaults to albedo for water, or 0.08.
#' @param Gsc solar constant (MJ/m^2/min). Defaults to 0.0820 MJ/m^2/min.
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#'
#' @return {Rnwb}{net radiation (MJ/m^2/day) at wet-bulb temperature}
#'
#' @export
McMahon_Rnwb <- function(loc, weather, albedo = 0.08, Gsc = 0.0820,
                         SBc = 4.903e-9) {

  Rns   <- FAO_Rns(weather$Rs, albedo)
  Ra    <- FAO_Ra(weather$dt, weather$datetimes, loc$phi, loc$Lm, loc$Lz, Gsc)
  Rso   <- FAO_Rso(Ra, loc$z)
  Cf    <- McJannet_Cf(weather$Rs, Rso)
  Ril   <- McMahon_Ril(Cf, weather$atmp, SBc)
  wbtmp <- McJannet_wbtmp(weather$atmp, weather$RH)
  Rolwb <- McMahon_Rolwb(weather$atmp, wbtmp, SBc)

  Rnwb  <- Rns + (Ril - Rolwb)

  return(Rnwb)
}

# ------------------------------------------------------------------------------
#' Outgoing Longwave Radiation at Wet Bulb Temperature
#'
#' Calculates the outgoing longwave radiation for a lake as a function of wet
#' bulb temperature based on McMahon et al. (2013) Equation S11.32.
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R.: Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363,
#'   https://doi.org/10.5194/hess-17-1331-2013, 2013.
#'
#' @param atmp air temperature (degrees C). When dt is "daily" or larger,
#'             argument should be a list with elements "min" and "max" for daily
#'             min and max temperatures. When dt is "hourly", argument should be
#'             a vector with hourly recorded temperatures.
#' @param wbtmp daily wet bulb temperature (degrees C).
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#'
#' @return {Rolwb}{outgoing longwaver radiation (MJ/m^2/day) at wet-bulb
#'                 temperature}
#'
#' @export
McMahon_Rolwb <- function(atmp, wbtmp, SBc = 4.903e-9) {
  if (class(atmp) == "list") { atmp <- (atmp$max + atmp$min)/2 }
  Rolwb <- SBc*(atmp + 273.15)^4 + 4*SBc*(atmp + 273.15)^3*(wbtmp - atmp)
  return(Rolwb)
}

# ------------------------------------------------------------------------------
#' Outgoing Longwave Radiation at Water Temperature
#'
#' Calculates the outgoing longwave radiation for a lake as a function of wet
#' bulb temperature based on McMahon et al. (2013) Equation S11.27.
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R.: Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363,
#'   https://doi.org/10.5194/hess-17-1331-2013, 2013.
#'
#' @param wtmp daily water temperature (degrees C).
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#'
#' @return {Rol}{outgoing longwaver radiation (MJ/m^2/day) at the temperature of
#'               the water}
#'
#' @export
McMahon_Rol <- function(wtmp, SBc = 4.903e-9) {
  Rol <- 0.97*SBc*(wtmp + 273.15)^4
  return(Rol)
}

# ------------------------------------------------------------------------------
#' Incoming Longwave Radiation
#'
#' Calculates the incoming longwave radiation using cloud cover fraction and
#' mean daily air temperature based on McMahon et al. (2013) Equation S11.27.
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R.: Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363,
#'   https://doi.org/10.5194/hess-17-1331-2013, 2013.
#'
#' @param Cf fraction of cloud cover (-)
#' @param atmp air temperature (degrees C). When dt is "daily" or larger,
#'             argument should be a list with elements "min" and "max" for daily
#'             min and max temperatures. When dt is "hourly", argument should be
#'             a vector with hourly recorded temperatures.
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#'
#' @return {Cf}{the fraction of cloud cover (-)}
#'
#' @export
McMahon_Ril <- function(Cf, atmp, SBc = 4.903e-9) {
  if (class(atmp) == "list") { atmp <- (atmp$max + atmp$min)/2 }
  Ril <- (Cf + (1 - Cf)*(1 - (0.261*exp(-7.77e-4*atmp^2))))*SBc*(atmp + 273.15)^4
  return(Ril)
}

# ------------------------------------------------------------------------------
#' Cloudiness Factor
#'
#' Estimates the fraction of cloud cover based on McJannet et al. (2008b,
#' Equations 14 and 15), as presented by McMahon et al. (2013), Equations S3.17
#' and S3.18.
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R.: Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363,
#'   https://doi.org/10.5194/hess-17-1331-2013, 2013.
#'
#' @param Rs incoming solar radiation (MJ/m^2/day)
#' @param Rso clear sky radiation (MJ/m^2/day).
#'
#' @return {Cf}{the fraction of cloud cover (-)}
#'
#' @export
McJannet_Cf <- function(Rs, Rso) {
  Kratio <- Rs/Rso
  Cf     <- NULL
  for (i in 1:length(Kratio)) {
    if (Kratio[i] <= 0.9) {
      Cf[i] <- 1.1 - Kratio[i]
    } else {
      Cf[i] <- 2*(1 - Kratio[i])
    }
  }
  return(Cf)
}
