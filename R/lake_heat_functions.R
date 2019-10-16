#lake_heat_functions.R
# Includes
# - McJannet_time_const
# - McJannet_eqtmp

# ------------------------------------------------------------------------------
#' Change in Heat Storage
#'
#' Calculates the change in heat storage based on water temperature following
#' McJannet et al. (2008) Equation 31 as presented in McMahon et al. S11.33
#' (2013).
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R.: Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363,
#'   https://doi.org/10.5194/hess-17-1331-2013, 2013.
#'
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
#' \item z0 - aerodynamic roughness of measurement site (m)
#' \item wtmp0 - initial water temperature for first day in datetimes.
#' }
#' @param rho_w density of water (kg/m^3), defaults to 997.9 kg/m^3 at 20 deg C
#' @param cw specific heat of water (MJ/kg/K), defaults to 0.00419
#'
#' @return {Gw}{change in heat storage (MJ/m^2/day)}
#'
#' @export
McJannet_Gw <- function(loc, lake, weather, rho_w = 997.9, cw = 0.00419) {
  wtmp0   <- weather$wtmp0
  wtmp    <- McJannet_wtmp(loc, lake, weather)
  Gw      <- NULL
  if (length(lake$depth_m) == 1) {
    lake$depth_m <- rep(lake$depth_m, length(wtmp))
  }
  for (i in 1:length(wtmp)){
    Gw[i]    <- rho_w*cw*lake$depth_m[i]*(wtmp[i] - wtmp0)
    wtmp0    <- wtmp[i]
  }

  return(Gw)
}

# ------------------------------------------------------------------------------
#' Water Temperature
#'
#' Calculates the water temperature based on the water temperature of the
#' previous day following McJannet et al. (2008) Equation 23 and de Bruin (1982)
#' Equation 10 as presented in McMahon et al. S11.28 (2013).
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R.: Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363,
#'   https://doi.org/10.5194/hess-17-1331-2013, 2013.
#'
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
#' \item z0 - aerodynamic roughness of measurement site (m)
#' \item wtmp0 - initial water temperature for first day in datetimes.
#' }
#'
#' @return {wtmp}{water temperature (degrees C)}
#'
#' @export
McJannet_wtmp <- function(loc, lake, weather) {
  lst   <- lake$lst
  eqtmp <- McJannet_eqtmp(loc, lake, weather)
  Ctime <- McJannet_time_const(loc, lake, weather)
  wtmp0 <- weather$wtmp0
  wtmp  <- NULL
  for (i in 1:length(weather$datetimes)) {
    today <- weather$datetimes[i]
    if (today %in% floor_date(lst$date, unit = "day")) {
      # Use input temp if have it
      wtmp[i] <- mean(lst$ltmp[which(floor_date(lst$date, unit = "day") == today)])
    } else {
      # Estimate from previous day temp if no input temp
      wtmp[i] <- eqtmp[i] + (wtmp0 - eqtmp[i])*exp(-1/Ctime[i])
    }
    wtmp0   <- wtmp[i]
  }
  return(wtmp)
}

# ------------------------------------------------------------------------------
#' Equilibrium Temperature
#'
#' Calculates the equilibrium temperature based on de Bruin (1982) Equation 3,
#' as presented in McMahon et al. (2013) S11.30.
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
#' \item wind - vector with mean windspead (m/s), corresponding with datetimes
#'              vector.
#' \item wind_elev - atomic number, height at which wind is measured (m)
#' \item z0 - aerodynamic roughness of measurement site (m)
#' }
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#'
#' @return {eqtmp}{equilibrium temperature (degrees C)}
#'
#' @export
McJannet_eqtmp <- function(loc, lake, weather, SBc = 4.903e-9) {
  u10     <- McMahon_u10(weather$wind, weather$wind_elev, weather$z0)
  u_fcn   <- McJannet_u_fcn(u10, lake$A)
  gamma   <- FAO_psychrometric_constant(loc$z)
  wbtmp   <- McJannet_wbtmp(weather$atmp, weather$RH)
  Deltawb <- FAO_slope_es_curve(wbtmp)
  Rnwb    <- McMahon_Rnwb(loc, weather)

  eqtmp   <- wbtmp + Rnwb/(4*SBc*(wbtmp + 273.15)^3 + u_fcn*(Deltawb + gamma))

  return(eqtmp)
}

# ------------------------------------------------------------------------------
#' Time Constant
#'
#' Calculates the time constant (day) based on McJannet et al. (2008) Equation 5
#' and de Bruin (1982) Equation 4, as presented in McMahon et al. (2013) S11.29.
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
#' \item wind - vector with mean windspead (m/s), corresponding with datetimes
#'              vector.
#' \item wind_elev - atomic number, height at which wind is measured (m)
#' \item z0 - aerodynamic roughness of measurement site (m)
#' }
#' @param rho_w density of water (kg/m^3), defaults to 997.9 kg/m^3 at 20 deg C
#' @param cw specific heat of water (MJ/kg/K), defaults to 0.00419
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#'
#' @return {Ctime}{time constant (day)}
#'
#' @export
McJannet_time_const <- function(loc, lake, weather, rho_w = 997.9, cw = 0.00419,
                                SBc = 4.903e-9) {
  u10     <- McMahon_u10(weather$wind, weather$wind_elev, weather$z0)
  u_fcn   <- McJannet_u_fcn(u10, lake$A)
  gamma   <- FAO_psychrometric_constant(loc$z)
  wbtmp   <- McJannet_wbtmp(weather$atmp, weather$RH)
  Deltawb <- FAO_slope_es_curve(wbtmp)

  Ctime   <- rho_w*cw*lake$depth_m/(4*SBc*(wbtmp + 273.15)^3 + u_fcn*(Deltawb + gamma))

  return(Ctime)
}
