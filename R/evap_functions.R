# evap_functions.R
# Includes:
# - FAO_evap
# - lake_evap

# ------------------------------------------------------------------------------
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
FAO_evap <- function(loc, weather, albedo = 0.23, Gsc = 0.0820, SBc = 4.903e-9){
  # Mean temperature as deg C and K
  if (weather$dt == "hourly"){
    tmp_C <- weather$atmp
  } else {
    tmp_C <- (weather$atmp$min + weather$atmp$max)/2
  }
  tmp_K <- NISTdegCtOk(tmp_C)

  # Vapour parameters
  Delta <- FAO_slope_es_curve(weather$atmp)
  gamma <- FAO_psychrometric_constant(loc$z)
  vpd   <- FAO_vpd(weather$atmp, weather$RH)

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

# ------------------------------------------------------------------------------
#' Lake Evaporation
#'
#' Calculates daily evaporation from a shallow (10m or less) lake based on the
#' Penman equation, as modified by Finch (2001), and as presented in Equations
#' S11.1 and S11.2 in McMahon et al. (2013).
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R.: Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363,
#'   https://doi.org/10.5194/hess-17-1331-2013, 2013.
#'
#' @param loc a list with location information that includes:
#' \itemize{
#'   \item z - elevation above mean sea level (m)
#'   \item phi - latitude (radians). Positive for northern hemisphere, negative
#'             for southern hemisphere.
#'   \item Lm - longitude of location (degrees west of Greenwich). Only needed if
#'            input data is hourly.
#'   \item Lz - longitude of location's timezone (degrees west of Greenwich).
#'              For example, Lz = 75, 90, 105 and 120° for the Eastern, Central,
#'              Rocky Mountain and Pacific time zones (United States) and Lz =
#'              0° for Greenwich, 330° for Cairo (Egypt), and 255° for Bangkok
#'              (Thailand). Only needed if input is hourly.
#' }
#' @param lake a list with lake data that includes:
#' \itemize{
#'   \item A - surface area of the lake (km^2).
#'   \item depth_m - depth of the lake (m).
#'   \item lst - data frame with date (datetime) and ltmp (lake temperature, degC)
#' }
#' @param weather a list with weather data that includes:
#' \itemize{
#'   \item dt - string indicating the timestep of input weather series. Expects
#'              "hourly", "daily", or "monthly".
#'   \item datetimes - datetimes of weather records [POSIXct]. If monthly
#'                     timestep, make sure date is the 15th of each month.
#'   \item atmp - If hourly timestep, vector of air temperature (degrees C)
#'                corresponding with datetimes vector. If daily or monthly
#'                timestep, list with two vectors, "min" and "max", with mean
#'                daily min and max air temperature (degrees C) corresponding
#'                with datetimes vector.
#'   \item RH - If hourly timestep, vector of relative humidity (percent)
#'              corresponding with datetimes vector. If daily or monthly
#'              timestep, list with two vectors, "min" and "max", with mean
#'              daily min and max relative humidity (percent) corresponding with
#'              datetimes vector.
#'   \item Rs - vector of incoming solar or shortwave radiation (MJ/m^2/timestep),
#'              corresponding with datetimes vector.
#'   \item wind - vector with mean windspead (m/s), corresponding with datetimes
#'                vector.
#'   \item wind_elev - atomic number, height at which wind is measured (m)
#'   \item z0 - aerodynamic roughness of measurement site (m)
#'   \item wtmp0 - initial water temperature for first day in datetimes.
#' }
#' @param albedo a list with albedos for different surfaces including:
#' \itemize{
#'   \item ref_crop - albedo of the hypothetical grass reference crop, 0.23
#'   \item lake - albedo of the lake, albedo for water is 0.08.
#' }
#' @param rho_a density of the air (kg/m^3), defaults to 1.20 kg/m^3 at 20 deg C
#' @param ca specific heat of the air (MJ/kg/K), defaults to 0.001013 MJ/kg/K
#' @param no_condensation defaults to TRUE to force negative evapotranspiration
#'                        values (i.e., condensation) to zero
#'
#' @return \item{ET}{Lake Evapotranspiration (mm/day)}
#'
#' @export

lake_evap <- function(loc, lake, weather, albedo, rho_a = 1.20, ca = 0.001013,
                      no_condensation = TRUE) {
  # Aerodynamic resistance
  ra   <- aero_resist(weather$wind, weather$wind_elev, weather$z0, lake$A, loc$z)

  # Lake water temperature
  wtmp <- lake_wtmp(loc, lake, weather)

  # Vapour parameters
  Delta_w <- FAO_slope_es_curve(wtmp)
  es_w    <- FAO_mean_es(wtmp)
  ea      <- FAO_mean_ea(weather$atmp, weather$RH)
  lambda  <- latent_heat_vapor(weather$atmp)
  gamma   <- FAO_psychrometric_constant(loc$z)

  # Radiation
  Rn <- lake_Rn(loc, lake, weather, albedo$lake)

  # Water Heat Flux
  Gw <- lake_Gw(loc, lake, weather)

  # Lake evaporation
  evap <- (Delta_w*(Rn - Gw) + 60*60*24*rho_a*ca*(es_w - ea)/ra)/
          (lambda*(Delta_w + gamma))
  if (no_condensation){
    evap[evap < 0] <- 0
  }


  return(evap)
}
