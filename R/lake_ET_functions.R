#' Evaporation
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
#' \item z - elevation above mean sea level (m)
#' \item phi - latitude (radians). Positive for northern hemisphere, negative
#'             for southern hemisphere.
#' \item Lm - longitude of location (degrees west of Greenwich). Only needed if
#'            input data is hourly.
#' \item Lz - longitude of location's timezone (degrees west of Greenwich). For
#'            example, Lz = 75, 90, 105 and 120° for the Eastern, Central, Rocky
#'            Mountain and Pacific time zones (United States) and Lz = 0° for
#'            Greenwich, 330° for Cairo (Egypt), and 255° for Bangkok
#'            (Thailand). Only needed if input is hourly.
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
#' @param albedo a list with albedos for different surfaces including:
#' \itemize{
#' \item ref_crop - albedo of the hypothetical grass reference crop, 0.23
#' \item lake - albedo of the lake, albedo for water is 0.08.
#' }
#' @param Gsc solar constant (MJ/m^2/min). Defaults to 0.0820 MJ/m^2/min.
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#' @param rho_a density of the air (kg/m^3), defaults to 1.20 kg/m^3 at 20 deg C
#' @param ca specific heat of the air (MJ/kg/K), defaults to 0.001013 MJ/kg/K
#'
#' @return \item{ET}{Lake Evapotranspiration (mm/day)}
#'
#' @export

McJannet_lake_ET <- function(loc, lake, weather, albedo, Gsc = 0.0820,
                             SBc = 4.903e-9, rho_a = 1.20, ca = 0.001013) {
  # Aerodynamic resistance
  ra <- McJannet_aero_resist(weather$wind, weather$wind_elev, weather$z0,
                             lake$A, loc$z)

  wtmp <- McJannet_wtmp(loc, lake, weather)

  # Vapour parameters
  Delta_w <- FAO_slope_es_curve(wtmp)
  es_w    <- FAO_mean_es(wtmp)
  ea      <- FAO_mean_ea(weather$atmp, weather$RH)
  lambda  <- latent_heat_vapor(weather$atmp)
  gamma   <- FAO_psychrometric_constant(loc$z)

  # Radiation
  Rn <- McMahon_Rn(loc, lake, weather, albedo$lake)

  # Water Heat Flux
  Gw <- McJannet_Gw(loc, lake, weather)

  ET <- (Delta_w*(Rn - Gw) + 60*60*24*rho_a*ca*(es_w - ea)/ra)/
        (lambda*(Delta_w + gamma))

  return(ET)
}
