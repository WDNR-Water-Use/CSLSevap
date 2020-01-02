#' Evaporation (All Methods)
#'
#' Calculates evaporation using input data and specified calculation method.
#' Options include "FAO" Penman-Montieth reference evapotranspiration or
#' "McJannet" lake evaporation. Returns evaporation in mm/day, unless specify
#' "FAO" and submt hourly imput data, which returns evaporation in mm/hour.
#'
#' "FAO" calculates FAO Penman-Monteith reference evapotranspiration (potential
#' evapotranspiration for a reference grass crop) based on Allen et al. (1998).
#'
#' "McJannet" calculates daily evaporation from a lake based on the method in
#' McJannet et al. (2008) as presented in Equation S11.22 in McMahon et al.
#' (2013).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @references McJannet, D. L., Webster, I. T., Stenson, M. P., and Sherman,
#'   B.S. (2008). Estimating open water evaporation for the Murray-Darling
#'   Basin. A report to the Australian Government from the CSIRO Murray-Darling
#'   Basin Sustainable Yields Project, CSIRO, Australia, 50 pp. Retrieved from
#'   http://www.clw.csiro.au/publications/waterforahealthycountry/mdbsy/technical/U-OpenWaterEvaporation.pdf.
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R. (2013). Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363.
#'   https://doi.org/10.5194/hess-17-1331-2013.
#'
#' @param method denotes which evaporation method to use ("FAO" or "McJannet").
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
#' \item z0 - aerodynamic roughness of weather measurement site (m)
#' }
#'
#' @param lake A list with lake data. Defaults to NULL, but for lake evaporation
#'             calculations should include:
#' \itemize{
#'   \item A - surface area of the lake (km^2).
#'   \item depth_m - depth of the lake (m).
#'   \item lst - data frame with date (datetime) and ltmp (lake temperature, degC)
#'   \item wtmp0 - initial water temperature for first day in datetimes.
#'
#' }
#' @param albedo a list with albedos for different surfaces including (defaults):
#' \itemize{
#'   \item ref_crop - albedo of the hypothetical grass reference crop, 0.23
#'   \item water - albedo of water, 0.08.
#' }
#' @param rho_a density of the air (kg/m^3), defaults to 1.20 kg/m^3 at 20 deg C
#' @param ca specific heat of the air (MJ/kg/K), defaults to 0.001013 MJ/kg/K
#' @param no_condensation defaults to TRUE to force negative evapotranspiration
#'                        values (i.e., condensation) to zero
#' @param Gsc solar constant (MJ/m^2/min). Defaults to 0.0820 MJ/m^2/min.
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#'
#' @return \item{evap}{Evapo(transpi)ration (mm/hour or mm/day)}
#'
#' @importFrom NISTunits NISTdegCtOk
#'
#' @export
evaporation <- function(method = "FAO", loc = NULL, weather, lake = NULL,
                        albedo = list(ref_crop = 0.23, water = 0.08),
                        no_condensation = TRUE, rho_a = 1.20, ca = 0.001013){
  if (method == "FAO") {
    # FAO PENMAN-MONTIETH REFERENCE EVAPOTRANSPIRATION -------------------------
    # Mean temperature as deg C and K
    if (weather$dt == "hourly"){
      tmp_C <- weather$atmp
    } else {
      tmp_C <- (weather$atmp$min + weather$atmp$max)/2
    }
    tmp_K <- NISTdegCtOk(tmp_C)

    # Vapour parameters
    Delta <- vp_sat_curve_slope(weather$atmp)
    gamma <- psychrometric_constant(loc$z)
    vpd   <- vpd(weather$atmp, weather$RH)

    # Wind at 2m
    if (weather$wind_elev != 2){
      weather$wind <- FAO_u2(weather$wind, weather$wind_elev)
    }

    # Radiation
    Rn <- R_n(method, loc, lake, weather, albedo$ref_crop)

    # Soil Heat
    G  <- FAO_G(weather$dt, weather$datetimes, loc$phi, loc$Lm, loc$Lz, Rn,
                weather$atmp)

    evap <- (0.408*Delta*(Rn - G) + gamma*(900/tmp_K)*weather$wind*vpd)/
            (Delta + gamma*(1 + 0.34*weather$wind))

  } else if (method == "McJannet") {
    # MCJANNET LAKE EVAPORATION ------------------------------------------------
    # Aerodynamic resistance
    ra   <- aero_resist(weather$wind, weather$wind_elev, weather$z0, lake$A, loc$z)

    # Lake water temperature
    wtmp <- lake_wtmp(loc, lake, weather)

    # Vapour parameters
    Delta_w <- vp_sat_curve_slope(wtmp)
    es_w    <- vp_sat_mean(wtmp)
    ea      <- vp_act_mean(weather$atmp, weather$RH)
    lambda  <- latent_heat_vapor(weather$atmp)
    gamma   <- psychrometric_constant(loc$z)

    # Radiation
    Rn <- R_n(method, loc, lake, weather, albedo$water)

    # Water Heat Flux
    Gw <- heat_flux(method, loc, lake, weather)

    # Lake evaporation
    evap <- (Delta_w*(Rn - Gw) + 60*60*24*rho_a*ca*(es_w - ea)/ra)/
            (lambda*(Delta_w + gamma))
  }

  # Adjust for no condensation, if needed
  if (no_condensation){
    evap[evap < 0] <- 0
  }

  return(evap)
}
