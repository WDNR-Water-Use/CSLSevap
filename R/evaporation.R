#' Evaporation (all methods)
#'
#' Calculates evaporation (mm/day given daily or monthly inputs, or mm/hour
#' given hourly inputs) using specified method. Options include:
#' \enumerate{
#'   \item "FAO" - FAO Penman-Montieth reference evapotranspiration (potential
#'                 evaporation for a reference grass crop) based on Allen et al.
#'                 (1998) for hourly, daily, or monthly timestep.
#'   \item "McJannet" - Daily lake evaporation based on the method in McJannet et al.
#'                      (2008) as presented in McMahon et al. (2013)
#'   \item "Hamon" - Daily lake evaporation based on the unmodified Hamon method
#'                   as presented in Harwell (2012)
#' }
#'
#' When assessing which parameters are required for a given method, use the
#' following guidelines:
#' \itemize{
#' \item "loc" - FAO (all), McJanet (all), Hamond (only "phi")
#' \item "weather" - FAO (all but "z0"), McJanet (all), Hamond (only "datetimes"
#'                   and "atmp" list with min and max daily atmp)
#' \item "lake" - FAO (not required, use default NULL value), McJanet (all but
#'                "lst"), Hamond (not required, use default NULL value).
#'  }
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
#' @references Harwell, G.R., 2012, Estimation of evaporation from open water—A
#'   review of selected studies, summary of U.S. Army Corps of Engineers data
#'   collection and methods, and evaluation of two methods for estimation of
#'   evaporation from five reservoirs in Texas: U.S. Geological Survey
#'   Scientific Investigations Report 2012–5202, 96 p.
#'
#' @param method denotes which evaporation method to use ("FAO", "McJannet", or
#'               "Hamon").
#' @param loc a list with location information that includes:
#' \itemize{
#' \item "z" - elevation above mean sea level (m)
#' \item"phi" - latitude of location (radians). Positive for northern
#'              hemisphere, negative for southern hemisphere.
#' \item"Lm" - longitude of location (degrees west of Greenwich).
#' \item"Lz" - longitude of location's measurement timezone (degrees west of
#'             Greenwich). For example, Lz = 75, 90, 105 and 120° for
#'             measurement times based on the Eastern, Central, Rocky Mountain
#'             and Pacific time zones (United States) and Lz = 0° for Greenwich,
#'             330° for Cairo (Egypt), and 255° for Bangkok (Thailand).
#' }
#' @param weather a list with weather data that includes:
#' \itemize{
#' \item "dt" - string indicating the timestep of input weather series. Expects
#'              "hourly", "daily", or "monthly".
#' \item "datetimes" - datetimes of weather records [POSIXct]. If monthly
#'                     timestep, make sure date is the 15th of each month.
#' \item "atmp" - If hourly timestep, vector of air temperature (degrees C)
#'                corresponding with datetimes vector. If daily or monthly
#'                timestep, list with two vectors, "min" and "max", with mean
#'                daily min and max air temperature (degrees C) corresponding
#'                with datetimes vector
#' \item "RH" - If hourly timestep, vector of relative humidity (percent)
#'              corresponding with datetimes vector. If daily or monthly
#'              timestep, list with two vectors, "min" and "max", with mean
#'              daily min and max relative humidity (percent) corresponding with
#'              datetimes vector.
#' \item "Rs" - vector of incoming solar or shortwave radiation (MJ/m^2/hr if
#'              hourly timestep, MG/m^2/day if daily or monthly), corresponding
#'              with datetimes vector.
#' \item "wind" - vector with mean wind speed (m/s), corresponding with
#'                datetimes vector.
#' \item "wind_elev" - atomic number, height at which wind is measured (m)
#' \item "z0" - aerodynamic roughness of weather measurement site (m)
#' }
#'
#' @param lake A list with lake data. Defaults to NULL, but for McJannet lake
#'             evaporation calculations, should include:
#' \itemize{
#'   \item "A" - surface area of the lake (km^2).
#'   \item "depth_m" - depth of the lake (m). Can be a static value or vector
#'                     corresponding with datetimes vector.
#'   \item "lst" - optional data frame with date (datetime) and ltmp (lake
#'                 temperature, degC).
#'   \item "wtmp0" - required initial water temperature for first day in
#'                   datetimes (degC).
#'
#' }
#' @param albedo a list with albedos for different surfaces, defaults to:
#' \itemize{
#'   \item "ref_crop" - albedo of the hypothetical grass reference crop, 0.23
#'   \item "water" - albedo of water, 0.08.
#' }
#' @param rho_a density of the air (kg/m^3), defaults to 1.20 kg/m^3 at 20 deg C
#' @param ca specific heat of the air (MJ/kg/K), defaults to 0.001013 MJ/kg/K
#' @param no_condensation defaults to TRUE to force negative evapotranspiration
#'                        values (i.e., condensation) to zero
#'
#' @return \item{evap}{Evapo(transpi)ration (mm/hour or mm/day)}
#'
#' @importFrom NISTunits NISTdegCtOk NISTinchTOmeter
#' @import lubridate
#'
#' @export
evaporation <- function(method = "FAO", loc, weather, lake = NULL,
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
    vpd   <- vp_deficit(weather$atmp, weather$RH)

    # Wind at 2m
    if (weather$wind_elev != 2){
      weather$wind <- u2_fcn(weather$wind, weather$wind_elev)
    }

    # Radiation
    Rn <- R_n(method, loc, lake, weather, albedo)

    # Soil Heat
    G  <- heat_flux(method, loc, lake, weather, albedo, Rn)

    evap <- (0.408*Delta*(Rn - G) + gamma*(900/tmp_K)*weather$wind*vpd)/
            (Delta + gamma*(1 + 0.34*weather$wind))

  } else if (method == "McJannet") {
    # MCJANNET LAKE EVAPORATION ------------------------------------------------
    # Aerodynamic resistance
    ra   <- aero_resist(weather$wind, weather$wind_elev, weather$z0, lake$A, loc$z)

    # Lake water temperature
    wtmp <- tmp_water(loc, lake, weather, albedo)

    # Vapour parameters
    Delta_w <- vp_sat_curve_slope(wtmp)
    es_w    <- vp_sat_mean(wtmp)
    ea      <- vp_act_mean(weather$atmp, weather$RH)
    lambda  <- latent_heat_vapor(weather$atmp)
    gamma   <- psychrometric_constant(loc$z)

    # Radiation
    Rn <- R_n(method, loc, lake, weather, albedo)

    # Water Heat Flux
    Gw <- heat_flux(method, loc, lake, weather, albedo)

    # Lake evaporation
    evap <- (Delta_w*(Rn - Gw) + 60*60*24*rho_a*ca*(es_w - ea)/ra)/
            (lambda*(Delta_w + gamma))
  } else if (method == "Hamon") {
    # HAMON UNMODIFIED LAKE EVAPORATION ----------------------------------------
    # Maximum possible daylight hours
    J       <- yday(weather$datetimes)
    delta   <- declination(J)
    omega_s <- hour_angle_sunset(loc$phi, delta)
    D       <- daylight_hours(omega_s)

    # Saturated vapor density
    svd <- sat_vapor_density(weather$atmp)

    # Lake evaporation
    evap <- 0.55*((D/12)^2)*(svd/100)
    evap <- NISTinchTOmeter(evap)*1000
  }

  # Adjust for no condensation, if needed
  if (no_condensation){
    evap[evap < 0] <- 0
  }

  return(evap)
}
