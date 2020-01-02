# vapour_functions.R
# Includes:
# - vp_sat
# - vp_sat_mean
# - vp_act_mean
# - vp_deficit
# - vp_sat_curve_slope
# - psychrometric_constant
# - latent_heat_vapor
# - tmp_dew
# - tmp_wet_bulb
#
# ------------------------------------------------------------------------------
#' Saturation Vapour Pressure
#'
#' Calculates the saturation vapour pressure at a given temperature based on
#' equation 11 of Allen et al. (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param tmp temperature of air or water (degrees C)
#'
#' @return eo - saturation vapour pressure at the temperature (kPa)
#'
#' @export

vp_sat <- function(tmp) {
  eo <- 0.6108*exp(17.27*tmp/(237.3 + tmp))
  return(eo)
}

# ------------------------------------------------------------------------------
#' Mean Saturation Vapour Pressure
#'
#' Calculates the mean saturation vapour pressure for an hourly, daily, or
#' larger time period based on Equation 12 (for daily or larger timestep) on pg.
#' 36 of Allen et al. (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param atmp air temperature (degrees C). When dt is "daily" or larger,
#'             argument should be a list with elements "min" and "max" for daily
#'             min and max temperatures. When dt is "hourly", argument should be
#'             a vector with hourly recorded temperatures.
#'
#' @return \item{es}{mean saturation vapour pressure (kPa) during time period of
#'                   interest}
#'
#' @export

vp_sat_mean <- function(atmp) {
  if (class(atmp) == "list"){
    es <- (vp_sat(atmp$max) + vp_sat(atmp$min))/2
  } else {
    es <- vp_sat(atmp)
  }
  return(es)
}

# ------------------------------------------------------------------------------
#' Mean Actual Vapour Pressure
#'
#' Calculates the mean actual vapour pressure for an hourly, daily, or larger
#' time period based on Equation 17 (for daily or larger timestep) on pg. 38 of
#' Allen et al. (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param atmp air temperature (degrees C). When dt is "daily" or larger,
#'             argument should be a list with elements "min" and "max" for daily
#'             min and max temperatures. When dt is "hourly", argument should be
#'             a vector with hourly recorded temperatures.
#' @param RH relative humidity (percent). When dt is "daily" or larger, argument
#'           should be a list with elements "min" and "max" for daily min and
#'           max relative humidities. When dt is "hourly", argument should be a
#'           vector with hourly recorded relative humidities.
#'
#' @return \item{ea}{mean actual vapour pressure (kPa) during time period of
#'                   interest (hourly or daily timestep)}
#'
#' @export

vp_act_mean <- function(atmp, RH) {
  if (class(atmp) == "list") {
    ea <- (vp_sat(atmp$min)*RH$max/100 + vp_sat(atmp$max)*RH$min/100)/2
  } else {
    ea <- vp_sat(atmp)*RH/100
  }
  return(ea)
}

# ------------------------------------------------------------------------------
#' Vapour Pressure Deficit
#'
#' Calculates the vapour pressure deficit for a given time period (daily timestep
#' or larger) using the mean daily minimum temperature, mean daily maximum
#' temperature, mean daily minimum relative humidity, and mean daily maxiumum
#' humidity during that time period. See p.39 in Allen et al. (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param atmp air temperature (degrees C). When dt is "daily" or larger,
#'             argument should be a list with elements "min" and "max" for daily
#'             min and max temperatures. When dt is "hourly", argument should be
#'             a vector with hourly recorded temperatures.
#' @param RH relative humidity (percent). When dt is "daily" or larger, argument
#'           should be a list with elements "min" and "max" for daily min and
#'           max relative humidities. When dt is "hourly", argument should be a
#'           vector with hourly recorded relative humidities.
#'
#' @return \item{vpd}{vapour pressure deficit (kPa) for the given time period}
#'
#' @export

vp_deficit <- function(atmp, RH) {
  ea  <- vp_act_mean(atmp, RH)
  es  <- vp_sat_mean(atmp)
  vpd <- es - ea
  return(vpd)
}

# ------------------------------------------------------------------------------
#' Slope of Saturation Vapour Pressure Curve
#'
#' Calculates the slope of the saturation vapour ressure curve (i.e., the slope
#' of the relationship between saturation vapour pressure and temperature) based
#' on Equation 13 in Allen et al. (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param atmp air temperature (degrees C). When dt is "daily" or larger,
#'             argument should be a list with elements "min" and "max" for daily
#'             min and max temperatures. When dt is "hourly", argument should be
#'             a vector with hourly recorded temperatures.
#'
#' @return \item{Delta}{slope of the saturation vapour pressure curve at the
#'                      given air temperature}
#'
#' @export

vp_sat_curve_slope <- function(atmp) {
  if (class(atmp) == "list"){
    atmp <- (atmp$min + atmp$max)/2
  }
  Delta <- (4098*(vp_sat(atmp)))/
           ((atmp + 237.3)^2)
  return(Delta)
}

# ------------------------------------------------------------------------------
#' Psychrometric Constant
#'
#' Calculates the psychrometric constant for a given elevation based on Equation
#' 7 and 8 in Allen et al. (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param z elevation above sea level (m)
#' @param atmp air temperature (degrees C). Optional value used to calculate the
#'             latent heat of vaporization. Defaults to NULL to use constant
#'             value of lambda.
#' @param lambda latent heat of vaporization (MJ/kg). Default is 2.45.
#' @param cp specific heat at constant pressure (MJ/kg/degC). Default is
#'           1.1013e-3.
#' @param epsilon ratio molecular weight of water vapour/dry air (-). Default is
#'                0.622
#'
#' @return \item{gamma}{the psychrometric constant for the given elevation
#'                      (kPa/degrees C)}
#'
#' @export

psychrometric_constant <- function(z, atmp = NULL, lambda = 2.45,
                                   cp = 1.013e-3, epsilon = 0.622) {
  #If air temp provided, use to calculate latent heat of vaporization
  if (is.null(atmp) == FALSE){lambda <- latent_heat_vapor(atmp)}

  #Atmospheric pressure (kPa)
  P     <- 101.3*(((293 - 0.0065*z)/293)^5.26)

  gamma <- (cp*P)/(epsilon*lambda)
  return(gamma)
}

# ------------------------------------------------------------------------------
#' Latent Heat of Vaporization
#'
#' Calculates the latent heat of vaporization as a function of air temperature.
#'
#' @references \url{https://cran.r-project.org/web/packages/bigleaf/bigleaf.pdf}
#'
#' @param atmp air temperature (degrees C)
#'
#' @return \item{lambda}{the latent heat of vaporization (MJ/kg)}
#'
#' @export

latent_heat_vapor <- function(atmp) {
  if (class(atmp) == "list"){ atmp <- (atmp$max + atmp$min)/2 }
  lambda <- 2.501 - 0.00237*atmp
  return(lambda)
}

# ------------------------------------------------------------------------------
#' Dew Point Temperature
#'
#' Calculates the dew point temperature as a function of actual vapour pressure,
#' based on McJannet et al. (2008) Equation 26 as presented in McMahon et al.
#' (2013) Equation S2.3
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R. (2013). Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363.
#'   https://doi.org/10.5194/hess-17-1331-2013.
#'
#' @param atmp air temperature (degrees C). When dt is "daily" or larger,
#'             argument should be a list with elements "min" and "max" for daily
#'             min and max temperatures. When dt is "hourly", argument should be
#'             a vector with hourly recorded temperatures.
#' @param RH relative humidity (percent). When dt is "daily" or larger, argument
#'           should be a list with elements "min" and "max" for daily min and
#'           max relative humidities. When dt is "hourly", argument should be a
#'           vector with hourly recorded relative humidities.
#'
#' @return \item{dewtmp}{the dew point temperature (degrees C)}
#'
#' @export

tmp_dew <- function(atmp, RH) {
  ea     <- vp_act_mean(atmp, RH)
  dewtmp <- (116.9 + 237.3*log(ea))/(16.78 - log(ea))
  return(dewtmp)
}

# ------------------------------------------------------------------------------
#' Wet Bulb Temperature
#'
#' Calculates the wet bulb temperature as a function of actual vapour pressure,
#' dew point temperature, and air temperature based on McJannet et al. (2008)
#' Equation 25 as presented in McMahon et al. (2013) Equation S2.2
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
#' @param RH relative humidity (percent). When dt is "daily" or larger, argument
#'           should be a list with elements "min" and "max" for daily min and
#'           max relative humidities. When dt is "hourly", argument should be a
#'           vector with hourly recorded relative humidities.
#'
#' @return \item{wbtmp}{the wet bulb temperature (degrees C)}
#'
#' @export

tmp_wet_bulb <- function(atmp, RH) {
  ea     <- vp_act_mean(atmp, RH)
  dewtmp <- dew_tmp(atmp, RH)

  if (class(atmp) == "list"){ atmp <- (atmp$max + atmp$min)/2 }

  wbtmp  <- (0.00066*100*atmp + 4098*ea*dewtmp/(dewtmp + 237.3)^2)/
            (0.00066*100 + 4098*ea/(dewtmp + 237.3)^2)

  return(wbtmp)
}
