# vapour_functions.R
# Includes:
# - FAO_eo
# - FAO_mean_es
# - FAO_mean_ea
# - FAO_vpd
# - FAO_slope_es_curve
# - FAO_psychrometric_constant
#
# ------------------------------------------------------------------------------
#' FAO Saturation Vapour Pressure
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

FAO_eo <- function(tmp) {
  eo <- 0.6108*exp(17.27*tmp/(237.3 + tmp))
  return(eo)
}

# ------------------------------------------------------------------------------
#' FAO Mean Saturation Vapour Pressure - Daily Timestep or Larger
#'
#' Calculates the mean saturation vapour pressure for a time period (daily
#' timestep or larger) using the mean daily minimum temperature and the mean
#' daily maximum temperature during that time period. Based on equation 12 of
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
#'
#' @return \item{es}{mean saturation vapour pressure (kPa) during time period of
#'                   interest}
#'
#' @export

FAO_mean_es <- function(atmp) {
  if (class(atmp) == "list"){
    es <- (FAO_eo(atmp$max) + FAO_eo(atmp$min))/2
  } else {
    es <- FAO_eo(atmp)
  }
  return(es)
}

# ------------------------------------------------------------------------------
#' FAO Mean Actual Vapour Pressure - Daily Timestep or Larger
#'
#' Calculates the mean saturation vapour pressure for an hourly, daily, or
#' larger time period based on Equation 17 (for daily or larger timestep) or pg. of Allen et al. (1998).
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

FAO_mean_ea <- function(atmp, RH) {
  if (class(atmp) == "list") {
    ea <- (FAO_eo(atmp$min)*RH$max/100 + FAO_eo(atmp$max)*RH$min/100)/2
  } else {
    ea <- FAO_eo(atmp)*RH/100
  }
  return(ea)
}

# ------------------------------------------------------------------------------
#' FAO Vapour Pressure Deficit
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

FAO_vpd <- function(atmp, RH) {
  ea  <- FAO_mean_ea(atmp, RH)
  es  <- FAO_mean_es(atmp)
  vpd <- es - ea
  return(vpd)
}

# ------------------------------------------------------------------------------
#' FAO Slope of Saturation Vapour Pressure Curve
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

FAO_slope_es_curve <- function(atmp) {
  if (class(atmp) == "list"){
    atmp <- (atmp$min + atmp$max)/2
  }
  Delta <- (4098*(FAO_eo(atmp)))/
           ((atmp + 237.3)^2)
  return(Delta)
}

# ------------------------------------------------------------------------------
#' FAO Psychrometric Constant
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

FAO_psychrometric_constant <- function(z, atmp = NULL, lambda = 2.45,
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
#' @return \item{dewtmp}{the dew point temperature (degrees C)}
#'
#' @export

McJannet_dewtmp <- function(atmp, RH) {
  ea     <- FAO_mean_ea(atmp, RH)
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

McJannet_wbtmp <- function(atmp, RH) {
  ea     <- FAO_mean_ea(atmp, RH)
  dewtmp <- McJannet_dewtmp(atmp, RH)

  if (class(atmp) == "list"){ atmp <- (atmp$max + atmp$min)/2 }

  wbtmp  <- (0.00066*100*atmp + 4098*ea*dewtmp/(dewtmp + 237.3)^2)/
            (0.00066*100 + 4098*ea/(dewtmp + 237.3)^2)

  return(wbtmp)
}
