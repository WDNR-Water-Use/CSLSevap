#soil_heat_functions.R
# Includes:
# - FAO_G_hour
# - FAO_G_day
# - FAO_G_month

# ------------------------------------------------------------------------------
#' Soil Heat Flux - Hourly
#'
#' Calculates soil heat flux for hourly, daily, or monthly calculations beneath a
#' dense cover of grass. Based on Equations 42-46 in Allen et al. (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param dt string indicating the timestep of input weather series. Expects
#'            "hourly", "daily", or "monthly".
#' @param datetimes datetimes of weather records [POSIXct]. If monthly
#'                   timestep, make sure date is the 15th of each month.
#' @param phi latitude of the location (radians). Positive for northern
#'            hemisphere, negative for southern hemisphere.
#' @param Lm longitude of the measurement site (degrees west of Greenwich)
#' @param Lz longitude of the center of the local time zone (degrees west of
#'            Greenwich). For example, Lz = 75, 90, 105 and 120° for the
#'            Eastern, Central, Rocky Mountain and Pacific time zones (United
#'            States) and Lz = 0° for Greenwich, 330° for Cairo (Egypt), and
#'            255° for Bangkok (Thailand).
#' @param Rn vector with net radiation (MJ/m^2/timestep)
#' @param atmp If hourly timestep, vector of air temperature (degrees C)
#'            corresponding with datetimes vector. If daily or monthly timestep,
#'            list with two vectors, "min" and "max", with mean daily min and
#'            max air temperature (degrees C) corresponding with datetimes
#'            vector
#'
#' @import lubridate
#' @importFrom utils tail
#'
#' @return \item{G}{soil heat flux (MJ/m^2/hr)}
#'
#' @export
FAO_G <- function(dt, datetimes, phi, Lm, Lz, Rn, atmp) {
  if (dt == "hourly"){
    for (i in 1:length(datetimes)-1){
      tL   <- int_length(datetimes[i] %--% datetimes[i+1])
      tmid <- datetimes[i] + seconds(tL/2)
      J    <- yday(tmid)

      delta        <- FAO_declination(J)
      omega_sunset <- FAO_hour_angle_sunset(phi, delta)
      omega        <- FAO_hour_angle(tmid, tL, Lm, Lz)
      omega_mid    <- (omega$omega1 + omega$omega2)/2
      if (omega_mid > omega_sunset | omega_mid < -omega_sunset) {
        night[i]   <- 1
      } else {
        night[i]   <- 0
      }
    }
    night <- c(night, tail(night,1))
    day   <- 1 - night

    G     <- day*0.1*Rn + night*0.5*Rn

  } else if (dt == "daily") {

    G <- 0*Rn

  } else if (dt == "monthly") {
    tmp_C  <- (atmp$min + atmp$max)/2
        G    <- NULL
        G[1] <- NA
    if (length(tmp_C) == 1) {
    } else if (length(tmp_C) == 2) {
        G[2] <- 0.14*(tmp_C[2] - tmp_C[1])
    } else {
      for (i in 2:(length(tmp_C)-1)){
        G[i] <- 0.07*(tmp_C[i+1] - tmp_C[i-1])
      }
        G[i+1] <- 0.14*(tmp_C[i+1] - tmp_C[i])
    }

  }

  return(G)
}
