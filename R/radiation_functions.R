# radiation_functions.R
# Includes:
# - FAO_Rns
# - FAO_Ra
# - FAO_inverse_dist
# - FAO_declination
# - FAO_hour_angle_sunset
# - FAO_hour_angle

# ------------------------------------------------------------------------------
#' Net Solar or Shortwave Radiation
#'
#' Calculates the net solar or shortwave radiation for a given location based on
#' incoming solar radiation timeseries. Based on Equation 38 in Allen et al.
#' (1998).
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
#' @return \item{Rn}{net radiation (MJ/m^2/timestep)}
#'
#' @export
FAO_Rn <- function(loc, weather, albedo = 0.23, Gsc = 0.0820, SBc = 4.903e-9){
  # Net longwave radiation
  if (weather$dt == "hourly") {
    Ra <- FAO_Ra_night_adjusted(weather$datetimes, loc$phi, loc$Lm, loc$Lz, Gsc)
  } else {
    Ra <- FAO_Ra(weather$dt, weather$datetimes, loc$phi, loc$Lm, loc$Lz, Gsc)
  }
  Rso  <- FAO_Rso(Ra, loc$z)
  ea   <- FAO_mean_ea(weather$dt, weather$atmp, weather$RH)
  Rnl  <- FAO_Rnl(weather$dt, weather$Rs, Rso, ea, weather$atmp, SBc)

  # Net shortwave radiation
  Rns  <- FAO_Rns(weather$Rs, albedo)

  # Net radiation
  Rn   <- Rns - Rnl
  return(Rn)
}

# ------------------------------------------------------------------------------
#' Net Longwave Radiation
#'
#' Calculates the net outgoing longwave radiation for a given location based on
#' incoming solar radiation, incoming clear-sky radiation, and air temperature
#' records. Based on Equation 39 and p.74-75 in Allen et al. (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param dt string indicating the timestep of input weather series. Expects
#'            "hourly", "daily", or "monthly".
#' @param Rs incoming solar radiation (MJ/m^2/timestep), vector or atomic number
#' @param Rso clear-sky solar radiation (MJ/m^2/timestep), vector or atomic
#'            number.
#' @param ea actual vapour pressure (kPa), vector or atomic number.
#' @param atmp air temperature (degrees C). When dt is "daily" or larger,
#'             argument should be a list with elements "min" and "max" for daily
#'             min and max temperatures. When dt is "hourly", argument should be
#'             a vector or atomic number with hourly recorded temperatures.
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#'
#' @importFrom NISTunits NISTdegCtOk
#'
#' @return \item{Rnl}{Net outgoing longwave radiation (MJ/m^2/hour or
#'                    MJ/m^2/day, depending on dt)}
#'
#' @export
FAO_Rnl <- function(dt, Rs, Rso, ea, atmp, SBc = 4.903e-9){
  if (dt == "hourly") {
    SBc   <- SBc/24
    temp4 <- SBc*(NISTdegCtOk(atmp)^4)
  } else {
    temp4 <- SBc*(NISTdegCtOk(atmp$min)^4 + NISTdegCtOk(atmp$max)^4)/2
  }
  ratio   <- Rs/Rso
  ratio[ratio > 1] <- 1
  Rnl     <- temp4*(0.34 - 0.14*sqrt(ea))*(1.35*ratio - 0.35)

  return(Rnl)
}

# ------------------------------------------------------------------------------
#' Net Solar or Shortwave Radiation
#'
#' Calculates the net solar or shortwave radiation for a given location based on
#' incoming solar radiation timeseries. Based on Equation 38 in Allen et al.
#' (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param Rs incoming solar radiation (MJ/m^2/timestep)
#' @param albedo albedo or canopy reflection coefficient (-). Defaults to 0.23
#'               for the hypothetical grass reference crop.
#'
#' @return \item{Rns}{Net solar or shortwave radiation (MJ/m^2/timestep)}
#'
#' @export
FAO_Rns <- function(Rs, albedo = 0.23){
  Rns <- (1 - albedo)*Rs
  return(Rns)
}

# ------------------------------------------------------------------------------
#' Clear-Sky Solar Radiation
#'
#' Calculates the clear-sky solar radiation based on Equation 37 in Allen et al.
#' (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param Ra extraterrestrial solar radiation (MJ/m^2/day)
#' @param z station elevation above sea level (m).
#'
#' @return \item{Rso}{clear-sky solar radiation (MJ/m^2/day)}
#'
#' @export
FAO_Rso <- function(Ra, z){
  Rso   <- (0.75 + 2e-5*z)*Ra
  return(Rso)
}

# ------------------------------------------------------------------------------
#' Nighttime Extraterrestrial Solar Radiation
#'
#' Calculates the nighttime extraterrestrial radiation for a given location and
#' timestep. Based on recomendations on p. 75 in Allen et al. (1998) for hourly
#' ET caclculations.
#'
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param datetimes hourly timeseries to use in calculations [POSIXct].
#' @param phi latitude of the location (radians). Positive for northern
#'            hemisphere, negative for southern hemisphere.
#' @param Lm longitude of the measurement site (degrees west of Greenwich)
#' @param Lz longitude of the center of the local time zone (degrees west of
#'   Greenwich). For example, Lz = 75, 90, 105 and 120° for the Eastern,
#'   Central, Rocky Mountain and Pacific time zones (United States) and Lz = 0°
#'   for Greenwich, 330° for Cairo (Egypt), and 255° for Bangkok (Thailand).
#' @param Gsc solar constant (MJ/m^2/min). Defaults to 0.0820 MJ/m^2/min.
#'
#' @return \item{Ra}{extraterrestrial radiation (MJ/m^2/hour) with nighttime Ra
#'                   values adjusted per Allen et al. (1998)}
#'
#' @import lubridate
#'
#' @export
FAO_Ra_night_adjusted <- function(datetimes, phi, Lm, Lz, Gsc = 0.0820){
  tL_night   <- NA
  tmid_night <- NA
  for (i in 1:(length(datetimes)-1)) {
    # Date and times
    tL    <- int_length(datetimes[i] %--% datetimes[i+1])
    tmid  <- datetimes[i] + seconds(tL/2)
    J     <- yday(tmid)

    # Sun position stuff
    dr    <- FAO_inverse_dist(J)
    delta <- FAO_declination(J)

    # Solar time angles
    omega_sunset <- FAO_hour_angle_sunset(phi, delta)
    omega        <- FAO_hour_angle(tmid, tL, Lm, Lz = 90)
    omega_mid    <- (omega$omega1 + omega$omega2)/2

    # Check if night or just before sunset, adjust omega values if needed
    if ((omega_mid >= omega_sunset-0.79) & (omega_mid <= omega_sunset-0.52)) {
      tL_night   <- tL
      tmid_night <- tmid
    } else if (omega_mid > omega_sunset | omega_mid < -omega_sunset) {
      omega      <- FAO_hour_angle(tmid_night, tL_night, Lm, Lz)
    }
    omega1       <- omega$omega1
    omega2       <- omega$omega2

    # Calculate Ra, with nightime substitutions made for nighttime values
    C            <- (24*60/pi)*Gsc*dr
    Ra[i]        <- C*((omega2-omega1)*sin(phi)*sin(delta) +
                         cos(phi)*cos(delta)*(sin(omega2) - sin(omega1)))
  }
  Ra <- c(Ra, tail(Ra,1))
  return(Ra)
}

# ------------------------------------------------------------------------------
#' Extraterrestrial Solar Radiation
#'
#' Calculates the extraterrestrial radiation for a given location and timestep.
#'
#' For daily periods, uses Equation 21 from Allen et al. (1998). For hourly or
#' shorter periods, uses Equation 28 from Allen et al. (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param dt a string indicating the timestep of calculation. If not "hourly",
#'           assumes a daily or larger timestep
#' @param datetimes timeseries of datetimes to use in calculations [POSIXct]
#' @param phi latitude of the location (radians). Positive for northern
#'            hemisphere, negative for southern hemisphere.
#' @param Lm longitude of the measurement site (degrees west of Greenwich).
#'           Defaults to NULL since is not needed for daily timesteps.
#' @param Lz longitude of the center of the local time zone (degrees west of
#'            Greenwich). For example, Lz = 75, 90, 105 and 120° for the
#'            Eastern, Central, Rocky Mountain and Pacific time zones (United
#'            States) and Lz = 0° for Greenwich, 330° for Cairo (Egypt), and
#'            255° for Bangkok (Thailand). Defaults to 90 for Wisconsin (US
#'            Central Time Zone). Defaults to NULL since is not needed for daily
#'            timesteps.
#' @param Gsc solar constant (MJ/m^2/min). Defaults to 0.0820 MJ/m^2/min.
#'
#' @return \item{Ra}{extraterrestrial radiation (MJ/m^2/day) or (MJ/m^2/hour),
#'              depending on length of input time period}
#'
#' @import lubridate
#'
#' @export
FAO_Ra <- function(dt, datetimes, phi, Lm, Lz, Gsc = 0.0820){
  # Day and time
  if (dt == "hourly") {
    for (i in 1:(length(datetimes)-1)){
      tL[i]   <- int_length(datetimes[i] %--% datetimes[i+1])
      tmid[i] <- datetimes[i] + seconds(tL[i]/2)
    }
    tmid   <- as_datetime(tmid)
  } else {
    tmid <- datetimes
  }
  J      <- yday(tmid)

  # Sun position stuff
  delta  <- FAO_declination(J)
  dr     <- FAO_inverse_dist(J)

  # Ra
  C  <- (24*60/pi)*Gsc*dr
  if (dt == "hour") {
    omega   <- FAO_hour_angle(tmid, tL, Lm, Lz = 90)
    omega1  <- omega$omega1
    omega2  <- omega$omega2
    Ra      <- C*((omega2-omega1)*sin(phi)*sin(delta) +
                    cos(phi)*cos(delta)*(sin(omega2) - sin(omega1)))
  } else {
    omega_s <- FAO_hour_angle_sunset(phi, delta)
    Ra      <- C*(omega_s*sin(phi)*sin(delta) + cos(phi)*cos(delta)*sin(omega_s))
  }
  return(Ra)
}

# ------------------------------------------------------------------------------
#' Inverse Relative Distance Earth-Sun
#'
#' Calculates the inverse relative distance Earth-Sun using Equation 23 from
#' Allen et al. (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param J Julian day, i.e., the number of the day in the year between 1 (1
#'          January) and 365 or 366 (31 December).
#'
#' @return \item{dr}{inverse relative distance Earth-Sun (unitless)}
#'
#' @export
FAO_inverse_dist <- function(J) {
  dr <- 1 + 0.033*cos(2*pi*J/365)
  return(dr)
}

# ------------------------------------------------------------------------------
#' Solar Declination
#'
#' Calculates the solar declination (radians) using Equation 24 from Allen et
#' al. (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param J Julian day, i.e., the number of the day in the year between 1 (1
#'          January) and 365 or 366 (31 December).
#'
#' @return \item{delta}{solar declination (radians)}
#'
#' @export
FAO_declination <- function(J) {
  delta <- 0.409*sin(2*pi*J/365 - 1.39)
  return(delta)
}

# ------------------------------------------------------------------------------
#' Sunset Hour Angle
#'
#' Calculates the sunset hour angle (radians) for a given location and day of
#' the year using Equation 25 from Allen et al. (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param phi the latitude of the location (radians). Positive for northern
#'            hemisphere, negative for southern hemisphere.
#' @param delta the solar declination based on the Julian day (radians)
#'
#' @return \item{omega_s}{sunset hour angle (radians)}
#'
#' @export
FAO_hour_angle_sunset <- function(phi, delta) {
  omega_s <- acos(-tan(phi)*tan(delta))
  return(omega_s)
}

# ------------------------------------------------------------------------------
#' Solar Time Angle
#'
#' Calculates the beginning, midpoint, and end solar time angle (radians) for a
#' given location and day of the year using Equations 29-33 from Allen et al.
#' (1998). Used for calculations at an hourly or shorter timestep.
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param tmid date and time of the midpoint of the time period [POSIXct].
#' @param tL length of the time period, from lubridate "interval_len".
#' @param Lm longitude of the measurement site (degrees west of Greenwich)
#' @param Lz longitude of the center of the local time zone (degrees west of
#'           Greenwich). For example, Lz = 75, 90, 105 and 120° for the Eastern,
#'           Central, Rocky Mountain and Pacific time zones (United States) and
#'           Lz = 0° for Greenwich, 330° for Cairo (Egypt), and 255° for Bangkok
#'           (Thailand).
#'
#' @return
#' \item{omega1}{solar time angle at the beginning of the period (radians)}
#' \item{omega2}{solar time angle at the end of the period (radians)}
#'
#' @import lubridate
#'
#' @export
FAO_hour_angle <- function(tmid, tL, Lm, Lz) {
  J    <- yday(tmid)
  tmid <- hour(tmid)
  tL   <- hour(seconds_to_period(tL))


  b         <- 2*pi*(J - 81)/364
  Sc        <- 0.1645*sin(2*b) - 0.1225*cos(b) - 0.025*sin(b)

  omega_mid <- (pi/12)*((tmid + 0.06667*(Lz - Lm) + Sc) - 12)
  omega1    <- omega_mid - pi*tL/24
  omega2    <- omega_mid + pi*tL/24

  return(list(omega1 = omega1, omega2 = omega2))
}
