# radiation_functions.R
# Includes:
# - R_n (FAO, McJannet, wet_bulb)
# - R_nl
# - R_ns
# - R_ol (McJannet, wet_bulb)
# - R_il
# - R_so
# - R_a
# - inverse_dist
# - declination
# - hour_angle_sunset
# - hour_angle
# - cloud_factor

# ------------------------------------------------------------------------------
#' Net Solar or Shortwave Radiation
#'
#' Calculates the net solar or shortwave radiation using the "FAO" Penman
#' Monteith approach for reference evapotranspiration (Allen et al., 1998,
#' Equation 38), the "McJannet" approach for lake temperature (McMahon et al.,
#' 2013, Equation S11.25), or the McJannet approach for "wet_bulb" temperature
#' (McMahon et al., 2013, Equation S11.31).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R. (2013). Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363.
#'   https://doi.org/10.5194/hess-17-1331-2013.
#'
#' @param type  type of net radiation function to use. Defaults to "FAO" for FAO
#'              Penman-Monteith reference evapotranspiration. Other options
#'              include "McJannet" for net radiation at lake temperature and
#'              "wet_bulb" for net radiation at wet-bulb temperature.
#' @param Gsc solar constant (MJ/m^2/min). Defaults to 0.0820 MJ/m^2/min.
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#' @inheritParams evaporation
#'
#' @return \item{Rn}{net radiation (MJ/m^2/timestep)}
#'
#' @export
R_n <- function(type = "FAO", loc, lake, weather, albedo){
  # Shortwave
  Rns  <- R_ns(weather$Rs, albedo)

  # Longwave
  Ra   <- R_a(weather$dt, weather$datetimes, loc$phi, loc$Lm, loc$Lz)
  Rso  <- R_so(Ra, loc$z)
  if (type == "FAO") {
    # FAO NET RADIATION --------------------------------------------------------
    ea   <- vp_act_mean(weather$atmp, weather$RH)
    Rnl  <- R_nl(weather$dt, weather$Rs, Rso, ea, weather$atmp)
  } else {
    # MCJANNET LAKE EVAPORATION ------------------------------------------------
    Cf   <- cloud_factor(weather$Rs, Rso)
    Ril  <- R_il(Cf, weather$atmp)
    if (type == "McJannet") {
      # At Lake Temperature ----------------------------------------------------
      wtmp <- tmp_water(loc, lake, weather, albedo)
      Rol  <- R_ol(type, wtmp)
    } else if (type == "wet_bulb") {
      # At Wet-Bulb Temperature ------------------------------------------------
      wbtmp <- tmp_wet_bulb(weather$atmp, weather$RH)
      Rol   <- R_ol(type, wbtmp, weather$atmp)
    }
    Rnl  <- Rol - Ril
  }

  # Net
  Rn   <- Rns - Rnl
  return(Rn)
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
R_ns <- function(Rs, albedo = 0.23){
  Rns <- (1 - albedo)*Rs
  return(Rns)
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
R_nl <- function(dt, Rs, Rso, ea, atmp, SBc = 4.903e-9){
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
R_so <- function(Ra, z){
  Rso   <- (0.75 + 2e-5*z)*Ra
  return(Rso)
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
R_a <- function(dt, datetimes, phi, Lm, Lz, Gsc = 0.0820){
  # Day and time
  if (dt == "hourly") {
    for (i in 1:(length(datetimes)-1)){
      tL[i]   <- int_length(datetimes[i] %--% datetimes[i+1])
      tmid[i] <- datetimes[i] + seconds(tL[i]/2)
    }
    tmid      <- as_datetime(tmid)
  } else {
    tmid <- datetimes
  }
  J      <- yday(tmid)

  # Sun position stuff
  delta  <- declination(J)
  dr     <- inverse_dist(J)
  C      <- (24*60/pi)*Gsc*dr

  if (dt == "hour") {
    # HOURLY, NIGHT-ADJUSTED ---------------------------------------------------
    for (i in 1:(length(datetimes)-1)){
      # Solar time angles
      omega_sunset <- hour_angle_sunset(phi, delta[i])
      omega        <- hour_angle(tmid[i], tL[i], Lm, Lz = 90)
      omega_mid    <- (omega$omega1 + omega$omega2)/2

      # Check if night or just before sunset, adjust omega values if needed
      if ((omega_mid >= omega_sunset-0.79) & (omega_mid <= omega_sunset-0.52)) {
        tL_night   <- tL[i]
        tmid_night <- tmid[i]
      } else if (omega_mid > omega_sunset | omega_mid < -omega_sunset) {
        omega      <- hour_angle(tmid_night, tL_night, Lm, Lz)
      }
      omega1       <- omega$omega1
      omega2       <- omega$omega2

      # Calculate Ra, with nightime substitutions made for nighttime values
      Ra[i]        <- C[i]*((omega2-omega1)*sin(phi)*sin(delta[i]) +
                           cos(phi)*cos(delta)*(sin(omega2) - sin(omega1)))
    }
    Ra <- c(Ra, tail(Ra,1))
  } else {
    # DAILY --------------------------------------------------------------------
    omega_s <- hour_angle_sunset(phi, delta)
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
inverse_dist <- function(J) {
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
declination <- function(J) {
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
hour_angle_sunset <- function(phi, delta) {
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
hour_angle <- function(tmid, tL, Lm, Lz) {
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

# ------------------------------------------------------------------------------
#' Outgoing Longwave Radiation at Wet Bulb Temperature
#'
#' Calculates the outgoing longwave radiation for a lake as a function of lake
#' temperature ("McJannet"; McMahon et al., 2013, Equation S11.27) or as a
#' function of wet bulb temperature  ("wet_bulb"; McMahon et al., 2013, Equation
#' S11.32).
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R.: Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363,
#'   https://doi.org/10.5194/hess-17-1331-2013, 2013.
#'
#' @param type type of outgoing longwave radiation function to use. Defaults to
#'             "McJannet" for outgoing longwave radiation at lake temperature.
#'             Set to "wet_bulb" for outgoing longwave radiation at wet-bulb
#'             temperature.
#' @param tmp1 daily lake water temperature (if "McJannet") or wet bulb
#'             temperature (if "wet_bulb") (degrees C).
#' @param tmp2 defaults to NULL. If "wet_bulb", set to air temperature (degrees
#'             C). When dt is "daily" or larger, argument should be a list with
#'             elements "min" and "max" for daily min and max temperatures. When
#'             dt is "hourly", argument should be a vector with hourly recorded
#'             temperatures.
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#'
#' @return {Rol}{outgoing longwave radiation (MJ/m^2/day)}
#'
#' @export
R_ol <- function(type = "McJannet", tmp1, tmp2 = NULL, SBc = 4.903e-9) {
  if (type == "McJannet") {
    # OUTGOING RADIATION AT LAKE TEMPERATURE -----------------------------------
    Rol <- 0.97*SBc*(tmp1 + 273.15)^4
  } else if (type == "wet_bulb") {
    # OUTGOING RADIATION AT WET-BULB TEMPERATURE -------------------------------
    if (class(tmp2) == "list") { tmp2 <- (tmp2$max + tmp2$min)/2 }
    Rol <- SBc*(tmp2 + 273.15)^4 + 4*SBc*(tmp2 + 273.15)^3*(tmp1 - tmp2)
  }
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
R_il <- function(Cf, atmp, SBc = 4.903e-9) {
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
cloud_factor <- function(Rs, Rso) {
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
