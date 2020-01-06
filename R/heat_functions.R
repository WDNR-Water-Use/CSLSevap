#heat_functions.R
# Includes:
# - heat_flux
# - tmp_water
# - tmp_equil
# - time_const

# ------------------------------------------------------------------------------
#' Heat Flux
#'
#' Calculates soil heat flux for hourly, daily, or monthly calculations beneath a
#' dense cover of grass. Based on Equations 42-46 in Allen et al. (1998).
#'
#' Calculates the change in heat storage based on water temperature following
#' McJannet et al. (2008) Equation 31 as presented in McMahon et al. S11.33
#' (2013).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R.: Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363,
#'   https://doi.org/10.5194/hess-17-1331-2013, 2013.
#'
#' @inheritParams evap
#' @param rho_w density of water (kg/m^3), defaults to 997.9 kg/m^3 at 20 deg C
#' @param cw specific heat of water (MJ/kg/K), defaults to 0.00419
#'
#' @import lubridate
#' @importFrom utils tail
#'
#' @return \item{G}{heat flux (MJ/m^2/hr)}
#'
#' @export
heat_flux <- function(method, loc, lake, weather, albedo,
                      Rn = NULL, rho_w = 997.9, cw = 0.00419){
  if (method == "FAO") {
    # FAO SOIL HEAT FLUX -------------------------------------------------------
    dt        <- weather$dt
    datetimes <- weather$datetimes
    atmp      <- weather$atmp
    phi       <- loc$phi
    Lm        <- loc$Lm
    Lz        <- loc$Lz
    if (dt == "hourly"){
      for (i in 1:length(datetimes)-1){
        tL   <- int_length(datetimes[i] %--% datetimes[i+1])
        tmid <- datetimes[i] + seconds(tL/2)
        J    <- yday(tmid)

        delta        <- declination(J)
        omega_sunset <- hour_angle_sunset(phi, delta)
        omega        <- hour_angle(tmid, tL, Lm, Lz)
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
  } else if (method == "McJannet"){ # end FAO soil heat flux
    # MCJANNET LAKE HEAT FLUX ----------------------------------------------------
    wtmp0   <- lake$wtmp0
    wtmp    <- tmp_water(loc, lake, weather, albedo)
    G       <- NULL
    if (length(lake$depth_m) == 1) {
      lake$depth_m <- rep(lake$depth_m, length(wtmp))
    }
    for (i in 1:length(wtmp)){
      G[i]     <- rho_w*cw*lake$depth_m[i]*(wtmp[i] - wtmp0)
      wtmp0    <- wtmp[i]
    }
  }

  return(G)
}

# ------------------------------------------------------------------------------
#' Lake Water Temperature
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
#' @inheritParams lake_evap
#'
#' @return {wtmp}{water temperature (degrees C)}
#'
#' @export
tmp_water <- function(loc, lake, weather, albedo) {
  lst   <- lake$lst
  eqtmp <- tmp_equil(loc, lake, weather, albedo)
  Ctime <- time_const(loc, lake, weather)
  wtmp0 <- lake$wtmp0
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
#' @inheritParams evap
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#'
#' @return {eqtmp}{equilibrium temperature (degrees C)}
#'
#' @export
tmp_equil <- function(loc, lake, weather, albedo, SBc = 4.903e-9) {
  u10     <- uz_to_u10(weather$wind, weather$wind_elev, weather$z0)
  ufcn    <- u_fcn(u10, lake$A)
  gamma   <- psychrometric_constant(loc$z)
  wbtmp   <- tmp_wet_bulb(weather$atmp, weather$RH)
  Deltawb <- vp_sat_curve_slope(wbtmp)
  Rnwb    <- R_n("wet_bulb", loc, lake, weather, albedo$water)

  eqtmp   <- wbtmp + Rnwb/(4*SBc*(wbtmp + 273.15)^3 + ufcn*(Deltawb + gamma))

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
#' @inheritParams lake_evap
#' @param rho_w density of water (kg/m^3), defaults to 997.9 kg/m^3 at 20 deg C
#' @param cw specific heat of water (MJ/kg/K), defaults to 0.00419
#' @param SBc Stefan-Boltzman constant (MJ/m^2/day). Defaults to 4.903e-9
#'            MJ/m^2/day.
#'
#' @return {Ctime}{time constant (day)}
#'
#' @export
time_const <- function(loc, lake, weather, rho_w = 997.9, cw = 0.00419,
                       SBc = 4.903e-9) {
  u10     <- uz_to_u10(weather$wind, weather$wind_elev, weather$z0)
  ufcn    <- u_fcn(u10, lake$A)
  gamma   <- psychrometric_constant(loc$z)
  wbtmp   <- tmp_wet_bulb(weather$atmp, weather$RH)
  Deltawb <- vp_sat_curve_slope(wbtmp)

  Ctime   <- rho_w*cw*lake$depth_m/(4*SBc*(wbtmp + 273.15)^3 + ufcn*(Deltawb + gamma))

  return(Ctime)
}
