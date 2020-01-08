context("McJannet Lake Evaporation")
library(CSLSevap)
library(lubridate)
library(NISTunits)

test_that("McMahon Example S19 - McJannet Lake Evap", {
  method  <- "McJannet"
  loc     <- list(z = 546,
                  phi = NISTdegTOradian(-23.7951),
                  Lm = NULL,
                  Lz = NULL)
  lake    <- list(A = 5,
                  depth_m = 10,
                  wtmp0 = 10.8734)
  weather <- list(dt = "daily",
                  datetimes = mdy("7/20/80"),
                  atmp = list(min = 2, max = 21),
                  RH   = list(min = 25, max = 71),
                  Rs   = 17.1940,
                  wind = 0.5903,
                  wind_elev = 2,
                  z0 = 0.0002)
  albedo  <- list(ref_crop = 0.23, water = 0.08)

  u10     <- uz_to_u10(weather$wind, weather$wind_elev, weather$z0)
  gamma   <- psychrometric_constant(loc$z)
  ea      <- vp_act_mean(weather$atmp, weather$RH)
  dewtmp  <- tmp_dew(weather$atmp, weather$RH)
  wbtmp   <- tmp_wet_bulb(weather$atmp, weather$RH)
  Deltawb <- vp_sat_curve_slope(wbtmp)
  ufcn    <- u_fcn(u10, lake$A)
  ra      <- aero_resist(weather$wind, weather$wind_elev, weather$z0, lake$A,
                         loc$z)
  Rns     <- R_ns(weather$Rs, albedo$water)
  Ra      <- R_a(weather$dt, weather$datetimes, loc$phi, loc$Lm, loc$Lz)
  Rso     <- R_so(Ra, loc$z)
  Cf      <- cloud_factor(weather$Rs, Rso)
  Ril     <- R_il(Cf, weather$atmp)
  Rolwb   <- R_ol("wet_bulb", wbtmp, weather$atmp)
  Rnwb    <- R_n("wet_bulb", loc, lake, weather, albedo)
  Ctime   <- time_const(loc, lake, weather)
  eqtmp   <- tmp_equil(loc, lake, weather, albedo)
  wtmp    <- tmp_water(loc, lake, weather, albedo)
  Gw      <- heat_flux(method, loc, lake, weather, albedo)
  Rol     <- R_ol(method, wtmp)
  Rn      <- R_n(method, loc, lake, weather, albedo)
  Delta_w <- vp_sat_curve_slope(wtmp)
  es_w    <- vp_sat_mean(wtmp)
  ET      <- evaporation(method, loc, weather, lake)

  expect_equal(u10, 0.6934505, tolerance=1e-5)
  expect_equal(dewtmp, -1.158446, tolerance=1e-5)
  expect_equal(wbtmp, 6.630961, tolerance=1e-5)
  expect_equal(ufcn, 4.888717, tolerance=1e-5)
  expect_equal(ra, 340.1621, tolerance=1e-5)
  expect_equal(Rns, 15.81848, tolerance=1e-5)
  expect_equal(Cf, 0.08653416, tolerance=1e-5)
  expect_equal(Ril, 25.26406, tolerance=1e-5)
  expect_equal(Rolwb, 29.98653, tolerance=1e-5)
  expect_equal(Rnwb, 11.09601, tolerance=1e-5)
  expect_equal(Ctime, 39.18152, tolerance=1e-5)
  expect_equal(eqtmp, 17.0289, tolerance=1e-5)
  expect_equal(wtmp, 11.02851, tolerance=1e-5)
  expect_equal(Gw, 6.485639, tolerance=1e-5)
  expect_equal(Rol, 31.01691, tolerance=1e-5)
  expect_equal(Rn, 10.06563, tolerance=1e-5)
  expect_equal(Delta_w, 0.08740011, tolerance=1e-5)
  expect_equal(es_w, 1.315204, tolerance=1e-5)
  expect_equal(ET, 1.465042, tolerance=1e-5)
})
