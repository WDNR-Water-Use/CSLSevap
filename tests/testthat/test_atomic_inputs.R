context("Atomic Inputs - FAO Examples")
library(CSLSevap)
library(lubridate)
library(NISTunits)

test_that("FAO Example #18 - Daily ET", {
  method  <- "FAO"
  loc     <- list(z = 100,
                  phi = NISTdegTOradian(50 + 48/60),
                  Lm = NULL,
                  Lz = NULL)
  weather <- list(dt = "daily",
                  datetimes = mdy("7/6/19"),
                  atmp = list(min = 12.3, max = 21.5),
                  RH   = list(min = 63, max = 84),
                  Rs   = 22.07,
                  wind = 2.78,
                  wind_elev = 10)
  albedo  <- list(ref_crop = 0.23, water = 0.08)

  u2    <- u2_fcn(weather$wind, weather$wind_elev)
  gamma <- psychrometric_constant(loc$z)
  Delta <- vp_sat_curve_slope(weather$atmp)
  ea    <- vp_act_mean(weather$atmp, weather$RH)
  vpd   <- vp_deficit(weather$atmp, weather$RH)
  Ra    <- R_a(weather$dt, weather$datetimes, loc$phi, loc$Lm, oc$Lz)
  Rso   <- R_so(Ra, loc$z)
  Rnl   <- R_nl(weather$dt, weather$Rs, Rso, ea, weather$atmp)
  Rn    <- R_n(method, loc, lake, weather, albedo$ref_crop)
  G     <- heat_flux(method, loc, lake, weather, Rn)
  ET    <- evaporation(method, loc, weather)

  expect_equal(u2, 2.079304, tolerance=1e-5)
  expect_equal(gamma, 0.06655628, tolerance=1e-5)
  expect_equal(Delta, 0.1221127, tolerance=1e-5)
  expect_equal(vpd, 0.5888618, tolerance=1e-5)
  expect_equal(ea, 1.408624, tolerance=1e-5)
  expect_equal(Ra, 41.08838, tolerance=1e-5)
  expect_equal(Rso, 30.89846, tolerance=1e-5)
  expect_equal(Rnl, 3.711241, tolerance=1e-5)
  expect_equal(Rn, 13.28266, tolerance=1e-5)
  expect_equal(G, 0)
  expect_equal(ET, 3.880142, tolerance=1e-5)
})

test_that("FAO Example #17 - Monthly ET", {
  method  <- "FAO"
  loc     <- list(z = 2,
                  phi = NISTdegTOradian(13 + 44/60),
                  Lm = NULL,
                  Lz = NULL)
  weather <- list(dt = "monthly",
                  datetimes = mdy(c("3/15/19","4/15/19")),
                  atmp = list(min = c(29.2, 25.6), max = c(29.2, 34.8)),
                  RH   = list(min = c(NA, 51), max = c(NA, 87)),
                  Rs   = c(1, 22.65),
                  wind = c(NA, 2),
                  wind_elev = 2)
  albedo  <- list(ref_crop = 0.23, water = 0.08)

  u2    <- u2_fcn(weather$wind, weather$wind_elev)
  gamma <- psychrometric_constant(loc$z)
  Delta <- vp_sat_curve_slope(weather$atmp)
  es    <- vp_sat_mean(weather$atmp)
  ea    <- vp_act_mean(weather$atmp, weather$RH)
  vpd   <- vp_deficit(weather$atmp, weather$RH)
  Ra    <- R_a(weather$dt, weather$datetimes, loc$phi, loc$Lm, oc$Lz)
  Rso   <- R_so(Ra, loc$z)
  Rnl   <- R_nl(weather$dt, weather$Rs, Rso, ea, weather$atmp)
  Rn    <- R_n(method, loc, lake, weather, albedo$ref_crop)
  G     <- heat_flux(method, loc, lake, weather, Rn)
  ET    <- evaporation(method, loc, weather)

  expect_equal(gamma, 0.06732263, tolerance=1e-5)
  expect_equal(Delta[2], 0.2458003, tolerance=1e-5)
  expect_equal(es[2], 4.421798, tolerance=1e-5)
  expect_equal(ea[2], 2.846016, tolerance=1e-5)
  expect_equal(vpd[2], 1.575782, tolerance=1e-5)
  expect_equal(Ra[2], 38.05769, tolerance=1e-5)
  expect_equal(Rso[2], 28.54479, tolerance=1e-5)
  expect_equal(Rnl[2], 3.11295, tolerance=1e-5)
  expect_equal(Rn[2], 14.32755, tolerance=1e-5)
  expect_equal(G[2], 0.14, tolerance=1e-5)
  expect_equal(ET[2], 5.718284, tolerance=1e-5)
})

test_that("McMahon S19 - McJannet Lake Evap", {
  method  <- "McJannet"
  loc     <- list(z = 546,
                  phi = NISTdegTOradian(-23.7951),
                  Lm = NULL,
                  Lz = NULL)
  lake    <- list(A = 5,
                  depth_m = 10)
  weather <- list(dt = "daily",
                  datetimes = mdy("7/20/80"),
                  atmp = list(min = 2, max = 21),
                  wtmp0 = 10.8734,
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
  Rnwb    <- R_n("wet_bulb", loc, lake, weather, albedo$water)
  Ctime   <- time_const(loc, lake, weather)
  eqtmp   <- tmp_equil(loc, lake, weather, albedo)
  wtmp    <- tmp_water(loc, lake, weather, albedo)
  Gw      <- heat_flux(method, loc, lake, weather, albedo)
  Rol     <- R_ol(method, wtmp)
  Rn      <- R_n(method, loc, lake, weather, albedo$water)
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
