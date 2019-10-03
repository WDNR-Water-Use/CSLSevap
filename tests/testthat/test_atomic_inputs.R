context("Atomic Inputs - FAO Examples")
library(faoET)
library(lubridate)
library(NISTunits)

test_that("FAO Example #18 - Daily ET", {
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

  u2    <- FAO_u2(weather$wind, weather$wind_elev)
  gamma <- FAO_psychrometric_constant(loc$z)
  Delta <- FAO_slope_es_curve(weather$atmp)
  ea    <- FAO_mean_ea(weather$atmp, weather$RH)
  vpd   <- FAO_vpd(weather$atmp, weather$RH)
  Ra    <- FAO_Ra(weather$dt, weather$datetimes, loc$phi, loc$Lm, oc$Lz)
  Rso   <- FAO_Rso(Ra, loc$z)
  Rnl   <- FAO_Rnl(weather$dt, weather$Rs, Rso, ea, weather$atmp)
  Rn    <- FAO_Rn(loc, weather)
  G     <- FAO_G(weather$dt, weather$datetimes, loc$phi, loc$Lm, loc$Lz, Rn,
                 weather$atmp)
  ET    <- FAO_ET(loc, weather)

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

  gamma <- FAO_psychrometric_constant(loc$z)
  Delta <- FAO_slope_es_curve(weather$atmp)
  es    <- FAO_mean_es(weather$atmp)
  ea    <- FAO_mean_ea(weather$atmp, weather$RH)
  vpd   <- FAO_vpd(weather$atmp, weather$RH)
  Ra    <- FAO_Ra(weather$dt, weather$datetimes, loc$phi, loc$Lm, oc$Lz)
  Rso   <- FAO_Rso(Ra, loc$z)
  Rnl   <- FAO_Rnl(weather$dt, weather$Rs, Rso, ea, weather$atmp)
  Rn    <- FAO_Rn(loc, weather)
  G     <- FAO_G(weather$dt, weather$datetimes, loc$phi, loc$Lm, loc$Lz, Rn, weather$atmp)
  ET    <- FAO_ET(loc, weather)

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
  albedo  <- list(lake = 0.08)
  u10     <- McMahon_u10(weather$wind, weather$wind_elev, weather$z0)
  gamma   <- FAO_psychrometric_constant(loc$z)
  ea      <- FAO_mean_ea(weather$atmp, weather$RH)
  dewtmp  <- McJannet_dewtmp(weather$atmp, weather$RH)
  wbtmp   <- McJannet_wbtmp(weather$atmp, weather$RH)
  Deltawb <- FAO_slope_es_curve(wbtmp)
  u_fcn   <- McJannet_u_fcn(u10, lake$A)
  ra      <- McJannet_aero_resist(weather$wind, weather$wind_elev, weather$z0, lake$A, loc$z)
  Rns     <- FAO_Rns(weather$Rs, albedo$lake)
  Ra      <- FAO_Ra(weather$dt, weather$datetimes, loc$phi, loc$Lm, loc$Lz)
  Rso     <- FAO_Rso(Ra, loc$z)
  Cf      <- McJannet_Cf(weather$Rs, Rso)
  Ril     <- McMahon_Ril(Cf, weather$atmp)
  Rolwb   <- McMahon_Rolwb(weather$atmp, wbtmp)
  Rnwb    <- McMahon_Rnwb(loc, weather, albedo$lake)
  Ctime   <- McJannet_time_const(loc, lake, weather)
  eqtmp   <- McJannet_eqtmp(loc, lake, weather)
  wtmp    <- McJannet_wtmp(loc, lake, weather)
  Gw      <- McJannet_Gw(loc, lake, weather)
  Rol     <- McMahon_Rol(wtmp, SBc = 4.903e-9)
  Rn      <- McMahon_Rn(loc, lake, weather, albedo$lake)
  Delta_w <- FAO_slope_es_curve(wtmp)
  es_w    <- FAO_mean_es(wtmp)
  ET      <- McJannet_lake_ET(loc, lake, weather, albedo)

  expect_equal(u10, 0.6934505, tolerance=1e-5)
  expect_equal(dewtmp, -1.158446, tolerance=1e-5)
  expect_equal(wbtmp, 6.630961, tolerance=1e-5)
  expect_equal(u_fcn, 4.888717, tolerance=1e-5)
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
