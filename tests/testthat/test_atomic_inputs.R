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
  Delta <- FAO_slope_es_curve(weather$dt, weather$atmp)
  ea    <- FAO_mean_ea(weather$dt, weather$atmp, weather$RH)
  vpd   <- FAO_vpd(weather$dt, weather$atmp, weather$RH)
  Ra    <- FAO_Ra(weather$dt, weather$datetimes, loc$phi, loc$Lm, oc$Lz)
  Rso   <- FAO_Rso(Ra, loc$z)
  Rnl   <- FAO_Rnl(weather$dt, weather$Rs, Rso, ea, weather$atmp)
  Rn    <- FAO_Rn(loc, weather)
  G     <- FAO_G(weather$dt, weather$datetimes, loc$phi, loc$Lm, loc$Lz, Rn,
                 weather$atmp)
  ET    <- FAO_ET(loc, weather)

  expect_equal(u2, 2.079304, tolerance=1e-5)
  expect_equal(gamma, 0.07235778, tolerance=1e-5)
  expect_equal(Delta, 0.1221127, tolerance=1e-5)
  expect_equal(vpd, 0.5888618, tolerance=1e-5)
  expect_equal(ea, 1.408624, tolerance=1e-5)
  expect_equal(Ra, 41.08838, tolerance=1e-5)
  expect_equal(Rso, 30.89846, tolerance=1e-5)
  expect_equal(Rnl, 3.711241, tolerance=1e-5)
  expect_equal(Rn, 13.28266, tolerance=1e-5)
  expect_equal(G, 0)
  expect_equal(ET, 3.813441, tolerance=1e-5)
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
  Delta <- FAO_slope_es_curve(weather$dt, weather$atmp)
  es    <- FAO_mean_es(weather$dt, weather$atmp)
  ea    <- FAO_mean_ea(weather$dt, weather$atmp, weather$RH)
  vpd   <- FAO_vpd(weather$dt, weather$atmp, weather$RH)
  Ra    <- FAO_Ra(weather$dt, weather$datetimes, loc$phi, loc$Lm, oc$Lz)
  Rso   <- FAO_Rso(Ra, loc$z)
  Rnl   <- FAO_Rnl(weather$dt, weather$Rs, Rso, ea, weather$atmp)
  Rn    <- FAO_Rn(loc, weather)
  G     <- FAO_G(weather$dt, weather$datetimes, loc$phi, loc$Lm, loc$Lz, Rn, weather$atmp)
  ET    <- FAO_ET(loc, weather)

  expect_equal(gamma, 0.07319093, tolerance=1e-5)
  expect_equal(Delta[2], 0.2458003, tolerance=1e-5)
  expect_equal(es[2], 4.421798, tolerance=1e-5)
  expect_equal(ea[2], 2.846016, tolerance=1e-5)
  expect_equal(vpd[2], 1.575782, tolerance=1e-5)
  expect_equal(Ra[2], 38.05769, tolerance=1e-5)
  expect_equal(Rso[2], 28.54479, tolerance=1e-5)
  expect_equal(Rnl[2], 3.11295, tolerance=1e-5)
  expect_equal(Rn[2], 14.32755, tolerance=1e-5)
  expect_equal(G[2], 0.14, tolerance=1e-5)
  expect_equal(ET[2], 5.714204, tolerance=1e-5)
})
