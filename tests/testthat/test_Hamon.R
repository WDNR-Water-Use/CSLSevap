context("Hamon Lake Evaporation")
library(CSLSevap)
library(lubridate)
library(NISTunits)

test_that("Harwell 2.3 - Daylight Hours", {
  datetime <- mdy("04-15-2020")
  phi      <- NISTdegTOradian(30)
  J        <- yday(datetime)
  delta    <- declination(J)
  omega_s  <- hour_angle_sunset(phi, delta)
  D        <- daylight_hours(omega_s)

  expect_equal(D, 12.7, tolerance=1e-1)
})
