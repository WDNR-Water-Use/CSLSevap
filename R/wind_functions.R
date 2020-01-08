#wind_functions.R
# Includes:
# - u2_fcn
# - uz_to_u10
# - u_fcn
# - aero_resist

# ------------------------------------------------------------------------------
#' Wind speed at 2m height
#'
#' Calculates the wind speed at 2 m height given a wind speed at z m height
#' based on Equation 47 in Allen et al. (1998).
#'
#' @references Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop
#'   evapotranspiration: Guidelines for computing crop water requirements. Rome:
#'   FAO. Retrieved from http://www.fao.org/docrep/X0490E/x0490e00.htm.
#'
#' @param uz wind speed at z m above ground surface (m/s)
#' @param z height of measurement above ground surface (m). Default is 10m, a
#'          common height for wind speed measurements in meteorology.
#'
#' @return \item{u2}{wind speed at 2 m above ground surface (m/s)}
#'
#' @export

u2_fcn <- function(uz, z = 10) {
  u2   <- uz*4.87/log(67.8*z - 5.42)
  return(u2)
}

# ------------------------------------------------------------------------------
#' Wind speed at 10m height
#'
#' Calculates the wind speed at 10 m height given a wind speed at z m height based
#' on Equation S4.4 in McMahon et al. (2013).
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R. (2013). Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363.
#'   https://doi.org/10.5194/hess-17-1331-2013.
#'
#' @param uz wind speed at z m above ground surface (m/s)
#' @param z height of measurement above ground surface (m). Default is 2m, a
#'          common height for wind speed measurements in agronomy
#' @param z0 roughness length for the measurement surface (m). Defaults to 0.02
#'           for a short grass.
#'
#' @return \item{u10}{wind speed at 10 m above ground surface (m/s)}
#'
#' @export

uz_to_u10 <- function(uz, z = 2, z0 = 0.02) {
  u10 <- uz*log(10/z0)/log(z/z0)
  return(u10)
}

# ------------------------------------------------------------------------------
#' Wind function for lakes
#'
#' Calculates the wind function given the wind speed at 10m and the surface area
#' of the lake, based on McJannet et al. (2008) Appendix B, Equation 10, as
#' presented in Equation S11.24 in McMahon et al. (2013).
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R. (2013). Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363.
#'   https://doi.org/10.5194/hess-17-1331-2013.
#'
#' @references McJannet, D. L., Webster, I. T., Stenson, M. P., and Sherman,
#'   B.S. (2008). Estimating open water evaporation for the Murray-Darling
#'   Basin. A report to the Australian Government from the CSIRO Murray-Darling
#'   Basin Sustainable Yields Project, CSIRO, Australia, 50 pp. Retrieved from
#'   http://www.clw.csiro.au/publications/waterforahealthycountry/mdbsy/technical/U-OpenWaterEvaporation.pdf.
#'
#' @param u10 wind speed at 10 m above ground surface (m/s)
#' @param A area of the water body (km^2)
#'
#' @return \item{ufcn}{wind function from McJannet et al., 2008}
#'
#' @export

u_fcn <- function(u10, A) {
  ufcn <- ((5/A)^0.05)*(3.80 + 1.57*u10)
  return(ufcn)
}

# ------------------------------------------------------------------------------
#' Aerodynamic resistance
#'
#' Calculates the aerodynamic resistance (s/m) over a lake as defined by Calder
#' and Neal (1984, pg. 93) and McJannet et al. (2008, Appendix B, Equation 10),
#' and as presented in McMahon et al. (2013) Equation S11.23.
#'
#' @references McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R., and
#'   McVicar, T. R. (2013). Estimating actual, potential, reference crop and pan
#'   evaporation using standard meteorological data: a pragmatic synthesis,
#'   Hydrol. Earth Syst. Sci., 17, 1331–1363.
#'   https://doi.org/10.5194/hess-17-1331-2013.
#'
#' @references McJannet, D. L., Webster, I. T., Stenson, M. P., and Sherman,
#'   B.S. (2008). Estimating open water evaporation for the Murray-Darling
#'   Basin. A report to the Australian Government from the CSIRO Murray-Darling
#'   Basin Sustainable Yields Project, CSIRO, Australia, 50 pp. Retrieved from
#'   http://www.clw.csiro.au/publications/waterforahealthycountry/mdbsy/technical/U-OpenWaterEvaporation.pdf.
#'
#' @param uz wind speed at wind_z height (m/s)
#' @param wind_z height at which uz is measured (m)
#' @param z0 aerodynamic roughness of land cover at measurement site (m)
#' @param A surface area of the lake (km^2)
#' @param lake_z elevation of lake above mean sea level (m)
#' @param rho_a density of the air (kg/m^3), defaults to 1.20 kg/m^3 at 20 deg C
#' @param ca specific heat of the air (MJ/kg/K), defaults to 0.001013 MJ/kg/K
#'
#' @return \item{ra}{aerodynamic resistance (s/m)}
#'
#' @export

aero_resist <- function(uz, wind_z, z0, A, lake_z, rho_a = 1.20, ca = 0.001013){
  u10   <- uz_to_u10(uz, wind_z, z0)
  ufcn  <- u_fcn(u10, A)
  gamma <- psychrometric_constant(lake_z)
  ra    <- rho_a*ca/(gamma*ufcn/(60*60*24))
  return(ra)
}
