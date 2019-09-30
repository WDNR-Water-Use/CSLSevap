#wind_functions.R
# Includes:
# - FAO_u2

# ------------------------------------------------------------------------------
#' FAO Wind Speed at 2m Height
#'
#' Calculates the wind speed at 2m height given a wind speed at zm height based
#' on Equation 47 in Allen et al. (1998).
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

FAO_u2 <- function(uz, z = 10) {
  u2 <- uz*4.87/log(67.8*z - 5.42)
  return(u2)
}
