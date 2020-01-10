#' Given CSLS data, format "lake" input parameter
#'
#' Extracts location information for the lake from the site dictionary and
#' reformats for input into lake evaporation functions.
#'
#' @param dictionary a data frame with the elev_m, lat_deg, and long_deg of
#'                   measurement sites, as in the dictionary dataset,
#'                   subset for a single lake.
#' @param Lz the longitude of the local timezone (degrees west of Greenwich,
#'           ranges from 0 to 360 degrees). Defaults to 90 for Central Time
#'           Zone, USA.
#'
#' @return loc, a list with the following location-specific parameters:
#' \item{A}{surface area of the lake (km^2)}
#' \item{depth_m}{depth of the lake (m). Can be a static value or vector
#'               corresponding with datetimes vector.}
#' \item{lst}{optional data frame with date (datetime) and ltmp (lake
#'            temperature, degC).}
#' \item{wtmp0}{required initial water temperature for first day in datetimes
#'              (degC)}
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise
#' @importFrom rlang .data
#' @importFrom NISTunits NISTdegTOradian
#'
#' @export
CSLS_format_lake <- function(elev_area_vol, lake_levels, lst = NULL, wtmp0) {
  lake <- list(A = 1e-6*(lake_levels$area_m2),
               depth_m = lake_levels$level_m - min(elev_area_vol$elev_m),
               lst = lst,
               wtmp0 = wtmp0)
  return(lake)
}
