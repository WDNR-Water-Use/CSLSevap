#' Given CSLS data, format "lake" input parameter
#'
#' Extracts location information for the lake from the site dictionary and
#' reformats for input into lake evaporation functions.
#'
#' @param elev_area_vol a data frame with the lake, stage_m, surf_area_m2, and
#'                  volume_m3 as in the elev_area_voldataset, subset
#'                  for a single lake.
#' @param lake_levels a data frame with daily water level measurements as
#'                    formatted in the lake_levels dataset,
#'                    subset to lake level records for the lake of interest.
#' @param lst a data frame with sub-monthly lake surface temperature
#'            measurements as formatted in the lst_HOBO dataset, subset
#'            for a single lake.
#' @param wtmp0 initial water temperature for first day in timeseries (degC)
#' @inheritParams CSLS_daily_met
#'
#' @return **lake**, a list with the following lake-specific parameters:
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
CSLS_format_lake <- function(elev_area_vol, lake_levels, lst = NULL, wtmp0, use_lst) {
  lake <- list(A = 1e-6*(lake_levels$area_m2),
               depth_m = lake_levels$level_m - min(elev_area_vol$elev_m),
               lst = lst,
               wtmp0 = wtmp0)
  if (!use_lst){
    lake$lst <- NULL
  }
  return(lake)
}
