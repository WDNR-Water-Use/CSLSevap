#' Given CSLS lake and method, calculate daily evaporation
#'
#' This function calculates daily evaporation for a CSLS lake given the desired
#' method and lake name. While this function can use any of the methods
#' availaible in CSLSevap::evaporation, it is currently hardwired to load input
#' data from the CSLS field campaign (as available in CSLSdata) and analysis is
#' restricted to dates with available weather, lake surface temperature, and
#' lake level data regardless of the requirements of the chosen method.
#'
#' @param method denotes which evaporation method to use ("FAO", "McJannet", or
#'               "Hamon").
#' @param use_lst logical defaults to TRUE to use available lake surface
#'                temperature data.
#' @param lakename name of lake currently evaluating
#' @param Lz longitude of location's measurement timezone (degrees west of
#'           Greenwich). For example, Lz = 75, 90, 105 and 120° for measurement
#'           times based on the Eastern, Central, Rocky Mountain and Pacific
#'           time zones (United States).
#' @param wind_elev height at which wind is measured (m), default: 3
#' @param z0 aerodynamic roughness of weather measurement site (m), default: 0.2
#' @param no_condensation defaults to FALSE. If TRUE, forces negative
#'                        evapotranspiration values (i.e., condensation) to zero
#'
#' @return weather, a data frame with daily weather information including:
#' \describe{
#' \item{date}{day of each weather observation}
#' \item{atmp_min}{minimum air temperature for the day (deg C)}
#' \item{atmp_max}{maximum air temperature for the day (deg C)}
#' \item{RH_min}{minimum relative humidity for the day (percent)}
#' \item{RH_max}{maximum relative humidity for the day (percent)}
#' \item{P}{total precipitation for the day (mm)}
#' \item{E}{total lake evaporation for the day (mm)}
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr group_by summarise
#'
#' @export
CSLS_daily_met <- function(method = "McJannet",
                           use_lst = TRUE,
                           lakename,
                           Lz = 90,
                           wind_elev = 3,
                           z0 = 0.02,
                           no_condensation = FALSE) {

  weather       <- CSLSdata::weather
  lst           <- CSLSdata::lst_HOBO[[lakename]]
  lake_levels   <- CSLSdata::lake_levels[[lakename]]
  elev_area_vol <- CSLSdata::elev_area_vol[[lakename]]
  dictionary    <- CSLSdata::dictionary[[lakename]]

  # Subset to common interval
  inputs <- CSLS_find_common_dates(weather, lake_levels, lst)

  # Reformat inputs
  loc     <- CSLS_format_loc(dictionary, Lz)
  lake    <- CSLS_format_lake(elev_area_vol, inputs$lake_levels, inputs$lst,
                              inputs$wtmp0, use_lst)
  weather <- CSLS_format_weather(inputs$weather, wind_elev, z0)

  # Lake evaporation
  E <- evaporation(method, loc, weather, lake,
                   no_condensation = no_condensation)

  # Combine daily weather values into new weather dataframe
  weather <- as.data.frame(cbind(weather$datetimes,
                                 weather$atmp$min,
                                 weather$atmp$max,
                                 weather$RH$min,
                                 weather$RH$max,
                                 weather$P,
                                 E))
  colnames(weather) <- c("date",
                         "atmp_min",
                         "atmp_max",
                         "RH_min",
                         "RH_max",
                         "P",
                         "E")
  weather$date <- as_datetime(weather$date)

  return(weather)
}
