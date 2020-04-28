#' Given CSLS lake and method, calculate daily evaporation
#'
#' This function calculates daily evaporation for a CSLS lake given the desired
#' method and lake name. While this function can use any of the methods
#' availaible in \code{\link{evaporation}}, it is currently hardwired to load
#' input data from the CSLS field campaign (as available in \pkg{CSLSdata}) and
#' analysis is restricted to dates with available weather, lake surface
#' temperature, and lake level data (via \code{\link{CSLS_find_common_dates}})
#' regardless of the requirements of the chosen method.
#'
#' @param method denotes which evaporation method to use ("FAO", "McJannet", or
#'               "Hamon").
#' @param use_lst logical defaults to TRUE to use available lake surface
#'                temperature data.
#' @param Lz longitude of location's measurement timezone (degrees west of
#'           Greenwich). For example, Lz = 75, 90, 105 and 120° for measurement
#'           times based on the Eastern, Central, Rocky Mountain and Pacific
#'           time zones (United States).
#' @param wind_elev height at which wind is measured (m), default: 3
#' @param z0 aerodynamic roughness of weather measurement site (m), default: 0.2
#' @param no_condensation defaults to FALSE. If TRUE, forces negative
#'                        evapotranspiration values (i.e., condensation) to zero
#'
#' @return **weather**, a data frame with daily weather information including:
#' \item{date}{day of each weather observation [POSIXct]}
#' \item{atmp_min}{minimum air temperature for the day (deg C)}
#' \item{atmp_max}{maximum air temperature for the day (deg C)}
#' \item{RH_min}{minimum relative humidity for the day (percent)}
#' \item{RH_max}{maximum relative humidity for the day (percent)}
#' \item{P}{total precipitation for the day (mm)}
#' \item{E}{total lake evaporation for the day (mm)}
#'
#' @examples
#' daily_met <- CSLS_daily_met("McJannet")
#' daily_met <- CSLS_daily_met("McJannet", use_lst = FALSE)
#' daily_met <- CSLS_daily_met("Hamon")
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @import dplyr
#'
#' @export
CSLS_daily_met <- function(method = "McJannet",
                           use_lst = TRUE,
                           Lz = 90,
                           wind_elev = 3,
                           z0 = 0.02,
                           no_condensation = FALSE) {

  weather       <- CSLSdata::weather
  water_chem    <- CSLSdata::water_chem
  water_chem    <- water_chem %>%
                   filter(.data$description == "TEMPERATURE HOBO") %>%
                   mutate(result = as.numeric(.data$result))  %>%
                   group_by(lake = .data$lake,
                            date = floor_date(.data$date, unit = "day"),
                            depth_m = .data$depth1_m) %>%
                   summarise(ltmp = mean(.data$result, na.rm = TRUE)) %>%
                   filter(!is.nan(.data$ltmp)) %>%
                   ungroup()
  lst           <- water_chem %>%
                   group_by(lake = .data$lake,
                            date = .data$date) %>%
                   mutate(min_depth = min(.data$depth_m)) %>%
                   filter(.data$depth_m == .data$min_depth) %>%
                   ungroup() %>%
                   select(.data$lake, .data$date, .data$ltmp)
  lake_levels   <- CSLSdata::lake_levels
  elev_area_vol <- CSLSdata::elev_area_vol
  dictionary    <- CSLSdata::dictionary

  # Subset to common interval
  inputs <- CSLS_find_common_dates(weather, lake_levels, lst)

  big_df <- NULL
  for (this_lake in unique(water_chem$lake)) {
    # Reformat inputs
    loc     <- CSLS_format_loc(filter(dictionary, .data$lake == this_lake), Lz)
    lake    <- CSLS_format_lake(filter(elev_area_vol, .data$lake == this_lake),
                                filter(inputs$lake_levels, .data$lake == this_lake),
                                filter(inputs$lst, .data$lake == this_lake),
                                filter(inputs$wtmp0, .data$lake == this_lake)$ltmp,
                                use_lst)
    weather <- CSLS_format_weather(inputs$weather, wind_elev, z0)

    # Lake evaporation
    E <- evaporation(method, loc, weather, lake,
                     no_condensation = no_condensation)

    # Combine daily weather values into new weather dataframe
    df <- as.data.frame(cbind(weather$datetimes,
                              weather$atmp$min,
                              weather$atmp$max,
                              weather$RH$min,
                              weather$RH$max,
                              weather$P,
                              E))
    colnames(df) <- c("date",
                           "atmp_min",
                           "atmp_max",
                           "RH_min",
                           "RH_max",
                           "P",
                           "E")
    df$date <- as_datetime(df$date)
    df$lake <- this_lake

    big_df <- rbind(big_df, df)
  }

  return(big_df)
}
