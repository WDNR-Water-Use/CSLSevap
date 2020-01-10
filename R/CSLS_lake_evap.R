#' Calculate Lake Evaporation
#'
#' Calculates Lake Evaporation
#'
#' @param
#'
#' @return weather, a list with daily weather information including:
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
#' @importFrom CSLSevap evaporation
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr group_by summarise
#'
#' @export
CSLS_daily_met <- function(lakename, Lz = 90, wind_elev = 3, z0 = 0.02,
                           no_condensation = FALSE) {

  weather       <- CSLSdata::weather
  lst           <- CSLSdata::lst_HOBO[[lakename]]
  lake_levels   <- CSLSdata::lake_levels[[lakename]]
  elev_area_vol <- CSLSdata::elev_area_vol[[lakename]]
  dictionary    <- CSLSdata::dictionary[[lakename]]

  # Determine first day lst, weather, and lake level data
  lst_day0         <- floor_date(min(lst$date), unit = "day")
  weather_day0     <- floor_date(min(weather$date), unit = "day")
  lake_levels_day0 <- floor_date(min(lake_levels$date), unit = "day")
  day0             <- max(c(lst_day0, weather_day0, lake_levels_day0))

  lst         <- lst %>% filter(date >= day0)
  weather     <- weather %>% filter(date > day0)
  lake_levels <- lake_levels %>% filter(date > day0)

  # Location info
  loc <- CSLS_format_loc(dictionary, Lz)

  # Lake information
  daily_lst  <- lst %>%
                group_by(date = floor_date(.data$date, unit = "day")) %>%
                summarise(ltmp = mean(.data$ltmp, na.rm = TRUE)) %>%
                ungroup()
  wtmp0      <- daily_lst$ltmp[daily_lst$date == day0]
  lst        <- daily_lst %>% filter(date > day0)
  lake       <- CSLS_format_lake(elev_area_vol, lake_levels, lst, wtmp0)

  # Weather information
  weather <- CSLS_format_weather(weather, wind_elev, z0)


  # Lake evaporation
  E <- evaporation("McJannet", loc, weather, lake,
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
