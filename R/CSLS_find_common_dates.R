#' Given CSLS data, find common dates of record
#'
#' Finds date interval with records for lst, weather, and lake levels.
#'
#' @param weather a data frame with hourly weather data incl. air temperature
#'                (atmp), relative humidity (RH), incoming solar radiation (Rs),
#'                precipitation (P), and wind speed (wind).
#' @param lake_levels a data frame with daily water level measurements as
#'                    formatted in the lake_levels dataset,
#'                    subset to lake level records for the lake of interest.
#' @param lst a data frame with sub-monthly lake surface temperature
#'            measurements as formatted in the lst_HOBO dataset, subset
#'            for a single lake.
#'
#' @return :
#' \item{z}{elevation above mean sea level (m)}
#' \item{phi}{latitude of location (radians). Positive for northern
#'            hemisphere, negative for southern hemisphere}
#' \item{Lm}{longitude of location (degrees west of Greenwich)}
#' \item{Lz}{longitude of location's measurement timezone (degrees west of
#'             Greenwich). For example, Lz = 75, 90, 105 and 120° for
#'             measurement times based on the Eastern, Central, Rocky Mountain
#'             and Pacific time zones (United States) and Lz = 0° for Greenwich,
#'             330° for Cairo (Egypt), and 255° for Bangkok (Thailand).}
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise ungroup arrange select
#' @importFrom rlang .data
#' @importFrom zoo read.zoo na.approx
#'
#' @export
CSLS_find_common_dates <- function(weather, lake_levels, lst) {

  # Determine first day lst, weather, and lake level data
  lst_day0            <- floor_date(min(lst$date), unit = "day")
  weather_day0        <- floor_date(min(weather$date), unit = "day")
  lake_levels_day0    <- floor_date(min(lake_levels$date), unit = "day")
  day0                <- max(c(lst_day0, weather_day0, lake_levels_day0))

  # Determine last day lst, weather, and lake level data
  lst_day_end         <- floor_date(max(lst$date), unit = "day")
  weather_day_end     <- floor_date(max(weather$date), unit = "day")
  lake_levels_day_end <- floor_date(max(lake_levels$date), unit = "day")
  day_end             <- min(c(lst_day_end, weather_day_end, lake_levels_day_end))

  # Determine intervals of overlap
  date_interval1 <- interval(day0, day_end) # inclusive of day0
  date_interval2 <- interval(day0 + days(1), day_end) # exclusive of day0

  # Create vecotor with all days in interval
  i           <- 1
  date_vec    <- NULL
  this_date   <- day0
  while (this_date < day_end){
    this_date   <- day0 + days(i)
    date_vec[i] <- this_date
    i           <- i + 1
  }
  date_vec <- as.data.frame(date_vec)
  colnames(date_vec) <- 'date'
  date_vec$date <- as_datetime(date_vec$date)

  # Subset data to date intervals
  lst         <- lst[lst$date %within% date_interval1,]
  daily_lst   <- lst %>%
                 group_by(date = floor_date(.data$date, unit = "day")) %>%
                 summarise(ltmp = mean(.data$ltmp, na.rm = TRUE)) %>%
                 ungroup()
  wtmp0       <- daily_lst$ltmp[daily_lst$date == day0]
  lst         <- daily_lst[daily_lst$date %within% date_interval2,]

  weather     <- weather[weather$date %within% date_interval2,]
  lake_levels <- lake_levels[lake_levels$date %within% date_interval2,]

  # Ensure all have all dates for lake levels
  lake_levels       <- merge(lake_levels, date_vec, all.y = TRUE)
  lake_levels       <- lake_levels %>%
                       arrange(.data$date) %>%
                       select(.data$date,
                              .data$level_m,
                              .data$area_m2,
                              .data$vol_m3)
  zoo.levels        <- read.zoo(lake_levels)
  zoo.levels        <- as.data.frame(na.approx(zoo.levels))
  lake_levels[,2:4] <- zoo.levels

  return(list(weather = weather,
              lst = lst,
              lake_levels = lake_levels,
              wtmp0 = wtmp0))
}