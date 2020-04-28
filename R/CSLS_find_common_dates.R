#' Given CSLS data, find common dates of record
#'
#' Finds date interval with records for lst, weather, and lake levels. Checks
#' for the latest "start" date in each of these timeseries, then checks for the
#' earliest "end" date in each timeseries, then restricts all input data to be
#' within the overlapping interval. Lastly, interpolates the \code{lake_levels}
#' dataset over \code{NA} values to ensure a continuous timeseries of lake area and
#' lake depth for evaporation calculations.
#'
#' Note that this function may cause undesired behavior when using
#' \code{\link{CSLS_daily_met}} with the Unmodified Hammon method or the
#' McJannet method with \code{use_lst = FALSE}. Neither uses the \code{lst}
#' dataset for calculations, yet other input data is restricted based on the
#' availability of the \code{lst} dataset. Unmodified Hamon also does not use
#' \code{lake_levels}, but is restricted by the avilability of that dataset as
#' well.
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
#' @return a list with the following items:
#' \item{weather}{same as input weather data frame, but filtered to only include
#'                data during the overlap time period}
#' \item{lst}{same as input lst data frame, but filtered to only include data
#'            during the overlap time period and summarized at a daily time
#'            step}
#' \item{lake_levels}{same as input lake_levels data frame, but filtered to only
#'                    include data during the overlap time period and with
#'                    missing (NA) lake levels interpolated}
#' \item{wtmp0}{initial lake surface temperatuer (degC) from the day before the
#'              weather, lst, and lake_levels timeseries begin}
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise ungroup arrange select
#' @importFrom rlang .data
#' @importFrom zoo read.zoo na.approx
#'
#' @export
CSLS_find_common_dates <- function(weather, lake_levels, lst) {

  # Determine first day lst, weather, and lake level data
  days_01 <- lst %>%
             group_by(.data$lake) %>%
             summarise(day0 = floor_date(min(.data$date), unit = "day"))
  days_02 <- weather %>%
             summarise(day0 = floor_date(min(.data$date), unit = "day"))
  days_03 <- lake_levels %>%
             group_by(.data$lake) %>%
             summarise(day0 = floor_date(min(.data$date), unit = "day"))
  day0    <- max(c(days_01$day0, days_02$day0, days_03$day0))

  # Determine last day lst, weather, and lake level data
  days_01 <- lst %>%
             group_by(.data$lake) %>%
             summarise(day0 = floor_date(max(.data$date), unit = "day"))
  days_02 <- weather %>%
             summarise(day0 = floor_date(max(.data$date), unit = "day"))
  days_03 <- lake_levels %>%
             group_by(.data$lake) %>%
             summarise(day0 = floor_date(max(.data$date), unit = "day"))
  day_end <- min(c(days_01$day0, days_02$day0, days_03$day0))

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
  daily_lst   <- lst[lst$date %within% date_interval1,]
  wtmp0       <- daily_lst[daily_lst$date == day0,]
  lst         <- daily_lst[daily_lst$date %within% date_interval2,]

  weather     <- weather[weather$date %within% date_interval2,]
  lake_levels <- lake_levels[lake_levels$date %within% date_interval2,]

  return(list(weather = weather,
              lst = lst,
              lake_levels = lake_levels,
              wtmp0 = wtmp0))
}
