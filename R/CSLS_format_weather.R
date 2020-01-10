#' Given CSLS data, format "weather" input parameter
#'
#' Takes hourly weather data, summarizes as daily weather data, and re-formats
#' for input to lake evaporation functions
#'
#' @param weather a data frame with hourly weather data incl. air temperature
#'                (atmp), relative humidity (RH), incoming solar radiation (Rs),
#'                precipitation (P), and wind speed (wind).
#' @param wind_elev height at which wind is measured (m), default: 3
#' @param z0 aerodynamic roughness of weather measurement site (m), default: 0.02
#'
#' @return weather a list with weather data that includes:
#' \item{dt}{string indicating the timestep of input weather series. Expects
#'              "hourly", "daily", or "monthly".}
#' \item{datetimes}{datetimes of weather records [POSIXct]. If monthly timestep,
#'                  make sure date is the 15th of each month.}
#' \item{atmp}{If hourly timestep, vector of air temperature (degrees C)
#'             corresponding with datetimes vector. If daily or monthly
#'             timestep, list with two vectors, "min" and "max", with mean daily
#'             min and max air temperature (degrees C) corresponding with
#'             datetimes vector}
#' \item{RH}{If hourly timestep, vector of relative humidity (percent)
#'           corresponding with datetimes vector. If daily or monthly timestep,
#'           list with two vectors, "min" and "max", with mean daily min and max
#'           relative humidity (percent) corresponding with datetimes vector.}
#' \item{Rs}{vector of incoming solar or shortwave radiation (MJ/m^2/hr if
#'           hourly timestep, MG/m^2/day if daily or monthly), corresponding
#'           with datetimes vector.}
#' \item{wind}{vector with mean wind speed (m/s), corresponding with datetimes
#'             vector.}
#' \item{wind_elev}{atomic number, height at which wind is measured (m)}
#' \item{z0}{aerodynamic roughness of weather measurement site (m)}
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise
#' @importFrom rlang .data
#' @importFrom NISTunits NISTdegTOradian
#'
#' @export
CSLS_format_weather <- function(weather, wind_elev = 3, z0 = 0.02) {
  # Summarize at daily weather
  daily_weather <- weather %>%
                   group_by(datetimes = floor_date(.data$date, unit = "day")) %>%
                   summarise(atmp_min = min(.data$atmp),
                             atmp_max = max(.data$atmp),
                             RH_min = min(.data$RH),
                             RH_max = max(.data$RH),
                             P = sum(.data$P),
                             Rs = sum(.data$Rs),
                             wind = mean(.data$wind)) %>%
                   ungroup()

  # Convert to list for input to lake evap function (minus lake temp info)
  weather           <- daily_weather %>%
                       select(.data$datetimes, .data$P, .data$Rs, .data$wind) %>%
                       as.list()
  weather$atmp      <- list(min = daily_weather$atmp_min,
                            max = daily_weather$atmp_max)
  weather$RH        <- list(min = daily_weather$RH_min,
                            max = daily_weather$RH_max)
  weather$wind_elev <- wind_elev
  weather$dt        <- "daily"
  weather$z0        <- z0

  return(weather)
}
