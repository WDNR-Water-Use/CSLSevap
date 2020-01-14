#' Given CSLS data, format "loc" input parameter
#'
#' Extracts location information for the lake from the dictionary and reformats
#' for input into lake evaporation functions.
#'
#' @param dictionary a data frame with the obs_type ("LK"), elev_m, lat_deg,
#'                   and long_deg for the lake of interest.
#' @param Lz the longitude of the local timezone (degrees west of Greenwich,
#'           ranges from 0 to 360 degrees). Defaults to 90 for Central Time
#'           Zone, USA.
#'
#' @return **loc*, a list with the following location-specific parameters:
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
#' @importFrom dplyr select filter
#' @importFrom rlang .data
#' @importFrom NISTunits NISTdegTOradian
#'
#' @export
CSLS_format_loc <- function(dictionary, Lz = 90) {
  loc     <- dictionary %>%
             filter(.data$obs_type == "LK") %>%
             select(z = .data$elev_m,
                    phi = .data$lat_deg,
                    Lm = .data$long_deg)
  loc$phi <- NISTdegTOradian(loc$phi)
  loc$Lm  <- loc$Lm
  loc$Lz  <- Lz
  loc     <- as.list(loc)
  return(loc)
}
