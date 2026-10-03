#' Daily bike trips from and to the Duomo district of Milan (BikeMi)
#'
#' A dataset containing the daily number of trips of the BikeMi bike sharing
#' service of Milan that started and ended in the Duomo district, from the
#' 25th of January to the 6th of March 2016 (the 25th of February is
#' excluded), together with daily meteorological covariates.
#'
#' @format A data frame with 41 rows and 6 variables:
#' \describe{
#'   \item{start}{number of trips started in Duomo on a given day}
#'   \item{end}{number of trips ended in Duomo on a given day}
#'   \item{we}{1 if the day is a weekend day, 0 otherwise}
#'   \item{rain}{mean amount of rain during the day}
#'   \item{dtemp}{difference between the average temperature of the day and
#'   that of the period}
#'   \item{we_rain}{interaction between weekend and rain}
#' }
#' @source Torti, Pini, Vantini (2021), "Modelling time-varying mobility flows
#'   using function-on-function regression: Analysis of a bike sharing system
#'   in the city of Milan", Journal of the Royal Statistical Society Series C
#'   70(1), 226-247.
"bikeMi"
