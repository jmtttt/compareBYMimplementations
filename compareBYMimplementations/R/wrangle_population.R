#' Bereite Populationsdaten vor
#'
#' @param file Dateipfad zu Datei mit Bevölkerungszahlen
#' @param years Bilde Mittelwert über die hier genannten Jahre
#' @param verbose Soll hin und wieder ein Status auf der Konsole ausgegeben werden?
#'
#' @return ein data-frame mit Bevölkerungszahlen nach Kreis, Geschlecht und Altersgruppe
#' @export
#'
#' @import checkmate
#' @import readr
#' @import dplyr
#'
wrangle_population <- function(file,
                               years   = 1900:2050,
                               verbose = TRUE){
  checkmate::assert_file(file)
  checkmate::assert_file_exists(file)
  checkmate::assert_vector(years)
  checkmate::assert_logical(verbose)

  if(verbose){
    print("Read population data")
  }
  suppressMessages(
    bevdata <- readr::read_csv2(file = file)
  )
  bevdata %<>%
    dplyr::filter(year %in%  years) %>%
    dplyr::ungroup() %>%
    dplyr::rename(agegroup = age_class) %>%
    dplyr::group_by(region_ID,agegroup,sex) %>%
    dplyr::summarise(bevoelkerung= sum(count)/length(unique(year))) %>%
    dplyr::ungroup()

  return(bevdata)
}
