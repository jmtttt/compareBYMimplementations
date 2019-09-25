#' Geographical data for Germany
#'
#' @format A SpatialPolygonsDataFrame with 402 polygons and 2 variables:
#' \describe{
#'   \item{region_ID}{Amtlicher Gemeindeschlüssel, see https://en.wikipedia.org/wiki/Community_Identification_Number#Germany}
#'   \item{n}{average population in the time interval from 2008 to 2015}
#' }
#'
#' @details Access via load(file.path(find.package("compareBYMimplementations"), "data", "shape_GER.Rdata"))
#'
#' @note This is an edited shapefile.
#'       The original can be found at https://gdz.bkg.bund.de/index.php/default/open-data/gebietseinheiten-1-250-000-ge250.html
#'       by Bundesamt für Kartographie und Geodäsie, 2019
#'       with data sources from Statistisches Bundesamt (Destatis), Bundesinstitut für Bau-, Stadt- und Raumforschung (BBSR)
#'       under license "Data licence Germany – attribution – version 2.0" (http://www.govdata.de/dl-de/by-2-0)
#'
#'
"shape_GER"
