#' Prepare GIS-Data
#'
#' @param dir Directory of Shapefiles
#' @param layer Name of Layer
#' @param verbose Report status on console?
#'
#' @return a SpatialPolygonsDataFrame
#' @export
#'
wrangle_shapefile <- function(dir,
                              layer,
                              verbose = TRUE){
  checkmate::assert_directory(dir)
  checkmate::assert_directory_exists(dir)
  checkmate::assert_character(layer)
  checkmate::assert_logical(verbose)

  if(verbose){
  print("Read geographic data ...")
  }

  Dshape <- rgdal::readOGR(dsn = dir,
                           layer = layer,
                           verbose = FALSE,
                           encoding = "UTF-8",
                           use_iconv = TRUE)
  # Kartendaten vereinfachen (warum auch immer das notwendig ist?)
  names(Dshape@data)[5] <- "GKZ"
  Dshape <- Dshape[Dshape@data$GF==4,]
  Dshape@data$GKZ <- substring(Dshape@data$GKZ,1,5)

  # Gewuenschte Spalten extrahieren
  Dshape@data <- Dshape@data[,c(5,7)]
  names(Dshape@data)[2] <- "GKZname"
  Dshape@data <- as.data.frame(Dshape@data[,c(1,2)])

  # Shapefile vereinfachen
  if(verbose){
    print("... and simplity them")
  }
  Dshape <- rmapshaper::ms_simplify(Dshape, keep = .05)
  Dshape@data$GKZ <- as.numeric(Dshape@data$GKZ)

  Dshape@data <- data.frame(region_ID = sprintf("%05d", Dshape@data$GKZ),
                            row.names = sprintf("%05d", Dshape@data$GKZ),
                            stringsAsFactors = FALSE)
  return(Dshape)
}
