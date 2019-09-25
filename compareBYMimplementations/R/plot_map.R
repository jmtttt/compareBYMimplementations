#' Plotting and colouring a map
#'
#' @param formula a formula in the form plot_value ~ region_ID. Values of plot_value less (or equal to) 0 are ignored in the calculation of colorings and are assigned lowest coloring.
#' @param shapefile the SpatialPolygonsDataFrame shapefile for plotting with color values in data slot
#' @param title caption of legend as character string
#'
#' @details creates a map with the information passed in the data slot.
#'     the coloring scheme is viridisLite::cividis.
#'     Coloring steps are geometrically distributed.
#'     if the (non-zero) minimum is less tha 1, the colors are centered at 1
#'
#' @return a tmap object with the coloured map
#' @export
#'
#' @import tmap
#' @import viridisLite
#' @import stats
#'
plot_map <- function(formula,
                     shapefile,
                     title = ""){


  shapefile@data <- stats::model.frame(formula = formula,
                                data    = shapefile@data)
  colnames(shapefile@data) <- c("plot_value",
                                "region_ID")



  col_min <- log(min(shapefile@data$plot_value[shapefile@data$plot_value > 0]))
  col_max <- log(max(shapefile@data$plot_value))

  if(col_min < 0){
    col_min <- min( col_min,
                    -col_max)
    col_max <- max( -col_min,
                    col_min)
  }

  res_map <- tmap::tm_shape(shapefile) +
    tmap::tm_fill(col = "plot_value",
                  title = title,
                  breaks = exp(seq(from       = col_min,
                                   to         = col_max,
                                   length.out = 28)),

                  palette = viridisLite::cividis(n = 27, direction = -1))+
    tmap::tm_borders() +
    tmap::tm_layout(fontface = 2,
                    legend.outside = TRUE,
                    frame = FALSE,
                    legend.position = c("left",
                                        "bottom"))

  return(res_map)
}
