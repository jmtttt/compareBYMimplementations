#' Smoothing map with R package CARBayrs and BYM model
#'
#' @param formula a formula in one of the form
##' \itemize{
##'  \item{}{observed ~ expected + region_ID}
##' }
#' @param data data.frame corresponding to formula
#' @param shapefile of class SpatialPolygonsDataFrame (see package sp),
#'                  where rownames of shapefile\@data equal region_ID in formula
#' @param ... more parameters of other methods (to unify call)
#'
#' @return A two element list
##' \itemize{
##'  \item{"SIR"}{A tibble (=data.fame) with region_ID and SIR (observed / expected)}
##'  \item{"run"}{the full CARBayes object for debugging}
##' }
#' @export
#'
#' @details For number of MCMC-iterations, the default is choosen as in the CARBayes-Paper (2013)
#'
#' @import CARBayes
#' @import checkmate
#' @import dplyr
#' @import lubridate
#' @import spdep
#'
smooth_rpkg <- function(formula,
                        data,
                        shapefile,
                        ...){
  time_start <- lubridate::now()
  ###########################################
  #
  # check input data for validity and recode
  #
  ###########################################
  #
  # check, if input as expected
  #
  checkmate::assert_formula(x = formula)
  checkmate::assert_data_frame(x = data,
                               min.cols = 3)
  checkmate::assert_class(x = shapefile,
                          classes = "SpatialPolygonsDataFrame")
  #
  # create working dataset for calculation
  #
  working_data <- model.frame(formula = formula,
                              data = data,
                              na.action = NULL)
  names(working_data) <- c("O", "E", "K")
  #
  # only continue, if data and neighbours can be mapped
  #
  checkmate::assert_set_equal(x = working_data$K,
                              y = row.names(shapefile@data))
  #
  ###########################################
  #
  # Glue together working data
  #
  ###########################################
  #
  # Prepare shapefile-data
  #
  shapefile@data$K <- row.names(shapefile@data)
  shapefile@data <- dplyr::left_join(x = shapefile@data %>% dplyr::select(K),
                                     y = working_data,
                                     by = "K")
  #
  # prepare parameters to pass to CARBayes
  #
  ls_neighbours  <- spdep::poly2nb(shapefile)
  mat_neighbours <- spdep::nb2mat(neighbours = ls_neighbours,
                                  style = "B")

  run_CARbayes <- CARBayes::S.CARbym(formula  = O ~ E,
                                     data     = shapefile@data,
                                     family   = "poisson",
                                     W        = mat_neighbours,
                                     burnin   = 20000,
                                     n.sample = 100000,
                                     thin = 10)

  A <- run_CARbayes$samples$fitted
  B <- data.frame(region_ID = shapefile@data$K,
                  mean = summary(A)$statistics[, 1],
                  CI_lower = summary(A)$quantiles[,1],
                  CI_upper = summary(A)$quantiles[,5])

  B <- B %>%
    dplyr::mutate(mean = mean / shapefile@data$E,
                  CI_lower = CI_lower / shapefile@data$E,
                  CI_upper = CI_upper / shapefile@data$E)
  #
  # return
  #
  time_stop <- lubridate::now()
  t_diff <- as.numeric(lubridate::seconds(time_start %--% time_stop))
  return(list(sRR = B,
              run = run_CARbayes,
              duration = t_diff))
}

