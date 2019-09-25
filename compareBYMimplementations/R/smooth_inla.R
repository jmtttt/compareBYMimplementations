#' Smoothing map with INLA and BYM model
#'
#' @param formula a formula in one of the form
##' \itemize{
##'  \item{}{observed ~ expected + region_ID}
##' }
#' @param data data.frame corresponding to formula
#' @param shapefile of class SpatialPolygonsDataFrame (see package sp),
#' @param INLA_model one of "bym" and "bym2" - defines, how INLA should parametrize the model
#' @param ... more parameters of other methods (to unify call)
#'                  where rownames of shapefile\@data equal region_ID in formula
#' @return A two element list
##' \itemize{
##'  \item{"SIR"}{A tibble (=data.fame) with region_ID and SIR (observed / expected)}
##'  \item{"run"}{the full INLA object for debugging}
##' }
#' @export
#'
#' @import checkmate
#' @import spdep
#' @import INLA
#' @import tidyverse
#' @import lubridate
#'
#'
smooth_inla <- function(formula,
                        data,
                        shapefile,
                        INLA_model,
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
  # checkmate::assert_character(x = covariate_method)
  # checkmate::assert_choice(x = covariate_method,
  #                          choices = c("ignore",
  #                                      "additive"))
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
  # create data for INLA
  #
  ###########################################
  #
  # map region keys to IDs 1:N
  #
  shapefile@data <- shapefile@data %>% dplyr::select()
  shapefile@data$K  <- row.names(shapefile@data)
  shapefile@data$ID <- 1:nrow(shapefile@data)
  LINK <- shapefile@data %>% dplyr::select(K, ID)
  working_data <- dplyr::left_join(x = working_data,
                                   y = LINK,
                                   by = "K")
  working_data <- cbind(working_data,
                        ID.local = working_data$ID)
  #
  # Create model graph for INLA
  #
  ls_neighbours <- spdep::poly2nb(shapefile)
  INLA_graph    <- INLA::inla.read.graph(ls_neighbours)
  #
  ###########################################
  #
  # Call INLA
  #
  ###########################################
  #
  # define model
  #
  INLA_formula <- formula(O ~ 1)
  INLA_formula <- update.formula(old = INLA_formula,
                                 new = ~ . + f(ID.local,
                                               model=INLA_model,
                                               graph = INLA_graph,
                                               constr  = TRUE))
  #
  # run model
  #
  INLA_run  <-  INLA::inla(formula = INLA_formula,
                           family  = "poisson",
                           data    = working_data,
                           E       = E,
                           verbose = FALSE)
  #
  # extract result
  #
  #
  INLA_result <- tibble::tibble(ID = working_data$ID)
  INLA_result_u <- INLA_run$summary.random$ID.local[, c("ID", "mean")] %>%
    dplyr::rename(mean.u = mean)
  INLA_result <- dplyr::left_join(x = INLA_result,
                                  y = INLA_result_u,
                                  by = "ID")
  INLA_result <- INLA_result %>%
    tibble::add_column(alpha = INLA_run$summary.fixed[1, "mean"])

  INLA_result$log_theta = rowSums(INLA_result %>% dplyr::select(-ID))

  INLA_result <- INLA_result %>%
    dplyr::mutate(theta = exp(log_theta))

  INLA_result <- dplyr::left_join(x = working_data,
                                  y = INLA_result,
                                  by = "ID") %>%
    dplyr::select(ID, K, theta) %>%
    dplyr::rename(region_ID = K)

  CI_borders <- tibble(ID = 1:nrow(INLA_result),
                       "2.5%"  = exp(INLA_run$summary.random$ID.local$`0.025quant`[1:nrow(INLA_result)] +
                                       INLA_run$summary.fixed$mean),
                       "97.5%" = exp(INLA_run$summary.random$ID.local$`0.975quant`[1:nrow(INLA_result)] +
                                       INLA_run$summary.fixed$mean))
  INLA_result <- dplyr::left_join(x = INLA_result,
                                  y = CI_borders,
                                  by = "ID")

  INLA_result <- INLA_result %>%
    ungroup() %>%
    dplyr::select(setdiff(x = colnames(INLA_result),
                          y = c("ID"))) %>%
    dplyr::arrange(region_ID) %>%
    tibble::as_tibble()
  colnames(INLA_result) <- c("region_ID", "mean", "CI_lower", "CI_upper")
  #
  # return
  #
  time_stop <- lubridate::now()
  t_diff <- as.numeric(lubridate::seconds(time_start %--% time_stop))
  return(list(sRR = INLA_result,
              run = INLA_run,
              duration = t_diff))
}
