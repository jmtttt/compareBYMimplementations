#' Smoothing map with MCMC (OpenBUGS) and BYM model
#'
#' @param formula a formula in one of the form
##' \itemize{
##'  \item{}{observed ~ expected + region_ID}
##' }
#' @param data data.frame corresponding to formula
#' @param shapefile of class SpatialPolygonsDataFrame (see package sp),
#'                  where rownames of shapefile\@data equal region_ID in formula
#' @param time_scaling Openbugs uses 10,000 * time_scaling iterations,
#'                     First half is burnin
#' @param BUGS_tempdir Some (empty) directory for OpenBUGS to store model information
#' @param ... more parameters of other methods (to unify call)
#'
#' @return A two element list
##' \itemize{
##'  \item{"SIR"}{A tibble (=data.fame) with region_ID and SIR (observed / expected)}
##'  \item{"run"}{the full OpenBUGS object for debugging}
##' }
#' @export
#'
#' @import checkmate
#' @import spdep
#' @import tidyverse
#' @import R2OpenBUGS
#' @import lubridate
#'
smooth_bugs <- function(formula,
                        data,
                        shapefile,
                        BUGS_tempdir,
                        time_scaling = 1,
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
  checkmate::assert_character(x = BUGS_tempdir)
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
  # prepare parameters to pass to OpenBUGS
  #
  N <- nrow(shapefile@data)
  ls_neighbours <- spdep::poly2nb(shapefile)
  ls_neighbours <- spdep::nb2WB(ls_neighbours)
  #
  MCMC_param <- list(N = N,
                     observed = shapefile@data$O,
                     expected = shapefile@data$E,
                     adj = ls_neighbours$adj,
                     weights = ls_neighbours$weights,
                     num = ls_neighbours$num)
  #
  inits    <- list(list(u = rep(0, N), v = rep(0, N), alpha = 1, precu = 0.001, precv = 0.001),
                   list(u = rep(0, N), v = rep(0, N), alpha = 1, precu = 0.001, precv = 0.001),
                   list(u = rep(0, N), v = rep(0, N), alpha = 1, precu = 0.001, precv = 0.001) )
  #
  ###########################################
  #
  # Run OpenBUGS
  #
  ###########################################
  #
  # find model
  #
  MCMC_model <- model_incidence()

  #
  # store model to BUGS_tempdir
  #
  dir.create(path = BUGS_tempdir,
             showWarnings = F,
             recursive = T)
  write(x = MCMC_model,
        file = file.path(BUGS_tempdir,
                         "mcmc_model.txt"))
  #
  # call OpenBUGS
  #
  MCMC_run <- R2OpenBUGS::bugs(data = MCMC_param,
                               inits = inits,
                               parameters.to.save = c("theta",
                                                      "alpha",
                                                      "u",
                                                      "v",
                                                      "sigmau",
                                                      "sigmav"),
                               model.file = "mcmc_model.txt",
                               n.chains = 3,
                               n.iter =   ceiling(10000 * time_scaling),    # n.iter =   9800,
                               n.burnin = ceiling( 5000 * time_scaling),    # 2  n.burnin = 9700,
                               n.thin = 10,                                 # 3  n.thin = 10,
                               saveExec = TRUE,
                               restart = FALSE,
                               working.directory = BUGS_tempdir,
                               debug = F)



  shapefile@data <- dplyr::bind_cols(shapefile@data,
                                     as.data.frame(MCMC_run$summary)[1:N,])

  MCMC_result <- shapefile@data %>% dplyr::select(K, mean, "2.5%", "97.5%") %>%
    dplyr::rename(region_ID = K,
                  theta = mean) %>%
    dplyr::select(region_ID, theta, "2.5%", "97.5%") %>%
    dplyr::arrange(region_ID) %>%
    tibble::as_tibble()
  colnames(MCMC_result) <- c("region_ID", "mean", "CI_lower", "CI_upper")
  #
  # return
  #
  time_stop <- lubridate::now()
  t_diff <- as.numeric(lubridate::seconds(time_start %--% time_stop))
  return(list(sRR = MCMC_result,
              run = MCMC_run,
              duration = t_diff))
}

#' OpenBUGS-model for BYM model
#'
#' @return the model description in a character element (multiple lines separated by "\ n")
#'
model_incidence = function(){
  paste(
    "model",
    "{",
    "  for(i in 1:N)",
    "  {",
    "    observed[i]   ~ dpois(mu[i])",
    "    mu[i]         <- expected[i]*theta[i]",
    "    log(theta[i]) <- alpha + u[i] + v[i]",
    "    u[i]          ~ dnorm(0, precu)",
    "  }",
    "",
    "  v[1:N] ~ car.normal(adj[], weights[], num[], precv)",
    "  alpha  ~ dflat()",
    "  precu  ~ dgamma(0.001, 0.001)",
    "  precv  ~ dgamma(0.1, 0.1)",
    "  sigmau <- 1/precu",
    "  sigmav <- 1/precv",
    "}",
    sep = "\n"
  )
}
