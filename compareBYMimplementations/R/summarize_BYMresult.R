#' Calculate Coverage from "true" Simulation and BYM smoothing
#'
#' @param estimation_tbl What should be estimated
#' @param simulation_tbl What was estimated
#' @param method_name the software used
#' @param region an identifier for the simulated region
#' @param n the number of counties in region
#' @param ID an ID of the try
#' @param seed_risk seed for risk distribution
#' @param seed_poisson seed for case distribution
#' @param n_convolution number of related neighbours
#' @param rate_per_100k determinant for case numbers per county
#' @param duration length of calculation time
#'
#' @return a data.frame with one row and many information
#' @export
#'
summarize_BYMresult <- function(estimation_tbl,
                                simulation_tbl,
                                method_name = NA,
                                region = NA,
                                n = NA,
                                ID = NA,
                                seed_risk = NA,
                                seed_poisson = NA,
                                n_convolution = NA,
                                rate_per_100k = NA,
                                duration = NA){
  merged_res <- dplyr::left_join(x    = estimation_tbl,
                                 y    = simulation_tbl,
                                 by   = "region_ID")

  summarized_res <- merged_res %>%
    dplyr::summarise(  coverage_lower = mean(theta < CI_lower),
                       coverage_upper = mean(theta > CI_upper),
                       coverage       = mean((theta > CI_lower) & theta < CI_upper),
                       average_diff   = mean(abs(theta - mean)),
                       squared_diff   = mean((theta - mean)^2),
                       CI_length      = mean(CI_upper - CI_lower)) %>%
    tibble::add_column(method_name    = method_name,
                       region         = region,
                       n              = n,
                       ID             = ID,
                       seed_risk      = seed_risk,
                       seed_poisson   = seed_poisson,
                       n_convolution  = n_convolution,
                       rate_per_100k  = rate_per_100k,
                       duration       = duration,
                       .before        = "coverage_lower")

  return(summarized_res)
}
