#' Simulate spatially dependent random data
#'
#' @param shapefile a SpatialPolygonsDataFrame with data columns region_ID (= AGS-Key) and n (= population)
#' @param n_convolution number, how often local effect should be leveled with its neighbours
#' @param seed_risk seed for simulating theta
#' @param rate_per_100k average risk per 100k people
#' @param seed_poisson seed for simulatijng case numbers
#'
#' @return a SpatialPolygonsDataFrame with additional data columns u, v, mu, theta, O, E
#' @export
#'
#' @details The BYM-model assumes data to be like O_i ~ Poi(E_i * exp(u_i + v_i + mu))
#'     with observed (O) and expected (E) cases and random effects mu (constant), v (unique) and u (local).
#'     Exactly this model is simulated.
#'     The autocorrelation of u is simulated via averaging over its neighbours n_convolution times
#'
#' @import spdep
#' @import dplyr
#'
artificial_risk <- function(shapefile,
                            n_convolution,
                            rate_per_100k,
                            seed_risk,
                            seed_poisson){

  # calculate some useful variables
  N <- nrow(shapefile@data)
  nb_mat <- spdep::poly2nb(shapefile)

  # Simulate risk distribution
  set.seed(seed = seed_risk)

  shapefile@data <- shapefile@data %>%
    dplyr::mutate(u = rnorm(n = N, mean = 0, sd = 100),
                  v = rnorm(n = N, mean = 0, sd = .07),
                  mu = -.02)

  # force u to be regionally dependent and scale it to a reasonable value
  if(n_convolution > 0){
    for(i in 1:n_convolution){
      shapefile@data$u <- unlist(lapply(X = 1:N,
                                        FUN = function(x){
                                          res = mean(shapefile@data$u[nb_mat[[x]]])
                                        }))
    }
  }

  shapefile@data <- shapefile@data %>%
    dplyr::mutate(u = u / sd(u) / 7) %>%
    dplyr::mutate(u = u - mean(u)) %>%
    dplyr::mutate(theta = exp(u + v + mu))

  # Simulate events
  set.seed(seed = seed_poisson)

  shapefile@data <- shapefile@data %>%
    dplyr::mutate(E         = n / 100000 * rate_per_100k,
                  E_spatial = n / 100000 * rate_per_100k * theta)

  shapefile@data$O <- unlist(lapply(X = shapefile@data$E_spatial,
                                    FUN = function(x){
                                      rpois(n = 1,
                                            lambda = x)
                                    }))

  shapefile@data <- shapefile@data %>%
    dplyr::select(-E_spatial)
  rownames(shapefile@data) <- shapefile@data$region_ID

  return(shapefile)
}
