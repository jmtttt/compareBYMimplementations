#' Smoothing map with Rstan and BYM model
#'
#' @param formula a formula in one of the form
##' \itemize{
##'  \item{}{observed ~ expected + region_ID}
##' }
#' @param data data.frame corresponding to formula
#' @param shapefile of class SpatialPolygonsDataFrame (see package sp),
#' @param STAN_dir Some (empty) directory for Rstan to store model information
#' @param STAN_cores Number of cores for Rstan for parallel computing
#' @param ... more parameters of other methods (to unify call)
#' @param STAN_bym_model precompiled stan bym2-model
#'
#' @return A two element list
##' \itemize{
##'  \item{"SIR"}{A tibble (=data.fame) with region_ID and SIR (observed / expected)}
##'  \item{"run"}{the full Rstan object for debugging}
##' }
#' @export
#'
#' @import checkmate
#' @import dplyr
#' @import tibble
#' @import spdep
#' @import lubridate
#'
smooth_stan <- function(formula,
                        data,
                        shapefile,
                        STAN_dir,
                        STAN_cores = NULL,
                        STAN_bym_model,
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
  checkmate::assert_character(x = STAN_dir)
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
  nb_shapefile = spdep::poly2nb(shapefile);
  y = shapefile@data$O
  E = shapefile@data$E;




  nbs=myrstan_nb2graph(nb_shapefile);
  N = nbs$N;
  node1 = nbs$node1;
  node2 = nbs$node2;
  N_edges = nbs$N_edges;
  scaling_factor = myrstan_scale_nb_components(nb_shapefile)[1];

  if(!dir.exists(STAN_dir)){
    dir.create(path = STAN_dir,
               showWarnings = FALSE,
               recursive = TRUE)
  }

  if(!is.null(STAN_cores)){
    options(mc.cores = STAN_cores);
  }

  bym2_fit  = sampling(object        = STAN_bym_model,
                       data          = list(N,
                                            N_edges,
                                            node1,
                                            node2,
                                            y,
                                            E,
                                            scaling_factor),
                       control       = list(adapt_delta = 0.97,
                                            max_treedepth = 15),
                       chains        = 2,
                       warmup        = 1000,
                       iter          = 2000,
                       save_warmup   = FALSE,
                       open_progress = FALSE,
                       show_messages = FALSE);


  rstan_res <- summary(bym2_fit)$summary
  rstan_res <- as.data.frame(rstan_res)
  rstan_res$element <- rownames(rstan_res)

  rstan_res <- rstan_res %>%
    tibble::as_tibble() %>%
    dplyr::filter(substr(element, 1, 3) == "mu[") %>%
    dplyr::rename(CI_lower = "2.5%",
                  CI_upper = "97.5%") %>%
    tibble::add_column(region_ID = shapefile@data$K,
                       E         = shapefile@data$E) %>%
    dplyr::mutate(mean     = mean     / E,
                  CI_lower = CI_lower / E,
                  CI_upper = CI_upper / E) %>%
    dplyr::select(region_ID, mean, CI_lower, CI_upper)

  time_stop <- lubridate::now()
  t_diff    <- as.numeric(lubridate::seconds(time_start %--% time_stop))
  return(list(sRR      = rstan_res,
              run      = bym2_fit,
              duration = t_diff))
}




#' Rstan-model for BYM-2 model
#'
#' @return the model description in a character element (multiple lines separated by "\ n")
#' @export
#'
model_STAN_bym2 <- function(){
  paste(c(
    "data {",
    "  int<lower=0> N;",
    "  int<lower=0> N_edges;",
    "  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]",
    "  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]",
    "  int<lower=0> y[N];              // count outcomes",
    "  vector<lower=0>[N] E;           // exposure",
    "  real<lower=0> scaling_factor; // scales the variance of the spatial effects",
    "}",
    "transformed data {",
    "  vector[N] log_E = log(E);",
    "}",
    "parameters {",
    "  real beta0;                // intercept",
    "  real<lower=0> sigma;        // overall standard deviation",
    "  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance",
    "  vector[N] theta;       // heterogeneous effects",
    "  vector[N] phi;  // spatial effects",
    "}",
    "transformed parameters {",
    "  vector[N] convolved_re;",
    "  // variance of each component should be approximately equal to 1",
    "  convolved_re =  sqrt(1 - rho) * theta + sqrt(rho / scaling_factor) * phi;",
    "}",
    "model {",
    "  y ~ poisson_log(log_E + beta0 + convolved_re * sigma);",
    "  target += -0.5 * dot_self(phi[node1] - phi[node2]);",
    "  beta0 ~ normal(0, 1);",
    "  theta ~ normal(0, 1);",
    "  sigma ~ normal(0, 1);",
    "  rho ~ beta(0.5, 0.5);",
    "  // soft sum-to-zero constraint on phi)",
    "  sum(phi) ~ normal(0, 0.001 * N);  // equivalent to mean(phi) ~ normal(0,0.001)",
    "}",
    "generated quantities {",
    "  real log_precision = -2.0 * log(sigma);",
    "  real logit_rho = log(rho / (1.0 - rho));",
    "  vector[N] eta = log_E + beta0 + convolved_re * sigma;",
    "  vector[N] mu = exp(eta);",
    "}",
    ""
  ),
  sep = "\n")
}

#' Rstan-model for BYM model
#'
#' @return the model description in a character element (multiple lines separated by "\ n")
#' @export
#'
model_STAN_bym <- function(){
  paste(c(
    "data {",
    "  int<lower=0> N;",
    "  int<lower=0> N_edges;",
    "  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]",
    "  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]",
    "  int<lower=0> y[N];              // count outcomes",
    "  vector<lower=0>[N] E;           // exposure",
    "}",
    "transformed data {",
    "  vector[N] log_E = log(E);",
    "}",
    "parameters {",
    "  real beta0;                // intercept",
    "  real<lower=0> tau_theta;   // precision of heterogeneous effects",
    "  real<lower=0> tau_phi;     // precision of spatial effects",
    "  vector[N] theta;       // heterogeneous effects",
    "  vector[N] phi;         // spatial effects",
    "}",
    "transformed parameters {",
    "  real<lower=0> sigma_theta = inv(sqrt(tau_theta));  // convert precision to sigma",
    "  real<lower=0> sigma_phi = inv(sqrt(tau_phi));      // convert precision to sigma",
    "}",
    "model {",
    "  y ~ poisson_log(log_E + beta0 + phi * sigma_phi + theta * sigma_theta);",
    "  // NOTE:  no prior on phi_raw, it is used to construct phi",
    "  // the following computes the prior on phi on the unit scale with sd = 1",
    "  target += -0.5 * dot_self(phi[node1] - phi[node2]);",
    "  // soft sum-to-zero constraint on phi)",
    "  sum(phi) ~ normal(0, 0.001 * N);  // equivalent to mean(phi) ~ normal(0,0.001)",
    "   beta0 ~ normal(0, 5);",
    "  theta ~ normal(0, 1);",
    "  tau_theta ~ gamma(3.2761, 1.81);  // Carlin WinBUGS priors",
    "  tau_phi ~ gamma(1, 1);            // Carlin WinBUGS priors",
    "}",
    "generated quantities {",
    "  vector[N] mu = exp(log_E + beta0 + phi * sigma_phi + theta * sigma_theta);",
    "}",
    ""
  ),
  sep = "\n")
}

#' myrstan_nb2graph
#'
#' @param x nb_object
#'
#' @return data.frame containing num nodes, num edges,
#          and a list of graph edges from node1 to node2.
#'
myrstan_nb2graph = function(x) {
  N = length(x);
  n_links = 0;
  for (i in 1:N) {
    if (x[[i]][1] != 0) {
      n_links = n_links + length(x[[i]]);
    }
  }
  N_edges = n_links / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  idx = 0;
  for (i in 1:N) {
    if (x[[i]][1] > 0) {
      for (j in 1:length(x[[i]])) {
        n2 = unlist(x[[i]][j]);
        if (i < n2) {
          idx = idx + 1;
          node1[idx] = i;
          node2[idx] = n2;
        }
      }
    }
  }
  return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2));
}







#' myrstan_scale_nb_components
#'
#' @param x nb_object
#'
#' @return vector of per-component scaling factor (for BYM2 model) scaling factor for singletons is 0
#'
#' @import Matrix
#' @import spdep
#' @import INLA
#'
myrstan_scale_nb_components = function(x) {
  N = length(x);
  comp_ids = n.comp.nb(x)[[2]];
  offsets = myrstan_indexByComponent(comp_ids);

  comps = as.matrix(table(comp_ids));
  num_comps = nrow(comps);
  scales = vector("numeric", length=num_comps);
  for (i in 1:num_comps) {
    N_subregions = comps[i,1];
    scales[i] = 0.0;
    if (N_subregions > 1) {
      # get adj matrix for this component
      drops = comp_ids != i;
      nb_tmp = spdep::droplinks(x, drops);
      nb_graph = myrstan_nb2subgraph(nb_tmp, i, comp_ids, offsets);
      adj.matrix = sparseMatrix( i=nb_graph$node1, j=nb_graph$node2, x=1, dims=c(N_subregions,N_subregions), symmetric=TRUE);
      # compute ICAR precision matrix
      Q =  Diagonal(N_subregions, Matrix::rowSums(adj.matrix)) - adj.matrix;
      # Add a small jitter to the diagonal for numerical stability (optional but recommended)
      Q_pert = Q + Diagonal(N_subregions) * max(diag(Q)) * sqrt(.Machine$double.eps)
      # Compute the diagonal elements of the covariance matrix subject to the
      # constraint that the entries of the ICAR sum to zero.
      Q_inv = INLA::inla.qinv(Q_pert, constr=list(A = matrix(1,1,N_subregions),e=0))
      # Compute the geometric mean of the variances, which are on the diagonal of Q.inv
      scaling_factor = exp(mean(log(diag(Q_inv))))
      scales[i] = scaling_factor;
    }
  }
  return(scales);
}


# myrstan_indexByComponent
#
# input: vector of component ids
# returns:
#
# details: subfunction of myrstan_scale_nb_components
#
#' myrstan_indexByComponent
#'
#' @param x vector of component ids
#'
#' @return vector of per-component consecutive node ids
#'
#' @details  subfunction of myrstan_scale_nb_components
#'
myrstan_indexByComponent = function(x) {
  y = x;
  comps = as.matrix(table(x));
  num_comps = nrow(comps);
  for (i in 1:nrow(comps)) {
    idx = 1;
    rel_idx = 1;
    while (idx <= length(x)) {
      if (x[idx] == i) {
        y[idx] = rel_idx;
        rel_idx = rel_idx + 1;
      }
      idx = idx + 1;
    }
  }
  return(y);
}



#' myrstan_nb2subgraph for a given subcomponent, return graph as lists of node1, node2 pairs
#'
#' @param x nb object
#' @param c_id subcomponent id
#' @param comp_ids vector of subcomponent ids
#' @param offsets vector of subcomponent node numberings
#'
#' @return list of node1, node2 ids
#' @details subfunction of myrstan_scale_nb_components
#'
myrstan_nb2subgraph = function(x, c_id, comp_ids, offsets) {
  N = length(x);
  n_links = 0;
  for (i in 1:N) {
    if (comp_ids[i] == c_id) {
      if (x[[i]][1] != 0) {
        n_links = n_links + length(x[i]);
      }
    }
  }
  N_edges = n_links / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  idx = 0;
  for (i in 1:N) {
    if (comp_ids[i] == c_id) {
      if (x[[i]][1] != 0) {
        for (j in 1:length(x[[i]])) {
          n2 = unlist(x[[i]][j]);
          if (i < n2) {
            idx = idx + 1;
            node1[idx] = offsets[i];
            node2[idx] = offsets[n2];
          }
        }
      }
    }
  }
  return (list("node1"=node1,"node2"=node2));
}






