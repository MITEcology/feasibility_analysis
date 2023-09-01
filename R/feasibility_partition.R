#' Calculate the feasibility partition of all possible compositions
#'
#' @description calculate the feasibility partition of all possible compositions \eqn{\mathcal{C} \subseteq \mathcal{S}}
#'
#' @param matA Numeric, interaction matrix defining the dynamics of the system
#' @param nt Numeric, number of replications to reduce numerical instabilities, default to 30
#'
#' @importFrom purrr map_dbl
#' @importFrom pracma eye
#' @return a vector \eqn{\{\mathbf{\Omega}\}} of length \eqn{2^{|\mathcal{S}|}} for all \eqn{\mathcal{C} \subseteq \mathcal{S}}, and essentially \eqn{\sum_{\mathcal{C}} \Omega(D(\mathcal{C})) = 1}
#'
#' @note the interaction matrix matA needs to be globally stable in order for a meaningful partition of the parameter space.
#' @export
#'
#' @examples
#' matA <- generate_inte_gs(4, 0.8, 1, "norm")  ## Generate a random interaction matrix
#' feasibility_partition(matA)
feasibility_partition <- function(matA, nt = 30) {
  # function that partitions the parameter space from a given interaction matrix
  # params: inte = interaction matrix (representing LV dynamics)
  # return: a list of all possible regions (as matrices)
  get_region <- function(inte) {
    num <- ncol(inte)
    l <- 0
    A <- inte
    B <- eye(nrow(inte), ncol(inte))
    record0 <- get_compo(num, 0)
    region <- list()
    for (l in 1 : 2^num){
      region[[l]] <- A %*% diag(record0[l, ]) + B %*% diag(1 - (record0[l, ]))
    }
    return(region)
  }
  # here we use the partition on the sphere other than simplex
  omega_vec <- map_dbl(get_region(matA), ~feasibility_community(., nt = nt, raw = TRUE))
  return(omega_vec)
}