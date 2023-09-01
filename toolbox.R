# nolint start
# R package for the feasibility analysis in ecology

if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(mvtnorm)) {install.packages("mvtnorm"); library(mvtnorm)}
if (!require(mgcv)) {install.packages("mgcv"); library(mgcv)}
if (!require(geometry)) {install.packages("geometry"); library(geometry)}
if (!require(uniformly)) {install.packages("uniformly"); library(uniformly)}
if (!require(furrr)) {install.packages("furrr"); library(furrr)}
if (!require(pracma)) {install.packages("pracma"); library(pracma)}
if (!require(combinat)) {install.packages("combinat"); library(combinat)}
if (!require(deSolve)) {install.packages("deSolve"); library(deSolve)}


# 1. Matrix generation. -----
#Functions to generate different types of matrices to which feasibility can be computed. 

# generate random interaction matrix
# input:
  # S: number of species
  # sigma: standard deviation of interaction strength
  # conne: connectivity of interaction matrix
  # mu: mean of interaction strength
  # dist: distribution of interaction strength
# output: an S times S interaction matrix
generate_inte_rand <- function(S, sigma, conne = 1, mu = 0, dist = "norm") {
  if (dist == "norm") {
    matA <- rnorm(S * S, mean = mu, sd = sigma)
  } else if (dist == "lnorm") {
    matA <- -rlnorm(S * S, mean = mu, sd = sigma)
  }
  zeroes <- sample(
    c(rep.int(1, floor(S * S * conne)), rep.int(0, (S * S - floor(S * S * conne))))
  )
  matA[which(zeroes == 0)] <- 0
  matA <- matrix(matA, ncol = S, nrow = S)
  diag(matA) <- -1
  return(matA)
}

# generate interaction matrix with fixed fraction of postive/negative interactions and fully connectivity
# input:
  # S: number of species
  # sigma: standard deviation of interaction strength
  # mu: mean of interaction strength
  # frac: fraction of negative(competitive) interactions out of all off-diagonal interactions
  # dist: distribution of interaction strength (lnrom or halfnorm)
# output: an S times S interaction matrix
generate_inte_frac <- function(S, sigma, frac, mu = 0, dist = "halfnorm") {
  Nb <- S * (S-1); Nn <- floor(Nb*frac); Np <- Nb-floor(Nb*frac)
  if (dist == "lnorm") {
    values <- rlnorm(Nb, mean = mu, sd = sigma)
  } else if (dist == "halfnorm") {
    values <- abs(rnorm(Nb, mean = mu, sd = sigma))
  }
  values <- values * sample(c(rep.int(-1, Nn), rep.int(1, Np)))
  matA <- matrix(0, ncol = S, nrow = S)
  matA[which(diag(rep(1,S)) == 0)] <- values
  diag(matA) <- -1
  return(matA)
}

check_negative_definite <- function(A){
  all(eigen(A+t(A))$values < 0)
}
# generate negative definite random interaction matrix
# input:
  # S: number of species
  # sigma: standard deviation of interaction strength
  # conne: connectivity of interaction matrix
  # mu: mean of interaction strength
  # dist: distribution of interaction strength
# output: an S times S interaction matrix
generate_inte_gs <- function(S, sigma, conne = 1, dist = "norm", mu = -0.1*sigma, trail = 0) {
  if (dist == "norm") {
    matA <- rnorm(S * S, mean = mu, sd = sigma)
  } else if (dist == "lnorm") {
    matA <- -rlnorm(S * S, mean = mu, sd = sigma)
  }
  zeroes <- sample(
    c(rep.int(1, floor(S * S * conne)), rep.int(0, (S * S - floor(S * S * conne))))
  )
  matA[which(zeroes == 0)] <- 0
  matA <- matrix(matA, ncol = S, nrow = S)
  diag(matA) <- -1
  if(check_negative_definite(matA)) {
    return(matA)
  } else if (trail < 100) {
    trail <- trail + 1
    generate_inte_gs(S, sigma, conne, mu, dist, trail)
  } else {
    stop("Error: cannot generate negative definite matrix within 100 trails")
  }
}

#Define symmetric interaction matrix
A_int <- matrix(rep(0,9), ncol=3)
A_int[,1] <- -c(1, 0, 0)
A_int[,2] <- -c(0.3406884710289558, 0.9328944522580684, 0.11678744219579273)
A_int[,3] <- -c(0.3406884710289558, 0.12410827817548943, 0.9319487652208257)

species_names <- paste0("sp_",1:ncol(A_int))
colnames(A_int) <- species_names
rownames(A_int) <- species_names
A_int

# 2. Estimations the size and shape of the feasibility domain. ----
#Functions to estimate the size and the shape of the feasibility domain of a community.

# 3. Comparison of the feasibility domain between communities. ----
#Functions to compare the feasibility domain among communities that either differ in the number of species
#or in the multitrophic structure.

# 4. Estimations of the persistence of individual species and entire communities. ----
#Functions that take into account the size and shape of the feasibility domain of a given community
#to estimate whether individual species can persist as well as the whole community.

# 5. Effects of environmental variability on species persistence. ---- 
#Functions that take a probabilistic approach to estimate the role of environmental variability
#and perturbations to estimate the likelihood in individual species and entire communities to persist 





# calculate the feasibility of a community
# input:
  # A: interaction matrix
  # nt: number of trails
  # raw: TRUE for community-level; FALSE for species-level
# equivalent to calculate_omega(A, method = "sphere")
feasibility_community <- function(A, nt = 30, raw = TRUE) {
  S <- nrow(A)
  omega <- function(S, Sigma) {
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    # out <- d[1]^(1 / S) # species level
    out <- d[1] # community level
    return(out)
  }
  f <- function(m) class(try(solve(t(m) %*% m), silent = T)) == "matrix"
  if (all(f(A) == FALSE)) {
    return(0)
  }
  else {
    Sigma <- solve(t(A) %*% A)
    return(replicate(nt, omega(S, Sigma)) %>% mean())
  }
}

# calculate the feasibility overlap between two communities
# input:
  # A: interaction matrix of community 1
  # B: interaction matrix of community 2 (the same dimension as A)
  # nsamples: number of sampling points
  # nt: number of trails
feasibility_overlap <- function(A, B, raw = TRUE, nsamples = 3000, nt = 10) {
  num <- nrow(A)

  overlap_vertex <- vertex_detection(A, B) %>%
    cbind(vertex_detection(B, A)) %>%
    unique(MARGIN = 2)

  if (qr(overlap_vertex)$rank < num) {
    volume_overlap <- 0
  } else {
    volume_overlap <- tryCatch(
      {
        replicate(nt, calculate_omega(overlap_vertex, raw, nsamples)) %>% mean()
      },
      error = function(cond) {
        0
      }
    )
  }
  return(volume_overlap)
}

# calculate the asymmetry of the feasibility domain of a community
# input:
  # A: interaction matrix
feasibility_asymmetry <- function(matA) {
  # standard deviation of all column lengths
  sd(apply(matA, 2, function(x) sqrt(sum(x^2))))
}

# calculate the feasibility partition of all possible compositions
# input:
  # matA: interaction matrix
  # nt: number of trails
# output: a vector of feasibility values for all possible compositions
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

  # there are two definitions of partitions that need to be compared
  omega_vec <- map_dbl(get_region(matA), ~feasibility_community(., nt = nt, raw = TRUE))
  return(omega_vec)
}

# function that computes the normalized feasibility from an interaction matrix
calculate_omega <- function(vertex, raw = FALSE, nsamples = 100,
                            method = "convex_hull") {
  num <- nrow(vertex)
  vertex <- norm2(vertex)
  
  if (method == "convex_hull") {
    set.seed(1010)
    vertex <- cbind(
      vertex,
      vertex %*% t(abs(runif_on_sphere(n = nsamples, d = ncol(vertex), r = 1)))
    )
    if (num < 5) {
      vertex <- norm2(vertex) %*% diag(
        runif(
          ncol(vertex),
          (1 - .05 * (num - 2)),
          (1 + .05 * (num - 2))
        )
      )
    } else {
      vertex <- norm2(vertex) %*% diag(
        runif(
          ncol(vertex),
          (1 - .05 * (num - 2)),
          (1 + .1 * (num - 2))
        )
      )
    }
    
    vertex <- cbind(vertex, rep(0, num))
    
    vol_ori <- (convhulln(t(vertex), output.options = TRUE)$vol)
    vol_ball <- (pi^(num / 2) / gamma(num / 2 + 1))
    # vol_ball <- calculate_omega(diag(num), nsamples = nsamples)
    
    omega <- ifelse(raw == FALSE,
                    (vol_ori / vol_ball)^(1 / num),
                    vol_ori / vol_ball
    )
  }
  if (method == "sphere") {
    m <- matrix(0, num, 1)
    a <- matrix(0, num, 1)
    b <- matrix(Inf, num, 1)
    d <- pmvnorm(
      lower = rep(0, num),
      upper = rep(Inf, num),
      mean = rep(0, num), sigma = solve(t(vertex) %*% vertex)
    )
    omega <- ifelse(raw == FALSE,
                    d[1]^(1 / num),
                    d[1]
    )
  }
  return(omega)
}

# calculate the feasibility of species survival
# input:
  # matA: interaction matrix
  # sp: which species (index in the matrix) to calculate
  # nt: number of trails
feasibility_species <- function(matA, sp, nt = 30) {
  # add up all the conditional probabilities
  S <- ncol(matA)
  if(sp > S) stop("sp is out of range")
  feas_regions <- which(get_compo(S, 0)[, sp] == 1)
  return(feasibility_partition(matA, nt = nt)[feas_regions] %>% sum())
}

# calculate the full resistance of community S to parameter perturbations
# input:
  # matA: interaction matrix
  # r: intrinsic growth rate vector
  # nt: number of sampling on surface
  # all: TRUE for all distance values; FALSE for the minimum distance
  # norm: "l1" for L-1 norm; "l2" for L-2 norm
# output: distance to border (full resistance)
feasibility_resistance_full <- function(matA, r, norm = "l2", nt = 100, all = FALSE) {
  simplex_sampling <- function(m, n) {
    r <- list()
    for (j in 1:m) {
      dist <- c(sort(runif(n-1, 0, 1)), 1)
      r[[j]] <- c(dist[1], diff(dist))
    }
    return(r)
  }
  euclidean_distance <- function(a, b) {
    sqrt(sum((a - b)^2))
  }
  arc_length <- function(a, b) {
    acos(sum(a * b))
  }

  vertices <- combn(1:nrow(matA), nrow(matA)-1, simplify = FALSE)
  distances_list <- list()
  for (i in 1:length(vertices)) {
    vertex <- vertices[[i]]
    if (norm == "l1") {
      r <- norm1(r); matA <- norm1(matA)
      distances_list[[i]] <- 1:nt %>% 
        map_dbl(function(x) {
          t <- unlist(simplex_sampling(1, nrow(matA)-1))
          border_point <- c(matA[, vertex] %*% t)
          euclidean_distance(r, border_point)
        })
    }
    if (norm == "l2") {
      r <- norm2(r); matA <- norm2(matA)
      distances_list[[i]] <- 1:nt %>% 
        map_dbl(function(x) {
          t <- unlist(simplex_sampling(1, nrow(matA)-1))
          border_point <- c(matA[, vertex] %*% t)
          border_point_norm <- border_point / sqrt(sum(border_point^2))
          arc_length(r, border_point_norm)
        })
    }
  } 
  names(distances_list) <- unlist(lapply(vertices, paste, collapse = "_"))
  if (all) {
    return(distances_list)
  } else {
    return(sapply(distances_list, min))
  }
}

# calculate the partial resistance of community S to parameter perturbations
# input:
  # matA: interaction matrix
  # r: intrinsic growth rate vector
  # nt: number of sampling on surface
  # all: TRUE for all distance values; FALSE for the minimum distance
  # norm: "l1" for L-1 norm; "l2" for L-2 norm
# output: distance to border (full resistance)
feasibility_resistance_partial <- function(matA, r, norm = "l2") {
  euclidean_distance <- function(a, b) {
    sqrt(sum((a - b)^2))
  }
  arc_length <- function(a, b) {
    acos(sum(a * b))
  }

  if (norm == "l1") {
    r <- norm1(r); matA <- norm1(matA)
    distances <- 1:nrow(matA) %>% 
      map_dbl(~euclidean_distance(r, matA[, .x]))
  }
  if (norm == "l2") {
    r <- norm2(r); matA <- norm2(matA)
    distances <- 1:nrow(matA) %>% 
      map_dbl(~arc_length(r, matA[, .x]))
  }
  names(distances) <- 1:nrow(matA)
  return(distances)
}

# calculate the recovery of community S to abundance perturbations
# input:
  # matA: interaction matrix
  # r: intrinsic growth rate vector
  # type: "full" for full recovery; "part" for partial recovery
# output: partial recovery; full recovery
feasibility_recovery <- function(matA, r, type) {
  N <- solve(matA) %*% r
  matJ <- matA %*% diag(c(N))
  if (type == "full") {
    return(Re(eigen(matJ)$values)[1])
  }
  if (type == "part") {
    return(Re(eigen(matJ)$values)[ncol(matA)-1])
  }
}

# calculate contribution (to the orginal community) from a new species
# input:
  # matA: interaction matrix
  # sp: which species (index in the matrix) to calculate
  # nt: number of trails
feasibility_contribution <- function(matA, sp, nt = 30) {
  S <- ncol(matA)
  remaining_sp <- setdiff(1:S, sp)
  record0 <- get_compo(S,0)

  sub_matA <- matA[-sp, -sp]
  f1 <- feasibility_community(sub_matA, nt = nt)
  f2 <- sum(feasibility_partition(matA, nt = nt)[c(2^S, 2^S-sp)])

  long_term_eff <- f2 / f1
  return(long_term_eff)
}

# calculate feasibility of species invasion
# input:
  # A: interaction matrix
  # invader: which species (index in the matrix) to calculate
  # ns: number of sampling points
feasibility_invasion <- function(matA, invader, nt = 5000) {
  S <- nrow(matA)
  res <- c(1:S)[-invader]
  
  # sample r uniformly from the unit sphere
  r <- runif_on_sphere(n = nt, d = S, r = 1) %>% split(. , row(.)) %>% unname()
  spanA <- norm2(matA)

  # sample r from the geometric region of igr > 0 (assume global stability)
  check_igr <- map_dbl(r, ~inside_detection(cbind(diag(S)[, invader], spanA[, res]), .))
  res_only <- map_dbl(r, ~inside_detection(cbind(-diag(S)[, invader], spanA[, res]), .))

  # list of r inside the region
  r_inside <- r[sort(c(which(check_igr == 1), which(res_only == 1)))]

  # number of r inside the region
  n_point_inside <- length(r_inside) 
  
  # calculate the colonization probability with positive igr 
  igr <- sum(check_igr)/n_point_inside

  # output
  return(igr)
}




##################################################
#dependent functions for implementation of above functions ####
###################################################
###################################################




# function that normalizes a vector in the L1 norm
norm1 <- function(a) {
  if(is.matrix(a)) {
    a <- apply(a, 2, function(x) x / sum(x))
  } else {
    a <- a / sum(a)
  }
  return(a)
}

# function that normalizes a vector in the L2 norm
norm2 <- function(a) {
  if(is.matrix(a)) {
    a <- apply(a, 2, function(x) x / sqrt(sum(x^2)))
  } else {
    a <- a / sqrt(sum(a^2))
  }
  return(a)
}

inside_detection <- function(Span, vector) {
  lambda <- solve(Span, vector)
  if (sum(lambda >= -1e-10) == length(lambda)) {
    return(1)
  } else {
    return(0)
  }
}

# function that computes all the extreme points that belong to original vertexes
inside_vertex_detection <- function(A, B) {
  SpanA <- norm2(A)
  SpanB <- norm2(B)
  # to determine whether a vertex of one cone is inside another cone or not.

  inside_vertex <- list()
  l <- 1
  for (i in 1:ncol(B)) {
    auxi <- inside_detection(SpanA, SpanB[, i])
    if (auxi == 1) {
      inside_vertex[[l]] <- SpanB[, i]
      l <- l + 1
    }
  }
  for (i in 1:ncol(A)) {
    auxi <- inside_detection(SpanB, SpanA[, i])
    if (auxi == 1) {
      inside_vertex[[l]] <- SpanA[, i]
      l <- l + 1
    }
  }
  return(inside_vertex)
}

# function that computes all the extreme points generated that are generated by the intersections of the cones
intersection_vertex_detection <- function(S, M) {
  num <- ncol(S)
  if (num == 2) {
    return(list())
  } else {
    combination_S <- combn(1:ncol(S), 2)
    combination_M <- combn(1:ncol(S), (num - 1))
    Span_S <- norm2(S)
    Span_M <- norm2(M)

    border_M <- list()
    extreme_point_M <- list()
    for (i in 1:ncol(M)) {
      coeff_matrix <- matrix(0, ncol = num, nrow = num - 1)
      for (j in 1:(num - 1)) {
        coeff_matrix[j, ] <- Span_M[, combination_M[j, i]]
      }
      if (Rank(coeff_matrix) == (num - 1)) {
        border_M[[i]] <- nullspace(coeff_matrix)
      } else {
        stop("border_M_i is denegerated. Check the dimension.")
      }
      extreme_point_M[[i]] <- t(coeff_matrix)[1:(num - 1), 1:(num - 1)]
    }

    inside_face_detection <- function(extreme_point, test_vector) {
      lambda <- solve(extreme_point, test_vector)
      if (sum(lambda >= -1e-10) == length(lambda)) {
        return(1)
      } else {
        return(0)
      }
    }

    l <- 1
    intersection_vertex <- list()
    side <- c()
    for (i in 1:ncol(combination_S)) {
      vertex_1 <- Span_S[, combination_S[1, i]]
      vertex_2 <- Span_S[, combination_S[2, i]]
      for (j in 1:length(border_M)) {
        n1 <- sum(vertex_1 * border_M[[j]])
        n2 <- sum(vertex_2 * border_M[[j]])

        auxi <- n1 * n2
        if (auxi < -1e-10) {
          lambda <- n2 / (n2 - n1)
          possible <- lambda * vertex_1 + (1 - lambda) * vertex_2

          if (det(extreme_point_M[[j]]) != 0) {
            auxi2 <- inside_face_detection(extreme_point_M[[j]], possible[1:(num - 1)])
            if (auxi2 == 1) {
              intersection_vertex[[l]] <- possible
              side[l] <- j
              l <- l + 1
            }
          }
        }
      }
    }

    if (length(intersection_vertex) > 0) {
      for (i in 1:length(intersection_vertex)) {
        intersection_vertex[[i]] <- norm2(intersection_vertex[[i]])
      }
    }
  }

  return(intersection_vertex)
}

# function that computes all the extreme points
vertex_detection <- function(A, B) {
  num <- ncol(A)
  inside_vertex <- inside_vertex_detection(A, B)
  intersection_vertex <- intersection_vertex_detection(A, B)

  # combine the two vertex lists
  if (length(inside_vertex) > 0) {
    vertex <- matrix(unlist(inside_vertex), nrow = num, byrow = FALSE)
  } else {
    vertex <- matrix(0, nrow = num, ncol = 2)
  }
  if (length(intersection_vertex) > 0) {
    vertex <- cbind(vertex, matrix(unlist(intersection_vertex), nrow = num, byrow = FALSE))
  }

  # delete the points that are nonzero due to numerical error
  delete_zeroes <- c()
  for (i in 1:ncol(vertex)) {
    if (near(sum(vertex[, i]^2), 0)) {
      delete_zeroes <- c(delete_zeroes, i)
    }
  }
  if (length(delete_zeroes) > 0) vertex <- vertex[, -delete_zeroes]


  # delete the same ones
  if (length(vertex) > num) {
    for (test in 1:ncol(vertex)) {
      vertex[, test] <- norm2(vertex[, test])
    }
    delete_duplicates <- c()
    for (i in 1:(ncol(vertex) - 1)) {
      for (j in (i + 1):ncol(vertex)) {
        if (sum(near(vertex[, i], vertex[, j])) == nrow(vertex)) {
          delete_duplicates <- c(delete_duplicates, j)
        }
      }
    }
    if (length(delete_duplicates) > 0) vertex <- vertex[, -unique(delete_duplicates)]
  }
  return(vertex)
}




# function that generates a table of presence/absence combinations
# params: N = number of species, nv = initial notation for absence
# return: the table of all possible communities (as a matrix)
get_compo <- function(N, nv = 0) {
  record <- matrix(nv, nrow = 2^N, ncol = N)
  k <- 2
  for (s in 1:N){
    for (i in 1:choose(N, s)){
      record[k, utils::combn(N, s)[, i]] <- 1
      k <- k + 1
    }
  }
  return(record)
}

# Solves the system of ordinary differential equations
# given by the generalized Lotka-Volterra dynamics and
# returns the state variables over time 
# inputs: N0 = vector of initial population sizes (initial condition)
# r = vector of intrinsic growth rates
# K = vector of carrying capacities
# A = square interaction matrix
# times: sequence of times points to numerically integrate ODE
lotka_volterra <- function(N0, r, K, A, times = seq(0, 200, 0.01), 
                           formalism = "r", final_abundance = FALSE) {
  if (formalism == "r") {
    # list of parameters
    pars <- list(r = r, A = A)
    # function that returns rate of change
    model <- function(t, N0, pars) {
      dN_dt <- N0 * (pars$r + c(pars$A %*% N0))
      return(list(dN_dt))
    }
    # numerical integration ode45 Runge-Kutta method
    out <- ode(y = N0, times = times, func = model, parms = pars, method = "ode45")
  }
  if (formalism == "K") {
    pars <- list(r = r, K = K, A = A)
    model <- function(t, N0, pars) {
      dN_dt <- N0 * (pars$r / pars$K) * (pars$K + c(pars$A %*% N0))
      return(list(dN_dt))
    }
    out <- ode(y = N0, times = times, func = model, parms = pars, method = "ode45")
  }
  if (formalism == "K_typeII") {
    pars <- list(r = r, K = K, A = A)
    model <- function(t, N0, pars) {
      dN_dt <- pars$r * N0 * (1 + (1 / pars$K) * c(pars$A %*% diag(1 / (1 + N0)) %*% N0))
      return(list(dN_dt))
    }
    out <- ode(y = N0, times = times, func = model, parms = pars, method = "ode45")
  }
  if (formalism == "K_stochastic") {
    pars <- list(r = r, K = K, A = A)
    # defining deterministic part
    f <- function(u, p, t) {
      deterministic <- u * (p$r / p$K) * (p$K + c(p$A %*% u))
      return(deterministic)
    }
    # defining stochastic part
    g <- function(u, p, t) {
      s <- rep(1 / sqrt(length(p$K)), length(p$K))
      stochastic <- s * u * (p$r / p$K) * (p$K - c(p$A %*% u))
      return(stochastic)
    }
    # integration time steps
    time_step <- times[2] - times[1]
    # numerical integration
    sol <- sde.solve(f = f, g = g, u0 = N0, tspan = range(times), p = pars, saveat = time_step)
    out <- as.data.frame(cbind(sol$t, sol$u))
    names(out) <- paste("time", 1:length(pars$K))
  }
  out <- out[complete.cases(out),] # delete the NA cases (last row)
  if(final_abundance){
    out <- as.numeric(out[nrow(out), -1])
  }
  return(out)
}

# Function that solves the Lotka-Volterra dynamics and returns the 
# indexes of surviving species (i.e., abundance larger than threshold)
# N0: vector of initial population sizes
# r: vector of intrinsic growth rates
# K: vector of carrying capacities
# A: square interaction matrix
# times: sequence of times points to numerically integrate ODE
# extinct_tol: species with a final population size smaller than
# this value are considered extinct
# formalism: whether to use r or K Lotka-Volterra formalism
lv_pruning <- function(N0, r, K, A, times = seq(0, 200, 0.01), 
                        formalism = "r", extinct_tol = 1e-6) {
  # solve Lotka-Volterra dynamics
  eq <- lotka_volterra(N0, r, K, A, times, formalism, final_abundance = TRUE)
  # return pruned system
  which_surv_sp <- as.numeric(which(eq > extinct_tol))
  return(which_surv_sp)
}


# nolint end