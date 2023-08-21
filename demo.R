# nolint start

source("toolbox.R")

set.seed(2022)
matA <- generate_inte_rand(4, 1, 1, 0, "norm")
set.seed(2023)
matB <- generate_inte_rand(4, 1, 1, 0, "norm")
set.seed(2024)
r <- matA %*% c(runif(4, 0, 1))


calculate_feas_commun(matA)

calculate_feas_overlap(matA, matB)

calculate_asymmetry(matA)

calculate_feas_partition(matA)

calculate_feas_sp_survival(matA, sp = 2)

calculate_resistance_full(matA, r)
calculate_resistance_partial(matA, r)

calculate_recovery(matA, r, type = "full")

calculate_contribution(matA, sp = 1)

calculate_feas_sp_invasion(matA, invader = 1)

# nolint end