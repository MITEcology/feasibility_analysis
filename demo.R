rm(list = ls())
source("toolbox.R")

set.seed(2022)
matA <- generate_inte_rand(4, 1, 1, 0, "norm")
set.seed(2023)
matB <- generate_inte_rand(4, 1, 1, 0, "norm")
set.seed(2024)
r <- matA %*% c(runif(4, 0, 1))

set.seed(2025)
matC <- generate_inte_neg_def(4, 0.8, 0.8, "norm")
set.seed(2026)
matD <- generate_inte_neg_def(4, 0.8, 0.8, "norm")
# save(matA, matB, matC, matD, r, file = "testdata.RData")

feasibility_community(matC)
# [1] 0.2682132

feasibility_overlap(matA, matB)
# [1] 0.0001663809

feasibility_asymmetry(matA)
# [1] 0.4558857

feasibility_partition(matD)
# [1] 0.062500000 0.090493483 0.130122337 0.111260296 0.085156470 0.061411663 0.148784899 0.080082983 0.034231896 0.032469443 0.059323957 0.017098322 0.017826845
# [14] 0.057302734 0.007891408 0.004046566

feasibility_species(matC, sp = 2)
# [1] 0.7020132

feasibility_resistance_full(matA, r)
# results might be random due to sampling of the minimum
#      1_2_3      1_2_4      1_3_4      2_3_4
# 0.13204785 0.26007696 0.08515493 0.25073912

feasibility_resistance_partial(matA, r)
#         1         2         3         4
# 0.2658463 0.9230088 1.2967255 0.7911352

feasibility_recovery(matA, r, type = "full")
#     feasibility_recovery(matA, r, type = "full")
# [1] -1.441499

feasibility_recovery(matA, r, type = "part")
# [1] -0.6845685

feasibility_contribution(matC, sp = 2)
# 1.241339

feasibility_invasion(matC, invader = 1)
# $colonization
# [1] 0.5295405

# $augmentation
# [1] 0.1123268

# $replacement
# [1] 0.4172137