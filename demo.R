rm(list = ls())
source("toolbox.R") ## all the main functions

### Feasibility: the size of the feasibility domain compatible with the long-term persistence of a given set of species $\mathcal{C}$ from a given community (pool) of species $\mathcal{S}$.
### The size is estimated by the normalized surface area of a unit sphere (i.e., also known as the normalized solid angle).
### Note that feasibility (size of the feasibility domain) is not directly a probability measure.
### To do so, one needs to normalize the feasibility such that values lie within $[0,1]$ according to the characteristics of the community (see section feasibility of a community $\mathcal{S}$).

matA <- generate_inte_rand(4, 1, 1, 0, "norm")  ## Generate a random matrix
matB <- generate_inte_gs(4, 1, 1, 0, "norm") ## Generate a globally stable random matrix
r <- matA %*% c(runif(4, 0, 1)) ## This is the direction of the r-vector (assuming it is known). This is only needed to calculate resistance and recovery

feasibility_community(matA) ## If the user is interested in the species-specif measure, then the user needs to calculate raise this outcome to the power 1/|S|, where |S| is the matrix dimension
feasibility_community(matB)

feasibility_overlap(matA, matB) ## Matrices need to have the same dimension

feasibility_asymmetry(matA) ## This is a standard deviation. The larger the outcome, the larger the assymetry

feasibility_partition(matB) ## The feasibility of all the 2^|S| possible combinations

feasibility_species(matB, sp = 2) ## The sum of feasibility regions whith species i

feasibility_resistance_full(matA, r) ## This is a distance from a a given r-vector to all possible |S-1| regions. The larger the outcome, the larger the resistance

feasibility_resistance_partial(matA, r) ## This is a distance from a a given r-vector to all possible vertices. The larger the outcome, the larger the resistance

feasibility_recovery(matA, r, type = "full") ## This is the largest eigenvalue (notice is negative)

feasibility_recovery(matA, r, type = "part") ## This is the second smallest eigenvalue (notice is negative)

feasibility_contribution(matB, sp = 2) ## This is a ratio between the feasibility with |S| species and with |S-1| species. Otcomes above (resp. below) 1 mean a positive (resp. negative) contribution

feasibility_invasion(matB, invader = 1) ## This is the feasibility region for invasion
