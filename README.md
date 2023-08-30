[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8289566.svg)](https://doi.org/10.5281/zenodo.8289566)


# feasibilityR

The goal of feasibilityR is to is to facilitate the analysis of feasibility within the largest framework of structural stability in ecology.

Feasibility (Omega): the size of the feasibility domain compatible with the long-term persistence of a given set of species C from a given community (pool) of species S. The size is estimated by the normalized surface area of a unit sphere (i.e., also known as the normalized solid angle). Note that feasibility (size of the feasibility domain) is not directly a probability measure. To do so, one needs to normalize the feasibility such that values lie within [0,1] according to the characteristics of the community (see section feasibility of a community S).

Formally, the feasibility domain of any set C_S under generalized Lotka-Volterra (gLV) dynamics can be expressed as Deng et al 2022.   

*D*(𝒞) = {∑<sub>*j* ∈ 𝒞</sub>*N*<sub>*j*</sub><sup>\*</sup>**a**<sub>*j*</sub>−∑<sub>*k* ∈ 𝒮/𝒞</sub>*N*<sub>*k*</sub><sup>\*</sup>**e**<sub>*k*</sub>:    

**A**=**a**<sub>1</sub> **a**<sub>2</sub> ⋯ **a**<sub>*S*</sub>,*N*<sub>*j*</sub><sup>\*</sup>,*N*<sub>*k*</sub><sup>\*</sup>∈ℝ<sub> \> 0</sub>} ∈ ℝ<sup>*S*</sup>.  

**Note: I am using this to convert latex to md files: https://alldocs.app/convert-latex-to-markdown**

## Installation

``` r
# install.packages("devtools")
devtools::install_github("MITEcology/feasibility_analysis")
```

## Example

This is a basic example:

``` r
library(feasibilityR)

#First we can generate some random matrices
matA <- generate_inte_rand(4, 1, 1, 0, "norm")  ## Generate a random matrix
matB <- generate_inte_gs(4, 1, 1, 0, "norm") ## Generate a globally stable random matrix
#And a known r vector 
r <- matA %*% c(runif(4, 0, 1)) ## This is the direction of the r-vector (assuming it is known). This is only needed to calculate resistance and recovery

#Now we are ready to calculate the feasabilty domain of this random communities.
feasibility_community(matA) ## If the user is interested in the species-specif measure, then the user needs to calculate raise this outcome to the power 1/|S|, where |S| is the matrix dimension.
feasibility_community(matB)

#etc...
```

## Citation

If using this package, please cite it:

``` r
citation("feasibilityR")
```

## Acknowledgements

TBA

## To Do:
- Create tests using testthat
- Create a vignette
