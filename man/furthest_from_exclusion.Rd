% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/furthest_from_exclusion.R
\name{furthest_from_exclusion}
\alias{furthest_from_exclusion}
\title{Parameters of the species furthest away from exclusion}
\usage{
furthest_from_exclusion(center, A_int)
}
\arguments{
\item{center}{a S vector that contains the initial intrinsic growth
rates that are used to computed the minimal isotropic variation that
excludes one or more species from the community.}

\item{A_int}{a SxS interaction matrix, where S is the number of species.}
}
\value{
A vector whose length varies between 1 and S. It element contains the
column index of the species that is/are furthest away from exclusion.
}
\description{
Given an interaction matrix and an initial vector of intrinsic growing rates,
denoted by center, this function returns the column number/s of the
interaction matrix with the species that is/are further away from exclusion.
}
\examples{
A_int <- -1*diag(c(3,5,7))
incenter_inradius_isoprob <- incenter_inradius_isoprob_calculation(A_int)
I <- incenter_inradius_isoprob[[1]]
center <- I + runif(nrow(A_int))
furthest_from_exclusion(center, A_int)

}
