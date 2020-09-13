\name{HuHuCAR.sim}
\alias{HuHuCAR.sim}
\title{Hu and Hu's General Covariate-Adaptive Randomization with Covariate Data Generating Mechanism}
\description{
Allocates patients to one of two treatments using general covariate-adaptive randomization proposed by Hu Y, Hu F (2012) <Doi:10.1214/12-AOS983>, by simulating covariate profiles based on the assumption of independence between covariates and levels within each covariate.
}
\usage{
HuHuCAR.sim(n = 1000, cov_num = 2, level_num = c(2, 2), 
            pr = rep(0.5, 4), omega = NULL, p = 0.85)
}
\arguments{
  \item{n}{the number of patients. The default is \code{1000}.}
  \item{cov_num}{the number of covariates. The default is \code{2}.}
  \item{level_num}{the vector of level numbers for each covariate. Hence the length of \code{level_num} should be equal to the number of covariates. The default is \code{c(2, 2)}.}
  \item{pr}{the vector of probabilities. Under the assumption of independence between covariates, \code{pr} is a vector containing probabilities for each level of each covariate. The length of \code{pr} should correspond to the number of all levels, and the vector sum of \code{pr} should be equal to \code{cov_num}. The default is \code{pr = rep(0.5, 4)}, which implies that \code{cov_num = 2}, and \code{level_num = c(2, 2)}.}
  \item{omega}{the vector of weights at the overall, within-stratum, and maginal levels. It is required that at least one element is larger than 0. If \code{omega = NULL} (default), it weights the overall, within-stratum as well as marginal levels with porportion \code{1/cov_num}.}
  \item{p}{the probability of assigning one patient to treatment 1. \code{p} should be larger than \code{1/2} to obtain balance. The default is \code{0.85}.}
}
\details{
See \code{\link{HuHuCAR}}.
}
\value{
See \code{\link{HuHuCAR}}. 
}
\seealso{
  See \code{\link{HuHuCAR}} for allocating patients with complete covariate data; See \code{\link{HuHuCAR.ui}} for the command-line user interface. 
}