\name{HuHuCAR.sim}
\alias{HuHuCAR.sim}
\title{Hu and Hu's General Covariate-Adaptive Randomization with Covariate Data Generating Mechanism}
\description{
Allocates patients to one of two treatments using general covariate-adaptive randomization proposed by Hu Y, Hu F (2012) <doi:10.1214/12-AOS983>, by simulating covariate profiles based on the assumption of independence between covariates and levels within each covariate.
}
\usage{
HuHuCAR.sim(n = 1000, cov_num = 2, level_num = c(2, 2), 
            pr = rep(0.5, 4), omega = NULL, p = 0.85)
}
\arguments{
  \item{n}{the number of patients. The default is \code{1000}.}
  \item{cov_num}{the number of covariates. The default is \code{2}.}
  \item{level_num}{a vector of level numbers for each covariate. Hence the length of \code{level_num} should be equal to the number of covariates. The default is \code{c(2, 2)}.}
  \item{pr}{a vector of probabilities. Under the assumption of independence between covariates, \code{pr} is a vector containing probabilities for each level of each covariate. The length of \code{pr} should correspond to the number of all levels, and the sum of the probabilities for each margin should be \code{1}. The default is \code{rep(0.5, 4)}, which corresponds to \code{cov_num = 2}, and \code{level_num = c(2, 2)}.}
  \item{omega}{a vector of weights at the overall, within-stratum, and within-covariate-margin levels. It is required that at least one element is larger than 0. If \code{omega = NULL} (default), the overall, within-stratum, and within-covariate-margin imbalances are weighted with porportions \code{0.2}, \code{0.3}, and \code{0.5/cov_num} for each covariate-margin, respectively, where \code{cov_num} is the number of covariates of interest.}
  \item{p}{the biased coin probability. \code{p} should be larger than \code{1/2} and less than \code{1}. The default is \code{0.85}.}
}
\details{
See \code{\link{HuHuCAR}}.
}
\value{
See \code{\link{HuHuCAR}}. 
}
\references{
Hu Y, Hu F. \emph{Asymptotic properties of covariate-adaptive randomization}[J]. The Annals of Statistics, 2012, 40(3): 1794-1815.
}
\seealso{
  See \code{\link{HuHuCAR}} for allocating patients with complete covariate data; See \code{\link{HuHuCAR.ui}} for the command-line user interface. 
}