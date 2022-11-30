\name{StrBCD.sim}
\alias{StrBCD.sim}
\title{
Shao's Method in the Two-Arms Case with Covariate Data Generating Mechanism
}
\description{
Allocates patients to one of two treatments using Shao's method proposed by Shao J, Yu X, Zhong B (2010) <doi:10.1093/biomet/asq014>, by simulating covariate profiles under the assumption of independence between covariates and levels within each covariate.
}
\usage{
StrBCD.sim(n = 1000, cov_num = 2, level_num = c(2, 2), 
           pr = rep(0.5, 4), p = 0.85)
}
\arguments{
  \item{n}{the number of patients. The default is \code{1000}.}
  \item{cov_num}{the number of covariates. The default is \code{2}.}
  \item{level_num}{a vector of level numbers for each covariate. Hence the length of \code{level_num} should be equal to the number of covariates. The default is \code{c(2, 2)}. }
  \item{pr}{a vector of probabilities. Under the assumption of independence between covariates, \code{pr} is a vector containing probabilities for each level of each covariate. The length of \code{pr} should correspond to the number of all levels, and the sum of the probabilities for each margin should be \code{1}. The default is \code{rep(0.5, 4)}, which corresponds to \code{cov_num = 2}, and \code{level_num = c(2, 2)}.}
  \item{p}{the biased coin probability. \code{p} should be larger than \code{1/2} and less than \code{1}. The default is \code{0.85}.}
}
\details{
See \code{\link{StrBCD}}.
}
\value{
See \code{\link{StrBCD}}.
}
\references{
Shao J, Yu X, Zhong B. \emph{A theory for testing hypotheses under covariate-adaptive randomization}[J]. Biometrika, 2010, 97(2): 347-360. 
}
\seealso{
  See \code{\link{StrBCD}} for allocating patients with complete covariate data; See \code{\link{StrBCD.ui}} for the command-line user interface. 
}