\name{PocSimMIN.sim}
\alias{PocSimMIN.sim}
\title{
Pocock and Simon's Method in the Two-Arms Case with Covariate Data Generating Mechanism
}
\description{
Allocates patients to one of two treatments using Pocock and Simon's method proposed by Pocock S J, Simon R (1975) <doi:10.2307/2529712>, by simulating covariate profiles under the assumption of independence between covariates and levels within each covariate.
}
\usage{
PocSimMIN.sim(n = 1000, cov_num = 2, level_num = c(2, 2), 
              pr = rep(0.5, 4), weight = NULL, p = 0.85)
}
\arguments{
  \item{n}{the number of patients. The default is \code{1000}.}
  \item{cov_num}{the number of covariates. The default is \code{2}.}
  \item{level_num}{a vector of level numbers for each covariate. Hence the length of \code{level_num} should be equal to the number of covariates. The default is \code{c(2, 2)}. }
  \item{pr}{a vector of probabilities. Under the assumption of independence between covariates, \code{pr} is a vector containing probabilities for each level of each covariate. The length of \code{pr} should correspond to the number of all levels, and the sum of the probabilities for each margin should be \code{1}. The default is \code{rep(0.5, 4)}, which corresponds to \code{cov_num = 2}, and \code{level_num = c(2, 2)}.}
  \item{weight}{a vector of weights for within-covariate-margin imbalances. It is required that at least one element is larger than 0. If \code{weight = NULL} (default), the within-covariate-margin imbalances are weighted with an equal proportion, \code{1/cov_num}, for each covariate-margin.}
  \item{p}{the biased coin probability. \code{p} should be larger than \code{1/2} and less than \code{1}. The default is \code{0.85}.}
}
\details{
See \code{\link{PocSimMIN}}.
}
\value{
See \code{\link{PocSimMIN}}. 
}
\references{
Pocock S J, Simon R. \emph{Sequential treatment assignment with balancing for prognostic factors in the controlled clinical trial}[J]. Biometrics, 1975: 103-115.
}
\seealso{
  See \code{\link{PocSimMIN}} for allocating patients with complete covariate data; See \code{\link{PocSimMIN.ui}} for the command-line user interface. 
}