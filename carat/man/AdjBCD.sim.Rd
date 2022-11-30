\name{AdjBCD.sim}
\alias{AdjBCD.sim}
\title{
Covariate-adjusted Biased Coin Design with Covariate Data Generating Mechanism
}
\description{
Allocates patients to one of two treatments based on the covariate-adjusted biased coin design as proposed by Baldi Antognini A, Zagoraiou M (2011) <doi:10.1093/biomet/asr021>, by simulating the covariates-profile under the assumption of independence between covariates and levels within each covariate.  
}
\usage{
AdjBCD.sim(n = 1000, cov_num = 2, level_num = c(2, 2), 
           pr = rep(0.5, 4), a = 3)
}
\arguments{
  \item{n}{the number of patients. The default is \code{1000}.}
  \item{cov_num}{the number of covariates. The default is \code{2}.}
  \item{level_num}{a vector of level numbers for each covariate. Hence the length of \code{level_num} should be equal to the number of covariates. The default is \code{c(2,2)}.}
  \item{pr}{a vector of probabilities. Under the assumption of independence between covariates, \code{pr} is a vector containing probabilities for each level of each covariate. The length of \code{pr} should correspond to the number of all levels, and the sum of the probabilities for each margin should be \code{1}. The default is \code{rep(0.5, 4)}, which corresponds to \code{cov_num = 2}, and \code{level_num = c(2, 2)}.}
  \item{a}{a design parameter governing the degree of randomness. The default is \code{3}.}
}
\details{
See \code{\link{AdjBCD}}.
}
\value{
See \code{\link{AdjBCD}}. 
}
\references{
Baldi Antognini A, Zagoraiou M. \emph{The covariate-adaptive biased coin design for balancing clinical trials in the presence of prognostic factors}[J]. Biometrika, 2011, 98(3): 519-535.
}
\seealso{
See \code{\link{AdjBCD}} for allocating patients with complete covariate data; See \code{\link{AdjBCD.ui}} for the command-line user interface. 
}