\name{AdjBCD.sim}
\alias{AdjBCD.sim}
\title{
Covariate-adjusted Biased Coin Design with Covariate Data Generating Mechanism
}
\description{
Allocates patients to one of two treatments based on the covariate-adjusted biased coin design as proposed by Baldi Antognini A, Zagoraiou M (2011) <Doi:10.1093/biomet/asr021>, by simulating the covariates-profile under the assumption of independence between covariates and levels within each covariate.  
}
\usage{
AdjBCD.sim(n = 1000, cov_num = 2, level_num = c(2, 2), 
           pr = rep(0.5, 4), a = 2)
}
\arguments{
  \item{n}{the number of patients. The default is \code{1000}.}
  \item{cov_num}{the number of covariates. The default if \code{2}.}
  \item{level_num}{a vector of level numbers for each covariate. Hence the length of \code{level_num} should be equal to the number of covariates. The default is \code{c(2,2)}.}
  \item{pr}{a vector of probabilities. Under the assumption of independence between covariates, \code{pr} is a vector containing probabilities for each levels of each covariates. The length of \code{pr} should correspond to the number of all levels, and the vector sum of \code{pr} should be equal to \code{cov_num}. The default is \code{pr = rep(0.5, 4)}, which implies that \code{cov_num = 2}, and \code{level_num = c(2, 2)}.}
  \item{a}{a design parameter. The default is \code{2}. As \code{a} goes to \eqn{\infty}, the design becomes more deterministic.}
}
\details{
See \code{\link{AdjBCD}}.
}
\value{
See \code{\link{AdjBCD}}. 
}
\seealso{
See \code{\link{AdjBCD}} for allocating patients with complete covariate data; See \code{\link{AdjBCD.ui}} for the command-line user interface. 
}