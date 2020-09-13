\name{StrPBR.sim}
\alias{StrPBR.sim}
\title{
Stratified Permuted Block Randomization with Covariate Data Generating Mechanism
}
\description{
Allocates patients to one of two treatments using stratified randomization proposed by Zelen M (1974) <Doi: 10.1016/0021-9681(74)90015-0>, by simulating covariates-profile on assumption of independence between covariates and levels within each covariate.
}
\usage{
StrPBR.sim(n = 1000, cov_num = 2, level_num = c(2, 2), 
           pr = rep(0.5, 4), bsize = 4)
}
\arguments{
  \item{n}{the number of patients. The default is \code{1000}.}
  \item{cov_num}{the number of covariates. The default is \code{2}.}
  \item{level_num}{the vector of level numbers for each covariate. Hence the length of \code{level_num} should be equal to the number of covariates. The default is \code{c(2, 2)}. }
  \item{pr}{the vector of probabilities. Under the assumption of independence between covariates, \code{pr} is a vector containing probabilities for each level of each covariate. The length of \code{pr} should correspond to  the number of all levels, and the vector sum of \code{pr} should be equal to \code{cov_num}. The default is \code{pr = rep(0.5, 4)} (default), which implies that \code{cov_num = 2} and \code{level_num = c(2, 2)}.}
  \item{bsize}{the block size for the stratified randomization. It is required to be a multiple of 2. The default is \code{4}.}
}
\details{
See \code{\link{StrPBR}}.
}
\value{
See \code{\link{StrPBR}}.
}
\seealso{
  See \code{\link{StrPBR}} for allocating patients with complete covariate data; See \code{\link{StrPBR.ui}} for the command-line user interface. 
}