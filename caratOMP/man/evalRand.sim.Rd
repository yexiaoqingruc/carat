\name{evalRand.sim}
\alias{evalRand.sim}
\alias{evalRand.sim.careval}
\title{
Evaluation Randomization Procedures with Covariate Data Generating Mechanism
}
\description{
Evaluates randomization procedure based on several different quantities of imbalances by simulating patients' covariate-profile on the assumption of independence between covariates and levels within each covariate.
}
\usage{
\method{evalRand.sim}{careval}(n = 1000, N = 500, Replace = FALSE, cov_num = 2, 
             level_num = c(2, 2), pr = rep(0.5, 4), method = "HuHuCAR", ...)
}
\arguments{
  \item{N}{iteration number.}
  \item{n}{number of patients. Default is \code{1000}.}
  \item{Replace}{bool. If \code{Replace = FALSE}, the function does clinical trial design for \code{N} iterations for one group of patients. If \code{Replace = TRUE}, the function dose clinical trial design for \code{N} iterations for \code{N} different groups of patients. }
  \item{cov_num}{number of covariates. Default is \code{2}.}
  \item{level_num}{vector of level numbers for each covariates. Hence the length of \code{level_num} should equal to number of covariates. Default is \code{level_num = c(2, 2)}. }
  \item{pr}{vector of probabilities. Under assumption of independence between covariates, \code{pr} is a vector containing probabilities for each levels of each covariates. The length of pr should correspond to number of all levels; and sum of pr should equal \code{cov_num}. If \code{pr = (0.5, 0.5, 0.5, 0.5)} (default), it means \code{cov_num = 2}, and \code{level_num = c(2, 2)}.}
  \item{method}{the randomization to be used in allocating patients. The default randomization \dQuote{\code{HuHuCAR}} uses Hu and Hu's general covariate-adaptive randomization; the alternatives are \dQuote{\code{PocSimMIN}}, \dQuote{\code{StrBCD}}, \dQuote{\code{StrPBR}}, \dQuote{\code{DoptBCD}} as well as \dQuote{\code{AdjBCD}}.}
  \item{\dots}{arguments to be passed to methods. It depends on the \code{method} and the following arguments are accepted:
  \describe{
  \item{omega}{vector of weights at the overall, within-stratum as well as maginal levels. It is required that at least one element is larger than 0. Notice that \code{omega} is only needed when \code{HuHuCAR} are to be assessed.}
  \item{weight}{vector of weights for all involved margins. It is required that at least one element is NOT 0 and \code{length(weight) = cov_num}. Notice that \code{weight} is only needed when \code{PocSimMIN} is to be assessed.}
  \item{p}{probabilities of assigning one patinet to treatment 1. \code{p} should be larger than \code{1/2} to obtain balance. Notice that \code{p} is only needed when \code{"HuHuCAR"},  \code{"PocSimMIN"} and \code{"StrBCD"} is to be assessed.}
  \item{a}{a design parameter. As \code{a} goes to \eqn{\infty}, the design becomes more deteministic. Notice that \code{a} is only needed when \code{"AdjBCD"} is to be assessed.}
  \item{bsize}{block size for stratified permuted block randomization. It should be multiple of \code{2}. Notice that \code{bsize} is only needed when \code{"StrPBR"} is to be assessed.}
  }}
}
\details{
See \code{\link{evalRand}}.
}
\value{
See \code{\link{evalRand}}.
}
\seealso{
  See \code{\link{evalRand}} to evaluate a randomization procedure with complete covariate data.
}