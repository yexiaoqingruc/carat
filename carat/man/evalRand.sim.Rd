\name{evalRand.sim}
\alias{evalRand.sim}
\title{
Evaluation Randomization Procedures with Covariate Data Generating Mechanism
}
\description{
Evaluates randomization procedure based on several different quantities of imbalances by simulating patients' covariate profiles under the assumption of independence between covariates and levels within each covariate.
}
\usage{
evalRand.sim(n = 1000, N = 500, Replace = FALSE, cov_num = 2, 
             level_num = c(2, 2), pr = rep(0.5, 4), method = "HuHuCAR", ...)
}
\arguments{
  \item{N}{the iteration number.}
  \item{n}{the number of patients. The default is \code{1000}.}
  \item{Replace}{bool. If \code{Replace = FALSE}, the function does clinical trial design for \code{N} iterations for one group of patients. If \code{Replace = TRUE}, the function dose clinical trial design for \code{N} iterations for \code{N} different groups of patients. }
  \item{cov_num}{the number of covariates. The default is \code{2}.}
  \item{level_num}{the vector of level numbers for each covariate. Hence the length of \code{level_num} should be equal to the number of covariates. The default is \code{c(2, 2)}. }
  \item{pr}{the vector of probabilities. Under the assumption of independence between covariates, \code{pr} is a vector containing probabilities for each level of each covariate. The length of \code{pr} should correspond to the number of all levels, and the vector sum of \code{pr} should be equal to \code{cov_num}. The default is \code{pr = (0.5, 0.5, 0.5, 0.5)}, which implies that \code{cov_num = 2}, and \code{level_num = c(2, 2)}.}
  \item{method}{the randomization method to be used in allocating patients. The default randomization \dQuote{\code{HuHuCAR}} uses Hu and Hu's general covariate-adaptive randomization; the alternatives are \dQuote{\code{PocSimMIN}}, \dQuote{\code{StrBCD}}, \dQuote{\code{StrPBR}}, \dQuote{\code{DoptBCD}} and \dQuote{\code{AdjBCD}}.}
  \item{\dots}{arguments to be passed to methods. These depends on \code{method}, and the following arguments are accepted:
  \describe{
  \item{omega}{the vector of weights at the overall, within-stratum, and marginal levels. It is required that at least one element is larger than 0. Note that \code{omega} is only needed when \code{HuHuCAR} are to be assessed.}
  \item{weight}{the vector of weights for marginal imbalances. It is required that at least one element is NOT 0 and \code{length(weight) = cov_num}. Note that \code{weight} is only needed when \code{PocSimMIN} is to be assessed.}
  \item{p}{the probability of assigning one patinet to treatment 1. \code{p} should be larger than \code{1/2} to obtain balance. Note that \code{p} is only needed when \code{"HuHuCAR"},  \code{"PocSimMIN"} and \code{"StrBCD"} is to be assessed.}
  \item{a}{a design parameter. As \code{a} goes to \eqn{\infty}, the design becomes more deteministic. Note that \code{a} is only needed when \code{"AdjBCD"} is to be assessed.}
  \item{bsize}{the block size for stratified permuted block randomization. It is required to be a multiple of \code{2}. Note that \code{bsize} is only needed when \code{"StrPBR"} is to be assessed.}
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