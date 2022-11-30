\name{getData}
\alias{getData}
\title{Data Generation}
\usage{
  getData(n, cov_num, level_num, pr, type, beta, 
          mu1, mu2, sigma = 1, method = "HuHuCAR", \dots)
}
\description{
  Generates continuous or binary outcomes given patients' covariates, the underlying model and the randomization procedure.
}
\arguments{
  \item{n}{the number of patients.}
  \item{cov_num}{the number of covariates.}
  \item{level_num}{a vector of level numbers for each covariate. Hence the length of \code{level_num} should be equal to the number of covariates.}
  \item{pr}{a vector of probabilities. Under the assumption of independence between covariates, \code{pr} is a vector containing probabilities for each level of each covariate. The length of \code{pr} should correspond to the number of all levels, and the sum of the probabilities for each margin should be \code{1}.}
  \item{type}{a data-generating method. Optional input: \code{"linear"} or \code{"logit"}.}
  \item{beta}{a vector of coefficients of covariates. The length of \code{beta} must correspond to the sum of all covariates' levels.}
  \item{mu1,mu2}{main effects of treatment \code{1} and treatment \code{2}.}
  \item{sigma}{the error variance for the linear model. The default is 1. This should be a positive value and is only used when \code{type = linear}.}
  \item{method}{the randomization procedure to be used for generating randomization sequences. This package provides data-generating function for \code{"HuHuCAR"}, \code{"PocSimMIN"}, \code{"StrBCD"}, \code{"StrPBR"}, \code{"AdjBCD"}, and \code{"DoptBCD"}.}
  \item{\dots}{arguments to be passed to \code{method}. These arguments depend on the randomization method used and the following arguments are accepted:
  \describe{
  \item{omega}{a vector of weights at the overall, within-stratum, and within-covariate-margin levels. It is required that at least one element is larger than 0. Note that \code{omega} is only needed when \code{HuHuCAR} is to be used.}
      \item{weight}{a vector of weights for within-covariate-margin imbalances. It is required that at least one element is larger than 0. Note that \code{weight} is only needed when \code{PocSimMIN} is to be used.}
      \item{p}{the biased coin probability. \code{p} should be larger than \code{1/2} and less than \code{1}. Note that \code{p} is only needed when \code{"HuHuCAR", "PocSimMIN"} and \code{"StrBCD"} are to be used.}
      \item{a}{a design parameter governing the degree of randomness. Note that \code{a} is only needed when \code{"AdjBCD"} is to be used.}
  \item{bsize}{the block size for stratified randomization. It is required to be a multiple of 2. Note that \code{bsize} is only needed when \code{"StrPBR"} is to be used.}
  }}
}

\details{
  To generate continuous outcomes, we use the linear model:\deqn{y_i = \mu_j+x_i^T\beta+\epsilon_i,}
  
  to generate binary outcomes, we use the logit link function:\deqn{P(y_i=1) = \frac{exp\{\mu_j+x_i^T\beta \}}{1+exp \{\mu_j+x_i^T\beta }},
  
  where \eqn{j} indicates patient \eqn{i} belongs to treatment \eqn{j}.
}

\value{
  \code{getData} returns a size \eqn{cov_num+2 \times n} dataframe. The first \code{cov_num} rows represent patients' profile. The next row consists of patients' assignments and the final row consists of generated outcomes.
}

\examples{
#Parameters' Setting
set.seed(100)
n = 1000
cov_num = 5
level_num = c(2,2,2,2,2)
beta = c(1,4,3,2,5,5,4,3,2,1)
mu1 = 0
mu2 = 0
sigma = 1
type = "linear"
p = 0.85
omega = c(0.1, 0.1, rep(0.8 / 5, times = 5))
pr = rep(0.5,10)

#Data Generation
dataH = getData(n, cov_num,level_num, pr, type, beta,
                mu1, mu2, sigma, "HuHuCAR", omega, p)
dataH[1:(cov_num+2),1:5]
}
