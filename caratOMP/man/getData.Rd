\name{getData}
\alias{getData}
\title{Data Generation}
\usage{
  getData(n, cov_num, level_num, pr, type, beta, 
          mu1, mu2, sigma = 1, method = HuHuCAR, \dots)
}
\description{
  Generates continuous or binary outcomes given patients' covariates, underlying model and randomization procedure.
}
\arguments{
  \item{n}{number of patients.}
  \item{cov_num}{number of covariates.}
  \item{level_num}{vector of level numbers for each covariates. Hence the length of level_num should equal to number of covariates.}
  \item{pr}{vector of probabilities. Under assumption of independence between covariates, \code{pr} is a vector containing probabilities for each levels of each covariates. The length of pr should correspond to number of all levels; and sum of pr should equal \code{cov_num}.}
  \item{type}{way of generating data. Optional input: \code{linear} or \code{logit}.}
  \item{beta}{vector of coefficients of covariates. The length of \code{beta} must correspond to \code{cov_num}.}
  \item{mu1,mu2}{main effects of treatment \code{1} and treatment \code{2}.}
  \item{sigma}{error variance for linear model. Default is 1. It is only used when \code{type} is \code{linear}.}
  \item{method}{the randomization to be used in allocating patients. The default randomization \code{HuHuCAR} uses Hu and Hu's general covariate-adaptive randomization; the alternatives are \code{PocSimMIN}, \code{StrBCD}, \code{StrPBR}, \code{DoptBCD} as well as \code{AdjBCD}.}
  \item{\dots}{arguments to be passed to methods. It depends on the \code{method} used and the following arguments are accepted:
  \describe{
  \item{omega}{vector of weights at the overall, within-stratum as well as maginal levels. It is required that at least one element is larger than 0. Notice that \code{omega} is only needed when \code{HuHuCAR} is to be used.}
      \item{weight}{vector of weights for maginal imbalances. It is required that at least one element is larger than 0. Notice that \code{weight} is only needed when \code{PocSimMIN} is to be used.}
      \item{p}{probabilities of assigning one patient to treatment 1. p should be larger than \code{1/2} to obtain balance. Notice that \code{p} is only needed when \code{"HuHuCAR", "PocSimMIN"} and \code{"StrBCD"} are to be used.}
      \item{a}{a design parameter. As \code{a} goes to \eqn{\infty}, the design becomes more deterministic.}
  \item{bsize}{block size for stratified randomization. It is required to be multiple of 2. Notice that \code{bsize} is only needed when \code{"StrPBR"} is to be used.}
  }}
}

\details{
  To generate continuous outcomes, we use the linear model:\deqn{y_i = \mu_j+x_i^T\beta+\epsilon_i,}
  
  to generate binary outcomes, we use the logit link function:\deqn{P(y_i=1) = \frac{exp\{\mu_j+x_i^T\beta \}}{1+exp \{\mu_j+x_i^T\beta }},
  
  where \eqn{j} indicates patient \eqn{i} belongs to treatment \eqn{j}.
}

\value{
  \code{getData} returns a size \eqn{cov_num+2 \times n} dataframe. The first \code{cov_num} rows represent patients' profile. The next row consists patients' assignments and the final row consists of generated outcomes.
}

\examples{
#Parameters' Setting
set.seed(100)
n = 1000
cov_num = 5
level_num = c(2,2,2,2,2)
beta = c(1,4,3,2,5)
mu1 = 0
mu2 = 0
sigma = 1
type = "linear"
method = HuHuCAR
p = 0.85
omega = c(0.1, 0.1, rep(0.8 / 5, times = 5))
pr = rep(0.5,10)

#Data Generation
dataH = getData(n, cov_num,level_num, pr, type, beta,
                mu1, mu2, sigma, HuHuCAR, omega, p)
dataH[1:(cov_num+2),1:5]
}
