\name{evalPower}
\alias{evalPower}
\title{Evaluation of Tests and Randomization Procedures through Power}
\usage{
  evalPower(n, cov_num, level_num, pr, type, beta, di = seq(0,0.5,0.1), sigma = 1,
            Iternum, sl = 0.05, method = HuHuCAR, test, plot = "TRUE", ...)
}
\description{
  Returns powers and a plot of chosen test and method under different treatment effects.
}
\arguments{
  \item{n}{number of patients.}
  \item{cov_num}{number of covariates.}
  \item{level_num}{vector of level numbers for each covariates. Hence the length of level_num should equal to number of covariates.}
  \item{pr}{vector of probabilities. Under assumption of independence between covariates, \code{pr} is a vector containing probabilities for each levels of each covariates. The length of pr should correspond to number of all levels; and sum of pr should equal \code{cov_num}.}
  \item{type}{way of generating data. Optional input: \code{linear} or \code{logit}.}
  \item{beta}{vector of coefficients of covariates. The length of \code{beta} must correspond to \code{cov_num}.}
  \item{di}{vector of values of difference in treatment effects. Default value is a sequence from \code{0} to \code{0.5} and the unit is \code{0.1}.}
  \item{sigma}{error variance for linear model. Default value is 1. It is only used when \code{type} is \code{linear}.}
  \item{Iternum}{integer, the number of iterations to calculate the average power.}
  \item{sl}{significance level, if the p-value returned by the test is less than \code{sl}, we will reject the null hypothesis. Default value is \code{0.05}.}
  \item{method}{the randomization to be used in allocating patients. The default randomization \code{HuHuCAR} uses Hu and Hu's general covariate-adaptive randomization; the alternatives are \code{PocSimMIN}, \code{StrBCD}, \code{StrPBR}, \code{DoptBCD} as well as \code{AdjBCD}.}
  \item{test}{test used to verify hypothesis. Optional input: \code{rand.test}, \code{boot.test} or \code{corr.test}. They are randomization test, bootstrap t-test and corrected t-test respectively.}
  \item{plot}{string, whether to plot or not. Optional input: \code{TRUE} or \code{FALSE}.}
  \item{...}{arguments to be passed to methods. It depends on the \code{method} and \code{test} used and the following arguments are accepted:
  \describe{
      \item{omega}{vector of weights at the overall, within-stratum as well as maginal levels. It is required that at least one element is larger than 0. Notice that \code{omega} is only needed when \code{HuHuCAR} is to be used.}
      \item{weight}{vector of weights for maginal imbalances. It is required that at least one element is larger than 0. Notice that \code{weight} is only needed when \code{PocSimMIN} is to be used.}
      \item{p}{probabilities of assigning one patient to treatment 1. p should be larger than \code{1/2} to obtain balance. Notice that \code{p} is only needed when \code{"HuHuCAR", "PocSimMIN"} and \code{"StrBCD"} are to be used.}
      \item{a}{a design parameter. As \code{a} goes to \eqn{\infty}, the design becomes more deterministic.}
  \item{bsize}{block size for stratified randomization. It is required to be multiple of 2. Notice that \code{bsize} is only needed when \code{"StrPBR"} is to be used.}
  \item{B}{integer,the number of bootstrap samplings. It is needed only when \code{test} is \code{boot.test}.}
  \item{Reps}{integer,the number of randomized replications. It is needed only when \code{test} is \code{rand.test}.}
  \item{nthreads}{number of threads used in parallel computation. It is needed only under \code{rand.test} or \code{boot.test}.}
  }}
}
\value{
  This function returns a list. The first element is a dataframe representing the powers of chosen test under different values of treatment effects. The second element of the list is the plot of power in which \code{di} forms the vertical axis. 
}

\examples{
##settings
set.seed(2019)
n = 100#<<for demonstration,it is suggested to be larger than 1000
cov_num = 5
level_num = c(2,2,2,2,2)
pr = rep(0.5,10)
beta = c(0.1,0.4,0.3,0.2,0.5)
omega = c(0.1, 0.1, rep(0.8 / 5, times = 5))
di = seq(0,0.5,0.1)
sigma = 1
type = "linear"
p = 0.85
Iternum = 10#<<for demonstration,it is suggested to be around 1000
sl = 0.05
Reps = 10#<<for demonstration,it is suggested to be 200
nthreads = 1

#Evaluation of Power
library("ggplot2")
Strtp=evalPower(n, cov_num, level_num, pr, type,
          beta, di, sigma, Iternum, sl, HuHuCAR,
          rand.test, TRUE, omega, p, Reps, nthreads)
Strtp
}
