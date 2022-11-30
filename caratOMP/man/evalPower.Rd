\name{evalPower}
\alias{evalPower}
\title{Evaluation of Tests and Randomization Procedures through Power}
\usage{
  evalPower(n, cov_num, level_num, pr, type, beta, di = seq(0,0.5,0.1), sigma = 1,
          Iternum, sl = 0.05, method = c("HuHuCAR", "PocSimMIN", "StrBCD", "StrPBR", 
                                         "DoptBCD","AdjBCD"), 
          test = c("boot.test", "corr.test", "rand.test"), plot = TRUE, ...)
}
\description{
  Returns powers and a plot of the chosen test and method under different treatment effects.
}
\arguments{
  \item{n}{the number of patients.}
  \item{cov_num}{the number of covariates.}
  \item{level_num}{a vector of level numbers for each covariate. Hence the length of \code{level_num} should be equal to the number of covariates.}
  \item{pr}{a vector of probabilities. Under the assumption of independence between covariates, \code{pr} is a vector containing probabilities for each level of each covariate. The length of \code{pr} should correspond to the number of all levels, and the sum of the probabilities for each margin should be \code{1}.}
  \item{type}{a data-generating method. Optional input: \code{"linear"} or \code{"logit"}.}
  \item{beta}{a vector of coefficients of covariates. The length of \code{beta} must correspond to the sum of all covariates' levels.}
  \item{di}{a value or a vector of values of difference in treatment effects. The default value is a sequence from \code{0} to \code{0.5} with increments of \code{0.1}. The value(s) forms the horizontal axis of the plot.}
  \item{sigma}{the error variance for the linear model. The default is 1. This should be a positive value and is only used when \code{type = linear}.}
  \item{Iternum}{an integer. It is the number of iterations required for power calculation.}
  \item{sl}{the significance level. If the \eqn{p} value returned by the test is less than \code{sl}, the null hypothesis will be rejected. The default value is \code{0.05}.}
  \item{method}{the randomization procedure to be used for power calculation. This package provides power calculation for \code{"HuHuCAR"}, \code{"PocSimMIN"}, \code{"StrBCD"}, \code{"StrPBR"}, \code{"AdjBCD"}, and \code{"DoptBCD"}.}
  \item{test}{a character string specifying the alternative tests used to verify hypothesis, must be one of \code{"boot.test"}, \code{"corr.test"} or \code{"rand.test"}, which are the bootstrap \eqn{t} test, the corrected \eqn{t} test, and the randomization test, respectively. The arguments associated with the testing function can be specified; otherwise, the default value will be used.}
  \item{plot}{a bool. It indicates whether to plot or not. Optional input: \code{TRUE} or \code{FALSE}.}
  \item{\dots}{arguments to be passed to \code{method}. These arguments depend on the randomization method used and the following arguments are accepted:
  \describe{
      \item{omega}{a vector of weights at the overall, within-stratum, and within-covariate-margin levels. It is required that at least one element is larger than 0. Note that \code{omega} is only needed when \code{HuHuCAR} is to be used.}
      \item{weight}{a vector of weights for within-covariate-margin imbalances. It is required that at least one element is larger than 0. Note that \code{weight} is only needed when \code{PocSimMIN} is to be used.}
      \item{p}{the biased coin probability. \code{p} should be larger than \code{1/2} and less than \code{1}. Note that \code{p} is only needed when \code{"HuHuCAR", "PocSimMIN"} and \code{"StrBCD"} are to be used.}
      \item{a}{a design parameter governing the degree of randomness. Note that \code{a} is only needed when \code{"AdjBCD"} is to be used.}
  \item{bsize}{the block size for the stratified randomization. It is required to be a multiple of 2. Note that \code{bsize} is only needed when \code{"StrPBR"} is to be used.}
  \item{B}{an integer. It is the number of bootstrap samples. It is needed only when \code{test} is \code{boot.test}.}
  \item{Reps}{an integer. It is the number of randomized replications used in the randomization test. It is needed only when \code{test} is \code{rand.test}.}
  \item{nthreads}{the number of threads to be used in parallel computation. This is needed only under \code{rand.test} and \code{boot.test}. The default is \code{1}. }
  }}
}
\value{
  This function returns a list. The first element is a dataframe representing the powers of the chosen test under different values of treatment effects. The second element is the execution time. An optional element is the plot of power in which \code{di} forms the vertical axis. 
}

\examples{
##settings
set.seed(2019)
n = 100#<<for demonstration,it is suggested to be larger than 1000
cov_num = 5
level_num = c(2,2,2,2,2)
pr = rep(0.5,10)
beta = c(0.1,0.4,0.3,0.2,0.5,0.5,0.4,0.3,0.2,0.1)
omega = c(0.1, 0.1, rep(0.8 / 5, times = 5))
di = seq(0,0.5,0.1)
sigma = 1
type = "linear"
p = 0.85
Iternum = 10#<<for demonstration,it is suggested to be around 1000
sl = 0.05
Reps = 10#<<for demonstration,it is suggested to be 200

#Evaluation of Power
library("ggplot2")
Strtp=evalPower(n,cov_num,level_num,pr,type,beta,di,sigma,
                Iternum,sl,"HuHuCAR","rand.test",TRUE,omega,p,Reps, nthreads = 1)
Strtp
}
