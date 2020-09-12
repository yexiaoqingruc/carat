\name{evalPower}
\alias{evalPower}
\title{Evaluation of Tests and Randomization Procedures through Power}
\usage{
  evalPower(n, cov_num, level_num, pr, type, beta, di = seq(0,0.5,0.1), sigma = 1,
          Iternum, sl = 0.05, method = c("HuHuCAR", "PocSimMIN", "StrBCD", "StrPBR", 
                                         "DoptBCD","AdjBCD"), 
          test = c("rand.test", "boot.test", "corr.test"), plot = TRUE, ...)
}
\description{
  Returns powers and a plot of the chosen test and method under different treatment effects.
}
\arguments{
  \item{n}{number of patients.}
  \item{cov_num}{number of covariates.}
  \item{level_num}{the vector of level numbers for each covariate. Hence the length of \code{level_num} should be equal to the number of covariates.}
  \item{pr}{the vector of probabilities. Under the assumption of independence between covariates, \code{pr} is a vector containing probabilities for each level of each covariate. The length of pr should correspond to the number of all levels, and the vector sum of \code{pr} should be equal to \code{cov_num}.}
  \item{type}{the type of models when generating data. Optional input: \code{linear} or \code{logit}.}
  \item{beta}{the vector of coefficients of covariates. The length of \code{beta} must correspond to \code{cov_num}.}
  \item{di}{the vector of values of difference in treatment effects. The default value is a sequence from \code{0} to \code{0.5} with increment being \code{0.1}.}
  \item{sigma}{the error variance for the linear model. The default value is 1. It is only used when \code{type} is \code{linear}.}
  \item{Iternum}{an integer. It is the number of iterations required for calculating the average power.}
  \item{sl}{the significance level. If the p-value returned by the test is less than \code{sl}, we will reject the null hypothesis. The default value is \code{0.05}.}
  \item{method}{a character string specifying the alternative randomization methods to be used in allocating patients, must be one of \code{"HuHuCAR"} (default), \code{"PocSimMIN"}, \code{"StrBCD"}, \code{"StrPBR"}, \code{"DoptBCD"} or \code{"AdjBCD"}.}
  \item{test}{a charater string specifying the alternative test used to verify hypothesis, must one of \code{"rand.test"}, \code{"boot.test"} or \code{"corr.test"}, which are the randomization test, the bootstrap t-test and the corrected t-test respectively.}
  \item{plot}{bool. It shows whether to plot or not. Optional input: \code{TRUE} or \code{FALSE}.}
  \item{...}{arguments to be passed to methods. These depend on the \code{method} and \code{test} used and the following arguments are accepted:
  \describe{
      \item{omega}{the vector of weights at the overall, within-stratum, and marginal levels. It is required that at least one element is larger than 0. Note that \code{omega} is only needed when \code{HuHuCAR} is to be used.}
      \item{weight}{the vector of weights for marginal imbalances. It is required that at least one element is larger than 0. Note that \code{weight} is only needed when \code{PocSimMIN} is to be used.}
      \item{p}{the probability of assigning one patient to treatment 1, where \code{p} should be larger than \code{1/2} to obtain balance. Note that \code{p} is only needed when \code{"HuHuCAR", "PocSimMIN"} and \code{"StrBCD"} are to be used.}
      \item{a}{a design parameter. As \code{a} goes to \eqn{\infty}, the design becomes more deterministic.}
  \item{bsize}{block size for the stratified randomization. It is required to be a multiple of 2. Note that \code{bsize} is only needed when \code{"StrPBR"} is to be used.}
  \item{B}{an integer. It is the number of bootstrap samplings. It is needed only when \code{test} is \code{boot.test}.}
  \item{Reps}{an integer. It represents the number of randomized replications. It is needed only when \code{test} is \code{rand.test}.}
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
beta = c(0.1,0.4,0.3,0.2,0.5)
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
