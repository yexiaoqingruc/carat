\name{rand.test}
\alias{rand.test}
\title{Randomization Test}
\usage{
  rand.test(data, Reps = 200, method = c("HuHuCAR", "PocSimMIN", "StrBCD", 
                                       "StrPBR", "DoptBCD", "AdjBCD"), 
          conf = 0.95, binwidth = 30, \dots)
}
\description{
  Performs randomization test on treatment effects.
}
\arguments{
  \item{data}{a dataframe. It consists of patients' profiles, treatment assignments and outputs. See \code{\link{getData}}.}
  \item{Reps}{an integer. It represents the number of randomized replications. It is suggested to be \code{200}.}
   \item{method}{a character string specifying the alternative randomization methods to be used in allocating patients, must be one of \code{"HuHuCAR"} (default), \code{"PocSimMIN"}, \code{"StrBCD"}, \code{"StrPBR"}, \code{"DoptBCD"} or \code{"AdjBCD"}.}
  \item{conf}{confidence level of the interval. Default is 0.95.}
  \item{binwidth}{the number of bins for each bar in histogram. The default is 30.}
  \item{\dots}{arguments to be passed to methods. These depends on the \code{method} used and the following arguments are accepted:
  \describe{
      \item{omega}{the vector of weights at the overall, within-stratum, and marginal levels. It is required that at least one element is larger than 0. Note that \code{omega} is only needed when \code{HuHuCAR} is to be used.}
      \item{weight}{the vector of weights for marginal imbalances. It is required that at least one element is larger than 0. Note that \code{weight} is only needed when \code{PocSimMIN} is to be used.}
      \item{p}{the probability of assigning one patient to treatment 1. \code{p} should be larger than \code{1/2} to obtain balance. Note that \code{p} is only needed when \code{"HuHuCAR", "PocSimMIN"} and \code{"StrBCD"} are to be used.}
      \item{a}{a design parameter. As \code{a} goes to \eqn{\infty}, the design becomes more deterministic.}
  \item{bsize}{the block size for stratified randomization. It is required to be a multiple of 2. Note that \code{bsize} is only needed when \code{"StrPBR"} is to be used.}
  }}
}

\details{
  The randomization test is described as follows:
  1) For the observed responses \eqn{Y_1,\dots,Y_n} and the treatment assignments \eqn{T_1,T_2,\dots,T_n}, compute the observed test statistic
\deqn{S_{obs} = \frac{-\sum_{i=1}^nY_i*(T_i-2)}{n_1}-\frac{\sum_{i=1}^n Y_i*(T_i-1)}{n_0}}
  where \eqn{n_1} is the number of patients assigned to treatment \eqn{1} and \eqn{n_0} is the number of patients assigned to treatment \eqn{2};

2) Perform the covariate-adaptive randomization procedure to obtain the new treatment assignments and calculate the corresponding test statistic \eqn{S_i}. And repeat this process \eqn{L} times;

3) Calculate the two-sided Monte Carlo p-value estimator
\deqn{p = \frac{\sum_{l=1}^L I(|S_l|\ge |S_{obs}|)}{L}}
}
\value{
  It returns an object of class \code{"htest"}.
  
  The function \code{print} is used to obtain results.The generic accessor functions \code{statistic}, \code{p.value} and others extract various useful features of the value returned by \code{rand.test}.
  
  An object of class \code{"htest"} is a list containing at least the following components:
  
  \item{data.name}{a character string giving the name(s) of the data.}
  \item{statistic}{the value of the t-statistic. As the randomization test is a nonparametric method, we cannot calculate the t-statistic, so it is hidden in this result.}
  \item{p.value}{p-value of the test, the null hypothesis is rejected if the p-value is less than \code{sl}.}
  \item{conf.int}{a confidence interval under the chosen level \code{conf} for the difference in treatment effect between treatment \code{1} and treatment \code{2}. As the randomization test is a nonparametric method, we cannot calculate the confidence interval, so it is hidden in this result.}
  \item{estimate}{the estimated difference in treatment effects between treatment \code{1} and treatment \code{2}.}
  \item{method}{a character string indicating what type of test was performed.}
}
\references{
Rosenberger W F, Lachin J M. Randomization in clinical trials: \emph{theory and practice}[M]. John Wiley & Sons, 2015.
}
\examples{
##generate data
set.seed(100)
n = 1000
cov_num = 5
level_num = c(2,2,2,2,2)
pr = rep(0.5,10)
beta = c(0.1,0.4,0.3,0.2,0.5)
mu1 = 0
mu2 = 0.01
sigma = 1
type = "linear"
p = 0.85

dataS = getData(n, cov_num, level_num, pr, type,
               beta, mu1, mu2, sigma, "StrBCD", p)

#run the randomization test
library("ggplot2")
Strt = rand.test(data = dataS, Reps = 200,method = "StrBCD", 
                conf = 0.95, binwidth = 30,
                p = 0.85)
Strt
}
