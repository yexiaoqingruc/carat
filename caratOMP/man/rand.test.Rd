\name{rand.test}
\alias{rand.test}
\title{Randomization Test}
\usage{
  rand.test(data, Reps = 200, method = HuHuCAR, conf = 0.95, 
                         plot = TRUE, binwidth = 30, \dots)
}
\description{
  Performs randomization test on treatment effects.
}
\arguments{
  \item{data}{a dataframe, consists of patients' profiles, treatment assignments and outputs. See \code{\link{getData}}.}
  \item{Reps}{integer, the number of randomized replications. It is suggested to be   200.}
   \item{method}{the randomization to be used in allocating patients. The default randomization \code{HuHuCAR} uses Hu and Hu's general covariate-adaptive randomization; the alternatives are \code{PocSimMIN}, \code{StrBCD}, \code{StrPBR}, \code{DoptBCD} as well as \code{AdjBCD}.}
  \item{conf}{confidence level of the interval. Default is 0.95.}
  \item{plot}{whether plot or not. Optional input: TRUE or FALSE.}
  \item{binwidth}{number of bins for each bar in histogram. Default is 30.}
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
  The randomization test is decribed as follows:
  1) For the observed responses \eqn{Y_1,\dots,Y_n} and the treatment assignments \eqn{T_1,T_2,\dots,T_n}, compute the observed test statistic
\deqn{S_{obs} = \frac{-\sum_{i=1}^nY_i*(T_i-2)}{n_1}-\frac{\sum_{i=1}^n Y_i*(T_i-1)}{n_0}}
  where \eqn{n_1} is the number of patients assigned to treatment \eqn{1} and \eqn{n_0} is the number of patients assigned to treatment \eqn{2};

2) Perform the covariate-adaptive randomization procedure, get new treatment assignments and calculate the corresponding test statistic \eqn{S_i}, repeat \eqn{L} times;

3) Calculate the two-sided Monte Carlo p-value estimator
\deqn{p = \frac{\sum_{l=1}^L I(|S_l|\ge |S_{obs}|)}{L}}
}
\value{
  It returns an object of class \code{"htest"}.
  
  The function \code{print} is used to obtain results.The generic accessor functions \code{statistic}, \code{p.value} and so on extract various useful features of the value returned by that function.
  
  An object of class \code{"htest"} is a list containing at least the following components:
  
  \item{data.name}{a character string giving the name(s) of the data.}
  \item{statistic}{value of the t-statistic. Since the randomization test is a nonparametric method, we can't calculate the t-statistic, so it is hidden in this result.}
  \item{p.value}{p-value of the test,the null hypothesis is rejected if p-value is less than \code{sl}.}
  \item{conf.int}{a confidence interval under chosen significance level \code{conf} for the treatment effect difference between treatment \code{1} and treatment \code{2}. Since the randomization test is a nonparametric method, we can't calculate the confidence interval, so it is hidden in this result.}
  \item{estimate}{estimated treatment effect difference between treatment \code{1} and treatment \code{2}.}
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
               beta, mu1, mu2, sigma, StrBCD, p)

#run the randomization test
library("ggplot2")
Strt = rand.test(data = dataS, Reps = 200,method = StrBCD, 
                conf = 0.95, plot = TRUE, binwidth = 30,
                p = 0.85)
Strt
}
