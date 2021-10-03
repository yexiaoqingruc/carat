\name{corr.test}
\alias{corr.test}
\title{Corrected t-test}
\usage{
  corr.test(data, conf = 0.95)
}
\description{
   Performs corrected t-test on treatment effects. This test follows the idea of Ma et al. (2015) <doi:10.1080/01621459.2014.922469>.
}
\arguments{
\item{data}{a data frame. It consists of patients' profiles, treatment assignments and outputs. See \code{\link{getData}}.}
\item{conf}{confidence level of the interval. The default is \code{0.95}.}
}
\details{
  When the working model is the true underlying linear model, and the chosen covariate-adaptive design achieves that the overall imbalance and marginal imbalances for all covariates are bounded in probability, we can derive the asymptotic distribution under the null distribution, where the treatment effect of each group is the same. Subsequently, we can replace the variance estimator in a simple two sample t-test with an adjusted variance estimator. Details can be found in Ma et al.(2015). 
}

\value{
  It returns an object of class \code{"htest"}.
  
  An object of class \code{"htest"} is a list containing the following components:
  
  \item{statistic}{the value of the t-statistic.}
  \item{p.value}{the p-value of the test,the null hypothesis is rejected if p-value is less than the pre-determined significance level.}
  \item{conf.int}{a confidence interval under the chosen level \code{conf} for the difference in treatment effect between treatment \code{1} and treatment \code{2}.}
  \item{estimate}{the estimated treatment effect difference between treatment \code{1} and treatment \code{2}.}
  \item{stderr}{the standard error of the mean (difference), used as denominator in the t-statistic formula.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
}

\references{
  Ma W, Hu F, Zhang L. \emph{Testing hypotheses of covariate-adaptive randomized clinical trials}[J]. Journal of the American Statistical Association, 2015, 110(510): 669-680.
}

\examples{
##generate data
set.seed(100)
n = 1000
cov_num = 5
level_num = c(2,2,2,2,2)
pr = rep(0.5,10)
beta = c(0.1,0.4,0.3,0.2,0.5)
omega = c(0.1, 0.1, rep(0.8 / 5, times = 5))
mu1 = 0
mu2 = 0.7
sigma = 1
type = "linear"
p = 0.85

dataH = getData(n,cov_num,level_num,pr,type,beta,
                mu1,mu2,sigma,"HuHuCAR",omega,p)

#run the corrected t-test
HHct=corr.test(dataH)
HHct
}
