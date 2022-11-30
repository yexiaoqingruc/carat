\name{boot.test}
\alias{boot.test}
\title{Bootstrap t-test}
\usage{
  boot.test(data, B = 200, method = c("HuHuCAR", "PocSimMIN", "StrBCD", 
                                  "StrPBR", "DoptBCD", "AdjBCD"), 
          conf = 0.95, \dots)
}
\description{
  Performs bootstrap t-test on treatment effects. This test is proposed by Shao et al. (2010) <doi:10.1093/biomet/asq014>.
}
\arguments{
  \item{data}{a data frame. It consists of patients' profiles, treatment assignments and outputs. See \code{\link{getData}}.}
  \item{B}{an integer. It is the number of bootstrap samples. The default is \code{200}.}
  \item{method}{the randomization procedure to be used for testing. This package provides tests for \code{"HuHuCAR"}, \code{"PocSimMIN"}, \code{"StrBCD"}, \code{"StrPBR"}, \code{"AdjBCD"}, and \code{"DoptBCD"}.}
  \item{conf}{confidence level of the interval. The default is \code{0.95}.}
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
  The bootstrap t-test is described as follows:

1) Generate bootstrap data (\eqn{Y_1^*,Z_1^*), \dots, (Y_n^*,Z_n^*)} as a simple random sample with replacement from the original data \eqn{(Y_1,Z_1), \dots,(Y_n,Z_n)}, where \eqn{Y_i} denotes the outcome and \eqn{Z_i} denotes the profile of the \eqn{i}th patient.

2) Perform covariate-adaptive procedures on the patients' profiles to obtain new treatment assignments \eqn{T_1^*,\dots,T_n^*}, and define
\deqn{\hat{\theta}^* = -\frac{1}{n_1^*}\sum\limits_{i=1}^n (T_i^*-2) \times Y_i^* - \frac{1}{n_0^*}\sum\limits_{i=1}^n (T_i^*-1) \times Y_i}
where \eqn{n_1^*} is the number of patients assigned to treatment \eqn{1} and \eqn{n_0^*} is the number of patients assigned to treatment \eqn{2}.

3) Repeat step 2 \eqn{B} times to generate \eqn{B} independent boostrap samples to obtain \eqn{\hat{\theta}^*_b}, \eqn{b = 1,\dots,B}. The variance of \eqn{\bar{Y}_1 - \bar{Y}_0} can then be approximated by the sample variance of \eqn{\hat{\theta}^*_b}. 
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
  Shao J, Yu X, Zhong B. \emph{A theory for testing hypotheses under covariate-adaptive randomization}[J]. Biometrika, 2010, 97(2): 347-360.
}
\examples{
#Suppose the data used is patients' profile from real world, 
#  while it is generated here. Data needs to be preprocessed 
#  and then get assignments following certain randomization.
set.seed(100)
df<- data.frame("gender" = sample(c("female", "male"), 100, TRUE, c(1 / 3, 2 / 3)),
                "age" = sample(c("0-30", "30-50", ">50"), 100, TRUE),
                "jobs" = sample(c("stu.", "teac.", "other"), 100, TRUE, c(0.4, 0.2, 0.4)), 
                stringsAsFactors = TRUE)
##data preprocessing
data.pd <- StrPBR(data = df, bsize = 4)$Cov_Assig

#Then we need to combine patients' profiles and outcomes after randomization and treatments.
outcome = runif(100)
data.combined = data.frame(rbind(data.pd,outcome), stringsAsFactors = TRUE)

#run the bootstrap t-test
B = 200
Strbt = boot.test(data.combined, B, "StrPBR", bsize = 4)
Strbt
}
