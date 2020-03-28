\name{boot.test}
\alias{boot.test}
\title{Bootstrap t-test}
\usage{
  boot.test(data, B = 200, method = HuHuCAR, conf = 0.95, \dots)
}
\description{
  Performs bootstrap t-test on treatment effects. This test is proposed by Shao et al. (2010) <doi:10.1093/biomet/asq014>.
}
\arguments{
  \item{data}{a dataframe, consisting of patients' profiles, treatment assignments and outputs. See \code{\link{getData}}.}
  \item{B}{integer, the number of bootstrap samples. Default is 200.}
  \item{method}{the randomization to be used in allocating patients. The default randomization \code{HuHuCAR} uses Hu and Hu's general covariate-adaptive randomization; the alternatives are \code{PocSimMIN}, \code{StrBCD}, \code{StrPBR}, \code{DoptBCD} as well as \code{AdjBCD}.}
  \item{conf}{confidence level of the interval. Default is 0.95.}
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
  The bootstrap t-test is described as follows:

1) Generate bootstrap data (\eqn{Y_1^*,Z_1^*), \dots, (Y_n^*,Z_n^*)} as a simple random sample with replacement from the original data \eqn{(Y_1,Z_1), \dots,(Y_n,Z_n)}, where \eqn{Y_i} denotes the outcome, \eqn{Z_i} denotes the profile of the \eqn{i}th patient.

2) Perform covariate-adaptive procedures on patients' profile, get new treatment assignments \eqn{T_1^*,\dots,T_n^*}, and define
\deqn{\hat{\theta}^* = -\frac{1}{n_1^*}\sum\limits_{i=1}^n (T_i^*-2) \times Y_i^* - \frac{1}{n_0^*}\sum\limits_{i=1}^n (T_i^*-1) \times Y_i}
where \eqn{n_1^*} is the number of patients assigned to treatment \eqn{1} and \eqn{n_0^*} is the number of patients assigned to treatment \eqn{2}.

3) Repeat step 2 \eqn{B} times, we can get \eqn{B} independent boostrap samples to obtain \eqn{\hat{\theta}^*_b}, \eqn{b = 1,\dots,B}. Then the variance of \eqn{\bar{Y}_1 - \bar{Y}_0} can be approximated by the sample variance of \eqn{\hat{\theta}^*_b}. 
}
\value{
  It returns an object of class \code{"htest"}.
  
  The function \code{print} is used to obtain results. The generic accessor functions \code{statistic}, \code{p.value}, \code{conf.int} and so on extract various useful features of the value returned by that function.
  
  An object of class \code{"htest"} is a list containing at least the following components:
  \item{data.name}{a character string giving the name(s) of the data.}
  \item{statistic}{value of the t-statistic.}
  \item{pval}{p-value of the test,the null hypothesis is rejected if p-value is less than \code{sl}.}
  \item{conf.int}{a confidence interval under chosen significance level \code{conf} for the treatment effect difference between treatment \code{1} and treatment \code{2}.}
  \item{estimate}{estimated treatment effect difference between treatment \code{1} and treatment \code{2}.}
  \item{method}{a character string indicating what type of test was performed.}
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
                "jobs" = sample(c("stu.", "teac.", "other"), 100, TRUE, c(0.4, 0.2, 0.4)))
##data preprocessing
data.pd <- StrPBR(data = df, bsize = 4)$Cov_Assig

#Then we need to combine patients' profiles and outcomes after randomization and treatments.
outcome = runif(100)
data.combined = data.frame(rbind(data.pd,outcome))

#run the bootstrap t-test
B = 200
Strbt = boot.test(data.combined, B, StrPBR, bsize = 4)
Strbt
}
