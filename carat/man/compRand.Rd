\name{compRand}
\alias{compRand}
\title{Compare Different Randomization Procedures via Tables and Plots}
\usage{
compRand(...)
}
\description{
Compares randomization procedures based on several different quantities of imbalances. Among all included randomization procedures of class \code{"careval"}, two or more procedures can be compared in this function. 
}
\arguments{
  \item{\dots}{objects of class \code{"careval"}.}
}
\details{
The primary goal of using covariate-adaptive randomization in practice is to achieve balance with respect to the key covariates and to the overall treatment assignments. We choose four rules to measure the absolute imbalances at overall, marginal and withinn-stratum levels, which are maximal, 95\%quantile, median and mean of the absolute imbalances at different aspects.

(1) Maximal
\deqn{\max_{i = 1, \dots, n}|D_n(\cdot)|.}
(2) 95\% quantile
\deqn{|D_{\lceil0.95n\rceil}(\cdot)|.}
(3) Median

\deqn{(|D_n(\cdot)|)= |D_{(n+1)/2}(\cdot)|}
for \eqn{n} is odd; 
\deqn{(|D_n(\cdot)|)=\frac{1}{2}(|D_{(n/2)}(\cdot)|+|D_{(n/2+1)}(\cdot)|)}
for \eqn{n} is even.

(4) Mean
\deqn{\frac{1}{n}\sum_{j = 1}^{n}|D_{j}(\cdot)|.}

The Monte Carlo method is used to calculate the four types of imbalances.
}
\value{
It returns an object of \code{\link{class}} \code{"carcomp"}. 

The function \code{\link{print}} is used to obtain results. The generic accessor functions \code{Assig}, \code{Diff}, \code{data}, \code{All strata} and others extract various useful features of the value returned by \code{compRand}. 

An object of class \code{"carcomp"} is a list containing at least the following components: 

\item{Overall Imbalances}{a matrix containing maximum, 95\%-quantile, median, mean, and loss of absolute overall imbalances for all the input methods.}
\item{Marginal Imbalances}{a matrix containing maximum, 95\%-quantile, median, mean, and loss of absolute marginal imbalances for all the input methods.}
\item{Within-stratum Imbalances}{a matrix containing maximum, 95\%-quantile, median, mean, loss of absolute imbalances, and also containing mean absolute imbalances of the strata with \eqn{i} patients falling in, where \eqn{i = 1, \dots, bsize} for all the input methods.}
}
\seealso{
  See \code{\link{evalRand}} or \code{\link{evalRand.sim}} to evaluate a specific randomization procedure.
}
\references{
Atkinson A C. \emph{Optimum biased coin designs for sequential clinical trials with prognostic factors}[J]. Biometrika, 1982, 69(1): 61-67.

Baldi Antognini A, Zagoraiou M. \emph{The covariate-adaptive biased coin design for balancing clinical trials in the presence of prognostic factors}[J]. Biometrika, 2011, 98(3): 519-535.

Hu Y, Hu F. \emph{Asymptotic properties of covariate-adaptive randomization}[J]. The Annals of Statistics, 2012, 40(3): 1794-1815.

Pocock S J, Simon R. \emph{Sequential treatment assignment with balancing for prognostic factors in the controlled clinical trial}[J]. Biometrics, 1975: 103-115.

Shao J, Yu X, Zhong B. \emph{A theory for testing hypotheses under covariate-adaptive randomization}[J]. Biometrika, 2010, 97(2): 347-360. 

Zelen M. \emph{The randomization and stratification of patients to clinical trials}[J]. Journal of chronic diseases, 1974, 27(7): 365-375.
}
\examples{
## Compare stratified permuted block randomization and Hu and Hu's general CAR
cov_num <- 2
level_num <- c(2, 2)
pr <- rep(0.5, 4)
n <- 500
N <- 20 # <<adjust according to CPU
bsize <- 4
# set weight for Hu and Hu's method, it satisfies
# (1)Length should equal to cov_num
omega <- c(1, 2, 1, 1)
# Assess Hu and Hu's general CAR
Obj1 <- evalRand.sim(n = n, N = N, Replace = FALSE, cov_num = cov_num, 
                     level_num = level_num, pr = pr, method = "HuHuCAR", 
                     omega, p = 0.85)
# Assess stratified permuted block randomization
Obj2 <- evalRand.sim(n = n, N = N, Replace = FALSE, cov_num = cov_num, 
                     level_num = level_num, pr = pr, method = "StrPBR", 
                     bsize)

RES <- compRand(Obj1, Obj2)

}