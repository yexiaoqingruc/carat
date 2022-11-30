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
The primary goal of using covariate-adaptive randomization in practice is to achieve balance with respect to the key covariates. We choose four rules to measure the absolute imbalances at overall, within-covariate-margin, and within-stratum levels, which are maximal, 95\%quantile, median and mean of the absolute imbalances at different aspects. The Monte Carlo method is used to calculate the four types of imbalances. Let \eqn{D_{n,i}(\cdot)} be the final difference at the corresponding level for \eqn{i}th iteration, \eqn{i=1,\ldots}, \code{N}, and \code{N} is the number of iterations.

(1) Maximal
\deqn{\max_{i = 1, \dots, N}|D_{n,i}(\cdot)|.}

(2) 95\% quantile
\deqn{|D_{n,\lceil0.95N\rceil}(\cdot)|.}

(3) Median
\deqn{|D_{n,(N+1)/2}(\cdot)|}
for \code{N} is odd, and
\deqn{\frac{1}{2}(|D_{(N/2)}(\cdot)|+|D_{(N/2+1)}(\cdot)|)}
for \code{N} is even.

(4) Mean
\deqn{\frac{1}{N}\sum_{i = 1}^{N}|D_{n, i}(\cdot)|.}
}
\value{
It returns an object of \code{\link{class}} \code{"carcomp"}. 

An object of class \code{"carcomp"} is a list containing the following components: 

\item{Overall Imbalances}{a matrix containing the maximum, 95\%-quantile, median, and mean of the absolute overall imbalances for the randomization method(s) to be evaluated.}
\item{Within-covariate-margin Imbalances Imbalances}{a matrix containing the maximum, 95\%-quantile, median, and mean of the absolute within-covariate-margin imbalances for the randomization method(s) to be evaluated.}
\item{Within-stratum Imbalances}{a matrix containing the maximum, 95\%-quantile, median, and mean of the absolute within-stratum imbalances for the randomization method(s) to be evaluated.}
\item{dfmm}{a data frame containing the mean absolute imbalances at the overall, within-stratum, and within-covariate-margin levels for the randomization method(s) to be evaluated.}
\item{df_abm}{a data frame containing the absolute imbalances at the overall, within-stratum, and within-covariate-margin levels.}
\item{mechanism}{a character string giving the randomization method(s) to be evaluated.}
\item{n}{the number of patients.}
\item{iteration}{the number of iterations.}
\item{cov_num}{the number of covariates.}
\item{level_num}{a vector of level numbers for each covariate.}
\item{Data Type}{a character string giving the data type, \code{Real} or \code{Simulated}. }
\item{DataGeneration}{a bool vector indicating whether the data used for all the iterations is the same for the randomization method(s) to be evaluated.}
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