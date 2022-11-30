\name{AdjBCD}
\alias{AdjBCD}
\title{Covariate-adjusted Biased Coin Design}
\usage{
AdjBCD(data, a = 3)
}
\description{
Allocates patients to one of two treatments based on covariate-adjusted biased coin design as proposed by Baldi Antognini A, Zagoraiou M (2011) <doi:10.1093/biomet/asr021>. 
}
\arguments{
  \item{data}{a data frame. A row of the dataframe corresponds to the covariate profile of a patient.}
  \item{a}{a design parameter governing the degree of randomness. The default is \code{3}.}
}
\details{
    Consider \eqn{I} covaraites and \eqn{m_i} levels for the \eqn{i}th covariate, \eqn{i=1,\ldots,I}. \eqn{T_j} is the assignment of the \eqn{j}th patient and \eqn{Z_j = (k_1,\dots,k_I)} indicates the covariate profile of the \eqn{j}th patient, \eqn{j=1,\ldots,n}. For convenience, \eqn{(k_1,\dots,k_I)} and \eqn{(i;k_i)} denote stratum and margin, respectively. \eqn{D_j(.)} is the difference between numbers of patients assigned to treatment \eqn{1} and treatment \eqn{2} at the corresponding level after \eqn{j} patients have been assigned.
    
    Let \eqn{F^a} be a decreasing and symmetric function of \eqn{D_j(.)}, which depends on a design parameter \eqn{a\ge 0}. Then, the probability of allocating the \eqn{(j+1)}th patient to treatment 1 is \eqn{F^a(D_j(.))}, where
    \deqn{F^a(x)=\frac{1}{x^a + 1},}  for \eqn{x\ge 1},
    \deqn{F^a(x)=1 / 2, } for \eqn{ x = 0}, and
    \deqn{F^a(x)=\frac{|x|^a}{|x|^a + 1}, } for \eqn{x\le -1.}
As \eqn{a} goes to \eqn{\infty}, the design becomes more deteministic.

    Details of the procedure can be found in Baldi Antognini and M. Zagoraiou (2011).
}
\value{
It returns an object of \code{\link{class}} \code{"carandom"}. 

An object of class \code{"carandom"} is a list containing the following components: 

\item{datanumeric}{a bool indicating whether the data is a numeric data frame.}
\item{covariates}{a character string giving the name(s) of the included covariates.}
\item{strt_num}{the number of strata.}
\item{cov_num}{the number of covariates.}
\item{level_num}{a vector of level numbers for each covariate.}
\item{n}{the number of patients.}
\item{Cov_Assig}{a \code{(cov_num + 1) * n} matrix containing covariate profiles for all patients and the corresponding assignments. The \eqn{i}th column represents the \eqn{i}th patient. The first \code{cov_num} rows include patients' covariate profiles, and the last row contains the assignments.}
\item{assignments}{the randomization sequence.}
\item{All strata}{a matrix containing all strata involved.}
\item{Diff}{a matrix with only one column. There are final differences at the overall, within-stratum, and within-covariate-margin levels.}
\item{method}{a character string describing the randomization procedure to be used.}
\item{Data Type}{a character string giving the data type, \code{Real} or \code{Simulated}. }
\item{framework}{the framework of the used randomization procedure: stratified randomization, or model-based method.}
\item{data}{the data frame.}
}
\seealso{
  See \code{\link{AdjBCD.sim}} for allocating patients with covariate data generating mechanism; 
  See \code{\link{AdjBCD.ui}} for the command-line user interface.
}
\references{
Baldi Antognini A, Zagoraiou M. \emph{The covariate-adaptive biased coin design for balancing clinical trials in the presence of prognostic factors}[J]. Biometrika, 2011, 98(3): 519-535.
}
\examples{
# a simple use
## Real Data
## create a dataframe
df <- data.frame("gender" = sample(c("female", "male"), 1000, TRUE, c(1 / 3, 2 / 3)), 
                 "age" = sample(c("0-30", "30-50", ">50"), 1000, TRUE), 
                 "jobs" = sample(c("stu.", "teac.", "others"), 1000, TRUE), 
                 stringsAsFactors = TRUE)
Res <- AdjBCD(df, a = 2)
## view the output
Res
\donttest{
  ## view all patients' profile and assignments
  Res$Cov_Assig
  }

## Simulated Data
n <- 1000
cov_num <- 3
level_num <- c(2, 3, 5) 
# Set pr to follow two tips:
#(1) length of pr should be sum(level_num);
#(2) sum of probabilities for each margin should be 1.
pr <- c(0.4, 0.6, 0.3, 0.4, 0.3, rep(0.2, times = 5))
# set the design parameter
a <- 1.8
# obtain result
Res.sim <- AdjBCD.sim(n, cov_num, level_num, pr, a)
\donttest{
  # view the assignments of patients
  Res.sim$Cov_Assig[cov_num + 1, ]
  # view the differences between treatment 1 and treatment 2 at all levels
  Res.sim$Diff
  }
}
\keyword{carandom}
\keyword{Covariate-ajusted biased coin design}
