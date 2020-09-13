\name{AdjBCD}
\alias{AdjBCD}
\title{Covariate-adjusted Biased Coin Design}
\usage{
AdjBCD(data, a = 2)
}
\description{
Allocates patients to one of two treatments based on covariate-adjusted biased coin design as proposed by Baldi Antognini A, Zagoraiou M (2011) <Doi:10.1093/biomet/asr021>. 
}
\arguments{
  \item{data}{a dataframe. A row of the dataframe contains the covariate profile of a certain patient.}
  \item{a}{a design parameter. The default is \code{2}. As \code{a} goes to \eqn{\infty}, the design becomes more deterministic.}
}
\details{
    Consider \eqn{I} covaraites and \eqn{m_i} levels for the \eqn{i}th covariate. \eqn{T_j} is the assignment of the \eqn{j}th patient and \eqn{Z_j = (k_1,\dots,k_I)} indicates the covariate profile of the \eqn{j}th patient. For convenience, \eqn{(k_1,\dots,k_I) and (i;k_i)} denote stratum and margin respectively. \eqn{D_n(.)} is the difference between numbers of assigned patients in treatment \eqn{1} and treatment \eqn{2} at the corresponding level after \eqn{n} patients have been assigned.
    
    Let \eqn{F^a} be a decreasing and symmetric function of \eqn{D_n(.)}, which depends on a design parameter \eqn{a\ge 0}. Then the probability of allocating the \eqn{(n+1)}th patient to treatment 1 is \eqn{F^a(D_n(.))}, where
    \deqn{F^a(x)=\frac{|x|^a}{|a|^a + 1}, } for \eqn{x\le -1},
    \deqn{F^a(x)=1 / 2, } for \eqn{ x = 0}, and
    \deqn{F^a(x)=\frac{1}{|x|^a + 1},}  for \eqn{x\ge 1.}
As \code{a} goes to \eqn{\infty}, the design becomes more deteministic.

    Details of the procedure can be found in Baldi Antognini and M. Zagoraiou (2011).
}
\value{
It returns an object of \code{\link{class}} \code{"carandom"}. 

The function \code{\link{print}} is used to obtain results. The generic accessor functions \code{Cov_Assig}, \code{Diff}, \code{data}, \code{All strata} and others extract various useful features of the value returned by \code{AdjBCD}.

An object of class \code{"carandom"} is a list containing at least the following components: 

\item{cov_num}{number of covariates.}
\item{n}{number of patients.}
\item{Cov_Assign}{a \code{(cov_num + 1) * n} matrix containing covariate profiles for all patients and the corresponding assignments. The \eqn{i}th column represents the \eqn{i}th patient. The first \code{cov_num} rows include patients' covariate profiles, and the last row contains the assignment.}
\item{All strata}{a matrix containing all strata involved.}
\item{Diff}{a matrix with only one column. There are final differences at the overall, within-stratum, and marginal levels.}
\item{Data Type}{data type. \code{Real} or \code{Simulated}. }
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
