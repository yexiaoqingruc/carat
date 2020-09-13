\name{StrBCD}
\alias{StrBCD}
\title{Shao's Method in the Two-Arms Case}
\usage{
StrBCD(data, p = 0.85)
}
\description{
Allocates patients to one of the two treatments using Shao's method proposed by Shao J, Yu X, Zhong B (2010) <Doi:10.1093/biomet/asq014>. 
}
\arguments{
  \item{data}{a dataframe. A row of the dataframe contains the covariate profile of a patient.}
  \item{p}{the probability of assigning one patient to treatment 1. \code{p} should be larger than \code{1/2} to obtain balance. The default is \code{0.85}.}
}
\details{
    Consider \eqn{I} covaraites and \eqn{m_i} levels for the \eqn{i}th covariate. \eqn{T_j} is the assignment of the \eqn{j}th patient and \eqn{Z_j = (k_1,\dots,k_I)} indicates the covariate profile of this patient. For convenience, \eqn{(k_1,\dots,k_I)} and \eqn{(i;k_i)} denote the stratum and margin respectively. \eqn{D_n(.)} is the difference between the numbers of assigned patients in treatment \eqn{1} and treatment \eqn{2} at the corresponding level after \eqn{n} patinets have been assigned.Then Shao's procedure is as follows:

(1) The first patient is assigned to treatment \eqn{1} with probability \eqn{1/2};

(2) Suppose \eqn{n-1} patients have each been assigned to a treatment (\eqn{n>1}) and the \eqn{n}th patinent falls within \eqn{(k_1^*,\dots,k_I^*)};

(3)If the \eqn{n}th patient was assigned to treatment \eqn{1}, then the potential within-stratum difference between the two groups is
\deqn{D_n^{(1)}(k_1^*,\dots,k_I^*)=D_n(k_1^*,\dots,k_I^*)+1.}
Similarly, the potential differences would be obtained in the same way if the \eqn{n}th patinent was assigned to treatment \eqn{2}.

(4) An imbalance measure is defined by
\deqn{Imb_n^{(l)}=[D_n^{(1)}(k_1^*,\dots,k_I^*)]^2,l=1,2;}

(5)Conditional on the assignments of the first (\eqn{n-1}) patients as well as the covariates'profiles of the first \eqn{n} patients, assign the nth patient to treatment \eqn{1} with probability
\deqn{P(T_n=1|Z_n,T_1,\dots,T_{n-1})=q,} for \eqn{Imb_n^{(1)}>Imb_n^{(2)},}
\deqn{P(T_n=1|Z_n,T_1,\dots,T_{n-1})=p,} for \eqn{Imb_n^{(1)}<Imb_n^{(2)}}, and
\deqn{P(T_n=1|Z_n,T_1,\dots,T_{n-1})=0.5,} for \eqn{Imb_n^{(1)}=Imb_n^{(2)}.}
}
\value{
It returns an object of \code{\link{class}} \code{"carandom"}. 

The function \code{\link{print}} is used to obtain results. The generic accessor functions \code{Cov_Assig}, \code{Diff}, \code{data}, \code{All strata} and others extract various useful features of the value returned by \code{StrBCD}. 

An object of class \code{"carandom"} is a list containing at least the following components: 

\item{cov_num}{the number of covariates.}
\item{n}{the number of patients.}
\item{Cov_Assign}{a \code{(cov_num + 1) * n} matrix containing covariate profiles for all patients and corresponding assignments. The \eqn{i}th column represents the \eqn{i}th patient. The first \code{cov_num} rows include patients' covariate profiles, and the last row contains the assignment.}
\item{All strata}{a matrix containing all the strata involved.}
\item{Diff}{a matrix with only one column. There are final differences at the overall, within-stratum, and marginal levels.}
\item{Data Type}{the data type. \code{Real} or \code{Simulated}. }
}
\seealso{
  See \code{\link{StrBCD.sim}} for allocating patients with covariate data generating mechanism. 
  See \code{\link{StrBCD.ui}} for command-line user interface.
}
\references{
Shao J, Yu X, Zhong B. \emph{A theory for testing hypotheses under covariate-adaptive randomization}[J]. Biometrika, 2010, 97(2): 347-360. 
}
\examples{
# a simple use
## Real Data
## creat a dataframe
df <- data.frame("gender" = sample(c("female", "male"), 1000, TRUE, c(1 / 3, 2 / 3)), 
                 "age" = sample(c("0-30", "30-50", ">50"), 1000, TRUE), 
                 "jobs" = sample(c("stu.", "teac.", "others"), 1000, TRUE), 
                 stringsAsFactors = TRUE)
Res <- StrBCD(data = df)
## view the output
Res
\donttest{
## view all patients' profile and assignments
Res$Cov_Assig}

## Simulated Data
cov_num = 3
level_num = c(2, 3, 3)
pr = c(0.4, 0.6, 0.3, 0.4, 0.3, 0.4, 0.3, 0.3)
Res.sim <- StrBCD.sim(n = 1000, cov_num, level_num, pr)
## view the output
Res.sim
\donttest{
## view the detials of difference
Res.sim$Diff}

\donttest{
N <- 5
n <- 1000
cov_num <- 3
level_num <- c(2, 3, 5) 
# Set pr to follow two tips:
# (1) length of pr should be sum(level_num);
# (2)sum of probabilities for each margin should be 1
pr <- c(0.4, 0.6, 0.3, 0.4, 0.3, rep(0.2, times = 5))
omega <- c(0.2, 0.2, rep(0.6 / cov_num, times = cov_num))

## generate a container to contain Diff
DH <- matrix(NA, ncol = N, nrow = 1 + prod(level_num) + sum(level_num))
DS <- matrix(NA, ncol = N, nrow = 1 + prod(level_num) + sum(level_num))
for(i in 1 : N){
  result <- HuHuCAR.sim(n, cov_num, level_num, pr, omega)
  resultS <- StrBCD.sim(n, cov_num, level_num, pr)
  DH[ , i] <- result$Diff; DS[ , i] <- resultS$Diff
}

## do some analysis
require(dplyr)

## analyze the overall imbalance
Ana_O <- matrix(NA, nrow = 2, ncol = 3)
rownames(Ana_O) <- c("NEW", "Shao")
colnames(Ana_O) <- c("mean", "median", "95\%quantile")
temp <- DH[1, ] \%>\% abs
tempS <- DS[1, ] \%>\% abs
Ana_O[1, ] <- c((temp \%>\% mean), (temp \%>\% median),
                (temp \%>\% quantile(0.95)))
Ana_O[2, ] <- c((tempS \%>\% mean), (tempS \%>\% median),
                (tempS \%>\% quantile(0.95)))

## analyze the within-stratum imbalances
tempW <- DH[2 : (1 + prod(level_num)), ] \%>\% abs
tempWS <- DS[2 : 1 + prod(level_num), ] \%>\% abs
Ana_W <- matrix(NA, nrow = 2, ncol = 3)
rownames(Ana_W) <- c("NEW", "Shao")
colnames(Ana_W) <- c("mean", "median", "95\%quantile")
Ana_W[1, ] = c((tempW \%>\% apply(1, mean) \%>\% mean),
               (tempW \%>\% apply(1, median) \%>\% mean),
               (tempW \%>\% apply(1, mean) \%>\% quantile(0.95)))
Ana_W[2, ] = c((tempWS \%>\% apply(1, mean) \%>\% mean),
               (tempWS \%>\% apply(1, median) \%>\% mean),
               (tempWS \%>\% apply(1, mean) \%>\% quantile(0.95)))

## analyze the marginal imbalance
tempM <- DH[(1 + prod(level_num) + 1) :
              (1 + prod(level_num) + sum(level_num)), ] \%>\% abs
tempMS <- DS[(1 + prod(level_num) + 1) :
               (1 + prod(level_num) + sum(level_num)), ] \%>\% abs
Ana_M <- matrix(NA, nrow = 2, ncol = 3)
rownames(Ana_M) <- c("NEW", "Shao")
colnames(Ana_M) <- c("mean", "median", "95\%quantile")
Ana_M[1, ] = c((tempM \%>\% apply(1, mean) \%>\% mean),
               (tempM \%>\% apply(1, median) \%>\% mean),
               (tempM \%>\% apply(1, mean) \%>\% quantile(0.95)))
Ana_M[2, ] = c((tempMS \%>\% apply(1, mean) \%>\% mean),
               (tempMS \%>\% apply(1, median) \%>\% mean),
               (tempMS \%>\% apply(1, mean) \%>\% quantile(0.95)))

AnaHP <- list(Ana_O, Ana_M, Ana_W)
names(AnaHP) <- c("Overall", "Marginal", "Within-stratum")

AnaHP
}
}

