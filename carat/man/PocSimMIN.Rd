\name{PocSimMIN}
\alias{PocSimMIN}
\title{Pocock and Simon's Method in the Two-Arms Case}
\usage{
PocSimMIN(data, weight = NULL, p = 0.85)
}
\description{
Allocates patients to one of two treatments using Pocock and Simon's method proposed by Pocock S J, Simon R (1975) <doi:10.2307/2529712>. 
}
\arguments{
  \item{data}{a data frame. A row of the dataframe corresponds to the covariate profile of a patient.}
  \item{weight}{a vector of weights for within-covariate-margin imbalances. It is required that at least one element is larger than 0. If \code{weight = NULL} (default), the within-covariate-margin imbalances are weighted with an equal proportion, \code{1/cov_num}, for each covariate-margin.}
  \item{p}{the biased coin probability. \code{p} should be larger than \code{1/2} and less than \code{1}. The default is \code{0.85}.}
}
\details{
     Consider \eqn{I} covariates and \eqn{m_i} levels for the \eqn{i}th covariate, \eqn{i=1,\ldots,I}. \eqn{T_j} is the assignment of the \eqn{j}th patient and \eqn{Z_j = (k_1,\dots,k_I)} indicates the covariate profile of this patient, \eqn{j=1,\ldots,n}. For convenience, \eqn{(k_1,\dots,k_I)} and \eqn{(i;k_i)} denote the stratum and margin, respectively. \eqn{D_j(.)} is the difference between the numbers of patients assigned to treatment \eqn{1} and treatment \eqn{2} at the corresponding levels after \eqn{j} patients have been assigned. The Pocock and Simon's minimization procedure is as follows:

(1) The first patient is assigned to treatment \eqn{1} with probability \eqn{1/2};

(2) Suppose that \eqn{j-1} patients have been assigned (\eqn{1<j\le n}) and the \eqn{j}th patient falls within \eqn{(k_1^*,\dots,k_I^*)};

(3) If the \eqn{j}th patient were assigned to treatment \eqn{1}, then the potential within-covariate-margin differences between the two treatments would be
\deqn{D_j^{(1)}(i;k_i^*)=D_{j-1}(i,k_i^*)+1}

for margin \eqn{(i;k_i^*)}. Similarly, the potential differences would be obtained in the same way if the \eqn{j}th patient were assigned to treatment \eqn{2};

(4) An imbalance measure is defined by
\deqn{Imb_j^{(l)}=\sum_{i=1}^{I}\omega_{m,i}[D_j^{(l)}(i;k_i^*)]^2,l=1,2;}

(5) Conditional on the assignments of the first (\eqn{j-1}) patients as well as the covariate profiles of the first \eqn{j} patients, assign the \eqn{j}th patient to treatment \eqn{1} with the probability
\deqn{P(T_j=1|Z_j,T_1,\dots,T_{j-1})=q} for \eqn{Imb_j^{(1)}>Imb_j^{(2)},}
\deqn{P(T_j=1|Z_j,T_1,\dots,T_{j-1})=p} for \eqn{Imb_j^{(1)}<Imb_j^{(2)}}, and
\deqn{P(T_j=1|Z_j,T_1,\dots,T_{j-1})=0.5} for \eqn{Imb_j^{(1)}=Imb_j^{(2)}.}

 Details of the procedure can be found in Pocock S J, Simon R (1975).
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
\item{weight}{a vector giving the weights imposed on each covariate.}
\item{framework}{the framework of the used randomization procedure: stratified randomization, or model-based method.}
\item{data}{the data frame.}
}
\seealso{
  See \code{\link{PocSimMIN.sim}} for allocating patients with covariate data generating mechanism. 
  See \code{\link{PocSimMIN.ui}} for the command-line user interface.
}
\references{
Pocock S J, Simon R. \emph{Sequential treatment assignment with balancing for prognostic factors in the controlled clinical trial}[J]. Biometrics, 1975: 103-115.
}
\examples{
# a simple use
## Real Data
## creat a dataframe
df <- data.frame("gender" = sample(c("female", "male"), 1000, TRUE, c(1 / 3, 2 / 3)), 
                 "age" = sample(c("0-30", "30-50", ">50"), 1000, TRUE), 
                 "jobs" = sample(c("stu.", "teac.", "others"), 1000, TRUE), 
                 stringsAsFactors = TRUE)
weight <- c(1, 2, 1)
Res <- PocSimMIN(data = df, weight)
## view the output
Res
\donttest{
## view all patients' profile and assignments
Res$Cov_Assig}

## Simulated Data
cov_num = 3
level_num = c(2, 3, 3)
pr = c(0.4, 0.6, 0.3, 0.3, 0.4, 0.4, 0.3, 0.3)
Res.sim <- PocSimMIN.sim(n = 1000, cov_num, level_num, pr)
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
# (2)sum of probabilities for each margin should be 1.
pr <- c(0.4, 0.6, 0.3, 0.4, 0.3, rep(0.2, times = 5))
omega <- c(0.2, 0.2, rep(0.6 / cov_num, times = cov_num))
weight <- c(2, rep(1, times = cov_num - 1))

## generate a container to contain Diff
DH <- matrix(NA, ncol = N, nrow = 1 + prod(level_num) + sum(level_num))
DP <- matrix(NA, ncol = N, nrow = 1 + prod(level_num) + sum(level_num))
for(i in 1 : N){
  result <- HuHuCAR.sim(n, cov_num, level_num, pr, omega)
  resultP <- PocSimMIN.sim(n, cov_num, level_num, pr, weight)
  DH[ , i] <- result$Diff; DP[ , i] <- resultP$Diff
}

## do some analysis
require(dplyr)

## analyze the overall imbalance
Ana_O <- matrix(NA, nrow = 2, ncol = 3)
rownames(Ana_O) <- c("NEW", "PS")
colnames(Ana_O) <- c("mean", "median", "95\%quantile")
temp <- DH[1, ] \%>\% abs
tempP <- DP[1, ] \%>\% abs
Ana_O[1, ] <- c((temp \%>\% mean), (temp \%>\% median),
                (temp \%>\% quantile(0.95)))
Ana_O[2, ] <- c((tempP \%>\% mean), (tempP \%>\% median),
                (tempP \%>\% quantile(0.95)))

## analyze the within-stratum imbalances
tempW <- DH[2 : (1 + prod(level_num)), ] \%>\% abs
tempWP <- DP[2 : 1 + prod(level_num), ] \%>\% abs
Ana_W <- matrix(NA, nrow = 2, ncol = 3)
rownames(Ana_W) <- c("NEW", "PS")
colnames(Ana_W) <- c("mean", "median", "95\%quantile")
Ana_W[1, ] = c((tempW \%>\% apply(1, mean) \%>\% mean),
               (tempW \%>\% apply(1, median) \%>\% mean),
               (tempW \%>\% apply(1, mean) \%>\% quantile(0.95)))
Ana_W[2, ] = c((tempWP \%>\% apply(1, mean) \%>\% mean),
               (tempWP \%>\% apply(1, median) \%>\% mean),
               (tempWP \%>\% apply(1, mean) \%>\% quantile(0.95)))

## analyze the marginal imbalance
tempM <- DH[(1 + prod(level_num) + 1) :
              (1 + prod(level_num) + sum(level_num)), ] \%>\% abs
tempMP <- DP[(1 + prod(level_num) + 1) :
               (1 + prod(level_num) + sum(level_num)), ] \%>\% abs
Ana_M <- matrix(NA, nrow = 2, ncol = 3)
rownames(Ana_M) <- c("NEW", "PS")
colnames(Ana_M) <- c("mean", "median", "95\%quantile")
Ana_M[1, ] = c((tempM \%>\% apply(1, mean) \%>\% mean),
               (tempM \%>\% apply(1, median) \%>\% mean),
               (tempM \%>\% apply(1, mean) \%>\% quantile(0.95)))
Ana_M[2, ] = c((tempMP \%>\% apply(1, mean) \%>\% mean),
               (tempMP \%>\% apply(1, median) \%>\% mean),
               (tempMP \%>\% apply(1, mean) \%>\% quantile(0.95)))

AnaHP <- list(Ana_O, Ana_M, Ana_W)
names(AnaHP) <- c("Overall", "Marginal", "Within-stratum")

AnaHP
}
}
