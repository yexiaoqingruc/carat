\name{DoptBCD}
\alias{DoptBCD}
\alias{DoptBCD.carandom}
\title{Atkinson's \eqn{D_A}-optimal Biased Coin Design}
\usage{
\method{DoptBCD}{carandom}(data)
}
\description{
Allocates patients to one of two treatments based on \eqn{D_A}-optimal biased coin design with in the presence of the prognostic factors proposed by Atkinson A C (1982) <Doi:10.2307/2335853>.
}
\arguments{
  \item{data}{a dataframe. A row of the dataframe contains the covariate profile of a patient.}
}
\details{
In order to minimize the loss associated with an experiment involving \eqn{n} patients, Atkinson's optimal applied \eqn{D_A}-optimality to the method, in which the probability in the presence of prognostic factors of assigning the \eqn{(n+1)}th patient to treatment 1 is
    \deqn{\frac{[1 - (1; \bold{x}^t_{n+1})(\bold{F}^t_n\bold{F}_n)^{-1}\bold{b}_n]^2}{[1-(1; \bold{x}^t_{n+1})(\bold{F}_n^t\bold{F}_n)^{-1}\bold{b}_n]^2+[1 + (1; \bold{x}_{n+1}^t)(\bold{F}_n^t\bold{F}_n)^{-1}\bold{b}_n]^2},}
where \eqn{\bold{X} = (\bold{x_i}, i = 1, \dots, n)} and \eqn{\bold{x}_i = (x_{i1}, \dots, x_{in})} denoting the covariate-profile of \eqn{i}th patient; and \eqn{\bold{F}_n = [\bold{1}_n; \bold{X}]} is the information matrix; \eqn{\bold{b}_n^T=(2\bold{T}_n-\bold{1}_n)^T\bold{F}_n}, \eqn{\bold{T}_n = (T_1, \dots, T_n)} which is a sequence containing the first \eqn{n} patients' allocations.

    Detials of the the procedure can be seen in A.C.Atkinson (1982).
}
\value{
It returns an object of \code{\link{class}} \code{"carandom"}. 

The functions \code{\link{print}} is used to obtain results. The generic accessor functions \code{Cov_Assig}, \code{Diff}, \code{data}, \code{All strata} and so on extract various useful features of the value returned by that function. 

An object of class \code{"carandom"} is a list containing at least the following components: 

\item{cov_num}{number of covariates.}
\item{n}{number of patients.}
\item{Cov_Assign}{a \code{(cov_num + 1) * n} matrix containing covariates' profile for all patients and corresponding assignments. The \eqn{i}th column represents the \eqn{i}th patient. The first \code{cov_num} rows include patient's covariate profile respectively; and the last row contains the assignment.}
\item{All strata}{a matrix containing all strata involved.}
\item{Diff}{a matrix with only \eqn{1} column. There are final differences at the overall, within-stratum, and marginal levels.}
\item{Data Type}{data type. \code{Real} or \code{Simulated}. }
}
\seealso{
  See \code{\link{DoptBCD.sim}} for allocating patients with covariate data generating mechanism.  
  See \code{\link{DoptBCD.ui}} for command-line user interface.
}
\references{
Atkinson A C. \emph{Optimum biased coin designs for sequential clinical trials with prognostic factors}[J]. Biometrika, 1982, 69(1): 61-67.
}
\examples{
# a simple use
## Real Data
df <- data.frame("gender" = sample(c("female", "male"), 100, TRUE, c(1 / 3, 2 / 3)), 
                 "age" = sample(c("0-30", "30-50", ">50"), 100, TRUE), 
                 "jobs" = sample(c("stu.", "teac.", "others"), 100, TRUE), 
                 stringsAsFactors = TRUE)
Res <- DoptBCD(df)
## view the output
Res
\donttest{
## view all patients' profile and assignments
## Res$Cov_Assig}

## Simulated Data
n <- 1000
cov_num <- 2 

level_num <- c(2, 5)
# Set pr to follow two tips:
#(1) length of pr should be sum(level_num);
#(2)sum of probabilities for each margin should be 1.
pr <- c(0.4, 0.6, rep(0.2, times = 5))
Res.sim <- DoptBCD.sim(n, cov_num, level_num, pr)
## view the output
Res.sim
\donttest{
## view the difference between treatment 1 and treatment 2 
##             at overall, within-strt. and overall levels
Res.sim$Diff}

\donttest{
N <- 5
n <- 100
cov_num <- 2
level_num <- c(3, 5) # << adjust to your CPU and the length should correspond to cov_num
## Set pr to follow two tips:
## (1) length of pr should be sum(level_num);
## (2)sum of probabilities for each margin should be 1
pr <- c(0.3, 0.4, 0.3, rep(0.2, times = 5))
omega <- c(0.2, 0.2, rep(0.6 / cov_num, times = cov_num))

## generate a container to contain Diff
DH <- matrix(NA, ncol = N, nrow = 1 + prod(level_num) + sum(level_num))
DA <- matrix(NA, ncol = N, nrow = 1 + prod(level_num) + sum(level_num))
for(i in 1 : N){
  result <- HuHuCAR.sim(n, cov_num, level_num, pr, omega)
  resultA <- StrBCD.sim(n, cov_num, level_num, pr)
  DH[ , i] <- result$Diff; DA[ , i] <- resultA$Diff
}
## do some analysis
require(dplyr)

## analyze the overall imbalance
Ana_O <- matrix(NA, nrow = 2, ncol = 3)
rownames(Ana_O) <- c("HuHuCAR", "DoptBCD")
colnames(Ana_O) <- c("mean", "median", "95\%quantile")
temp <- DH[1, ] \%>\% abs
tempA <- DA[1, ] \%>\% abs
Ana_O[1, ] <- c((temp \%>\% mean), (temp \%>\% median),
                (temp \%>\% quantile(0.95)))
Ana_O[2, ] <- c((tempA \%>\% mean), (tempA \%>\% median),
                (tempA \%>\% quantile(0.95)))

## analyze the within-stratum imbalances
tempW <- DH[2 : (1 + prod(level_num)), ] \%>\% abs
tempWA <- DA[2 : 1 + prod(level_num), ] \%>\% abs
Ana_W <- matrix(NA, nrow = 2, ncol = 3)
rownames(Ana_W) <- c("HuHuCAR", "DoptBCD")
colnames(Ana_W) <- c("mean", "median", "95\%quantile")
Ana_W[1, ] = c((tempW \%>\% apply(1, mean) \%>\% mean),
               (tempW \%>\% apply(1, median) \%>\% mean),
               (tempW \%>\% apply(1, mean) \%>\% quantile(0.95)))
Ana_W[2, ] = c((tempWA \%>\% apply(1, mean) \%>\% mean),
               (tempWA \%>\% apply(1, median) \%>\% mean),
               (tempWA \%>\% apply(1, mean) \%>\% quantile(0.95)))

## analyze the marginal imbalance
tempM <- DH[(1 + prod(level_num) + 1) :
              (1 + prod(level_num) + sum(level_num)), ] \%>\% abs
tempMA <- DA[(1 + prod(level_num) + 1) :
               (1 + prod(level_num) + sum(level_num)), ] \%>\% abs
Ana_M <- matrix(NA, nrow = 2, ncol = 3)
rownames(Ana_M) <- c("HuHuCAR", "DoptBCD")
colnames(Ana_M) <- c("mean", "median", "95\%quantile")
Ana_M[1, ] = c((tempM \%>\% apply(1, mean) \%>\% mean),
               (tempM \%>\% apply(1, median) \%>\% mean),
               (tempM \%>\% apply(1, mean) \%>\% quantile(0.95)))
Ana_M[2, ] = c((tempMA \%>\% apply(1, mean) \%>\% mean),
               (tempMA \%>\% apply(1, median) \%>\% mean),
               (tempMA \%>\% apply(1, mean) \%>\% quantile(0.95)))

AnaHP <- list(Ana_O, Ana_M, Ana_W)
names(AnaHP) <- c("Overall", "Marginal", "Within-stratum")

AnaHP}
}
