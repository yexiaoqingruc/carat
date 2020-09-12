\name{StrPBR}
\alias{StrPBR}
\title{Stratified Permuted Block Randomization}
\usage{
StrPBR(data, bsize = 4)
}
\description{
Allocates patients to one of two treatments using stratified permuted block randomization proposed by Zelen M (1974) <Doi: 10.1016/0021-9681(74)90015-0>. 
}
\arguments{
  \item{data}{a dataframe. A row of the dataframe contains the covariate profile of a patient.}
  \item{bsize}{the block size for stratified randomization. It is required to be a multiple of 2. The default is \code{4}.}
}
\details{
   Different covariate profiles are defined to be strata, and then permuted block randomization is applied to each stratum. It works efficiently when the number of strata is small, but when the number of strata increases, the stratified permuted block randomization fails to obtain balance between two treatments.
    
    \code{Permuted-block randomization}, or blocking, is used to balance treatment arms within a block so that there are the same number of subjects in each treatment arm. A block contains the same number of each treatment and blocks of different sizes are combined to make up the randomization list. 
}
\value{
It returns an object of \code{\link{class}} \code{"carandom"}. 

The functions \code{\link{print}} is used to obtain results. The generic accessor functions \code{Cov_Assig}, \code{Diff}, \code{data}, \code{All strata} and others extract various useful features of the value returned by \code{StrPBR}. 

An object of class \code{"carandom"} is a list containing at least the following components: 

\item{cov_num}{the number of covariates.}
\item{n}{the number of patients.}
\item{Cov_Assign}{a \code{(cov_num + 1) * n} matrix containing covariate profiles for all patients and corresponding assignments. The \eqn{i}th column represents the \eqn{i}th patient. The first \code{cov_num} rows include patients' covariate profiles, and the last row contains the assignments.}
\item{All strata}{a matrix containing all strata involved.}
\item{Diff}{a matrix with only one column. There are final differences at the overall, within-stratum, and marginal levels.}
\item{Data Type}{the data type. \code{Real} or \code{Simulated}. }
}
\seealso{
  See \code{\link{StrPBR.sim}} for allocating patients with covariate data generating mechanism. 
  See \code{\link{StrPBR.ui}} for the command-line user interface.
}
\references{
Zelen M. \emph{The randomization and stratification of patients to clinical trials}[J]. Journal of chronic diseases, 1974, 27(7): 365-375.
}
\examples{
# a simple use
## Real Data
## creat a dataframe
df <- data.frame("gender" = sample(c("female", "male"), 100, TRUE, c(1 / 3, 2 / 3)), 
                 "age" = sample(c("0-30", "30-50", ">50"), 100, TRUE), 
                 "jobs" = sample(c("stu.", "teac.", "others"), 100, TRUE), 
                 stringsAsFactors = TRUE)
Res <- StrPBR(data = df, bsize = 4)
## view the output
Res
\donttest{
## view all patients' profile and assignments
Res$Cov_Assig}

## Simulated data
cov_num <- 3
level_num <- c(2, 3, 3)
pr <- c(0.4, 0.6, 0.3, 0.4, 0.3, 0.4, 0.3, 0.3)
Res.sim <- StrPBR.sim(n = 100, cov_num, level_num, pr)
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
#(1) length of pr should be sum(level_num);
#(2)sum of probabilities for each margin should be 1.
pr <- c(0.4, 0.6, 0.3, 0.4, 0.3, rep(0.2, times = 5))
omega <- c(0.2, 0.2, rep(0.6 / cov_num, times = cov_num))
# Set block size for stratified randomization
bsize <- 4

## generate a container to contain Diff
DS <- matrix(NA, ncol = N, nrow = 1 + prod(level_num) + sum(level_num))
for(i in 1 : N){
  rtS <- StrPBR.sim(n, cov_num, level_num, pr, bsize)
  DS[ , i] <- rtS$Diff
}

## do some analysis
require(dplyr)

## analyze the overall imbalance
Ana_O <- matrix(NA, nrow = 1, ncol = 3)
rownames(Ana_O) <- c("Str.R")
colnames(Ana_O) <- c("mean", "median", "95\%quantile")
tempS <- DS[1, ] \%>\% abs
Ana_O[1, ] <- c((tempS \%>\% mean), (tempS \%>\% median),
                (tempS \%>\% quantile(0.95)))
## analyze the within-stratum imbalances
tempWS <- DS[2 : 1 + prod(level_num), ] \%>\% abs
Ana_W <- matrix(NA, nrow = 1, ncol = 3)
rownames(Ana_W) <- c("Str.R")
colnames(Ana_W) <- c("mean", "median", "95\%quantile")
Ana_W[1, ] = c((tempWS \%>\% apply(1, mean) \%>\% mean),
               (tempWS \%>\% apply(1, median) \%>\% mean),
               (tempWS \%>\% apply(1, mean) \%>\% quantile(0.95)))

## analyze the marginal imbalance
tempMS <- DS[(1 + prod(level_num) + 1) : (1 + prod(level_num) + sum(level_num)), ] \%>\% abs
Ana_M <- matrix(NA, nrow = 1, ncol = 3)
rownames(Ana_M) <- c("Str.R");
colnames(Ana_M) <- c("mean", "median", "95\%quantile")
Ana_M[1, ] = c((tempMS \%>\% apply(1, mean) \%>\% mean),
               (tempMS \%>\% apply(1, median) \%>\% mean),
               (tempMS \%>\% apply(1, mean) \%>\% quantile(0.95)))

AnaHP <- list(Ana_O, Ana_M, Ana_W)
names(AnaHP) <- c("Overall", "Marginal", "Within-stratum")

AnaHP
}
}
