\name{evalRand}
\alias{evalRand}
\title{Evaluation of Randomization Procedures}
\usage{
evalRand(data, method = "HuHuCAR", N = 500, ...)
}
\description{
Evaluates a specific randomization procedure based on several different quantities of imbalances. 
}
\arguments{
  \item{data}{a dataframe. A row of the dataframe contains the covariate profile of a patient.}
  \item{N}{the iteration number.}
  \item{method}{the randomization method to be used in allocating patients. The default randomization \dQuote{\code{HuHuCAR}} uses Hu and Hu's general covariate-adaptive randomization; the alternatives are \dQuote{\code{PocSimMIN}}, \dQuote{\code{StrBCD}}, \dQuote{\code{StrPBR}}, \dQuote{\code{DoptBCD}}, and \dQuote{\code{AdjBCD}}.}
  \item{\dots}{arguments to be passed to methods. These depend on the \code{method} and the following arguments are accepted:
  \describe{
  \item{omega}{the vector of weights at the overall, within-stratum, and marginal levels. It is required that at least one element is larger than 0. Note that \code{omega} is only needed when \code{HuHuCAR} is to be assessed.}
  \item{weight}{the vector of weights for all involved margins. It is required that at least one element is NOT 0 and \code{length(weight) = cov_num}. Note that \code{weight} is only needed when \code{PocSimMIN} is to be assessed.}
  \item{p}{the probability of assigning one patinet to treatment 1. \code{p} should be larger than \code{1/2} to obtain balance. Note that \code{p} is only needed when \code{"HuHuCAR"},  \code{"PocSimMIN"} and \code{"StrBCD"} are to be assessed.}
  \item{a}{a design parameter. As \code{a} goes to \eqn{\infty}, the design becomes more deteministic. Note that \code{a} is only needed when \code{"AdjBCD"} is to be assessed.}
  \item{bsize}{the block size for stratified permuted block randomization. It is required to be a multiple of \code{2}. Note that \code{bsize} is only needed when \code{"StrPBR"} is to be assessed.}
  }}
}
\details{
    The \code{data} is designed for \code{N} times using \code{method}.
}
\value{
It returns an object of \code{\link{class}} \code{"careval"}. 

The function \code{\link{print}} is used to obtain results. The generic accessor functions \code{Assig}, \code{Diff}, \code{data}, \code{All strata} and others extract various useful features of the value returned by \code{evalRand}. 

An object of class \code{"careval"} is a list containing at least the following components: 

\item{N}{the number of patients.}
\item{Assig}{a \code{n*N} matrix containing assignments for each patient for N iterations. }
\item{Imb}{a matrix containing maximum, 95\%-quantile, median, and mean of absolute imbalances at overall, each within-stratum and each marginal levels.}
\item{SNUM}{a matrix with \code{N} colunms containing number of patients at overall, each marginal and each within-stratum levels for each iteration.}
\item{Data Type}{the data type. \code{Real} or \code{Simulated}. }
}
\seealso{
  See \code{\link{evalRand.sim}} to evaluate a randomization procedure with covariate data generating mechanism.
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
# a simple use
## Access by real data
## create a dataframe
df <- data.frame("gender" = sample(c("female", "male"), 1000, TRUE, c(1 / 3, 2 / 3)), 
                 "age" = sample(c("0-30", "30-50", ">50"), 1000, TRUE), 
                 "jobs" = sample(c("stu.", "teac.", "others"), 1000, TRUE), 
                 stringsAsFactors = TRUE)
Res <- evalRand(data = df, method = "HuHuCAR", N = 500, 
                omega = c(1, 2, rep(1, ncol(df))), p = 0.85)
## view the output
Res
\donttest{
  ## view all patients' assignments
  Res$Assig}

## Assess by simulated data
cov_num <- 3
level_num <- c(2, 3, 5)
pr <- c(0.35, 0.65, 0.25, 0.35, 0.4, 0.25, 0.15, 0.2, 0.15, 0.25)
n <- 1000
N <- 50
omega = c(1, 2, 1, 1, 2)
# assess Hu and Hu's procedure with the same group of patients
Res.sim <- evalRand.sim(n = n, N = N, Replace = FALSE, cov_num = cov_num, 
                        level_num = level_num, pr = pr, method = "HuHuCAR", 
                        omega, p = 0.85)
\donttest{ 
  ## Compare four procedures
  cov_num <- 3
  level_num <- c(2, 10, 2)
  pr <- c(rep(0.5, times = 2), rep(0.1, times = 10), rep(0.5, times = 2))
  n <- 100
  N <- 200 # <<adjust according to CPU
  bsize <- 4
  ## set weights for HuHuCAR
  omega <- c(1, 2, rep(1, cov_num)); 
  ## set weights for PocSimMIN
  weight = rep(1, cov_num); 
  ## set biased probability
  p = 0.80
  # assess Hu and Hu's procedure
  RH <- evalRand.sim(n = n, N = N, Replace = FALSE, cov_num = cov_num, 
                     level_num = level_num, pr = pr, method = "HuHuCAR", 
                     omega = omega, p = p)
  # assess Pocock and Simon's method
  RPS <- evalRand.sim(n = n, N = N, Replace = FALSE, cov_num = cov_num, 
                      level_num = level_num, pr = pr, method = "PocSimMIN", 
                      weight, p = p)
  # assess Shao's procedure
  RS <- evalRand.sim(n = n, N = N, Replace = FALSE, cov_num = cov_num, 
                     level_num = level_num, pr = pr, method = "StrBCD", 
                     p = p)
  # assess stratified randomization
  RSR <- evalRand.sim(n = n, N = N, Replace = FALSE, cov_num = cov_num, 
                      level_num = level_num, pr = pr, method = "StrPBR", 
                      bsize)
  
  # create containers
  C_M = C_O = C_WS = matrix(NA, nrow = 4, ncol = 4)
  colnames(C_M) = colnames(C_O) = colnames(C_WS) = 
    c("max", "95\%quan", "med", "mean")
  rownames(C_M) = rownames(C_O) = rownames(C_WS) = 
    c("HH", "PocSim", "Shao", "StraRand")
  
  # assess the overall imbalance
  C_O[1, ] = RH$Imb[1, ]
  C_O[2, ] = RPS$Imb[1, ]
  C_O[3, ] = RS$Imb[1, ]
  C_O[4, ] = RSR$Imb[1, ]
  # view the result
  C_O
  
  # assess the marginal imbalances
  C_M[1, ] = apply(RH$Imb[(1 + RH$strt_num) : (1 + RH$strt_num + sum(level_num)), ], 2, mean)
  C_M[2, ] = apply(RPS$Imb[(1 + RPS$strt_num) : (1 + RPS$strt_num + sum(level_num)), ], 2, mean)
  C_M[3, ] = apply(RS$Imb[(1 + RS$strt_num) : (1 + RS$strt_num + sum(level_num)), ], 2, mean)
  C_M[4, ] = apply(RSR$Imb[(1 + RSR$strt_num) : (1 + RSR$strt_num + sum(level_num)), ], 2, mean)
  # view the result
  C_M
  
  # assess the within-stratum imbalances
  C_WS[1, ] = apply(RH$Imb[2 : (1 + RH$strt_num), ], 2, mean)
  C_WS[2, ] = apply(RPS$Imb[2 : (1 + RPS$strt_num), ], 2, mean)
  C_WS[3, ] = apply(RS$Imb[2 : (1 + RS$strt_num), ], 2, mean)
  C_WS[4, ] = apply(RSR$Imb[2 : (1 + RSR$strt_num), ], 2, mean)
  # view the result
  C_WS
  
  # Compare the four procedures through plots
  meth = rep(c("Hu", "PS", "Shao", "STR"), times = 3)
  shape <- rep(1 : 4, times = 3)
  crt <- rep(1 : 3, each = 4)
  crt_c <- rep(c("O", "M", "WS"), each = 4)
  mean <- c(C_O[, 4], C_M[, 4], C_WS[, 4])
  df_1 <- data.frame(meth, shape, crt, crt_c, mean, 
                     stringsAsFactors = TRUE)
  
  require(ggplot2)
  p1 <- ggplot(df_1, aes(x = meth, y = mean, color = crt_c, group = crt,
                         linetype = crt_c, shape = crt_c)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    xlab("method") +
    ylab("absolute mean") +
    theme(plot.title = element_text(hjust = 0.5))
  p1
}
}
