\name{AdjBCD.ui}
\alias{AdjBCD.ui}
\title{
Command-line User Interface Using Covariate-adjusted Biased Coin Design
}
\description{
A call to the user-interface function for allocation of patients to one of two treatments, using covariate-adjusted biased coin design, as proposed by Baldi Antognini A, Zagoraiou M (2011) <doi:10.1093/biomet/asr021>. 
}
\usage{
AdjBCD.ui(path, folder = "AdjBCD")
}
\arguments{
  \item{path}{the path in which a folder used to store variables will be created.}
  \item{folder}{name of the folder. If it is the default, a folder named \code{"AdjBCD"} will be created.}
}
\details{
    See \code{\link{AdjBCD}}.
}
\seealso{
  See \code{\link{AdjBCD}} for allocating patients with complete covariate data; See \code{\link{AdjBCD.sim}} for allocating patients with covariate data generating mechanism.
  
}
\note{
This function provides a command-line user interface, and users should follow the prompts to enter data including covariates as well as levels for each covariate, design parameter \code{a} and the covariate profile of the new patient.
}
\references{
Baldi Antognini A, Zagoraiou M. \emph{The covariate-adaptive biased coin design for balancing clinical trials in the presence of prognostic factors}[J]. Biometrika, 2011, 98(3): 519-535.
}
\value{
It returns an object of \code{\link{class}} \code{"carseq"}.

The function \code{\link{print}} is used to obtain results. The generic accessor functions \code{assignment}, \code{covariate}, \code{cov_num}, \code{cov_profile} and others extract various useful features of the value returned by \code{AdjBCD.ui}. 
}
\keyword{CAR}
\keyword{user-interface}