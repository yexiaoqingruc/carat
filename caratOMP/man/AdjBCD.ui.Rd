\name{AdjBCD.ui}
\alias{AdjBCD.ui}
\alias{AdjBCD.ui.carseq}
\title{
Command-line User Interface Using Covariate-adjusted Biased Coin Design
}
\description{
A call to user-interface function used to allocate patients to one of two treatments using covariate-adjusted biased coin design proposed by Baldi Antognini A, Zagoraiou M (2011) <Doi:10.1093/biomet/asr021>. 
}
\usage{
\method{AdjBCD.ui}{carseq}(path, folder = "AdjBCD")
}
\arguments{
  \item{path}{path in which a folder used to store variables would be created.}
  \item{folder}{name of the folder. If default, a folder names \code{"AdjBCD"} will be created.}
}
\details{
    See \code{\link{AdjBCD}}.
}
\seealso{
  See \code{\link{AdjBCD}} for allocating patients with complete covariate data; See \code{\link{AdjBCD.sim}} for allocating patients with covariate data generating mechanism.
  
}
\note{
This function provides command-line user interface so that users should follow instructions to enter data including covariates as well as levels for each covariate, design parameter \code{a} and covariate-profile of the new patient.
}
\value{
It returns an object of \code{\link{class}} \code{"carseq"}.

The functions \code{\link{print}} is used to obtain results. The generic accessor functions \code{assignment}, \code{covariate}, \code{cov_num}, \code{cov_profile} and so on extract various useful features of the value returned by that function. 
}
\keyword{CAR}
\keyword{user-interface}