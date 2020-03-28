\name{DoptBCD.ui}
\alias{DoptBCD.ui}
\alias{DoptBCD.ui.carseq}
\title{
Command-line User Interface Using Atkinson's \eqn{D_A}-optimal Biased Coin Design
}
\description{
A call to user-interface function used to allocate patients to one of two treatments using Atkinson's \eqn{D_A}-optimal biased coin design proposed by Atkinson A C (1982) <Doi:10.2307/2335853>. 
}
\usage{
\method{DoptBCD.ui}{carseq}(path, folder = "DoptBCD")
}
\arguments{
  \item{path}{path in which a folder used to storage variables would be created. }
  \item{folder}{name of the folder. If default, a folder names \code{"DoptBCD"} will be created.}
}
\details{
See \code{\link{DoptBCD}}.
}
\seealso{
  See \code{\link{DoptBCD}} for allocating patients with complete covariate data; See \code{\link{DoptBCD.sim}} for allocating patients with covariate data generating mechanism.
}
\note{
This function provides command-line interface so that users should follow instructions to enter data including covariates as well as levels for each covariate and covariate-profile of the new patient.
}
\value{
It returns an object of \code{\link{class}} \code{"carseq"}.

The functions \code{\link{print}} is used to obtain results. The generic accessor functions \code{assignment}, \code{covariate}, \code{cov_num}, \code{cov_profile} and so on extract various useful features of the value returned by that function. 
}
\keyword{ CAR }
\keyword{ user-interface }
