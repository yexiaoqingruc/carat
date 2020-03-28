\name{PocSimMIN.ui}
\alias{PocSimMIN.ui}
\alias{PocSimMIN.ui.carseq}
\title{
Command-line User Interface Using Pocock and Simon's Procedure with Two-arm Case
}
\description{
A call to user-iterface function used to allocate patients to one of two treatments using Pocock and Simon's method proposed by Pocock S J, Simon R (1975) <Doi:10.2307/2529712>. 
}
\usage{
\method{PocSimMIN.ui}{carseq}(path, folder = "PocSimMIN")
}
\arguments{
  \item{path}{path in which a folder used to storage variables would be created. }
  \item{folder}{name of the folder. If default, a folder names \code{"PocSimMIN"} will be created.}
}
\details{
See \code{\link{PocSimMIN}}.
}
\note{
This function provides command-line interface so that users should follow instructions to enter data including covariates as well as levels for each covariate, \code{weight}, biased probability \code{p} and covariate-profile of the new patient. 
}
\seealso{
  See \code{\link{PocSimMIN}} for allocating a given completely collected data; See \code{\link{PocSimMIN.sim}} for allocating patients with covariate data generating mechanism.
}
\value{
It returns an object of \code{\link{class}} \code{"carseq"}.

The functions \code{\link{print}} is used to obtain results. The generic accessor functions \code{assignment}, \code{covariate}, \code{cov_num}, \code{cov_profile} and so on extract various useful features of the value returned by that function. 
}
\keyword{ CAR }% use one of  RShowDoc("KEYWORDS")
\keyword{ user-interface }% __ONLY ONE__ keyword per line

