\name{PocSimMIN.ui}
\alias{PocSimMIN.ui}
\title{
Command-line User Interface Using Pocock and Simon's Procedure with Two-Arms Case
}
\description{
A call to the user-iterface function used to allocate patients to one of two treatments using Pocock and Simon's method proposed by Pocock S J, Simon R (1975) <Doi:10.2307/2529712>. 
}
\usage{
PocSimMIN.ui(path, folder = "PocSimMIN")
}
\arguments{
  \item{path}{the path in which a folder used to storage variables will be created. }
  \item{folder}{name of the folder. If default, a folder named \code{"PocSimMIN"} will be created.}
}
\details{
See \code{\link{PocSimMIN}}.
}
\note{
This function provides a command-line interface and users should follow the prompts to enter data including covariates as well as levels for each covariate, \code{weight}, biased probability \code{p} and the covariate profile of the new patient. 
}
\seealso{
  See \code{\link{PocSimMIN}} for allocating a given completely collected data; See \code{\link{PocSimMIN.sim}} for allocating patients with covariate data generating mechanism.
}
\value{
It returns an object of \code{\link{class}} \code{"carseq"}.

The function \code{\link{print}} is used to obtain results. The generic accessor functions \code{assignment}, \code{covariate}, \code{cov_num}, \code{cov_profile} and others extract various useful features of the value returned by \code{PocSimMIN.ui}. 
}
\keyword{ CAR }% use one of  RShowDoc("KEYWORDS")
\keyword{ user-interface }% __ONLY ONE__ keyword per line

