\name{StrBCD.ui}
\alias{StrBCD.ui}
\alias{StrBCD.ui.carseq}
\title{
Command-line User Interface Using Shao's Method 
}
\description{
A call to user-iterface function used to allocate patients to one of two treatments using Shao's method proposed by Shao J, Yu X, Zhong B (2010) <Doi:10.1093/biomet/asq014>. 
}
\usage{
\method{StrBCD.ui}{carseq}(path, folder = "StrBCD")
}
\arguments{
  \item{path}{path in which a folder used to storage variables would be created.}
  \item{folder}{name of the folder. If default, a folder names \code{"StrBCD"} will be created.}
}
\details{
See \code{\link{StrBCD}}.
}
\seealso{
  See \code{\link{StrBCD}} for allocating patients with complete covariate data; See \code{\link{StrBCD.sim}} for allocating patients with covariate data generating mechanism.
}
\value{
It returns an object of \code{\link{class}} \code{"carseq"}.

The functions \code{\link{print}} is used to obtain results. The generic accessor functions \code{assignment}, \code{covariate}, \code{cov_num}, \code{cov_profile} and so on extract various useful features of the value returned by that function. 
}
\note{
This function provides command-line interface so that users should follow instructions to enter data including covariates as well as levels for each covariate, biased probability \code{p} and covariate-profile of the new patient. 
}
\keyword{ CAR }
\keyword{ user-interface}