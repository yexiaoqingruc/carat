\name{StrBCD.ui}
\alias{StrBCD.ui}
\title{
Command-line User Interface Using Shao's Method 
}
\description{
A call to the user-interface function used to allocate patients to one of two treatments using Shao's method proposed by Shao J, Yu X, Zhong B (2010) <Doi:10.1093/biomet/asq014>. 
}
\usage{
StrBCD.ui(path, folder = "StrBCD")
}
\arguments{
  \item{path}{the path in which a folder used to storage variables will be created. }
  \item{folder}{name of the folder. If default, a folder named \code{"StrBCD"} will be created.}
}
\details{
See \code{\link{StrBCD}}.
}
\seealso{
  See \code{\link{StrBCD}} for allocating patients with complete covariate data; See \code{\link{StrBCD.sim}} for allocating patients with covariate data generating mechanism.
}
\value{
It returns an object of \code{\link{class}} \code{"carseq"}.

The function \code{\link{print}} is used to obtain results. The generic accessor functions \code{assignment}, \code{covariate}, \code{cov_num}, \code{cov_profile} and others extract various useful features of the value returned by \code{StrBCD.ui}. 
}
\note{
This function provides a command-line interface and users should follow the prompts to enter data including covariates as well as levels for each covariate, biased probability \code{p} and the covariate profile of the new patient. 
}
\keyword{ CAR }
\keyword{ user-interface}