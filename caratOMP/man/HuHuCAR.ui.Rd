\name{HuHuCAR.ui}
\alias{HuHuCAR.ui}
\alias{HuHuCAR.ui.carseq}
\title{
Command-line User Interface Using Hu and Hu's General Covariate-adaptive Randomization
}
\description{
A call to user-iterface function used to allocate patients to one of two treatments using Hu and Hu's general covariate-adaptive randomization proposed by Hu Y, Hu F (2012) <Doi:10.1214/12-AOS983>. 
}
\usage{
\method{HuHuCAR.ui}{carseq}(path, folder = "HuHuCAR")
}
\arguments{
  \item{path}{path in which a folder used to storage variables would be created. }
  \item{folder}{name of the folder. If default, a folder names \code{"HuHuCAR"} will be created.}
}
\details{
See \code{\link{HuHuCAR}}
}
\seealso{
  See \code{\link{HuHuCAR}} for allocating patients with complete covariate data; See \code{\link{HuHuCAR.sim}} for allocating patients with covariate data generating mechanism.
}
\value{
It returns an object of \code{\link{class}} \code{"carseq"}.

The functions \code{\link{print}} is used to obtain results. The generic accessor functions \code{assignment}, \code{covariate}, \code{cov_num}, \code{cov_profile} and so on extract various useful features of the value returned by that function. 
}
\note{This function provides command-line interface so that users should follow instructions to enter data including covariates as well as levels for each covariate, weights \code{omega}, biased probability \code{p} and covariate-profile of the new patient. 
}
\keyword{CAR}% use one of  RShowDoc("KEYWORDS")
\keyword{ user-interface}% __ONLY ONE__ keyword per line
