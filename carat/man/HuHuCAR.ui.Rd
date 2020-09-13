\name{HuHuCAR.ui}
\alias{HuHuCAR.ui}
\title{
Command-line User Interface Using Hu and Hu's General Covariate-adaptive Randomization
}
\description{
A call to the user-iterface function used to allocate patients to one of two treatments using Hu and Hu's general covariate-adaptive randomization method as proposed by Hu Y, Hu F (2012) <Doi:10.1214/12-AOS983>. 
}
\usage{
HuHuCAR.ui(path, folder = "HuHuCAR")
}
\arguments{
  \item{path}{the path in which a folder used to store variables will be created. }
  \item{folder}{name of the folder. If default, a folder named \code{"HuHuCAR"} will be created.}
}
\details{
See \code{\link{HuHuCAR}}
}
\seealso{
  See \code{\link{HuHuCAR}} for allocating patients with complete covariate data; See \code{\link{HuHuCAR.sim}} for allocating patients with covariate data generating mechanism.
}
\value{
It returns an object of \code{\link{class}} \code{"carseq"}.

The function \code{\link{print}} is used to obtain results. The generic accessor functions \code{assignment}, \code{covariate}, \code{cov_num}, \code{cov_profile} and others extract various useful features of the value returned by \code{HuHuCAR.ui}. 
}
\note{This function provides a command-line interface so that users should follow the prompts to enter data, including covariates as well as levels for each covariate, weights \code{omega}, biased probability \code{p} and the covariate profile of the new patient. 
}
\keyword{CAR}% use one of  RShowDoc("KEYWORDS")
\keyword{ user-interface}% __ONLY ONE__ keyword per line
