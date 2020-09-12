\name{StrPBR.ui}
\alias{StrPBR.ui}
\title{
Command-line User Interface Using Stratified Permuted Block Randomization with Two-Arms Case
}
\description{
A call to the user-iterface function used to allocate patients to one of two treatments using stratified permuted block randomization proposed by Zelen M (1974) <Doi: 10.1016/0021-9681(74)90015-0>. 
}
\usage{
StrPBR.ui(path, folder = "StrPBR")
}
\arguments{
  \item{path}{the path in which a folder used to storage variables will be created. }
  \item{folder}{name of the folder. If default, a folder named \code{"StrPBR"} will be created.}
}
\details{
See \code{\link{StrPBR}}.
}
\seealso{
 See \code{\link{StrPBR}} for allocating patients with complete covariate data; See \code{\link{StrPBR.sim}} for allocating patients with covariate data generating mechanism.
}
\note{
This function provides a command-line interface and users should follow the prompts to enter data including covariates as well as levels for each covariate, block size \code{bsize} and the covariate profile of the new patient.
}
\value{
It returns an object of \code{\link{class}} \code{"carseq"}.

The function \code{\link{print}} is used to obtain results. The generic accessor functions \code{assignment}, \code{covariate}, \code{cov_num}, \code{cov_profile} and others extract various useful features of the value returned by \code{StrPBR.ui}. 
}
\keyword{ CAR }
\keyword{user-interface}
