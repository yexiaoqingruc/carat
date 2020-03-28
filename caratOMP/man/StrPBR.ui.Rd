\name{StrPBR.ui}
\alias{StrPBR.ui}
\alias{StrPBR.ui.carseq}
\title{
Command-line User Interface Using Stratified Permuted Block Randomization with Two-arm Case
}
\description{
A call to user-iterface function used to allocate patients to one of two treatments using stratified permuted block randomization proposed by Zelen M (1974) <Doi: 10.1016/0021-9681(74)90015-0>. 
}
\usage{
\method{StrPBR.ui}{carseq}(path, folder = "StrPBR")
}
\arguments{
  \item{path}{path in which a folder used to storage variables would be created. }
  \item{folder}{name of the folder. If default, a folder names \code{"StrPBR"} will be created.}
}
\details{
See \code{\link{StrPBR}}.
}
\seealso{
 See \code{\link{StrPBR}} for allocating patients with complete covariate data; See \code{\link{StrPBR.sim}} for allocating patients with covariate data generating mechanism.
}
\note{
This function provides command-line interface so that users should follow instructions to enter data including covariates as well as levels for each covariate, block size \code{bsize} and covariate-profile of the new patient.
}
\value{
It returns an object of \code{\link{class}} \code{"carseq"}.

The functions \code{\link{print}} is used to obtain results. The generic accessor functions \code{assignment}, \code{covariate}, \code{cov_num}, \code{cov_profile} and so on extract various useful features of the value returned by that function. 
}
\keyword{ CAR }
\keyword{user-interface}
