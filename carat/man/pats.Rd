\name{pats}
\alias{pats}
\docType{data}
\title{
Data of Covariate Profile of Patients
}
\description{
gives the simulated covariate profile of patients for clincal trials. 
}
\usage{
data(pats)
}
\arguments{
  \item{pats}{a dataframe. Each row contains an individual's covariate profile and each column corresponds to a covariate. It contains the following columns
  \describe{
  \item{gender}{Options are \code{male} and \code{female}. }
  \item{employment status}{Options are "unemployment" (\code{unemp}), "part time" (\code{part.}), "full time" (\code{full.}).}
  \item{income}{Options are \code{>= 1w}, \code{<= 0.5w}, \code{0.5~1w}.}
  \item{marriage status}{Options are \code{unmarried}, \code{married}, \code{divorced}}
  }
}
}
\keyword{datasets }