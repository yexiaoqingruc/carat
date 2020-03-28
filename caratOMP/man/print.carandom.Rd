\name{print.carandom}
\alias{print.carandom}
\alias{print.careval}
\alias{print.carcomp}
\alias{print.carseq}
\title{
Print Methods for Sequences, Evaluation and Compariason for Covariate-adaptive Randomization Procedures
}
\description{
Printing objects of class \code{"carandom"}, \code{"careval"}, \code{"carcomp"} or \code{carseq} respectively for different functions, by simple \code{\link{print}} method.
}
\usage{
\method{print}{carandom}(x, digits = getOption("digits"), prefix = "\t", ...)

\method{print}{careval}(x, digits = getOption("digits"), prefix = "\t", ...)

\method{print}{carcomp}(x, digits = getOption("digits"), prefix = "\t", ...)

\method{print}{carseq}(x, digits = getOption("digits"), prefix = "\t", ...)
}
\arguments{
  \item{x}{objects of class \code{"carandom"}, \code{"careval"}, \code{"carcomp"}, or \code{carseq}. }
  \item{digits}{number of significant digits to be used.}
  \item{prefix}{string, passed to \code{\link{strwrap}} for displaying
    the \code{method} component of the \code{carandom} object, \code{careval} object, \code{carcomp} object and \code{carseq} object.}
    \item{\dots}{further arguments to be passed to or from methods.}
}
\value{
the argument \code{x}, invisibly, as for all \code{\link{print}}
  methods.
}
\seealso{
\code{\link{HuHuCAR}}, \code{\link{evalRand}}, \code{\link{compRand}}, 
\code{\link{HuHuCAR.ui}}
}
\keyword{ carandom }
\keyword{ careval }
\keyword{ carcomp }
\keyword{ carseq }
