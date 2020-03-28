#' Print Method for Assignment Sequence of Covariate-Adaptive Randomization
#' 
#' Prints the parameters of a covariate-adaptive randomization procedures
#' 
#' @export
#' @rdname print
#' @method print carseq
#' @param x objects of class\code{carseq}.
#' @param digits number of significant digits to be used.
#' @param prefix string, passed to \code{\link{strwrap}} for displaying the \code{method} component of the \code{carandom} object.
#' @param ... further arguments to be passed to or from methods.
#' @seealso \code{\link{compRand}}.


print.carseq = function(x, digits = getOption("digits"), prefix = "\t", ...)
{
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat("group", "=",  LETTERS[1 : 2], "\n", sep = " ")
  cat("Stamps for covariates: ", " ", paste(x$covariate, x$covr, sep = "--"), "\n"); 
  cat("Stamps for levels of each covariate: \n"); 
  for(i in 1 : x$cov_num){
    cat(" ", x$covariate[i], "--", x$covr[i], "\n")
    cat("\t",paste(x[[i]][, 3], x[[i]][, 4], sep = "--"), "\n")
  }
  cat("covariate profile: ", x$cov_profile, "\n"); 
  brid = c("A", "B"); 
  cat("\n");
  cat("assignment: ", brid[x$assignment]);
  cat("\n")
  invisible(x)
}

