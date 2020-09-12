#' Print Method for the Plot of Randomization Test
#' 
#' Prints the plot of Randomization Test
#' 
#' @export
#' @rdname print
#' @method print plotrandtest
#' @param x objects of class\code{plotrandtest}.
#' @param digits number of significant digits to be used.
#' @param prefix string, passed to \code{\link{strwrap}} for displaying the \code{method} component of the \code{plottest} object.
#' @param ... further arguments to be passed to or from methods.


print.plotrandtest <- function(x, digits = getOption("digits"), prefix = "\t", ...)
{
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat("data:  ", x$data.name, "\n", sep = "")
  out <- character()
  if(!is.null(x$statistic))
    out <- c(out, paste(names(x$statistic), "=",
                        format(x$statistic, digits = max(1L, digits - 2L))))
  if(!is.null(x$parameter))
    out <- c(out, paste(names(x$parameter), "=",
                        format(x$parameter, digits = max(1L, digits - 2L))))
  if(!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = max(1L, digits - 3L))
    out <- c(out, paste("p-value",
                        if(startsWith(fp, "<")) fp else paste("=",fp)))
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  if(!is.null(x$alternative)) {
    cat("alternative hypothesis: ")
    if(!is.null(x$null.value)) {
      if(length(x$null.value) == 1L) {
        alt.char <-
          switch(x$alternative,
                 two.sided = "not equal to",
                 less = "less than",
                 greater = "greater than")
        cat("true ", names(x$null.value), " is ", alt.char, " ",
            x$null.value, "\n", sep = "")
      }
      else {
        cat(x$alternative, "\nnull values:\n", sep = "")
        print(x$null.value, digits=digits, ...)
      }
    }
    else cat(x$alternative, "\n", sep = "")
  }
  if(!is.null(x$conf.int)) {
    cat(format(100 * attr(x$conf.int, "conf.level")),
        " percent confidence interval:\n", " ",
        paste(format(x$conf.int[1:2], digits=digits), collapse = " "),
        "\n", sep = "")
  }
  if(!is.null(x$estimate)) {
    cat("sample estimates:\n")
    print(x$estimate, digits=digits, ...)
  }
  if(!is.null(x$data)){
    pic<-ggplot2::ggplot(data = data.frame(x$data, stringsAsFactors = TRUE), 
                         ggplot2::aes(x = x$data))+
      ggplot2::geom_histogram(bins = x$binwidth)+
      ggplot2::geom_vline(ggplot2::aes(xintercept = x$estimate),colour = "#990000",linetype = "dashed")+
      ggplot2::xlab("Difference in means")+
      ggplot2::ylab("Frequency")
    print(pic)
  }
  cat("\n")
  invisible(x)
}
