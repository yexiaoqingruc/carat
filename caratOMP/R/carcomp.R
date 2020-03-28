#' Print Method for Comparison of Covariate-Adaptive Randomization
#' 
#' Prints the parameters of a covariate-adaptive randomization procedures
#' 
#' @export
#' @rdname print
#' @method print carcomp
#' @param x objects of class\code{carcomp}.
#' @param digits number of significant digits to be used.
#' @param prefix string, passed to \code{\link{strwrap}} for displaying the \code{method} component of the \code{carandom} object.
#' @param ... further arguments to be passed to or from methods.
#' @seealso \code{\link{compRand}}.



print.carcomp = function(x, digits = getOption("digits"), prefix = "\t", ...){
  cat("\n")
  # abb = c("HH", "PocSim", "Shao", "StraRand", "AtkinBCD", "BayesBCD", "AdBCD")
  # com = c("Hu and Hu's General CAR", "Pocock and Simon's Procedure with Two Arms", 
  #         "Shao's Procedure", "Stratified Randomization with Two Arms", 
  #         "Atkinson's Optimum Biased Coin Design", "Bayesian Biased Coin Design", 
  #         "Covariate-adjusted Biased Coin Design")
  # 
  # for(k in 1 : length(x$mechanism)){
  #   meth[k] = com[which(abb == x$mechanism[k], arr.ind = T)]
  # }
  cat("Comparison: ", sep = "\n")
  meth = character()
  meth[1] = x$mechanism[1]
  for(k in 2 : length(x$mechanism)){
    meth = paste(meth, x$mechanism[k], sep = ", ")
  }
  cat("Randomization", "=", meth, "\n", sep = " ")
  if(length(unique(x$DataType)) != 1){
    cat("Data Type: ", x$DataType, "\n"); 
  }else{
    cat("Data Type: ", x$DataType[1], "\n"); 
    if(!"TRUE" %in% is.na(x$DataGeneration)){
      cat("Data generation: ", x$DataGeneration[1], "\n"); 
    }
  }
  cat("group", "=",  LETTERS[1 : 2], "\n", sep = " ")
  if(is.na(x$N)){
    cat("iteration: no consistency between methods.", "\n")
  }else{
    cat("N", "=", x$N, "\n", sep = " ")
  }
  if(is.na(x$iteration)){
    cat("iteration: no consistency between methods.", "\n")
  }else{
    cat("iteration", "=", x$iteration, "\n", sep = " ")
  }
  cat("cov_num", "=", x$cov_num, "\n", sep = " ")
  if(length(is.na(x$level_num)) == 1 && is.na(x$level_num)){
    cat("level_num: ", "no consistency between methods.", "\n"); 
  }else{
    cat("level_num", "=", as.character(x$level_num), "\n", sep = " ")
  }
  
  cat("\n")
  
  cat("absolute mean at overall, marginal and within-strt. levels:\n")
  cat("Overall:", "\n")
  print(x$`Overall Imbalances`, digits = 3)
  cat("\n")
  cat("Marginal:\n")
  print(x$`Marginal Imbalances`, digits = 3)
  cat("\n")
  cat("Within-strt.:\n")
  print(x$`Within-stratum Imbalances`, digits = 3)
  cat("\n")
  invisible(x)
}
