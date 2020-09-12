#' Print Method for Evaluation of Covariate-Adaptive Randomization
#' 
#' Prints the parameters of a covariate-adaptive randomization procedures
#' 
#' @export
#' @rdname print
#' @method print careval
#' @param x objects of class\code{careval}.
#' @param digits number of significant digits to be used.
#' @param prefix string, passed to \code{\link{strwrap}} for displaying the \code{method} component of the \code{carandom} object.
#' @param ... further arguments to be passed to or from methods.
#' @seealso \code{\link{evalRand}}, \code{\link{evalRand.sim}}.

print.careval = function(x, digits = getOption("digits"), prefix = "\t", ...){
  cat("\n")
  abb = c("HuHuCAR", "PocSimMIN", "StrBCD", "StrPBR", "DoptBCD", "BayesBCD", "AdjBCD")
  com = c("Hu and Hu's General CAR", "Pocock and Simon's Procedure with Two Arms", 
          "Shao's Procedure", "Stratified Randomization with Two Arms", 
          "Atkinson's Optimum Biased Coin Design", "Bayesian Biased Coin Design", 
          "Covariate-adjusted Biased Coin Design")
  ind = which(abb == x$method, arr.ind = T)
  meth = com[ind]
  cat(strwrap(meth, prefix = prefix), sep = "\n")
  cat("\n")
  cat("call:\n", 
      paste("evalRand.sim(", "method = ", x$method, ")\n", sep = ""));
  cat("\n"); 
  cat("group", "=",  LETTERS[1 : 2], "\n", sep = " ")
  cat("N", "=", x$N, "\n", sep = " ")
  cat("iteration", "=", x$iteration, "\n", sep = " ")
  cat("cov_num", "=", x$cov_num, "\n", sep = " ")
  cat("level_num", "=", as.character(x$level_num), "\n", sep = " ")
  if(x$method == "BayesBCD"){
    cat("Categor class numbers", "=", x$J, "\n", sep = " ")
  }
  if(x$method == "StrPBR"){
    cat("block", "=", x$bsize, "\n", sep = " ")
  }
  cat("Data type: ", x$`Data Type`, "\n"); 
  
  if(x$`Data Type` == "Simulated"){
    cat("Data generation mode: ", x$DataGeneration, "\n", sep = " ")
  }
  
  cat("\n")
  if(x$N <= 7){K = x$N}else{K = 7}
  if(x$iteration <= 3){I = x$iteration}else{I = 3}
  cat("assignments of the first", I, "iterations for the first", K, 
      "patients", ":", "\n", sep = " ")
  ass = as.data.frame(t(x$Assig[1 : K, 1 : I]))
  for(l in 1 : I){
    ass[l, ] = LETTERS[as.numeric(ass[l, ])]
  }
  ass$' ' = rep("...", times = I)
  print(ass)
  cat("\n")
  cat("evaluation by imbalances: \n"); 
  cat("absolute overall imbalances:\n")
  print(x$Imb[1, ], digits = 3);
  cat("\n"); 
  if(x$strt_num <= 3){s = x$strt_num}else{s = 3}
  cat("absolute within-strt. imbalances for the first", s, "strata:", "\n", sep = " "); 
  print(x$Imb[2 : (s + 1), ], digits = 3)
  cat("\n"); 
  cat("absolute marginal imbalances for", x$cov_num, "margins:", "\n", sep = " "); 
  v = vector(); 
  r = 1 + x$strt_num + 1; 
  for(i in 1 :x$cov_num){
    v[i] = r; 
    r = r + x$level_num[i]; 
  }
  print(x$Imb[v, ], digits = 3); 
  cat("\n")
  invisible(x)
}
