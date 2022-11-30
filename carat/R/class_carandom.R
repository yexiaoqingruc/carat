#' Print Method for Covariate-Adaptive Randomization
#' 
#' Prints the parameters of a covariate-adaptive randomization procedures
#' 
#' @export
#' @rdname print
#' @method print carandom
#' @param x objects of class\code{carandom}.
#' @param digits number of significant digits to be used.
#' @param prefix string, passed to \code{\link{strwrap}} for displaying the \code{method} component of the \code{carandom} object.
#' @param ... further arguments to be passed to or from methods.
#' @seealso \code{\link{HuHuCAR}}, \code{\link{HuHuCAR.sim}}.

print.carandom = function(x, digits = getOption("digits"), prefix = "\t", ...)
{
  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat("Data: ", x$`Data Type`, "\n", sep = " ")
  cat("group", "=",  LETTERS[1 : 2], "\n", sep = " ")
  if(!is.null(x$bsize)){
    cat("block", "=", x$bsize, "\n", sep = " ")
  }
  cat("Sample size", "=", x$n, "\n", sep = " ")
  cat("cov_num", "=", x$cov_num, "\n", sep = " ")
  if(x$`Data Type` == "Real"){
    cat("considered covariates: ", x$covariates, "\n", sep = "  ")
  }
  cat("level_num", "=", as.character(x$level_num), "\n", sep = " ")
  if(!x$framework %in% c("Model-based approach", "Stratified randomization")){
    if(is.null(x$weight)){
      cat("set priority to (margins): ", "excluded from consideration", "\n", sep = " " );
    }else{
      impot = character(); 
      impot = which(x$weight == max(x$weight), arr.ind = T); 
      if(length(impot) == x$cov_num){
        cat("set priority to (margins): ", "equally important", "\n", sep = " " )
      }else{
        if(x$`Data Type` == "Real"){
          cat("set priority to (margins): ",
              "Covariate", x$covariates[impot], "\n", sep = " "); 
        }else{
          cat("set priority to (margins): ", "Covariate", impot, "\n", sep = " "); 
        }
      }
    }
  }
  cat("\n")
  CA = as.data.frame(t(x$Cov_Assig[, 1 : 3]))
  CA$assignment = LETTERS[CA$assignment]
  cat("the first three patients' covariate-profiles and assignments:\n")
  print(CA); #cat("...\n")
  cat("\n")
  out = vector()
  om = abs(t(x$Diff)[1, 1]); sm = mean(abs(t(x$Diff)[1, 2 : (1 + x$strt_num)]))
  mm = mean(abs(t(x$Diff)[(2 + x$strt_num) : (1 + x$strt_num + sum(x$level_num))]))
  if(x$method == "Stratified Permuted Block Randomization"){
    out[1] = om; out[2] = mm;
    par = x$`numbers of pats for each stratum` %% x$bsize
    str = t(x$Diff)[1, 2 : (1 + x$strt_num)]
    index = c(x$bsize, 1 : (x$bsize - 1))
    for(i in 1 : x$bsize){
      ind = which(par == i - 1, arr.ind = T)
      out[index[i] + 2] = mean(abs(str[ind[, 2]]))
    }
    names(out) = c("overall", "within-cov.-margin", BBCDname(x$bsize, "pnum = "))
  }else{
    out[1 : 3] = c(om, sm, mm)
    names(out) = c("overall", "within-strt.", "within-cov.-margin")
  }
  cat("Mean absolute imbalances at overall, within-strt., and within-cov.-margin levels:\n")
  print(out, digits = 4)
  cat("\n")
  if(x$`Data Type` == "Real"){
    Rlist = apply(x$data, 2, unique); 
    cat("Remark-Index: \n"); 
    if(!x$datanumeric){
      for(i in 1 : x$cov_num){
        cat(i, "--", x$covariates[i], "\n"); 
        if(length(unique(x$level_num)) > 1){
          cat("\t", paste(paste(1 : x$level_num[i], as.factor(Rlist[[i]]), 
                          sep = " <--> "), collapse = "; "), 
              sep = "  ", "\n"); 
        }else{
          cat("\t", paste(paste(1 : x$level_num[i], as.factor(Rlist[, i]), 
                          sep = " <--> "), collapse = "; "), 
              sep = "  ", "\n"); 
        }
      }
    }else{
      for(i in 1 : x$cov_num){
        cat(i, "--", x$covariates[i], "\n"); 
        if(length(unique(x$level_num)) > 1){
          cat("\t", paste(paste(1 : x$level_num[i], 
                          as.factor(Rlist[[i]])[match(1 : x$level_num[i], as.numeric(as.factor(Rlist[[i]])))], 
                          sep = " <--> "), collapse = "; "), 
              sep = "  ", "\n"); 
        }else{
          cat("\t", paste(paste(1 : x$level_num[i], 
                          as.factor(Rlist[, i])[match(1 : x$level_num[i], as.numeric(as.factor(Rlist[, i])))], 
                          sep = " <--> "), collapse = "; "), 
              sep = "  ", "\n"); 
        }
      }
    }
  }
  cat("\n");
  invisible(x)
}
