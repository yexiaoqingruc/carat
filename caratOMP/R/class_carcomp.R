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
  
  if(!is.null(x$dfmm) && !is.null(x$df_abm)){
    dfmm = x$dfmm; 
    df_abm = x$df_abm; 
    Randomization = df_abm$Randomization; 
    overall = df_abm$overall; 
    pos = dfmm$pos; 
    overallmm = dfmm$overallmm; 
    p1 = ggplot2::ggplot(df_abm) + 
      ggplot2::geom_boxplot(ggplot2::aes(x = Randomization, y = overall, fill = Randomization)) + 
      ggplot2::scale_fill_brewer(palette="Pastel2") + 
      ggplot2::xlab("") + 
      ggplot2::ylab("Mean abs. overall imbalance") + 
      ggplot2::geom_point(data = dfmm, ggplot2::aes(x = pos, y = overallmm), 
                          size = 1, color = "red", shape = 15) + 
      ggplot2::geom_line(data = dfmm, ggplot2::aes(x = pos, y = overallmm), 
                         size = 1, linetype = "dotted", color = "red") + 
      ggplot2::scale_x_discrete(limits = x$mechanism)
    
    within.strt. = df_abm$within.strt.; 
    within.strt.mm = dfmm$within.strt.mm; 
    p3 = ggplot2::ggplot(df_abm) + 
      ggplot2::geom_boxplot(ggplot2::aes(x = Randomization, y = within.strt., fill = Randomization)) + 
      ggplot2::scale_fill_brewer(palette="Pastel2") + 
      ggplot2::xlab("") + 
      ggplot2::ylab("Mean abs. within-strt. imbalance") + 
      ggplot2::geom_point(data = dfmm, ggplot2::aes(x = pos, y = within.strt.mm), 
                          size = 1, color = "red", shape = 15) + 
      ggplot2::geom_line(data = dfmm, ggplot2::aes(x = pos, y = within.strt.mm), 
                         size = 1, linetype = "dotted", color = "red") + 
      ggplot2::scale_x_discrete(limits = x$mechanism) 
    
    marginal = df_abm$marginal; 
    marginalmm = dfmm$marginalmm; 
    p2 = ggplot2::ggplot(df_abm) + 
      ggplot2::geom_boxplot(ggplot2::aes(x = Randomization, y = marginal, fill = Randomization)) + 
      ggplot2::scale_fill_brewer(palette="Pastel2") +
      ggplot2::xlab("") + 
      ggplot2::ylab("Mean abs. within-cov-margin imbalance") + 
      ggplot2::geom_point(data = dfmm, ggplot2::aes(x = pos, y = marginalmm), 
                          size = 1, color = "red", shape = 15) + 
      ggplot2::geom_line(data = dfmm, ggplot2::aes(x = pos, y = marginalmm), 
                         size = 1, linetype = "dotted", color = "red") + 
      ggplot2::scale_x_discrete(limits = x$mechanism) 
    
    position = "bottom";
    ncol = 3; nrow = 1; 
    
    plots = list(p1, p2, p3); 
    g = ggplot2::ggplotGrob(plots[[1]] + 
                              ggplot2::theme(legend.position = position))$grobs; 
    legend = g[[which(sapply(g, function(x) x$name) == "guide-box")]]; 
    lheight = sum(legend$height); 
    lwidth = sum(legend$width); 
    gl = lapply(plots, function(x) x + 
                  ggplot2::theme(legend.position="none")); 
    gl = c(gl, ncol = ncol, nrow = nrow); 
    
    combined = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl),
                                      legend,
                                      ncol = 1,
                                      heights = grid::unit.c(grid::unit(1, "npc") - lheight, lheight));
    
    
    grid::grid.newpage()
    grid::grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
  }
  
  cat("\n")
  invisible(x)
}
