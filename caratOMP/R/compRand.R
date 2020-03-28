grid_arrange_shared_legend = function(..., ncol = length(list(...)), 
                                      nrow = 1, 
                                      position = c("bottom", "right")){
  
  plots = list(...)
  g = ggplot2::ggplotGrob(plots[[1]] + 
                            ggplot2::theme(legend.position = position))$grobs
  legend = g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight = sum(legend$height)
  lwidth = sum(legend$width)
  gl = lapply(plots, function(x) x + 
                ggplot2::theme(legend.position="none"))
  gl = c(gl, ncol = ncol, nrow = nrow)
  
  combined = switch(position,
                    "bottom" = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl),
                                           legend,
                                           ncol = 1,
                                           heights = grid::unit.c(grid::unit(1, "npc") - lheight, lheight)),
                    "right" = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl),
                                          legend,
                                          ncol = 2,
                                          widths = grid::unit.c(grid::unit(1, "npc") - lwidth, lwidth)))
  
  grid::grid.newpage()
  grid::grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

compRand.carcomp = function(...) UseMethod("carcomp")

compRand = function(...){
  Objects = list(...); 
  clch = as.character(lapply(Objects, class)); 
  if(length(which(clch != "careval")) >= 1){
    stop("Inputs must be of class 'careval'!")
  }
  DataG = character(); 
  DataType = character(); 
  mechanism = character(); 
  leng = length(Objects); 
  nvec = vector();
  Nvec = vector(); 
  cnvec = vector();
  for(j in 1 : leng){
    R = Objects[[j]]; 
    mechanism[j] = R$method; 
    nvec[j] = R$N; 
    Nvec[j] = R$iteration; 
    DataType[j] = R$`Data Type`;
    cnvec[j] = R$cov_num; 
    if(is.null(R$DataGeneration)){
      DataG[j] = NA;
    }else{
      DataG[j] = R$DataGeneration;
    }
  }
  if(length(unique(cnvec)) == 1){
    cov_num = cnvec[1]; 
    level_num = Objects[[1]]$level_num; 
  }else{
    warning("Results don't make sense: different dataframes are used in different methods."); 
    cov_num = cnvec; 
    level_num = NA; 
  }
  if(length(unique(DataType)) > 1){
    warning("Results don't make sense: comparison between simulated data and real data.")
  }
  if(length(unique(nvec)) > 1){
    n = NA; 
    warning("Results don't make sense: different sample sizes for different methods.")
  }else{
    n = nvec[1];
  }
  if(length(unique(Nvec)) > 1){
    N = NA; 
    N = min(Nvec); 
    warning("Minimum number of iterations is adopted for different methods.")
  }else{
    N = Nvec[1]; 
  }
  #bmax = max(bsize); 
  # cname = vector()
  # cname[1 : 4] = c("max", "95%-quan", "median", "mean");
  # for(j in 1 : bmax){
  #   cname[4 + j] = paste("num", "=", j, seq = " ");
  # }
  
  C_O = C_M = matrix(NA, nrow = leng, ncol = 4); 
  C_S = matrix(NA, nrow = leng, ncol = 4);
  rownames(C_O) = rownames(C_M) = rownames(C_S) = mechanism; 
  colnames(C_O) = colnames(C_M) = colnames(C_S) = c("max", "95%-quan", "median", "mean");
  
  AbM = matrix(NA, nrow = 3, ncol = N * leng); 
  rownames(AbM) = c("overall", "within-strt.", "marginal"); 
  
  for(i in 1 : leng){
    R = Objects[[i]];
    C_O[i, ] = R$Imb[1, ]; 
    C_M[i, ] = apply(R$Imb[(2 + R$strt_num) : (sum(R$level_num) + 1 + R$strt_num), ], 2, mean); 
    C_S[i, 1 : 4] = apply(R$Imb[2 : (1 + R$strt_num), ], 2, mean); 
    #C_S[i, 5 : (4 + R$bsize)] = R$`Within-strt. by num of pats`; 
    DIF = abs(R$DIF); 
    AbM[1, ((i - 1) * N + 1) : (i * N)] = DIF[1, ]; 
    AbM[2, ((i - 1) * N + 1) : (i * N)] = apply(DIF[2 : (1 + R$strt_num), ], 2, mean); 
    AbM[3, ((i - 1) * N + 1) : (i * N)] = apply(DIF[(2 + R$strt_num) :(sum(R$level_num) + 1 + R$strt_num), ], 2, mean);
  }
  
  # meth = rep(mechanism, times = 4);
  # crt = rep(1 : 4, each = leng);
  # mean = c(C_O[, 4], C_M[, 4], C_S[, 4]);
  # df1 = data.frame(meth, levels, crt, mean); 
  # 
  # p1 = ggplot2::ggplot(df1, ggplot2::aes(x = meth, y = mean, color = levels, group = crt,
  #                                        linetype = levels, shape = levels)) +
  #   ggplot2::geom_line(size = 1) +
  #   ggplot2::geom_point(size = 2) +
  #   ggplot2::xlab("method") +
  #   ggplot2::ylab("absolute mean") +
  #   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
  #                  legend.position = "bottom")
  # 
  # val = vector()
  # for(k in 1 : bmax){
  #   val = c(val, C_S[, 4 + k]);
  # }
  # num = rep(cname[5 : (4 + bmax)], each = leng);
  # meths = rep(mechanism, times = bmax)
  # df2 = data.frame("meth" = meths, "num" = num, "val" = val);
  # 
  # p2 = ggplot2::ggplot(df2, ggplot2::aes(x = num, y = val, color = meth, group = meth,
  #                                        linetype = meth, shape = meth)) +
  #   ggplot2::geom_line(size = 1) +
  #   ggplot2::geom_point(size = 2) +
  #   ggplot2::xlab("numbers of patients for each strata") +
  #   ggplot2::ylab("absolute within-stratum mean") +
  #   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), 
  #                  legend.position = "bottom")
  # 
  # p = gridExtra::grid.arrange(p1, p2, ncol = 1);
  # 
  RR = list("Overall Imbalances" = C_O, 
            "Marginal Imbalances" = C_M, "Within-stratum Imbalances" = C_S);
  
  df_abm = data.frame(t(AbM)); 
  Randomization = rep(mechanism, each = N); 
  label = rep(1 :leng, each = N); 
  df_abm$"Randomization" = Randomization; 
  df_abm$"label" = label; 
  
  dfmm = data.frame("Randomization" = mechanism); 
  for(i in 1 : leng){
    ind = which(df_abm$Randomization == mechanism[i]); 
    dfmm[i, 2] = mean(df_abm$overall[ind]); 
    dfmm[i, 3] = mean(df_abm$within.strt.[ind]); 
    dfmm[i, 4] = mean(df_abm$marginal[ind]); 
  }
  colnames(dfmm) = c("Randomizationmm", "overallmm", "within.strt.mm", "marginalmm");
  dfmm$pos = 1 : leng
  p1 = ggplot2::ggplot(df_abm) + 
    ggplot2::geom_boxplot(ggplot2::aes(x = df_abm$Randomization, y = df_abm$overall, fill = df_abm$Randomization)) + 
    ggplot2::scale_fill_brewer(palette="Pastel2") + 
    ggplot2::xlab("") + 
    ggplot2::ylab("Mean abs. overall imbalance") + 
    ggplot2::geom_point(data = dfmm, ggplot2::aes(x = dfmm$pos, y = dfmm$overallmm), 
                        size = 1, color = "red", shape = 15) + 
    ggplot2::geom_line(data = dfmm, ggplot2::aes(x = dfmm$pos, y = dfmm$overallmm), 
                       size = 1, linetype = "dotted", color = "red") + 
    #ggplot2::theme(legend.position = "bottom")
    ggplot2::scale_x_discrete(limits = mechanism)
    
  p3 = ggplot2::ggplot(df_abm) + 
    ggplot2::geom_boxplot(ggplot2::aes(x = df_abm$Randomization, y = df_abm$within.strt., fill = df_abm$Randomization)) + 
    ggplot2::scale_fill_brewer(palette="Pastel2") + 
    ggplot2::xlab("") + 
    ggplot2::ylab("mean abs. within-strt. imbalance") + 
    ggplot2::geom_point(data = dfmm, ggplot2::aes(x = dfmm$pos, y = dfmm$within.strt.mm), 
                        size = 1, color = "red", shape = 15) + 
    ggplot2::geom_line(data = dfmm, ggplot2::aes(x = dfmm$pos, y = dfmm$within.strt.mm), 
                       size = 1, linetype = "dotted", color = "red") + 
    #ggplot2::theme(legend.position = "bottom")
    ggplot2::scale_x_discrete(limits = mechanism) 
  
  p2 = ggplot2::ggplot(df_abm) + 
    ggplot2::geom_boxplot(ggplot2::aes(x = df_abm$Randomization, y = df_abm$marginal, fill = df_abm$Randomization)) + 
    ggplot2::scale_fill_brewer(palette="Pastel2") +
    ggplot2::xlab("") + 
    ggplot2::ylab("mean abs. marginal imbalance") + 
    ggplot2::geom_point(data = dfmm, ggplot2::aes(x = dfmm$pos, y = dfmm$marginalmm), 
                        size = 1, color = "red", shape = 15) + 
    ggplot2::geom_line(data = dfmm, ggplot2::aes(x = dfmm$pos, y = dfmm$marginalmm), 
                       size = 1, linetype = "dotted", color = "red") + 
    #ggplot2::theme(legend.position = "bottom")
    ggplot2::scale_x_discrete(limits = mechanism) 
  
  position = "bottom";
  grid_arrange_shared_legend(p1, p2, p3, ncol = 3, nrow = 1, 
                             position = position)
  
  RR$mechanism = mechanism; 
  RR$N = n;
  RR$iteration = N;
  RR$cov_num = cov_num; 
  RR$level_num = level_num; 
  RR$DataType = DataType; 
  RR$DataGeneration = DataG; 
  
  class(RR) = "carcomp";
  return(RR)
}
