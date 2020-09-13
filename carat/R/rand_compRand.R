
#compRand.carcomp = function(...) UseMethod("carcomp")

compRand = function(...){
  Objects = list(...); 
  clch = as.numeric(lapply(Objects, methods::is, class2 = "careval")); 
  if(0 %in% clch){
    stop("Inputs must be of class 'careval'!"); 
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
  
  
  RR = list("Overall Imbalances" = C_O, 
            "Marginal Imbalances" = C_M, "Within-stratum Imbalances" = C_S);
  
  df_abm = data.frame(t(AbM), stringsAsFactors = TRUE); 
  Randomization = rep(mechanism, each = N); 
  label = rep(1 :leng, each = N); 
  df_abm$"Randomization" = Randomization; 
  df_abm$"label" = label; 
  
  dfmm = data.frame("Randomization" = mechanism, stringsAsFactors = TRUE); 
  for(i in 1 : leng){
    ind = which(df_abm$Randomization == mechanism[i]); 
    dfmm[i, 2] = mean(df_abm$overall[ind]); 
    dfmm[i, 3] = mean(df_abm$within.strt.[ind]); 
    dfmm[i, 4] = mean(df_abm$marginal[ind]); 
  }
  colnames(dfmm) = c("Randomizationmm", "overallmm", "within.strt.mm", "marginalmm");
  dfmm$pos = 1 : leng; 
  
  overall = df_abm$overall; 
  pos = dfmm$pos; 
  overallmm = dfmm$overallmm; 
  within.strt. = df_abm$within.strt.; 
  within.strt.mm = dfmm$within.strt.mm; 
  marginal = df_abm$marginal; 
  marginalmm = dfmm$marginalmm; 
  
  RR$dfmm = dfmm; 
  RR$df_abm = df_abm; 
  
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
