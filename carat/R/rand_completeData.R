###############################################################################
#############################   Real Data   ###################################
###############################################################################
## Hu and Hu's general covariate-adaptive randomization(CAR) ########
################################################################################
#HuHuCAR = function(data, omega = NULL, p = 0.85) UseMethod("HuHuCAR")

#HuHuCAR.carandom = function(data, omega = NULL, p = 0.85) UseMethod("carandom")

HuHuCAR = function(data, omega = NULL, p = 0.85){
  
  pdoub = as.double(p); 
  if(is.na(pdoub)){
    stop("p must be a positive number!");
  }else if(p <= 0.5){
    stop("set p larger than 0.5 to achieve balance!");
  }else if(p > 1){
    stop("p must be a positive number between 0 and 1!");
  }
  
  if(length(data[is.na(data)]) == 0){
    datap = data;
  }else{
    data[is.na(data)] = "HFFMWYXQTFY<= 1.0";
    datap = data;
  }
  
  rdata = Preprocess(datap); 
  data_proc = rdata$data; 
  cov_num = rdata$cov_num; level_num = rdata$level_num; 
  
  if(length(omega) != (2 + cov_num) && !is.null(omega)){
    stop("Length of omega must equal to ncols(data) + 2 !")
  }else if(is.null(omega)){
    omega = rep(1.0 / (cov_num + 2), times = cov_num + 2);
  }else{
    omega = abs(omega) / sum(abs(omega)); 
  }
  
  RES = C_RHPS(data_proc, cov_num, level_num, omega, p); 
  
  R = NULL;
  
  covn = colnames(data);
  R$covariates = covn;
  
  strt_num = ncol(RES[2, 1][[1]]);
  R$strt_num = strt_num; 
  
  R$cov_num = cov_num; 
  R$level_num = level_num; 
  
  CA = RES[3, 1][[1]]; 
  n = ncol(CA); 
  R$N = n;
  colnames(CA) = BBCDname(n, "pat"); 
  rownames(CA) = c(BBCDname(cov_num, "covariate"), "assignment"); 
  R$Cov_Assig = CA;
  
  assig_temp = CA[dim(CA)[1], ]; 
  R$assignments = LETTERS[assig_temp]; 
  
  AS = RES[2, 1][[1]];
  colnames(AS) = BBCDname(strt_num, "strt.");
  rownames(AS) = BBCDname(cov_num, "covariate"); 
  R$'All strata' = AS;
  
  Df = RES[4, 1][[1]]; 
  rownames(Df) = nameString(cov_num, level_num, strt_num, "All", "Real"); 
  R$Diff = t(Df);
  
  R$method = "Hu and Hu's General CAR";
  R$'Data Type' = "Real";
  R$weight = omega[3 : (2 + cov_num)];
  R$framework = "Minimization";
  R$data = data; 
  
  class(R) = "carandom";
  
  return(R);
}

###############################################################################
## Pocock and Simon's procedure ########
################################################################################
#PocSimMIN = function(data, weight = NULL, p = 0.85) UseMethod("PocSimMIN")

#PocSimMIN.carandom = function(data, weight = NULL, p = 0.85) UseMethod("carandom")

PocSimMIN = function(data, weight = NULL, p = 0.85){
  
  pdoub = as.double(p); 
  if(is.na(pdoub)){
    stop("p must be a positive number!");
  }else if(p <= 0.5){
    stop("set p larger than 0.5 to achieve balance!");
  }else if(p > 1){
    stop("p must be a positive number between 0 and 1!");
  }
  
  if(length(data[is.na(data)]) == 0){
    datap = data;
  }else{
    data[is.na(data)] = "NA";
    datap = data;
  }
  
  rdata = Preprocess(datap); 
  data_proc = rdata$data; 
  cov_num = rdata$cov_num; level_num = rdata$level_num; 
  
  if(length(weight) != cov_num && !is.null(weight)){
    stop("Length of weight must equal to ncols(data)!")
  }else if(is.null(weight)){
    omega = c(0, 0, rep(1.0 / cov_num, times = cov_num)); 
  }else{
    omega = c(0, 0, abs(weight) / sum(abs(weight)));
  }
  
  RES = C_RHPS(data_proc, cov_num, level_num, omega, p); 
  
  R = NULL;
  
  covn = colnames(data);
  R$covariates = covn;
  
  strt_num = ncol(RES[2, 1][[1]]);
  R$strt_num = strt_num; 
  
  R$cov_num = cov_num; 
  
  R$level_num = level_num; 
  
  CA = RES[3, 1][[1]]; 
  n = ncol(CA); 
  R$N = n;
  colnames(CA) = BBCDname(n, "pat"); 
  rownames(CA) = c(BBCDname(cov_num, "covariate"), "assignment"); 
  R$Cov_Assig = CA;
  
  assig_temp = CA[dim(CA)[1], ]; 
  R$assignments = LETTERS[assig_temp]; 
  
  AS = RES[2, 1][[1]];
  colnames(AS) = BBCDname(strt_num, "strt.");
  rownames(AS) = BBCDname(cov_num, "covariate"); 
  R$'All strata' = AS;
  
  Df = RES[4, 1][[1]]; 
  rownames(Df) = nameString(cov_num, level_num, strt_num, "All", "Real"); 
  R$Diff = t(Df);
  
  R$method = "Pocock and Simon's Procedure with Two Arms";
  R$'Data Type' = "Real";
  R$weight = weight;
  R$framework = "Minimization";
  R$data = data; 
  
  class(R) = "carandom";
  
  return(R);
}

###############################################################################
## Shao's randomization ########
################################################################################
#StrBCD = function(data, p = 0.85) UseMethod("StrBCD")

#StrBCD.carandom = function(data, p = 0.85) UseMethod("carandom")

StrBCD = function(data, p = 0.85){
  
  pdoub = as.double(p); 
  if(is.na(pdoub)){
    stop("p must be a positive number!");
  }else if(p <= 0.5){
    stop("set p larger than 0.5 to achieve balance!");
  }else if(p > 1){
    stop("p must be a positive number between 0 and 1!");
  }
  
  if(length(data[is.na(data)]) == 0){
    datap = data;
  }else{
    data[is.na(data)] = "NA";
    datap = data;
  }
  
  rdata = Preprocess(datap); 
  data_proc = rdata$data; 
  cov_num = rdata$cov_num; level_num = rdata$level_num; 
  omega = c(0, 1, rep(0, times = cov_num)); 
  
  RES = C_RHPS(data_proc, cov_num, level_num, omega, p); 
  
  R = NULL;
  
  covn = colnames(data);
  R$covariates = covn;
  
  strt_num = ncol(RES[2, 1][[1]]);
  R$strt_num = strt_num; 
  
  R$cov_num = cov_num; 
  R$level_num = level_num; 
  
  CA = RES[3, 1][[1]]; 
  n = ncol(CA); 
  R$N = n;
  colnames(CA) = BBCDname(n, "pat"); 
  rownames(CA) = c(BBCDname(cov_num, "covariate"), "assignment"); 
  R$Cov_Assig = CA; 
  
  assig_temp = CA[dim(CA)[1], ]; 
  R$assignments = LETTERS[assig_temp]; 
  
  AS = RES[2, 1][[1]];
  colnames(AS) = BBCDname(strt_num, "strt.");
  rownames(AS) = BBCDname(cov_num, "covariate"); 
  R$'All strata' = AS;
  
  Df = RES[4, 1][[1]]; 
  rownames(Df) = nameString(cov_num, level_num, strt_num, "All", "Real"); 
  R$Diff = t(Df);
  
  R$method = "Shao's Procedure";
  R$'Data Type' = "Real";
  R$framework = "Minimization/Stratified randomization";
  R$data = data; 
  
  class(R) = "carandom";
  
  return(R);
}

###############################################################################
## Stratified permuted block randomization ########
################################################################################
#StrPBR = function(data, bsize = 4) UseMethod("StrPBR")

#StrPBR.carandom = function(data, bsize = 4) UseMethod("carandom")

StrPBR = function(data, bsize = 4){
  
  bsint = as.integer(bsize);
  if(is.na(bsint)){
    stop("bsize must be a positive integer!")
  }else if(bsint %% 2 != 0){
    stop("bsize must be a multiple of 2!")
  }else if(bsize %% 2 != 0){
    warning("bsize is mandated to be an integer!", call. = FALSE); 
  }
  
  if(length(data[is.na(data)]) == 0){
    datap = data;
  }else{
    data[is.na(data)] = "NA"; 
    datap = data; 
  }
  
  rdata = Preprocess(datap); 
  data_proc = rdata$data; 
  cov_num = rdata$cov_num; level_num = rdata$level_num; 
  
  RES = C_RStrR(data_proc, cov_num, level_num, bsize, tr_num = 2); 
  
  R = NULL;
  
  covn = colnames(data);
  R$covariates = covn;
  
  strt_num = ncol(RES[2, 1][[1]]);
  R$strt_num = strt_num; 
  
  R$cov_num = cov_num; 
  R$level_num = level_num; 
  
  CA = RES[3, 1][[1]]; 
  n = ncol(CA); 
  R$N = n;
  colnames(CA) = BBCDname(n, "pat"); 
  rownames(CA) = c(BBCDname(cov_num, "covariate"), "assignment"); 
  R$Cov_Assig = CA; 
  
  assig_temp = CA[dim(CA)[1], ]; 
  R$assignments = LETTERS[assig_temp]; 
  
  AS = RES[2, 1][[1]];
  colnames(AS) = BBCDname(strt_num, "strt.");
  rownames(AS) = BBCDname(cov_num, "covariate"); 
  R$'All strata' = AS;
  
  Df = RES[4, 1][[1]]; 
  rownames(Df) = nameString(cov_num, level_num, strt_num, "All", "Real"); 
  R$Diff = t(Df);
  
  R$method = "Stratified Permuted Block Randomization";
  R$'Data Type' = "Real";
  R$framework = "Stratified randomization"; 
  R$data = data; 
  R$bsize = bsint; 
  
  st_num = RES[1, 1][[1]];
  colnames(st_num) = BBCDname(ncol(AS), "level-");
  
  R$`numbers of pats for each strata` = st_num;
  
  class(R) = "carandom";
  
  return(R);
}

###############################################################################
## Atkinson's Optimum Biased Coin Design ########
################################################################################
#DoptBCD = function(data) UseMethod("DoptBCD")

#DoptBCD.carandom = function(data) UseMethod("carandom")

DoptBCD = function(data){
  
  if(length(data[is.na(data)]) == 0){
    datap = data;
  }else{
    data[is.na(data)] = "NA"; 
    datap = data; 
  }
  
  rdata = Preprocess(datap); 
  data_proc = rdata$data; 
  cov_num = rdata$cov_num; level_num = rdata$level_num; 
  
  RES = C_RAtkinBCD(data_proc, cov_num, level_num); 
  
  R = NULL;
  
  covn = colnames(data);
  R$covariates = covn;
  
  strt_num = ncol(RES[2, 1][[1]]);
  R$strt_num = strt_num; 
  
  R$cov_num = cov_num; 
  R$level_num = level_num; 
  
  CA = RES[3, 1][[1]]; 
  n = ncol(CA); 
  R$N = n;
  colnames(CA) = BBCDname(n, "pat"); 
  rownames(CA) = c(BBCDname(cov_num, "covariate"), "assignment"); 
  R$Cov_Assig = CA; 
  
  assig_temp = CA[dim(CA)[1], ]; 
  R$assignments = LETTERS[assig_temp]; 
  
  AS = RES[2, 1][[1]];
  colnames(AS) = BBCDname(strt_num, "strt.");
  rownames(AS) = BBCDname(cov_num, "covariate"); 
  R$'All strata' = AS;
  
  Df = RES[4, 1][[1]]; 
  rownames(Df) = nameString(cov_num, level_num, strt_num, "All", "Real"); 
  R$Diff = t(Df);
  
  R$method = "Atkinson's Optimum Biased Coin Design";
  R$'Data Type' = "Real";
  R$framework = "Model-based approach"; 
  R$data = data; 
  
  class(R) = "carandom";
  
  return(R);
}

###############################################################################
## Covariate-adaptive Biased Coin Design ########
###############################################################################
#AdjBCD = function(data, a = 2.0) UseMethod("AdjBCD")

#AdjBCD.carandom = function(data, a = 2.0) UseMethod("carandom")

AdjBCD = function(data, a = 2.0){
  
  adoub = as.double(a); 
  if(is.na(adoub)){
    stop("a must be a positive number!");
  }else if(a == 0){
    stop("a must be a positive number!");
  }else if(a < 0){
    warning("a is mandated to be positive", call. = FALSE);
  }
  
  if(length(data[is.na(data)]) == 0){
    datap = data;
  }else{
    data[is.na(data)] = "NA"; 
    datap = data; 
  }
  
  rdata = Preprocess(datap); 
  data_proc = rdata$data; 
  cov_num = rdata$cov_num; level_num = rdata$level_num; 
  
  RES = C_RAdjustBCD(data_proc, cov_num, level_num, adoub); 
  
  R = NULL;
  
  covn = colnames(data);
  R$covariates = covn;
  
  strt_num = ncol(RES[2, 1][[1]]);
  R$strt_num = strt_num; 
  
  R$cov_num = cov_num; 
  R$level_num = level_num; 
  
  CA = RES[3, 1][[1]]; 
  n = ncol(CA); 
  R$N = n;
  colnames(CA) = BBCDname(n, "pat"); 
  rownames(CA) = c(BBCDname(cov_num, "covariate"), "assignment"); 
  R$Cov_Assig = CA;
  
  assig_temp = CA[dim(CA)[1], ]; 
  R$assignments = LETTERS[assig_temp]; 
  
  AS = RES[2, 1][[1]]; 
  colnames(AS) = BBCDname(strt_num, "strt.");
  rownames(AS) = BBCDname(cov_num, "covariate"); 
  R$'All strata' = AS;
  
  Df = RES[4, 1][[1]]; 
  rownames(Df) = nameString(cov_num, level_num, strt_num, "All", "Real"); 
  R$Diff = t(Df);
  
  R$method = "Covariate-adaptive Biased Coin Design";
  R$'Data Type' = "Real";
  R$framework = "Stratified randomization"; 
  R$data = data; 
  
  class(R) = "carandom";
  
  return(R);
}

# ###############################################################################
# ## Biased Coin Design with a Bayesian Bias ########
# ###############################################################################
# #BayesBCD = function(data, J = 2) UseMethod("BayesBCD")
# 
# BayesBCD.carandom = function(data, J = 2) UseMethod("carandom")
# 
# BayesBCD = function(data, J = 2){
#   Jint = as.integer(J); 
#   if(is.na(Jint)){
#     stop("J must be a positive integer!")
#   }else if(J %% Jint != 0){
#     warning("J is mandated to be an integer", call. = FALSE)
#   }
#   
#   if(length(data[is.na(data)]) == 0){
#     datap = data;
#   }else{
#     data[is.na(data)] = "NA"; 
#     datap = data; 
#   }
#   
#   rdata = Preprocess(datap); 
#   data_proc = rdata$data; 
#   cov_num = rdata$cov_num; level_num = rdata$level_num; 
#   
#   RES = C_RBayesBCD(data_proc, cov_num, level_num, Jint); 
#   
#   R = NULL;
#   
#   covn = colnames(data);
#   R$covariates = covn;
#   
#   strt_num = ncol(RES[4, 1][[1]]);
#   R$strt_num = strt_num; 
#   
#   R$cov_num = cov_num; 
#   R$level_num = level_num; 
#   
#   CA = RES[6, 1][[1]]; 
#   n = ncol(CA); 
#   R$N = n;
#   colnames(CA) = BBCDname(n, "pat"); 
#   rownames(CA) = c(BBCDname(cov_num, "covariate"), "assignment"); 
#   R$Cov_Assig = CA;
#   
#   AS = RES[4, 1][[1]]; 
#   colnames(AS) = BBCDname(strt_num, "strt.");
#   rownames(AS) = BBCDname(cov_num, "covariate"); 
#   R$'All strata' = AS;
#   
#   Df = RES[7, 1][[1]]; 
#   rownames(Df) = nameString(cov_num, level_num, strt_num, "All", "Real"); 
#   R$Diff = t(Df);
#   
#   R$method = "Covariate-adaptive Biased Coin Design";
#   R$'Data Type' = "Real";
#   R$framework = "Stratified randomization"; 
#   R$data = data; 
#   
#   numJ = RES[5, 1][[1]]; 
#   rownames(numJ) = BBCDname(2, "Treat."); 
#   colnames(numJ) = BBCDname(J, "Category-"); 
#   R$numJ = numJ; 
#   
#   RR = t(RES[2, 1][[1]]); 
#   colnames(RR) = c("no._of_t1", "no._of_t2", "diff"); 
#   R$num_diff = RR; 
#   
#   class(R) = "carandom";
#   
#   return(R);
# }
