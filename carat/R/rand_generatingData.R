###############################################################################
#########################  ##   Simulation   ##  ##############################
################################################################################
## Hu and Hu's general covariate-adaptive randomization(CAR) ########
################################################################################
# HuHuCAR.sim.carandom = function(n = 1000, cov_num = 2, level_num = c(2, 2), 
#                                 pr = rep(0.5, 4), omega = NULL, p = 0.85) UseMethod("carandom")

HuHuCAR.sim = function(n = 1000, cov_num = 2, level_num = c(2, 2), 
                       pr = rep(0.5, 4), omega = NULL, p = 0.85){
  
  if(length(level_num) != cov_num){
    stop("Length of level_num must be equal to cov_num !")
  }
  if(length(which(level_num <= 1.9)) > 0.1){
    stop("number of levels for each covariates must be larger than 2(including 2)!")
  }
  pdoub = as.double(p); 
  if(is.na(pdoub)){
    stop("p must be a positive number!");
  }else if(p <= 0.5){
    stop("set p larger than 0.5 to achieve balance!");
  }else if(p > 1){
    stop("p must be a positive number between 0 and 1!");
  }
  if(length(pr) != sum(level_num)){
    stop("Length of pr should be equal to number of all levels, i.e. sum(level_num)!")
  }
  pmat = Prob_S(cov_num, level_num, pr);
  prcheck = apply(pmat, 2, sum);
  if(length(which(prcheck != 1.0)) > 0){
    stop("probabilities of each margin must sum up to 1 !")
  }
  if(length(omega) != cov_num + 2 && !is.null(omega)){
    stop("Length of omega must be equal to (cov_num + 2)!")
  }else if(is.null(omega)){
    omega = rep(1.0 / (cov_num + 2), times = 2 + cov_num);
  }else{
    omega = abs(omega) / sum(abs(omega));
  }
  
  dat = genData_sim(n, cov_num, level_num, pmat); 
  RES = C_RHPS(dat, cov_num, level_num, omega, pdoub); 
  
  CA = RES[3, 1][[1]]; 
  rownames(CA) = c(BBCDname(cov_num, "covariate"), "assignment");
  colnames(CA) = BBCDname(n, "pat");
  
  R = NULL;
  R$Cov_Assig = CA;
  
  assig_temp = CA[dim(CA)[1], ]; 
  R$assignments = LETTERS[assig_temp]; 
  
  AS = RES[2, 1][[1]]; 
  rownames(AS) = c(BBCDname(cov_num, "covariate"));
  colnames(AS) = c(BBCDname(ncol(AS), "strt."));
  
  R$'All strata' = AS;
  
  Df = RES[4, 1][[1]];
  rownames(Df) = nameString(cov_num, level_num, ncol(AS), "All", "Real");
  
  R$Diff = t(Df);
  
  R$method = "Hu and Hu's General CAR";
  R$cov_num = cov_num;
  R$level_num = level_num;
  R$strt_num = ncol(AS);
  R$N = n;
  R$weight = omega[3 : (2 + cov_num)];
  R$'Data Type' = "Simulated";
  R$framework = "Minimization"; 
  
  class(R) <- "carandom";
  
  return(R);
}

################################################################################
## Pocock and Simon's procedure with two arms ########
################################################################################
# PocSimMIN.sim.carandom = function(n = 1000, cov_num = 2, level_num = c(2, 2), 
#                                   pr = rep(0.5, 4), weight = NULL, 
#                                   p = 0.85) UseMethod("carandom")

PocSimMIN.sim = function(n = 1000, cov_num = 2, level_num = c(2, 2), 
                         pr = rep(0.5, 4), weight = NULL, 
                         p = 0.85){
  
  if(length(level_num) != cov_num){
    stop("Length of level_num must be equal to cov_num !")
  }
  if(length(which(level_num <= 1.9)) > 0.1){
    stop("number of levels for each covariates must be larger than 2(including 2)!")
  }
  pdoub = as.double(p); 
  if(is.na(pdoub)){
    stop("p must be a positive number!");
  }else if(p <= 0.5){
    stop("set p larger than 0.5 to achieve balance!");
  }else if(p > 1){
    stop("p must be a positive number between 0 and 1!");
  }
  if(length(pr) != sum(level_num)){
    stop("Length of pr should be equal to number of all levels, i.e. sum(level_num)!")
  }
  pmat = Prob_S(cov_num, level_num, pr);
  prcheck = apply(pmat, 2, sum);
  if(length(which(prcheck != 1.0)) > 0){
    stop("probabilities of each margin must sum up to 1 !")
  }
  if(length(weight) != cov_num && !is.null(weight)){
    stop("Length of weight must be equal to cov_num!")
  }else if(is.null(weight)){
    omega = c(0, 0, rep(1.0 / cov_num, times = cov_num));
  }else{
    omega = c(0, 0, abs(weight) / sum(abs(weight)));
  }
  
  dat = genData_sim(n, cov_num, level_num, pmat); 
  RES = C_RHPS(dat, cov_num, level_num, omega, pdoub); 
  
  CA = RES[3, 1][[1]]; 
  rownames(CA) = c(BBCDname(cov_num, "covariate"), "assignment");
  colnames(CA) = BBCDname(n, "pat");
  
  R = NULL;
  R$Cov_Assig = CA;
  
  assig_temp = CA[dim(CA)[1], ]; 
  R$assignments = LETTERS[assig_temp]; 
  
  AS = RES[2, 1][[1]]; 
  rownames(AS) = c(BBCDname(cov_num, "covariate"));
  colnames(AS) = c(BBCDname(ncol(AS), "strt."));
  
  R$'All strata' = AS;
  
  Df = RES[4, 1][[1]];
  rownames(Df) = nameString(cov_num, level_num, ncol(AS), "All", "Real");
  
  R$Diff = t(Df);
  
  R$method = "Pocock and Simon's Procedure with Two Arms";
  R$cov_num = cov_num;
  R$level_num = level_num;
  R$strt_num = ncol(AS);
  R$N = n;
  R$weight = omega[3 : (2 + cov_num)];
  R$'Data Type' = "Simulated";
  R$framework = "Minimization"; 
  
  class(R) <- "carandom";
  
  return(R);
}

################################################################################
## Shao's randomization ########
################################################################################
# StrBCD.sim.carandom = function(n = 1000, cov_num = 2, level_num = c(2, 2),
#                                pr = rep(0.5, 4), p = 0.85) UseMethod("carandom")

StrBCD.sim = function(n = 1000, cov_num = 2, level_num = c(2, 2),
                      pr = rep(0.5, 4), p = 0.85){
  
  if(length(level_num) != cov_num){
    stop("Length of level_num must be equal to cov_num !")
  }
  if(length(which(level_num <= 1.9)) > 0.1){
    stop("number of levels for each covariates must be larger than 2(including 2)!")
  }
  pdoub = as.double(p); 
  if(is.na(pdoub)){
    stop("p must be a positive number!");
  }else if(p <= 0.5){
    stop("set p larger than 0.5 to achieve balance!");
  }else if(p > 1){
    stop("p must be a positive number between 0 and 1!");
  }
  if(length(pr) != sum(level_num)){
    stop("Length of pr should be equal to number of all levels, i.e. sum(level_num)!")
  }
  pmat = Prob_S(cov_num, level_num, pr);
  prcheck = apply(pmat, 2, sum);
  if(length(which(prcheck != 1.0)) > 0){
    stop("probabilities of each margin must sum up to 1 !")
  }
  omega = c(0, 1, rep(0, times = cov_num)); 
  
  dat = genData_sim(n, cov_num, level_num, pmat); 
  RES = C_RHPS(dat, cov_num, level_num, omega, pdoub); 
  
  CA = RES[3, 1][[1]]; 
  rownames(CA) = c(BBCDname(cov_num, "covariate"), "assignment");
  colnames(CA) = BBCDname(n, "pat");
  
  R = NULL;
  R$Cov_Assig = CA;
  
  assig_temp = CA[dim(CA)[1], ]; 
  R$assignments = LETTERS[assig_temp]; 
  
  AS = RES[2, 1][[1]]; 
  rownames(AS) = c(BBCDname(cov_num, "covariate"));
  colnames(AS) = c(BBCDname(ncol(AS), "strt."));
  
  R$'All strata' = AS;
  
  Df = RES[4, 1][[1]];
  rownames(Df) = nameString(cov_num, level_num, ncol(AS), "All", "Real");
  
  R$Diff = t(Df);
  
  R$method = "Shao's Procedure";
  R$cov_num = cov_num;
  R$level_num = level_num;
  R$strt_num = ncol(AS);
  R$N = n;
  R$'Data Type' = "Simulated";
  R$framework = "Minimization/Stratified randomization"; 
  
  class(R) <- "carandom";
  
  return(R);
}

################################################################################
## Stratified randomization (STR) with two arms ########
################################################################################
# StrPBR.sim.carandom = function(n = 1000, cov_num = 2, level_num = c(2, 2),
#                                pr = rep(0.5, 4), bsize = 4) UseMethod("carandom")

StrPBR.sim = function(n = 1000, cov_num = 2, level_num = c(2, 2), 
                      pr = rep(0.5, 4), bsize = 4){
  
  if(length(level_num) != cov_num){
    stop("Length of level_num must be equal to cov_num !")
  }
  if(length(which(level_num <= 1.9)) > 0.1){
    stop("number of levels for each covariates must be larger than 2(including 2)!")
  }
  if(length(pr) != sum(level_num)){
    stop("Length of pr should be equal to number of all levels, i.e. sum(level_num)!")
  }
  pmat = Prob_S(cov_num, level_num, pr);
  prcheck = apply(pmat, 2, sum);
  if(length(which(prcheck != 1.0)) > 0){
    stop("probabilities of each margin must sum up to 1 !")
  }
  if(bsize %% 2 != 0){
    stop("block size (bsize) is required to be a multiple of 2!")
  }
  
  dat = genData_sim(n, cov_num, level_num, pmat); 
  RES = C_RStrR(dat, cov_num, level_num, bsize, tr_num = 2); 
  
  CA = RES[3, 1][[1]]; 
  rownames(CA) = c(BBCDname(cov_num, "covariate"), "assignment");
  colnames(CA) = BBCDname(n, "pat");
  
  R = NULL;
  R$Cov_Assig = CA;
  
  assig_temp = CA[dim(CA)[1], ]; 
  R$assignments = LETTERS[assig_temp]; 
  
  AS = RES[2, 1][[1]]; 
  rownames(AS) = c(BBCDname(cov_num, "covariate"));
  colnames(AS) = c(BBCDname(ncol(AS), "strt."));
  
  R$'All strata' = AS;
  
  Df = RES[4, 1][[1]];
  rownames(Df) = nameString(cov_num, level_num, ncol(AS), "All", "Real");
  
  R$Diff = t(Df);
  
  st_num = RES[1, 1][[1]];
  colnames(st_num) = BBCDname(ncol(AS), "level-");
  
  R$`numbers of pats for each strata` = st_num;
  
  R$method = "Stratified Permuted Block Randomization";
  R$cov_num = cov_num;
  R$level_num = level_num;
  R$strt_num = ncol(AS);
  R$N = n;
  R$bsize = bsize;
  R$`Data Type` = "Simulated";
  R$framework = "Stratified randomization"; 
  
  class(R) <- "carandom";
  
  return(R);
}

################################################################################
## Atkinson's Optimum Biased Coin Design with two arms ########
################################################################################
# DoptBCD.sim.carandom = function(n = 1000, cov_num = 2, level_num = c(2, 2),
#                                 pr = rep(0.5, 4)) UseMethod("carandom")

DoptBCD.sim = function(n = 1000, cov_num = 2, level_num = c(2, 2), 
                       pr = rep(0.5, 4)){
  
  if(length(level_num) != cov_num){
    stop("Length of level_num must be equal to cov_num !")
  }
  if(length(which(level_num <= 1.9)) > 0.1){
    stop("number of levels for each covariates must be larger than 2(including 2)!")
  }
  if(length(pr) != sum(level_num)){
    stop("Length of pr should be equal to number of all levels, i.e. sum(level_num)!")
  }
  pmat = Prob_S(cov_num, level_num, pr);
  prcheck = apply(pmat, 2, sum);
  if(length(which(prcheck != 1.0)) > 0){
    stop("probabilities of each margin must sum up to 1 !")
  }
  
  dat = genData_sim(n, cov_num, level_num, pmat); 
  RES = C_RAtkinBCD(dat, cov_num, level_num); 
  
  CA = RES[3, 1][[1]]; 
  rownames(CA) = c(BBCDname(cov_num, "covariate"), "assignment");
  colnames(CA) = BBCDname(n, "pat");
  
  R = NULL;
  R$Cov_Assig = CA; 
  
  assig_temp = CA[dim(CA)[1], ]; 
  R$assignments = LETTERS[assig_temp]; 
  
  AS = RES[2, 1][[1]]; 
  rownames(AS) = c(BBCDname(cov_num, "covariate"));
  colnames(AS) = c(BBCDname(ncol(AS), "strt."));
  
  R$'All strata' = AS;
  
  Df = RES[4, 1][[1]];
  rownames(Df) = nameString(cov_num, level_num, ncol(AS), "All", "Real"); 
  
  R$Diff = t(Df);
  
  R$method = "Atkinson's Optimum Biased Coin Design with Two Arms";
  R$cov_num = cov_num;
  R$level_num = level_num;
  R$strt_num = ncol(AS);
  R$N = n;
  R$'Data Type' = "Simulated";
  R$framework = "Model-based approach"; 
  
  class(R) <- "carandom";
  
  return(R);
}

# ################################################################################
# ## Biased Coin Design with a Bayesian Bias with two arms ########
# ################################################################################
# BayesBCD.sim.carandom = function(n = 1000, cov_num = 2, level_num = c(2, 2),
#                                  pr = rep(0.5, 4), J = 2) UseMethod("carandom")
# 
# BayesBCD.sim = function(n = 1000, cov_num = 2, level_num = c(2, 2), 
#                         pr = rep(0.5, 4), J = 2){
#   if(length(level_num) != cov_num){
#     stop("Length of level_num must be equal to cov_num !")
#   }
#   if(length(which(level_num <= 1.9)) > 0.1){
#     stop("number of levels for each covariates must be larger than 2(including 2)!")
#   }
#   if(length(pr) != sum(level_num)){
#     stop("Length of pr should be equal to number of all levels, i.e. sum(level_num)!")
#   }
#   pmat = Prob_S(cov_num, level_num, pr);
#   prcheck = apply(pmat, 2, sum);
#   if(length(which(prcheck != 1.0)) > 0){
#     stop("probabilities of each margin must sum up to 1 !")
#   }
#   Jint = as.integer(J); 
#   if(is.na(Jint)){
#     stop("J must be a positive integer!")
#   }else if(J %% Jint != 0){
#     warning("J is mandated to be an integer", call. = FALSE)
#   }
#   
#   dat = genData_sim(n, cov_num, level_num, pmat); 
#   RES = C_RBayesBCD(dat, cov_num, level_num, Jint); 
#   
#   CA = RES[6, 1][[1]]; 
#   rownames(CA) = c(BBCDname(cov_num, "covariate"), "assignment");
#   colnames(CA) = BBCDname(n, "pat");
#   
#   R = NULL;
#   R$Cov_Assig = CA;
#   
#   AS = RES[4, 1][[1]]; 
#   rownames(AS) = c(BBCDname(cov_num, "covariate"));
#   colnames(AS) = c(BBCDname(ncol(AS), "strt."));
#   
#   R$'All strata' = AS;
#   
#   Df = RES[7, 1][[1]];
#   rownames(Df) = nameString(cov_num, level_num, ncol(AS), "All", "Real"); 
#   R$Diff = t(Df);
#   
#   R$method = "Biased Coin Design with a Bayesian Bias with Two Arms";
#   R$cov_num = cov_num;
#   R$level_num = level_num;
#   R$strt_num = ncol(AS);
#   R$N = n;
#   R$'Data Type' = "Simulated";
#   R$framework = "Model-based approach"; 
#   
#   RR = t(RES[2, 1][[1]]); 
#   colnames(RR) = c("no._of_t1", "no._of_t2", "diff"); 
#   R$num_diff = RR; 
#   
#   numJ = RES[5, 1][[1]]; 
#   rownames(numJ) = BBCDname(2, "Treat."); 
#   colnames(numJ) = BBCDname(J, "Category"); 
#   R$numJ = numJ; 
#   
#   class(R) <- "carandom";
#   
#   return(R);
# }


################################################################################
## Covariate-adaptive Biased Coin Design with two arms ########
################################################################################
# AdjBCD.sim.carandom = function(n = 1000, cov_num = 2, level_num = c(2, 2),
#                                pr = rep(0.5, 4), a = 2.0) UseMethod("carandom")

AdjBCD.sim = function(n = 1000, cov_num = 2, level_num = c(2, 2), 
                      pr = rep(0.5, 4), a = 2.0){
  if(length(level_num) != cov_num){
    stop("Length of level_num must be equal to cov_num !")
  }
  if(length(which(level_num <= 1.9)) > 0.1){
    stop("number of levels for each covariates must be larger than 2(including 2)!")
  }
  if(length(pr) != sum(level_num)){
    stop("Length of pr should be equal to number of all levels, i.e. sum(level_num)!")
  }
  pmat = Prob_S(cov_num, level_num, pr);
  prcheck = apply(pmat, 2, sum);
  if(length(which(prcheck != 1.0)) > 0){
    stop("probabilities of each margin must sum up to 1 !")
  }
  adoub = as.double(a); 
  if(is.na(adoub)){
    stop("a must be a positive number!");
  }else if(a == 0){
    stop("a must be a positive number!");
  }else if(a < 0){
    warning("a is mandated to be positive", call. = FALSE);
  }
  
  dat = genData_sim(n, cov_num, level_num, pmat); 
  RES = C_RAdjustBCD(dat, cov_num, level_num, adoub); 
  
  CA = RES[3, 1][[1]]; 
  rownames(CA) = c(BBCDname(cov_num, "covariate"), "assignment");
  colnames(CA) = BBCDname(n, "pat");
  
  R = NULL;
  R$Cov_Assig = CA;
  
  assig_temp = CA[dim(CA)[1], ]; 
  R$assignments = LETTERS[assig_temp]; 
  
  AS = RES[2, 1][[1]]; 
  rownames(AS) = c(BBCDname(cov_num, "covariate"));
  colnames(AS) = c(BBCDname(ncol(AS), "strt."));
  
  R$'All strata' = AS;
  
  Df = RES[4, 1][[1]];
  rownames(Df) = nameString(cov_num, level_num, ncol(AS), "All", "Real");
  R$Diff = t(Df);
  
  R$method = "Covariate-adaptive Biased Coin Design with Two Arms";
  R$cov_num = cov_num;
  R$level_num = level_num;
  R$strt_num = ncol(AS);
  R$N = n;
  R$'Data Type' = "Simulated";
  R$framework = "Stratified randomization"; 
  
  class(R) <- "carandom";
  
  return(R);
}
