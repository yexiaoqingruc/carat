
##  Foundation functions  ##
##############################################################################################################
## paths generation
##############################################################################################################
pathout = function(path, folder = "HHCAR", method = "HuHuCAR"){
  pathvec = character(); 
  pathf = paste(path, folder, sep = "/");
  pathff = paste(pathf, "CovLevInd", sep = "/"); 
  pathg = paste(pathf, "static", sep = "/");
  pathgg = paste(pathf, "dynamic", sep = "/");
  
  path1 = paste(pathff, "Reference.RData", sep = "/"); 
  path2 = paste(pathg, "cov_num.RData", sep = "/");
  path3 = paste(pathg, "level_num.RData", sep = "/");
  path4 = paste(pathff, "CovIndex.RData", sep = "/");
  
  path5 = paste(pathgg, "D.RData", sep = "/");
  
  ### paths to create folders
  pathvec[1] = pathf; 
  pathvec[2] = pathff; 
  pathvec[3] = pathg; 
  pathvec[4] = pathgg; 
  
  ### paths to save R-files
  pathvec[5] = path1; 
  pathvec[6] = path2;
  pathvec[7] = path3;
  pathvec[8] = path4;
  pathvec[9] = path5;
  
  if(method == "HuHuCAR" || method == "PocSimMIN"){
    path7 = paste(pathg, "omega.RData", sep = "/"); 
    path8 = paste(pathg, "p.RData", sep = "/");
    pathvec[10] = path7;
    pathvec[11] = path8;
    path6 = paste(pathgg, "strp.RData", sep = "/"); 
    pathvec[12] = path6; 
    pathvec[13] = paste(pathg, "pathlog", sep = "/")
    pathvec[14] = paste(pathvec[13], "pathvec.RData", sep = "/")
  }
  if(method == "StrBCD"){
    path8 = paste(pathg, "p.RData", sep = "/");
    pathvec[10] = path8;
    path6 = paste(pathgg, "strp.RData", sep = "/"); 
    pathvec[11] = path6; 
    pathvec[12] = paste(pathg, "pathlog", sep = "/")
    pathvec[13] = paste(pathvec[12], "pathvec.RData", sep = "/")
  }
  if(method == "StrPBR"){
    path9 = paste(pathgg, "BG.RData", sep = "/"); 
    path10 = paste(pathg, "B.RData", sep = "/");
    path11 = paste(pathgg, "strp.RData", sep = "/"); 
    path12 = paste(pathg, "bsize.RData", sep = "/");
    pathvec[10] = path9; 
    pathvec[11] = path10; 
    pathvec[12] = path11
    pathvec[13] = path12; 
    pathvec[14] = paste(pathg, "pathlog", sep = "/")
    pathvec[15] = paste(pathvec[14], "pathvec.RData", sep = "/")
  }
  if(method == "DoptBCD"){
    patha1 = paste(pathgg, "Fmatrix.RData", sep = "/");
    patha2 = paste(pathgg, "No.RData", sep = "/"); 
    patha3 = paste(pathgg, "b.RData", sep = "/"); 
    pathvec[10 : 12] = c(patha1, patha2, patha3); 
    path6 = paste(pathgg, "strp.RData", sep = "/"); 
    pathvec[13] = path6; 
    pathvec[14] = paste(pathg, "pathlog", sep = "/")
    pathvec[15] = paste(pathvec[14], "pathvec.RData", sep = "/")
  }
  if(method == "AdjBCD"){
    pathaj = paste(pathg, "a.RData", sep = "/"); 
    pathvec[10] = pathaj; 
    path6 = paste(pathgg, "strp.RData", sep = "/"); 
    pathvec[11] = path6; 
    pathvec[12] = paste(pathg, "pathlog", sep = "/")
    pathvec[13] = paste(pathvec[12], "pathvec.RData", sep = "/")
  }
  return(pathvec);
}

##############################################################################################################
## Function used to enter involved covariates
##############################################################################################################
# A simple version of function for Entering all covariates
##############################################################################################################
EnterSigCov = function(pathvec){
  covariate = vector();
  cov = "Begin"; 
  message("Please enter the involved covariates: "); 
  while (cov != "") {
    message("  Please enter a new covariate: "); 
    message("    Notice: If no more covariates to be entered, please PRESS Enter directly"); 
    cov = readline(prompt = "New Covariate: "); 
    while(cov %in% covariate){
      message("INVALID INPUT !", "\n", 
              "Please DO NOT enter covariate repeatedly!", "\n", 
              "Notice: If no more covariates to be entered, please PRESS Enter directly"); 
      cov = readline(prompt = "Reenter a new covariate: "); 
    }
    covariate = c(covariate, cov);
  }
  cov_num = length(covariate) - 1;
  covariate = covariate[1 : cov_num];
  #covr = as.numeric(as.factor(covariate));
  covr = 1 : cov_num;
  CovIndex = data.frame(covariate, covr, stringsAsFactors = TRUE);
  colnames(CovIndex) = c("Real Names", "factors"); 
  rownames(CovIndex) = BBCDname(cov_num, "covariate");
  message("According to your input, all covariates are stamped to be");
  message("\n");
  message(paste("\t", covariate, "--", covr, "\n", sep = " "))
  return(list("cov_num" = cov_num, "CovIndex" = CovIndex, 
              "covariate" = covariate, "covr" = covr));
}

##############################################################################################################
# A comprehensive version of function for Entering all covariates
##############################################################################################################
EnterCov = function(pathvec){
  RSigC = EnterSigCov(pathvec); 
  message("\n"); 
  message("Continue or not?", "\n", 
          " 'n' -- stop running", "  ", "input 'y' or PRESS Enter -- reenter or save\n"); 
  w1 = readline("Enter y or n: "); 
  if(w1 == "n" | w1 == "No" | w1 == "no" | w1 == "exit" | w1 == "Exit"){
    return(NULL);
  }else{
    message("Reenter involved covariates?"); 
    w2 = readline("Enter y or n: "); 
    while(w2 == "y" || w2 == "Y" || w2 == "yes" || w2 == "YES" || w2 == "Yes"){
      RSigC = EnterSigCov(pathvec); 
      message("\n")
      message("Reenter involved covariates?"); 
      w2 = readline("Enter y or n: "); 
    }
    return(list("cov_num" = RSigC$cov_num, "covariate" = RSigC$covariate, 
                "covr" = RSigC$covr, "CovIndex" = RSigC$CovIndex));
  }
}

##############################################################################################################
## Functions used to Enter levels for each covariztes
##############################################################################################################
# Function of entering a levels for a single covariate
##############################################################################################################
EnterSigLev = function(i, covariate, covr, pth, level_num){
  level = vector();
  lel = "Begin";
  while(lel != ""){
    message("  Enter the new LEVEL for covariate -- ", covariate[i], ": ");
    message("    Notice: If no more level to be entered for ", covariate[i],  ", please PRESS Enter directly"); 
    lel = readline(prompt = "New level: "); 
    while(lel %in% level){
      message("INVALID INPUT !", "\n", 
              "Please DO NOT enter level for ---", covariate[i], "repeatedly.", "\n", 
              "Notice: If no more level to be entered for", covariate[i],  ", please PRESS Enter directly", "\n");
      lel = readline(prompt = "New level: ");
    }
    level = c(level, lel);
  }
  level_num[covr[i]] = length(level) - 1;
  level = level[1 : level_num[covr[i]]];
  #lev = as.numeric(as.factor(level));
  lev = 1 : level_num[covr[i]]; 
  Ccha = rep(covariate[i], times = level_num[covr[i]]); 
  Cfac = rep(covr[i], times = level_num[covr[i]]); 
  LevIn = data.frame(Ccha, Cfac, level, lev, stringsAsFactors = TRUE); 
  save(LevIn, file = pth[i]); 
  message("According to you input, levels for covariate -- ", covariate[i], "\n\n", 
          covr[i], "--", covariate[i], sep = " ", "\n");
  for(k in 1 : level_num[covr[i]]){
    message("\t", level[k], "--", lev[k], sep = " ");
  }
  return(list("LevIn" = LevIn, "level_num" = level_num, 
              "lev" = lev, "level" = level));
}

##############################################################################################################
# A comprehensive function of entering levels for all covariate
##############################################################################################################
EnterLev = function(pathvec, cov_num, covariate, covr, pth){
  level_num = vector();
  Reference = data.frame(stringsAsFactors = TRUE);
  message("\n");
  message("Please enter LEVELs for each covariate: "); 
  for(i in 1 : cov_num){
    RSigL = EnterSigLev(i, covariate, covr, pth, level_num = level_num); 
    message("\n");
    message("Reenter LEVELs for covariate -- ", covariate[i], " or not?"); 
    w2 = readline(prompt = "Enter y or n: "); 
    while(w2 == "y" || w2 == "Y" || w2 == "yes" || w2 == "YES" || w2 == "Yes"){
      RSigL = EnterSigLev(i, covariate, covr, pth, level_num = level_num); 
      message("\n"); 
      message("Reenter LEVELs for covariate -- ", covariate[i], " or not?"); 
      w2 = readline(prompt = "Enter y or n: "); 
    }
    LevIn = RSigL$LevIn; level_num = RSigL$level_num;
    Reference = rbind(Reference, LevIn);
  }
  colnames(Reference) = c("covariate", "factor_cov", "level", "factor_lev");
  save(Reference, file = pathvec[5]);
  
  return(list("level_num" = level_num));
}

##############################################################################################################
## Function used to Enter omega for "HuHuCAR" and weight for "PocSimMIN"
##############################################################################################################
EnterWeig = function(cov_num, covariate, covr, method = "HuHuCAR"){
  warnings("off");
  aspect = reaspect = character();
  atp = "Enter the weight for the"; btp = "aspect: ";
  reatp = "Reenter the weight for the";
  aspect[1 : 2] = c(paste(atp, "OVERALL", btp, sep = " "),  paste(atp, "WITHIN-STRATUM", btp, sep = " ")); 
  reaspect[1 : 2] = c(paste(reatp, "OVERALL", btp, sep = " "),  paste(reatp, "WITHIN-STRATUM", btp, sep = " "));
  for(k in  1 : cov_num){
    aspect[2 + k] = paste(atp, paste("MARGIN", "--", covariate[which(covr == k)], ": "), sep = " ");
    reaspect[2 + k] = paste(reatp, paste("MARGIN", "--", covariate[which(covr == k)], ": "), sep = " ");
  }
  omega = vector();
  if(method == "HuHuCAR"){
    message("Please allocate WEIGHTs to each aspects: ", "\n", 
            "  Notice: larger the absolute value you enter, stronger tendency to obtain balance on the corresponding aspect.")
  }else if(method == "PocSimMIN"){
    message("Please allocate WEIGHTs to each margins: ", "\n", 
            "  Notice: larger the absolute value you enter, stronger tendency to obtain balance on corresponding margin.")
  }
  
  for(i in 1 : (2 + cov_num)){
    if(method == "PocSimMIN" && (i %in% 1 : 2)){
      omega[i] = 0; 
    }else{
      omov = readline(prompt = aspect[i]); 
      omov = as.integer(omov); 
      warnings("off"); 
      while(is.na(omov)){
        message("INVALID INPUT !", "\n", "Please enter a nonzero number: "); 
        omov = readline(prompt = reaspect[i]); 
        omov = as.integer(omov); 
        warnings("off"); 
      }
      omega[i] = omov; 
    }
  }
  omega = abs(omega) / sum(abs(omega)); 
  
  message("Weights for each aspects are: "); 
  message("\n")
  message("\t", "OVERALL", "--", omega[1], sep = " "); 
  message("\t", "WITHIN-STRT.", "--", omega[2], sep = " ", "\n"); 
  for(kc in 1 : (cov_num)){
    message("\t", covariate[which(covr == kc)], "--", omega[2 + kc], sep = " ");
  }
  return(omega); 
}

##############################################################################################################
## Function used to Enter new covariate profile
##############################################################################################################
EnterCovPrf = function(cov_num, pth, CovIndex, covariate, R){
  cov_profile = lecomvec = vector();
  message(" ", "Please enter COVARIATE PROFILE of the coming patients: ");
  for(l in 1 : cov_num){
    LevIn = 0; 
    load(pth[l]); 
    R[[l]] = LevIn;
    message("Please enter the level of covariate", "---", covariate[l], ": "); 
    lecoming = readline(prompt = "Enter the level: "); 
    while(!lecoming %in% LevIn$level){
      message("INVALID INPUT!", "\n", "Notice: Enter one of \n");
      print(t(LevIn[3]));
      lecoming = readline(prompt = "Reenter the level: "); 
    }
    lecomvec[l] = lecoming;
    cov_profile[CovIndex[l, 2]] = LevIn[which(LevIn[, 3] == lecoming, arr.ind = T), 4]; 
    # equal to ## cov_profile[l] = match(lecoming, LevIn$level); 
  }
  message("COVARIATE PROFILE of the coming patient is: \n");
  message("\n"); 
  for(ll in 1 : cov_num){
    message("\t", covariate[ll], " -- ", lecomvec[ll], "\n");
  }
  R$cov_profile = cov_profile; 
  return(R); 
}


##  Main Functions ##
##############################################################################################################
## Hu and Hu's general covariate-adaptive randomization
##############################################################################################################

#HuHuCAR.ui.carseq = function(path, folder = "HuHuCAR") UseMethod("carseq")

HuHuCAR.ui = function(path, folder = "HuHuCAR"){
  
  message("Is this the first patient? ");
  first = readline(prompt = "Enter T or F: ");
  if(first == "T" || first == "True" || first == "TRUE"){
    pathvec = pathout(path, folder, method = "HuHuCAR"); 
    while(dir.exists(pathvec[1])){
      folder = paste(folder, "tymh", sep = "_"); 
      warning(paste("Homonym exists! A file named '",  
                    folder, 
                    "' is created. ", sep = "")); 
      pathvec = pathout(path, folder, method = "HuHuCAR");
    }
    for(i in c(1 : 4, length(pathvec) - 1)){
      dir.create(pathvec[i]); 
    }
    warnings("off");
    save(pathvec, file = pathvec[length(pathvec)])
    
    Rcov = EnterCov(pathvec);
    if(is.null(Rcov)){
      message("Forced cessation of execution.\n")
      unlink(pathvec[1], recursive = TRUE); 
      return(NULL);
    }else{
      cov_num = Rcov$cov_num; 
      save(cov_num, file = pathvec[6])
      covariate = Rcov$covariate; 
      covr = Rcov$covr; 
      CovIndex = Rcov$CovIndex; 
      save(CovIndex, file = pathvec[8]); 
      
      pth = character();
      for(i in 1 : cov_num){
        a = paste("LevIndex", i, sep = "");
        b = paste(a, "RData", sep = ".");
        pth[i] = paste(pathvec[2], b, sep = "/");
      }
      
      Rlev = EnterLev(pathvec, cov_num, covariate, covr, pth);
      
      if(is.null(Rlev)){
        message("Forced cessation of execution.\n")
        unlink(pathvec[1], recursive = TRUE); 
        return(NULL); 
      }else{
        level_num = Rlev$level_num;
        save(level_num, file = pathvec[7]);
        
        D = matrix(0, ncol = 1 + prod(level_num) + sum(level_num), nrow = 1);
        strp = rep(0, times = prod(level_num)); 
        omega = EnterWeig(cov_num, covariate, covr); 
        message("\n"); 
        message("Reenter weights or not?"); 
        ww = readline(prompt = "Enter y or n: ");
        while (ww == "y" || ww == "Y" || ww == "yes" || ww == "YES" || ww == "Yes"){
          omega = EnterWeig(cov_num, covariate, covr); 
          message("\n"); 
          message("Reenter weights or not?"); 
          ww = readline(prompt = "Enter y or n: ");
        }
        save(omega, file = pathvec[10]); 
        
        message("Please enter the biased coin probability (0-1): "); 
        p = readline(prompt = "Enter the probability: ");
        p = as.double(p); 
        warnings("off"); 
        while (is.na(p)){
          message("INVALID INPUT !", "\n", 
                  "Please enter a positive number between 0 and 1: \n"); 
          p = readline(prompt = "Reenter the probability: "); 
          p = as.double(p); 
          warnings("off"); 
        }
        save(p, file = pathvec[11]); 
      }
    }
  }else{
    
    pathverify = paste(path, folder, sep = "/"); 
    while(dir.exists(pathverify)){
      folder = paste(folder, "tymh", sep = "_")
      pathverify = paste(path, folder, sep = "/"); 
    }
    folderpathlog = sub('.....$', '', folder)
    
    load(paste(path, folderpathlog, 
               "static", "pathlog", "pathvec.RData", 
               sep = "/")); 
    
    for(m in 6 : 12){
      load(pathvec[m]);
    }
    pth = character(); 
    for(i in 1 : cov_num){
      a = paste("LevIndex", i, sep = "");
      b = paste(a, "RData", sep = ".");
      pth[i] = paste(pathvec[2], b, sep = "/");
    }
    covariate = as.character(CovIndex[, 1]);
    covr = CovIndex[, 2];
  }
  R = NULL; 
  CPR = EnterCovPrf(cov_num, pth, CovIndex, covariate, R); 
  message("\n"); 
  message("Reenter COVARIATE PROFILE or not?"); 
  wc = readline("Enter y or n: "); 
  while(wc == "y" || wc == "Y" || wc == "yes" || wc == "YES" || wc == "Yes"){
    CPR = EnterCovPrf(cov_num, pth, CovIndex, covariate, R); 
    message("\n"); 
    message("Reenter COVARIATE PROFILE or not?"); 
    wc = readline("Enter y or n: "); 
  }
  
  PS = PStrGen(cov_num, level_num); 
  RES = HPSOne(t(D), PS, CPR$cov_profile, cov_num, 
               level_num, omega, strp, p);
  
  strp = RES[1, 1][[1]]; 
  rownames(strp) = BBCDname(prod(level_num), "level-")
  save(strp, file = pathvec[12]); 
  
  CPR$covariate = covariate; CPR$covr = covr; CPR$cov_num = cov_num; 
  CPR$method = "Hu and Hu's General CAR";
  
  D = t(RES[4, 1][[1]]); 
  colnames(D) = nameString(cov_num, level_num, prod(level_num), "All", PS); 
  save(D, file = pathvec[9]);
  
  ass = RES[3, 1][[1]]
  CPR$assignment = ass[1, 1]; 
  
  class(CPR) = "carseq";
  return(CPR);
}

##############################################################################################################
## Pocock and Simon's procedure
##############################################################################################################
#PocSimMIN.ui.carseq = function(path, folder = "PocSimMIN") UseMethod("carseq")

PocSimMIN.ui = function(path, folder = "PocSimMIN"){
  
  message("Is this the first patient? ");
  first = readline(prompt = "Enter T or F: ");
  if(first == "T" || first == "True" || first == "TRUE"){
    pathvec = pathout(path, folder, method = "PocSimMIN"); 
    while(dir.exists(pathvec[1])){
      folder = paste(folder, "tymh", sep = "_"); 
      warning(paste("Homonym exists! A file named '",  
                    folder, 
                    "' is created. ", sep = "")); 
      pathvec = pathout(path, folder, method = "PocSimMIN");
    }
    for(i in c(1 : 4, length(pathvec) - 1)){
      dir.create(pathvec[i]); 
    }
    warnings("off");
    save(pathvec, file = pathvec[length(pathvec)]); 
    
    Rcov = EnterCov(pathvec);
    if(is.null(Rcov)){
      message("Forced cessation of execution.\n")
      unlink(pathvec[1], recursive = TRUE); 
      return(NULL);
    }else{
      cov_num = Rcov$cov_num; 
      save(cov_num, file = pathvec[6])
      covariate = Rcov$covariate; 
      covr = Rcov$covr; 
      CovIndex = Rcov$CovIndex; 
      save(CovIndex, file = pathvec[8]); 
      
      pth = character();
      for(i in 1 : cov_num){
        a = paste("LevIndex", i, sep = "");
        b = paste(a, "RData", sep = ".");
        pth[i] = paste(pathvec[2], b, sep = "/");
      }
      
      Rlev = EnterLev(pathvec, cov_num, covariate, covr, pth);
      
      if(is.null(Rlev)){
        message("Forced cessation of execution.\n")
        unlink(pathvec[1], recursive = TRUE); 
        return(NULL); 
      }else{
        level_num = Rlev$level_num;
        save(level_num, file = pathvec[7]);
        
        D = matrix(0, ncol = 1 + prod(level_num) + sum(level_num), nrow = 1);
        strp = matrix(0, nrow = prod(level_num), ncol = 1); 
        
        omega = EnterWeig(cov_num, covariate, covr, method = "PocSimMIN"); 
        message("\n"); 
        message("Reenter weights or not?"); 
        ww = readline(prompt = "Enter y or n: ");
        while (ww == "y" || ww == "Y" || ww == "yes" || ww == "YES" || ww == "Yes"){
          omega = EnterWeig(cov_num, covariate, covr, method = "PocSimMIN"); 
          message("\n"); 
          message("Reenter weights or not?"); 
          ww = readline(prompt = "Enter y or n: ");
        }
        save(omega, file = pathvec[10]); 
        
        message("Please enter the biased coin probability (0-1):"); 
        p = readline(prompt = "Enter the probability: ");
        p = as.double(p); 
        warnings("off"); 
        while (is.na(p)){
          message("INVALID INPUT !", "\n", 
                  "Please enter a positive number between 0 and 1: \n"); 
          p = readline(prompt = "Reenter the probability: "); 
          p = as.double(p); 
          warnings("off"); 
        }
        save(p, file = pathvec[11]); 
      }
    }
  }else{
    
    pathverify = paste(path, folder, sep = "/"); 
    while(dir.exists(pathverify)){
      folder = paste(folder, "tymh", sep = "_")
      pathverify = paste(path, folder, sep = "/"); 
    }
    folderpathlog = sub('.....$', '', folder)
    
    load(paste(path, folderpathlog, 
               "static", "pathlog", "pathvec.RData", 
               sep = "/")); 
    
    for(m in 6 : 12){
      load(pathvec[m]);
    }
    pth = character(); 
    for(i in 1 : cov_num){
      a = paste("LevIndex", i, sep = "");
      b = paste(a, "RData", sep = ".");
      pth[i] = paste(pathvec[2], b, sep = "/");
    }
    covariate = as.character(CovIndex[, 1]);
    covr = CovIndex[, 2];
  }
  R = list(); 
  CPR = EnterCovPrf(cov_num, pth, CovIndex, covariate, R); 
  message("\n"); 
  message("Reenter COVARIATE PROFILE or not?"); 
  wc = readline("Enter y or n: "); 
  while(wc == "y" || wc == "Y" || wc == "yes" || wc == "YES" || wc == "Yes"){
    CPR = EnterCovPrf(cov_num, pth, CovIndex, covariate, R); 
    message("\n"); 
    message("Reenter COVARIATE PROFILE or not?"); 
    wc = readline("Enter y or n: "); 
  }
  
  PS = PStrGen(cov_num, level_num); 
  RES = HPSOne(t(D), PS, CPR$cov_profile, cov_num, 
               level_num, omega, strp, p);
  
  CPR$covariate = covariate; CPR$covr = covr; CPR$cov_num = cov_num; 
  CPR$method = "Pocock and Simon's Procedure with Two Arms";
  
  D = t(RES[4, 1][[1]]); 
  colnames(D) = nameString(cov_num, level_num, prod(level_num), "All", PS); 
  save(D, file = pathvec[9]);
  
  strp = RES[1, 1][[1]]; 
  rownames(strp) = BBCDname(prod(level_num), "level-"); 
  save(strp, file = pathvec[12]); 
  
  ass = RES[3, 1][[1]]; 
  CPR$assignment = ass[1, 1]; 
  
  class(CPR) = "carseq";
  return(CPR);
}

##############################################################################################################
## Shao's method
##############################################################################################################
#StrBCD.ui.carseq = function(path, folder = "StrBCD") UseMethod("carseq")

StrBCD.ui = function(path, folder = "StrBCD"){
  
  message("Is this the first patient? ");
  first = readline(prompt = "Enter T or F: ");
  if(first == "T" || first == "True" || first == "TRUE"){
    pathvec = pathout(path, folder, method = "StrBCD"); 
    while(dir.exists(pathvec[1])){
      folder = paste(folder, "tymh", sep = "_"); 
      warning(paste("Homonym exists! A file named '",  
                    folder, 
                    "' is created. ", sep = "")); 
      pathvec = pathout(path, folder, method = "StrBCD");
    }
    for(i in c(1 : 4, length(pathvec) - 1)){
      dir.create(pathvec[i]); 
    }
    warnings("off");
    save(pathvec, file = pathvec[length(pathvec)]); 
    
    Rcov = EnterCov(pathvec);
    if(is.null(Rcov)){
      message("Forced cessation of execution.\n")
      unlink(pathvec[1], recursive = TRUE); 
      return(NULL);
    }else{
      cov_num = Rcov$cov_num; 
      save(cov_num, file = pathvec[6])
      covariate = Rcov$covariate; 
      covr = Rcov$covr; 
      CovIndex = Rcov$CovIndex; 
      save(CovIndex, file = pathvec[8]); 
      
      pth = character();
      for(i in 1 : cov_num){
        a = paste("LevIndex", i, sep = "");
        b = paste(a, "RData", sep = ".");
        pth[i] = paste(pathvec[2], b, sep = "/");
      }
      
      Rlev = EnterLev(pathvec, cov_num, covariate, covr, pth);
      
      if(is.null(Rlev)){
        message("Forced cessation of execution.\n")
        unlink(pathvec[1], recursive = TRUE); 
        return(NULL); 
      }else{
        level_num = Rlev$level_num;
        save(level_num, file = pathvec[7]);
        
        D = matrix(0, ncol = 1 + prod(level_num) + sum(level_num), nrow = 1);
        strp = matrix(0, nrow = prod(level_num), ncol = 1); 
        
        message("Please enter the biased coin probability (0-1): "); 
        p = readline(prompt = "Enter the probability: ");
        p = as.double(p); 
        warnings("off"); 
        while (is.na(p)){
          message("INVALID INPUT !", "\n", 
                  "Please enter a positive number between 0 and 1: \n"); 
          p = readline(prompt = "Reenter the probability: "); 
          p = as.double(p); 
          warnings("off"); 
        }
        save(p, file = pathvec[10]); 
      }
    }
  }else{
    
    pathverify = paste(path, folder, sep = "/"); 
    while(dir.exists(pathverify)){
      folder = paste(folder, "tymh", sep = "_")
      pathverify = paste(path, folder, sep = "/"); 
    }
    folderpathlog = sub('.....$', '', folder)
    
    load(paste(path, folderpathlog, 
               "static", "pathlog", "pathvec.RData", 
               sep = "/")); 
    
    for(m in 6 : 11){
      load(pathvec[m]);
    }
    pth = character(); 
    for(i in 1 : cov_num){
      a = paste("LevIndex", i, sep = "");
      b = paste(a, "RData", sep = ".");
      pth[i] = paste(pathvec[2], b, sep = "/");
    }
    covariate = as.character(CovIndex[, 1]);
    covr = CovIndex[, 2]; 
  }
  R = list(); 
  CPR = EnterCovPrf(cov_num, pth, CovIndex, covariate, R); 
  message("\n"); 
  message("Reenter COVARIATE PROFILE or not?"); 
  wc = readline("Enter y or n: "); 
  while(wc == "y" || wc == "Y" || wc == "yes" || wc == "YES" || wc == "Yes"){
    CPR = EnterCovPrf(cov_num, pth, CovIndex, covariate, R); 
    message("\n"); 
    message("Reenter COVARIATE PROFILE or not?"); 
    wc = readline("Enter y or n: "); 
  }
  
  PS = PStrGen(cov_num, level_num); 
  RES = HPSOne(t(D), PS, CPR$cov_profile, cov_num, 
               level_num, omega = c(0, 1, rep(0, times = cov_num)), strp, p);
  
  CPR$covariate = covariate; CPR$covr = covr; CPR$cov_num = cov_num; 
  CPR$method = "Shao's Procedure";
  
  D = t(RES[4, 1][[1]]); 
  colnames(D) = nameString(cov_num, level_num, prod(level_num), "All", PS); 
  save(D, file = pathvec[9]);
  
  strp = RES[1, 1][[1]]; 
  rownames(strp) = BBCDname(prod(level_num), "level-"); 
  save(strp, file = pathvec[11]); 
  
  ass = RES[3, 1][[1]]; 
  CPR$assignment = ass[1, 1]; 
  
  class(CPR) = "carseq";
  return(CPR);
}

##############################################################################################################
## Stratified permuted block randomization
##############################################################################################################
#StrPBR.ui.carseq = function(path, folder = "StrPBR") UseMethod("carseq")

StrPBR.ui = function(path, folder = "StrPBR"){
  
  message("Is this the first patient? ");
  first = readline(prompt = "Enter T or F: ");
  if(first == "T" || first == "True" || first == "TRUE"){
    pathvec = pathout(path, folder, method = "StrPBR"); 
    while(dir.exists(pathvec[1])){
      folder = paste(folder, "tymh", sep = "_"); 
      warning(paste("Homonym exists! A file named '",  
                    folder, 
                    "' is created. ", sep = "")); 
      pathvec = pathout(path, folder, method = "StrPBR");
    }
    for(i in c(1 : 4, length(pathvec) - 1)){
      dir.create(pathvec[i]); 
    }
    warnings("off");
    save(pathvec, file = pathvec[length(pathvec)])
    
    Rcov = EnterCov(pathvec);
    if(is.null(Rcov)){
      message("Forced cessation of execution.\n")
      unlink(pathvec[1], recursive = TRUE); 
      return(NULL);
    }else{
      cov_num = Rcov$cov_num; 
      save(cov_num, file = pathvec[6])
      covariate = Rcov$covariate; 
      covr = Rcov$covr; 
      CovIndex = Rcov$CovIndex; 
      save(CovIndex, file = pathvec[8]); 
      
      pth = character();
      for(i in 1 : cov_num){
        a = paste("LevIndex", i, sep = "");
        b = paste(a, "RData", sep = ".");
        pth[i] = paste(pathvec[2], b, sep = "/");
      }
      
      Rlev = EnterLev(pathvec, cov_num, covariate, covr, pth);
      
      if(is.null(Rlev)){
        message("Forced cessation of execution.\n")
        unlink(pathvec[1], recursive = TRUE); 
        return(NULL); 
      }else{
        level_num = Rlev$level_num;
        save(level_num, file = pathvec[7]);
        
        D = matrix(0, ncol = 1 + prod(level_num) + sum(level_num), nrow = 1);
        
        message("Please enter the block size: \n"); 
        bsize = readline(prompt = "Enter the block size: ");
        bsize = as.integer(bsize); 
        warnings("off"); 
        while (is.na(bsize) || bsize %% 2 != 0){
          message("INVALID INPUT !", "\n", 
                  "Please enter a positive number which is also a multiple of 2. : \n"); 
          bsize = readline(prompt = "Reenter the block size "); 
          bsize = as.integer(bsize); 
          warnings("off"); 
        }
        save(bsize, file = pathvec[13]); 
        B = Bpert(bsize, 2); 
        save(B, file = pathvec[11]); 
        
        Bsize = ncol(B); 
        Psize = prod(level_num); 
        
        BG = B[, sample(1 : Bsize, Psize, replace = T)]; 
        save(BG, file = pathvec[10]);
        
        strp = rep(0, times = Psize); 
      }
    }
  }else{
    
    pathverify = paste(path, folder, sep = "/"); 
    while(dir.exists(pathverify)){
      folder = paste(folder, "tymh", sep = "_")
      pathverify = paste(path, folder, sep = "/"); 
    }
    folderpathlog = sub('.....$', '', folder)
    
    load(paste(path, folderpathlog, 
               "static", "pathlog", "pathvec.RData", 
               sep = "/")); 
    
    for(m in 6 : 13){
      load(pathvec[m]);
    }
    pth = character(); 
    for(i in 1 : cov_num){
      a = paste("LevIndex", i, sep = "");
      b = paste(a, "RData", sep = ".");
      pth[i] = paste(pathvec[2], b, sep = "/");
    }
    covariate = as.character(CovIndex[, 1]);
    covr = CovIndex[, 2]; 
  }
  R = list(); 
  CPR = EnterCovPrf(cov_num, pth, CovIndex, covariate, R); 
  message("\n"); 
  message("Reenter COVARIATE PROFILE or not?"); 
  wc = readline("Enter y or n: "); 
  while(wc == "y" || wc == "Y" || wc == "yes" || wc == "YES" || wc == "Yes"){
    CPR = EnterCovPrf(cov_num, pth, CovIndex, covariate, R); 
    message("\n"); 
    message("Reenter COVARIATE PROFILE or not?"); 
    wc = readline("Enter y or n: "); 
  }
  
  PS = PStrGen(cov_num, level_num);
  RES = StrROne(t(D), PS, CPR$cov_profile, cov_num, level_num, 
                bsize, B, BG, strp);
  
  CPR$covariate = covariate; CPR$covr = covr; CPR$cov_num = cov_num; 
  CPR$method = "Stratified Permuted Block Randomization";
  
  D = t(RES[4, 1][[1]]); 
  colnames(D) = nameString(cov_num, level_num, prod(level_num), "All", PS); 
  save(D, file = pathvec[9]);
  
  strp = RES[1, 1][[1]]; 
  rownames(strp) = BBCDname(prod(level_num), "level-"); 
  save(strp, file = pathvec[12]); 
  
  ass = RES[3, 1][[1]]; 
  CPR$assignment = ass[1, 1]; 
  
  class(CPR) = "carseq";
  return(CPR);
}

##############################################################################################################
## Atkinson's Optimum Biased Coin Design
##############################################################################################################
#DoptBCD.ui.carseq = function(path, folder = "DoptBCD") UseMethod("carseq")

DoptBCD.ui = function(path, folder = "DoptBCD"){
  
  message("Is this the first patient? ");
  first = readline(prompt = "Enter T or F: ");
  if(first == "T" || first == "True" || first == "TRUE"){
    pathvec = pathout(path, folder, method = "DoptBCD"); 
    while(dir.exists(pathvec[1])){
      folder = paste(folder, "tymh", sep = "_"); 
      warning(paste("Homonym exists! A file named '",  
                    folder, 
                    "' is created. ", sep = "")); 
      pathvec = pathout(path, folder, method = "DoptBCD");
    }
    for(i in c(1 : 4, length(pathvec) - 1)){
      dir.create(pathvec[i]); 
    }
    warnings("off");
    save(pathvec, file = pathvec[length(pathvec)]); 
    
    Rcov = EnterCov(pathvec);
    if(is.null(Rcov)){
      message("Forced cessation of execution.\n")
      unlink(pathvec[1], recursive = TRUE); 
      return(NULL);
    }else{
      cov_num = Rcov$cov_num; 
      save(cov_num, file = pathvec[6])
      covariate = Rcov$covariate; 
      covr = Rcov$covr; 
      CovIndex = Rcov$CovIndex; 
      save(CovIndex, file = pathvec[8]); 
      
      pth = character();
      for(i in 1 : cov_num){
        a = paste("LevIndex", i, sep = "");
        b = paste(a, "RData", sep = ".");
        pth[i] = paste(pathvec[2], b, sep = "/");
      }
      
      Rlev = EnterLev(pathvec, cov_num, covariate, covr, pth);
      
      if(is.null(Rlev)){
        message("Forced cessation of execution.\n")
        unlink(pathvec[1], recursive = TRUE); 
        return(NULL); 
      }else{
        level_num = Rlev$level_num;
        save(level_num, file = pathvec[7]);
        
        D = matrix(0, ncol = 1 + prod(level_num) + sum(level_num), nrow = 1);
        strp = matrix(0, nrow = prod(level_num), ncol = 1); 
        No = 0; 
        Fmatrix = matrix(NA, nrow = 1 + cov_num, ncol = 50000);
        Fmatrix[1, ] = 1; 
        b = rep(0, times = 50000); 
      }
    }
  }else{
    
    pathverify = paste(path, folder, sep = "/"); 
    while(dir.exists(pathverify)){
      folder = paste(folder, "tymh", sep = "_")
      pathverify = paste(path, folder, sep = "/"); 
    }
    folderpathlog = sub('.....$', '', folder)
    
    load(paste(path, folderpathlog, 
               "static", "pathlog", "pathvec.RData", 
               sep = "/")); 
    
    for(m in 6 : 13){
      load(pathvec[m]);
    }
    pth = character(); 
    for(i in 1 : cov_num){
      a = paste("LevIndex", i, sep = "");
      bc = paste(a, "RData", sep = ".");
      pth[i] = paste(pathvec[2], bc, sep = "/");
    }
    covariate = as.character(CovIndex[, 1]);
    covr = CovIndex[, 2]; 
  }
  R = list(); 
  CPR = EnterCovPrf(cov_num, pth, CovIndex, covariate, R); 
  message("\n"); 
  message("Reenter COVARIATE PROFILE or not?"); 
  wc = readline("Enter y or n: "); 
  while(wc == "y" || wc == "Y" || wc == "yes" || wc == "YES" || wc == "Yes"){
    CPR = EnterCovPrf(cov_num, pth, CovIndex, covariate, R); 
    message("\n"); 
    message("Reenter COVARIATE PROFILE or not?"); 
    wc = readline("Enter y or n: "); 
  }
  
  PS = PStrGen(cov_num, level_num); 
  RES = AtBCDOne(t(D), PS, CPR$cov_profile, cov_num,
                 level_num, Fmatrix, b, strp, No);
  
  Fmatrix = RES[3, 1][[1]]; 
  save(Fmatrix, file = pathvec[10]);
  No = RES[2, 1][[1]]; 
  save(No, file = pathvec[11]); 
  b = RES[4, 1][[1]]; 
  save(b, file = pathvec[12]); 
  CPR$covariate = covariate; CPR$covr = covr; CPR$cov_num = cov_num; 
  CPR$method = "Atkinson's Optimum Biased Coin Design";
  
  D = t(RES[6, 1][[1]]); 
  colnames(D) = nameString(cov_num, level_num, prod(level_num), "All", PS); 
  save(D, file = pathvec[9]);
  
  strp = RES[1, 1][[1]]; 
  rownames(strp) = BBCDname(prod(level_num), "level-"); 
  save(strp, file = pathvec[13]); 
  
  ass = RES[5, 1][[1]]; 
  CPR$assignment = ass[1, 1]; 
  
  class(CPR) = "carseq";
  return(CPR);
}

##############################################################################################################
## Covariate-adaptive Biased Coin Design
##############################################################################################################
#AdjBCD.ui.carseq = function(path, folder = "AdjBCD") UseMethod("carseq")

AdjBCD.ui = function(path, folder = "AdjBCD"){
  
  message("Is this the first patient? ");
  first = readline(prompt = "Enter T or F: ");
  if(first == "T" || first == "True" || first == "TRUE"){
    pathvec = pathout(path, folder, method = "AdjBCD"); 
    while(dir.exists(pathvec[1])){
      folder = paste(folder, "tymh", sep = "_"); 
      warning(paste("Homonym exists! A file named '",  
                    folder, 
                    "' is created. ", sep = "")); 
      pathvec = pathout(path, folder, method = "AdjBCD");
    }
    for(i in c(1 : 4, length(pathvec) - 1)){
      dir.create(pathvec[i]); 
    }
    warnings("off");
    save(pathvec, file = pathvec[length(pathvec)])
    
    Rcov = EnterCov(pathvec);
    if(is.null(Rcov)){
      message("Forced cessation of execution.\n")
      unlink(pathvec[1], recursive = TRUE); 
      return(NULL);
    }else{
      cov_num = Rcov$cov_num; 
      save(cov_num, file = pathvec[6])
      covariate = Rcov$covariate; 
      covr = Rcov$covr; 
      CovIndex = Rcov$CovIndex; 
      save(CovIndex, file = pathvec[8]); 
      
      pth = character();
      for(i in 1 : cov_num){
        ac = paste("LevIndex", i, sep = "");
        b = paste(ac, "RData", sep = ".");
        pth[i] = paste(pathvec[2], b, sep = "/");
      }
      
      Rlev = EnterLev(pathvec, cov_num, covariate, covr, pth);
      
      if(is.null(Rlev)){
        message("Forced cessation of execution.\n")
        unlink(pathvec[1], recursive = TRUE); 
        return(NULL); 
      }else{
        level_num = Rlev$level_num;
        save(level_num, file = pathvec[7]);
        
        D = matrix(0, ncol = 1 + prod(level_num) + sum(level_num), nrow = 1);
        strp = matrix(0, nrow = prod(level_num), ncol = 1); 
        
        message("Please Enter the degree of randomness parameter a: \n");
        a = readline(prompt = "Enter a: "); 
        a = as.double(a); 
        warnings("off"); 
        while(is.na(a)){
          message("  a must be a positive numeric number!")
          a = readline(prompt = "Reenter a: "); 
          a = as.integer(a);
          warnings("off"); 
        }
        save(a, file = pathvec[10])
      }
    }
  }else{
    
    pathverify = paste(path, folder, sep = "/"); 
    while(dir.exists(pathverify)){
      folder = paste(folder, "tymh", sep = "_")
      pathverify = paste(path, folder, sep = "/"); 
    }
    folderpathlog = sub('.....$', '', folder)
    
    load(paste(path, folderpathlog, 
               "static", "pathlog", "pathvec.RData", 
               sep = "/")); 
    
    for(m in 6 : 11){
      load(pathvec[m]);
    }
    pth = character(); 
    for(i in 1 : cov_num){
      ac = paste("LevIndex", i, sep = "");
      b = paste(ac, "RData", sep = ".");
      pth[i] = paste(pathvec[2], b, sep = "/");
    }
    covariate = as.character(CovIndex[, 1]);
    covr = CovIndex[, 2]; 
  }
  R = NULL; 
  CPR = EnterCovPrf(cov_num, pth, CovIndex, covariate, R); 
  message("\n"); 
  message("Reenter COVARIATE PROFILE or not?"); 
  wc = readline("Enter y or n: "); 
  while(wc == "y" || wc == "Y" || wc == "yes" || wc == "YES" || wc == "Yes"){
    CPR = EnterCovPrf(cov_num, pth, CovIndex, covariate, R); 
    message("\n"); 
    message("Reenter COVARIATE PROFILE or not?"); 
    wc = readline("Enter y or n: "); 
  }
  
  PS = PStrGen(cov_num, level_num); 
  RES = AdBCDOne(t(D), PS, CPR$cov_profile, 
                 cov_num, level_num, strp, a);
  
  CPR$covariate = covariate; CPR$covr = covr; CPR$cov_num = cov_num; 
  CPR$method = "Atkinson's Optimum Biased Coin Design";
  
  D = t(RES[3, 1][[1]]); 
  colnames(D) = nameString(cov_num, level_num, prod(level_num), "All", PS); 
  save(D, file = pathvec[9]);
  
  strp = RES[1, 1][[1]]; 
  rownames(strp) = BBCDname(prod(level_num), "level-"); 
  save(strp, file = pathvec[11]); 
  
  ass = RES[2, 1][[1]]; 
  CPR$assignment = ass[1, 1]; 
  
  class(CPR) = "carseq";
  return(CPR);
}

# ##############################################################################################################
# ## Biased Coin Design with a Bayesian Bias
# ##############################################################################################################
# BayesBCD.ui.carseq = function(path, folder = "BayesBCD") UseMethod("carseq")
# 
# BayesBCD.ui = function(path, folder = "BayesBCD"){
#   
#   
#   pathvec = pathout(path, folder, method = "BayesBCD"); 
#   
#   message("Is this the first patient? ");
#   first = readline(prompt = "Enter T or F: ");
#   if(first == "T" || first == "True" || first == "TRUE"){
#     for(i in 1 : 4){
#       dir.create(pathvec[i]); 
#     }
#     warnings("off");
#     
#     Rcov = EnterCov(pathvec);
#     if(is.null(Rcov)){
#       message("Forced cessation of execution.\n")
#       unlink(pathvec[1], recursive = TRUE); 
#       return(NULL);
#     }else{
#       cov_num = Rcov$cov_num; 
#       save(cov_num, file = pathvec[6])
#       covariate = Rcov$covariate; 
#       covr = Rcov$covr; 
#       CovIndex = Rcov$CovIndex; 
#       save(CovIndex, file = pathvec[8]); 
#       
#       pth = character();
#       for(i in 1 : cov_num){
#         a = paste("LevIndex", i, sep = "");
#         b = paste(a, "RData", sep = ".");
#         pth[i] = paste(pathvec[2], b, sep = "/");
#       }
#       
#       Rlev = EnterLev(pathvec, cov_num, covariate, covr, pth);
#       
#       if(is.null(Rlev)){
#         message("Forced cessation of execution.\n")
#         unlink(pathvec[1], recursive = TRUE); 
#         return(NULL); 
#       }else{
#         level_num = Rlev$level_num;
#         save(level_num, file = pathvec[7]);
#         
#         D = matrix(0, ncol = 1 + prod(level_num) + sum(level_num), nrow = 1);
#         strp = matrix(0, nrow = prod(level_num), ncol = 1); 
#         
#         message("Please Enter the number of category of interest: \n");
#         J = readline(prompt = "Enter J: "); 
#         J = as.integer(J); 
#         warnings("off"); 
#         while(is.na(a)){
#           message("  J must be a positive integer!")
#           J = readline(prompt = "Reenter J: "); 
#           J = as.integer(J);
#           warnings("off"); 
#         }
#         J = abs(J); 
#         save(J, file = pathvec[10])
#         
#         numJ = matrix(0, nrow = 2, ncol = J); 
#         No = 0; 
#       }
#     }
#   }else{
#     for(m in 6 : 13){
#       load(pathvec[m]);
#     }
#     pth = character(); 
#     for(i in 1 : cov_num){
#       ac = paste("LevIndex", i, sep = "");
#       b = paste(ac, "RData", sep = ".");
#       pth[i] = paste(pathvec[2], b, sep = "/");
#     }
#     covariate = as.character(CovIndex[, 1]);
#     covr = CovIndex[, 2]; 
#   }
#   R = list(); 
#   CPR = EnterCovPrf(cov_num, pth, CovIndex, covariate, R); 
#   message("\n"); 
#   message("Reenter COVARIATE PROFILE or not?"); 
#   wc = readline("Enter y or n: "); 
#   while(wc == "y" || wc == "Y" || wc == "yes" || wc == "YES" || wc == "Yes"){
#     CPR = EnterCovPrf(cov_num, pth, CovIndex, covariate, R); 
#     message("\n"); 
#     message("Reenter COVARIATE PROFILE or not?"); 
#     wc = readline("Enter y or n: "); 
#   }
#   
#   RES = BBCDOne(t(D), PStrGen(cov_num, level_num), CPR$cov_profile, 
#                 cov_num, level_num, numJ, strp, J, No);
#   
#   CPR$covariate = covariate; CPR$covr = covr; CPR$cov_num = cov_num; 
#   CPR$method = "Atkinson's Optimum Biased Coin Design";
#   
#   D = t(RES[5, 1][[1]]); 
#   colnames(D) = nameString(cov_num, level_num, prod(level_num), "All", "Real"); 
#   save(D, file = pathvec[9]);
#   
#   strp = RES[1, 1][[1]]; 
#   rownames(strp) = BBCDname(prod(level_num), "level-"); 
#   save(strp, file = pathvec[13]); 
#   
#   No_t = RES[2, 1][[1]];
#   No = No_t[1, 1]; 
#   save(No, file = pathvec[12]); 
#   
#   numJ = RES[3, 1][[1]]; 
#   rownames(numJ) = c("No._of_t1", "No._of_t2"); 
#   colnames(numJ) = BBCDname(J, "Category-"); 
#   save(numJ, file = pathvec[11]); 
#   
#   ass = RES[4, 1][[1]]; 
#   CPR$assignment = ass[1, 1]; 
#   
#   class(CPR) = "carseq";
#   return(CPR);
# }
# 
