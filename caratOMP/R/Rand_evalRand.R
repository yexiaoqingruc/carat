###############################################################################
#############################   Summary   ##################################
################################################################################
## Assess the procedure for real data
#evalRand.careval = function(data, method = "HuHuCAR", N = 500, ...) UseMethod("careval")

evalRand = function(data, method = "HuHuCAR", N = 500, ...){
  
  R = NULL; 
  meth = c("HuHuCAR", "PocSimMIN", "StrBCD", "StrPBR", 
           "DoptBCD", "AdjBCD", "BayesBCD"); 
  if(!method %in% meth){
    stop("Invalid input of method!")
  }
  var = list(...); 
  clasv = as.character(lapply(var, class)); 
  clasind = which(clasv == "numeric"); 
  varl = var[clasind]; 
  varnam = names(varl); 
  leng = length(varl); 
  lengv = as.numeric(lapply(varl, length)); 
  
  if(length(data[is.na(data)]) == 0){
    datap = data;
  }else{
    data[is.na(data)] = "HFFMWYXQTFY<= 1.0";
    datap = data;
  }
  
  dataverify = unique(apply(data[1, ], 2, class)); 
  if(length(dataverify) == 1 && dataverify[1] == "integer"||length(dataverify) == 1 && dataverify[1] == "numeric"){
    data_proc = t(data); 
    cov_num = nrow(data_proc); 
    level_num = vector(); 
    for(i in 1 : cov_num){
      level_num[i] = length(unique(data_proc[i, ])); 
    }
    R$datanumeric = TRUE; 
  }else{
    rdata = Preprocess(datap); 
    data_proc = rdata$data; 
    cov_num = rdata$cov_num; level_num = rdata$level_num; 
    R$datanumeric = FALSE; 
  }
  
  n = ncol(data_proc); 
  if("bsize" %in% varnam){
    bsize = varl[[which(varnam == "bsize")]]; 
  }else{
    bsize = 4; 
  }
  if(method == "HuHuCAR"){
    if(leng < 2){
      stop("omega and p are all needed for method == 'HuHuCAR'!"); 
    }else if(is.null(varnam)){
      indp1 = which(lengv == 1);
      pchoice = as.numeric(varl[indp1]); 
      indp2 = which((pchoice > 0.5) == (pchoice < 1)); 
      if(length(indp2) == 0){
        stop("please set biased probability p to be between 0.5 and 1!"); 
      }else{
        p = pchoice[indp2[1]];
      }
      indw = which(lengv == (cov_num + 2)); 
      if(length(indw) == 0){
        omega = rep(1.0 / (cov_num + 2), times = cov_num + 2); 
        warning("omega is mandated to be equal vector!")
      }else{
        otemp = varl[[indw[1]]]; 
        omega = abs(otemp) / sum(otemp); 
      }
    }else{
      if("p" %in% varnam){
        pind = which(varnam == "p"); 
        p = varl[[pind]]; 
        if(p <= 0.5 || p >= 1){
          stop("Set p to be between 0.5 and 1 to obtain balance!"); 
        }else{
          p = p; 
        }
      }else if(length(which(lengv[which(varnam == "")] == 1)) > 0){
        potherind1 = which(lengv[which(varnam == "")] == 1); 
        pother = as.numeric(varl[potherind1]); 
        potherind2 = which((pother > 0.5) == (pother < 1)); 
        if(length(potherind2) == 0){
          stop("please set biased probability p to be between 0.5 and 1!"); 
        }else{
          p = pother[potherind2[1]]; 
        }
      }else{
        stop("p is needed for method = 'HuHuCAR'!"); 
      }
      if("omega" %in% varnam){
        indweight = which(varnam == "omega");
        omeg = varl[[indweight]]; 
        if(length(omeg) != (cov_num + 2)){
          stop("Length of omega must equal to (cov_num + 2) !"); 
        }
        omega = abs(omeg) / sum(abs(omeg)); 
      }else if(length(which(lengv[which(varnam == "")] == cov_num + 2)) > 0){
        indweight = which(lengv[which(varnam == "")] == (cov_num + 2)); 
        omeg = varl[[indweight]]; 
        omega = abs(omeg) / sum(abs(omeg)); 
      }else{
        stop("omega is needed for method = 'HuHuCAR'!"); 
      }
    }
    
    RES = C_RSummarize(data_proc, cov_num, level_num, method, omega, p, 
                       bsize = bsize, a = 2, N = N); 
    R$weight = omega[2 : (2 + cov_num)]; 
    R$bsize = bsize; 
  }
  if(method == "PocSimMIN"){
    if(leng < 2){
      stop("weight and p are all needed for method = 'PocSimMIN'!"); 
    }else if(is.null(varnam)){
      indp1 = which(lengv == 1); 
      pchoice = as.numeric(varl[indp1]); 
      indp2 = which((pchoice > 0.5) == (pchoice < 1)); 
      if(length(indp2) == 0){
        stop("please set biased probability p to be between 0.5 and 1!"); 
      }else{
        p = pchoice[indp2[1]];
      }
      indw = which(lengv == cov_num); 
      if(length(indw) == 0){
        omega = c(0, 0, rep(1.0 / cov_num, times = cov_num)); 
        warning("weight is mandated to be equal vector!")
      }else{
        wtemp = varl[[indw[1]]]; 
        omega = c(0, 0, abs(wtemp) / sum(wtemp)); 
      }
    }else{
      if("p" %in% varnam){
        pind = which(varnam == "p"); 
        p = varl[[pind]]; 
        if(p <= 0.5 || p > 1){
          stop("Set p to be between 0.5 and 1 to obtain balance!"); 
        }else{
          p = p; 
        }
      }else if(length(which(lengv[which(varnam == "")] == 1)) > 0){
        potherind1 = which(lengv[which(varnam == "")] == 1); 
        pother = as.numeric(varl[potherind1]); 
        potherind2 = which((pother > 0.5) == (pother < 1)); 
        if(length(potherind2) == 0){
          stop("please set biased probability p to be between 0.5 and 1!"); 
        }else{
          p = pother[potherind2[1]]; 
        }
      }else{
        stop("p is needed for method = 'PocSimMIN'!"); 
      }
      if("weight" %in% varnam){
        indweight = which(varnam == "weight");
        weig = varl[[indweight]]; 
        if(length(weig) != cov_num){
          stop("Length of weight must equal to cov_num!"); 
        }
        omega = c(0, 0, abs(weig) / sum(abs(weig))); 
      }else if(length(which(lengv[which(varnam == "")] == cov_num)) > 0){
        indweight = which(lengv[which(varnam == "")] == cov_num); 
        weig = varl[[indweight]]; 
        omega = c(0, 0, abs(weig) / sum(abs(weig))); 
      }else{
        stop("weight is needed for method = 'PocSimMIN'!"); 
      }
    }
    
    RES = C_RSummarize(data_proc, cov_num, level_num, method, omega, p, 
                       bsize = bsize, a = 2, N = N); 
    R$weight = omega[2 : (2 + cov_num)]; 
    R$bsize = bsize; 
  }
  if(method == "StrBCD"){
    if(leng < 1){
      stop("p is needed for method = 'StrBCD'!"); 
    }else if(is.null(varnam)){
      indp1 = which(lengv == 1); 
      pchoice = as.numeric(varl[indp1]); 
      indp2 = which((pchoice > 0.5) == (pchoice < 1)); 
      if(length(indp2) == 0){
        stop("please set biased probability p to be between 0.5 and 1!"); 
      }else{
        p = pchoice[indp2[1]];
      }
    }else{
      if("p" %in% varnam){
        pind = which(varnam == "p"); 
        p = varl[[pind]]; 
        if(p <= 0.5 || p > 1){
          stop("Set p to be between 0.5 and 1 to obtain balance!"); 
        }else{
          p = p; 
        }
      }else if(length(which(lengv[which(varnam == "")] == 1)) > 0){
        potherind1 = which(lengv[which(varnam == "")] == 1); 
        pother = as.numeric(varl[potherind1]); 
        potherind2 = which((pother > 0.5) == (pother < 1)); 
        if(length(potherind2) == 0){
          stop("please set biased probability p to be between 0.5 and 1!"); 
        }else{
          p = pother[potherind2[1]]; 
        }
      }else{
        stop("p is needed for method = 'StrBCD'!"); 
      }
    }
    omega = c(0, 1, rep(0, times = cov_num)); 
    RES = C_RSummarize(data_proc, cov_num, level_num, method, omega, p, 
                       bsize = bsize, a = 2, N = N); 
    R$bsize = bsize; 
  }
  if(method == "StrPBR"){
    if(leng < 1){
      stop("bsize is needed for method = 'StrPBR' !"); 
    }else if (is.null(varnam)){
      indb1 = which(lengv == 1); 
      bchoice = as.numeric(varl[indb1]); 
      bremind = which((bchoice %% 2) == 0); 
      if(length(bremind) > 0){
        bsize = bchoice[bremind[1]]; 
      }else{
        bsint = as.integer(bchoice); 
        bsintremind = which((bsint %% 2) == 0); 
        if(length(bsintremind) > 0){
          bsize = bsint[bsintremind[1]]; 
          warning("bsize is mandated to be an integer!"); 
        }else{
          stop("bsize must be a multiple of 2!");
        }
      }
    }else{
      if("bsize" %in% varnam){
        bsind = which(varnam == "bsize"); 
        bsize = varl[[bsind]]; 
        if((bsize %% 2) == 0){
          bsize = bsize; 
        }else{
          if((as.integer(bsize) %% 2) == 0){
            bsize = as.integer(bsize); 
            warning("bsize is mandated to be a multiple of 2!"); 
          }else{
            stop("bsize must be a multiple of 2!")
          }
        }
      }else if(length(which(lengv[which(varnam == "")] == 1)) > 0){
        bsnullind1 = which(lengv[which(varnam == "")] == 1); 
        bsnull = as.numeric(varl[bsnullind1]); 
        bsnullind2 = which((bsnull %% 2) == 0); 
        bsnullind3 = which((as.integer(bsnull) %% 2) == 0); 
        if(length(bsnullind2) > 0){
          bsize = bsnull[bsnullind2[1]]; 
        }else if(length(bsnullind2) == 0 && length(bsnullind3) > 0){
          bsize = as.integer(bsnull[bsnullind3[1]]); 
          warning("bsize is mandated to be a multiple of 2!")
        }else{
          stop("bsize must be a multiple of 2!"); 
        }
      }else{
        stop("bsize is needed for method = 'StrPBR'!"); 
      }
    }
    RES = C_RSummarize(data_proc, cov_num, level_num, method, omega = c(1), p = 0, 
                       bsize = bsize, a = 2, N = N); 
    R$bsize = bsize; 
  }
  if(method == "DoptBCD"){
    RES = C_RSummarize(data_proc, cov_num, level_num, method, omega = c(1), p = 0, 
                       bsize = bsize, a = 2, N = N); 
    R$bsize = bsize; 
  }
  if(method == "AdjBCD"){
    if(leng < 1){
      stop("design parameter is needed for method = 'AdjBCD'!"); 
    }else if (is.null(varnam)){
      inda1 = which(lengv == 1); 
      achoice = as.numeric(varl[inda1]); 
      inda2 = which(achoice != 0);
      indaabs = which(achoice > 0); 
      if(length(inda2) > 0 && length(indaabs) > 0){
        a = varl[[indaabs[1]]]; 
      }else if(length(inda2) > 0 && length(indaabs) == 0){
        a = abs(varl[[inda2[1]]]); 
        warning("a is mandated to be positive!"); 
      }else{
        stop("a must be a positive numeric number!");
      }
    }else{
      if("a" %in% varnam){
        aind = which(varnam == "a"); 
        a = varl[[aind]]; 
      }else if(length(which(lengv[which(varnam == "")] == 1)) > 0){
        aotherind = which(lengv[which(varnam == "")] == 1); 
        a = varl[[aotherind[1]]]; 
      }else{
        stop("a is needed for method = 'AdjBCD'!"); 
      }
    }
    RES = C_RSummarize(data_proc, cov_num, level_num, method, omega = c(0), p = 0, 
                       bsize = 4, a = a, N = N); 
    R$bsize = bsize; 
  }
  
  covn = colnames(data); 
  R$covariates = covn; 
  
  A = RES[1, 1][[1]];
  colnames(A) = BBCDname(N, "iter"); 
  rownames(A) = BBCDname(n, "pat");
  R$Assig = A;
  
  PS0 = RES[4, 1][[1]]; 
  PS0ordinPS = MVReturnM(PStrGen(cov_num, level_num), PS0); 
  PS = PS0[, order(PS0ordinPS)]; 
  strt_num = ncol(PS); 
  R$strt_num = strt_num; 
  colnames(PS) = BBCDname(strt_num, "stratum"); 
  rownames(PS) = paste("covariate", 1 : cov_num,"(", covn, ")", sep = ""); 
  R$`All strata` = PS;
  
  Imbmat0 = RES[2, 1][[1]]; 
  Imbmat = Imbmat0[c(1, 1 + order(PS0ordinPS), (1 + strt_num + 1) : (1 + strt_num + sum(level_num))), ]
  colnames(Imbmat) = c("max", "95% quan", "median", "mean");
  rownames(Imbmat) = nameString(cov_num, level_num, strt_num, "All", PS); 
  R$Imb = Imbmat; 
  
  SNUM0 = RES[3, 1][[1]]; 
  SNUM = SNUM0[order(PS0ordinPS), , drop = FALSE]; 
  rownames(SNUM) = nameString(cov_num, level_num, strt_num, "All", PS)[2 : (1 + strt_num)];  
  colnames(SNUM) = BBCDname(N, "iter")
  R$SNUM = SNUM; 
  
  R$method = method; 
  R$cov_num = cov_num; 
  R$level_num = level_num; 
  R$n = ncol(data_proc);
  R$iteration = N;
  R$'Data Type' = "Real"; 
  
  DIF0 = RES[5, 1][[1]]; 
  DIF = DIF0[c(1, 1 + order(PS0ordinPS), (1 + strt_num + 1) : (1 + strt_num + sum(level_num))), , drop = FALSE]
  colnames(DIF) = BBCDname(N, "iter"); 
  rownames(DIF) = nameString(cov_num, level_num, strt_num, "All", PS); 
  R$DIF = DIF; 
  R$data = data; 
  
  
  class(R) = "careval"; 
  return(R);
}

## Assess the procedure for simulated data
# evalRand.sim.careval = function(n = 1000,  N = 500, Replace = FALSE, cov_num = 2, level_num = c(2, 2),
#                                 pr = rep(0.5, 4), method = "HuHuCAR", ...) UseMethod("careval")

evalRand.sim = function(n = 1000, N = 500, Replace = FALSE, cov_num = 2, 
                        level_num = c(2, 2), pr = rep(0.5, 4), 
                        method = "HuHuCAR", ...){
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
  R = NULL; 
  meth = c("HuHuCAR", "PocSimMIN", "StrBCD", "StrPBR", 
           "DoptBCD", "AdjBCD", "BayesBCD"); 
  if(!method %in% meth){
    stop("Invalid input of method!")
  }
  var = list(...); 
  clasv = as.character(lapply(var, class)); 
  clasind = which(clasv == "numeric"); 
  varl = var[clasind]; 
  varnam = names(varl); 
  leng = length(varl); 
  lengv = as.numeric(lapply(varl, length)); 
  if("bsize" %in% varnam){
    bsize = varl[[which(varnam == "bsize")]]; 
  }else{
    bsize = 4; 
  }
  if(method == "HuHuCAR"){
    if(leng < 2){
      stop("omega and p are all needed for method == 'HuHuCAR'!"); 
    }else if(is.null(varnam)){
      indp1 = which(lengv == 1);
      pchoice = as.numeric(varl[indp1]); 
      indp2 = which((pchoice > 0.5) == (pchoice < 1)); 
      if(length(indp2) == 0){
        stop("please set biased probability p to be between 0.5 and 1!"); 
      }else{
        p = pchoice[indp2[1]];
      }
      indw = which(lengv == (cov_num + 2)); 
      if(length(indw) == 0){
        omega = rep(1.0 / (cov_num + 2), times = cov_num + 2); 
        warning("omega is mandated to be equal vector!")
      }else{
        otemp = varl[[indw[1]]]; 
        omega = abs(otemp) / sum(otemp); 
      }
    }else{
      if("p" %in% varnam){
        pind = which(varnam == "p"); 
        p = varl[[pind]]; 
        if(p <= 0.5 || p >= 1){
          stop("Set p to be between 0.5 and 1 to obtain balance!"); 
        }else{
          p = p; 
        }
      }else if(length(which(lengv[which(varnam == "")] == 1)) > 0){
        potherind1 = which(lengv[which(varnam == "")] == 1); 
        pother = as.numeric(varl[potherind1]); 
        potherind2 = which((pother > 0.5) == (pother < 1)); 
        if(length(potherind2) == 0){
          stop("please set biased probability p to be between 0.5 and 1!"); 
        }else{
          p = pother[potherind2[1]]; 
        }
      }else{
        stop("p is needed for method = 'HuHuCAR'!"); 
      }
      if("omega" %in% varnam){
        indweight = which(varnam == "omega");
        omeg = varl[[indweight]]; 
        if(length(omeg) != (cov_num + 2)){
          stop("Length of omega must equal to (cov_num + 2) !"); 
        }
        omega = abs(omeg) / sum(abs(omeg)); 
      }else if(length(which(lengv[which(varnam == "")] == cov_num + 2)) > 0){
        indweight = which(lengv[which(varnam == "")] == cov_num + 2); 
        omeg = varl[[indweight]]; 
        omega = abs(omeg) / sum(abs(omeg)); 
      }else{
        stop("omega is needed for method = 'HuHuCAR'!"); 
      }
    }
    
    RES = C_Summarize(Replace, cov_num, level_num, pr, method, omega, p, 
                      bsize = bsize, a = 2, n = n, N = N); 
    R$weight = omega[2 : (2 + cov_num)]; 
    R$bsize = bsize; 
  }
  if(method == "PocSimMIN"){
    if(leng < 2){
      stop("weight and p are all needed for method = 'PocSimMIN'!"); 
    }else if(is.null(varnam)){
      indp1 = which(lengv == 1); 
      pchoice = as.numeric(varl[indp1]); 
      indp2 = which((pchoice > 0.5) == (pchoice < 1)); 
      if(length(indp2) == 0){
        stop("please set biased probability p to be between 0.5 and 1!"); 
      }else{
        p = pchoice[indp2[1]];
      }
      indw = which(lengv == cov_num); 
      if(length(indw) == 0){
        omega = c(0, 0, rep(1.0 / cov_num, times = cov_num)); 
        warning("weight is mandated to be equal vector!")
      }else{
        wtemp = varl[[indw[1]]]; 
        omega = c(0, 0, abs(wtemp) / sum(wtemp)); 
      }
    }else{
      if("p" %in% varnam){
        pind = which(varnam == "p"); 
        p = varl[[pind]]; 
        if(p <= 0.5 || p > 1){
          stop("Set p to be between 0.5 and 1 to obtain balance!"); 
        }else{
          p = p; 
        }
      }else if(length(which(lengv[which(varnam == "")] == 1)) > 0){
        potherind1 = which(lengv[which(varnam == "")] == 1); 
        pother = as.numeric(varl[potherind1]); 
        potherind2 = which((pother > 0.5) == (pother < 1)); 
        if(length(potherind2) == 0){
          stop("please set biased probability p to be between 0.5 and 1!"); 
        }else{
          p = pother[potherind2[1]]; 
        }
      }else{
        stop("p is needed for method = 'PocSimMIN'!"); 
      }
      if("weight" %in% varnam){
        indweight = which(varnam == "weight");
        weig = varl[[indweight]]; 
        if(length(weig) != cov_num){
          stop("Length of weight must equal to cov_num!"); 
        }
        omega = c(0, 0, abs(weig) / sum(abs(weig))); 
      }else if(length(which(lengv[which(varnam == "")] == cov_num)) > 0){
        indweight = which(lengv[which(varnam == "")] == cov_num); 
        weig = varl[[indweight]]; 
        omega = c(0, 0, abs(weig) / sum(abs(weig))); 
      }else{
        stop("weight is needed for method = 'PocSimMIN'!"); 
      }
    }
    
    RES = C_Summarize(Replace, cov_num, level_num, pr, method, omega, p, 
                      bsize = bsize, a = 2, n = n, N = N); 
    R$weight = omega[2 : (2 + cov_num)]; 
    R$bsize = bsize; 
  }
  if(method == "StrBCD"){
    if(leng < 1){
      stop("p is needed for method = 'StrBCD'!"); 
    }else if(is.null(varnam)){
      indp1 = which(lengv == 1); 
      pchoice = as.numeric(varl[indp1]); 
      indp2 = which((pchoice > 0.5) == (pchoice < 1)); 
      if(length(indp2) == 0){
        stop("please set biased probability p to be between 0.5 and 1!"); 
      }else{
        p = pchoice[indp2[1]];
      }
    }else{
      if("p" %in% varnam){
        pind = which(varnam == "p"); 
        p = varl[[pind]]; 
        if(p <= 0.5 || p > 1){
          stop("Set p to be between 0.5 and 1 to obtain balance!"); 
        }else{
          p = p; 
        }
      }else if(length(which(lengv[which(varnam == "")] == 1)) > 0){
        potherind1 = which(lengv[which(varnam == "")] == 1); 
        pother = as.numeric(varl[potherind1]); 
        potherind2 = which((pother > 0.5) == (pother < 1)); 
        if(length(potherind2) == 0){
          stop("please set biased probability p to be between 0.5 and 1!"); 
        }else{
          p = pother[potherind2[1]]; 
        }
      }else{
        stop("p is needed for method = 'StrBCD'!"); 
      }
    }
    omega = c(0, 1, rep(0, times = cov_num)); 
    RES = C_Summarize(Replace, cov_num, level_num, pr, method, omega, p, 
                      bsize = bsize, a = 2, n = n, N = N); 
    R$bsize = bsize; 
  }
  if(method == "StrPBR"){
    if(leng < 1){
      stop("bsize is needed for method = 'StrPBR' !"); 
    }else if (is.null(varnam)){
      indb1 = which(lengv == 1); 
      bchoice = as.numeric(varl[indb1]); 
      bremind = which((bchoice %% 2) == 0); 
      if(length(bremind) > 0){
        bsize = bchoice[bremind[1]]; 
      }else{
        bsint = as.integer(bchoice); 
        bsintremind = which((bsint %% 2) == 0); 
        if(length(bsintremind) > 0){
          bsize = bsint[bsintremind[1]]; 
          warning("bsize is mandated to be an integer!"); 
        }else{
          stop("bsize must be a multiple of 2!");
        }
      }
    }else{
      if("bsize" %in% varnam){
        bsind = which(varnam == "bsize"); 
        bsize = varl[[bsind]]; 
        if((bsize %% 2) == 0){
          bsize = bsize; 
        }else{
          if((as.integer(bsize) %% 2) == 0){
            bsize = as.integer(bsize); 
            warning("bsize is mandated to be a multiple of 2!"); 
          }else{
            stop("bsize must be a multiple of 2!")
          }
        }
      }else if(length(which(lengv[which(varnam == "")] == 1)) > 0){
        bsnullind1 = which(lengv[which(varnam == "")] == 1); 
        bsnull = as.numeric(varl[bsnullind1]); 
        bsnullind2 = which((bsnull %% 2) == 0); 
        bsnullind3 = which((as.integer(bsnull) %% 2) == 0); 
        if(length(bsnullind2) > 0){
          bsize = bsnull[bsnullind2[1]]; 
        }else if(length(bsnullind2) == 0 && length(bsnullind3) > 0){
          bsize = as.integer(bsnull[bsnullind3[1]]); 
          warning("bsize is mandated to be a multiple of 2!")
        }else{
          stop("bsize must be a multiple of 2!"); 
        }
      }else{
        stop("bsize is needed for method = 'StrPBR'!"); 
      }
    }
    RES = C_Summarize(Replace, cov_num, level_num, pr, method, omega = c(0), p = 0, 
                      bsize = bsize, a = 2, n = n, N = N); 
    R$bsize = bsize; 
  }
  if(method == "DoptBCD"){
    RES = C_Summarize(Replace, cov_num, level_num, pr, method, omega = c(0), p = 0, 
                      bsize = bsize, a = 2, n = n, N = N); 
    R$bsize = bsize; 
  }
  if(method == "AdjBCD"){
    if(leng < 1){
      stop("design parameter is needed for method = 'AdjBCD'!"); 
    }else if (is.null(varnam)){
      inda1 = which(lengv == 1); 
      achoice = as.numeric(varl[inda1]); 
      inda2 = which(achoice != 0);
      indaabs = which(achoice > 0); 
      if(length(inda2) > 0 && length(indaabs) > 0){
        a = varl[[indaabs[1]]]; 
      }else if(length(inda2) > 0 && length(indaabs) == 0){
        a = abs(varl[[inda2[1]]]); 
        warning("a is mandated to be positive!"); 
      }else{
        stop("a must be a positive numeric number!");
      }
    }else{
      if("a" %in% varnam){
        aind = which(varnam == "a"); 
        a = varl[[aind]]; 
      }else if(length(which(lengv[which(varnam == "")] == 1)) > 0){
        aotherind = which(lengv[which(varnam == "")] == 1); 
        a = varl[[aotherind[1]]]; 
      }else{
        stop("a is needed for method = 'AdjBCD'!"); 
      }
    }
    RES = C_Summarize(Replace, cov_num, level_num, pr, method, omega = c(0), p = 0, 
                      bsize = bsize, a = a, n = n, N = N); 
    R$bsize = bsize; 
  }
  
  A = RES[1, 1][[1]];
  colnames(A) = BBCDname(N, "iter"); 
  rownames(A) = BBCDname(n, "pat");
  R$Assig = A;
  
  PS0 = RES[4, 1][[1]]; 
  PS0ordinPS = MVReturnM(PStrGen(cov_num, level_num), PS0); 
  PS = PS0[, order(PS0ordinPS)]; 
  strt_num = ncol(PS); 
  R$strt_num = strt_num; 
  colnames(PS) = BBCDname(strt_num, "stratum"); 
  rownames(PS) = BBCDname(cov_num, "covariate"); 
  R$`All strata` = PS;
  
  Imbmat0 = RES[2, 1][[1]]; 
  Imbmat = Imbmat0[c(1, 1 + order(PS0ordinPS), (1 + strt_num + 1) : (1 + strt_num + sum(level_num))), , drop = FALSE]
  colnames(Imbmat) = c("max", "95% quan", "median", "mean");
  rownames(Imbmat) = nameString(cov_num, level_num, strt_num, "All", PS); 
  R$Imb = Imbmat; 
  
  SNUM0 = RES[3, 1][[1]]; 
  SNUM = SNUM0[order(PS0ordinPS), , drop = FALSE]; 
  rownames(SNUM) = nameString(cov_num, level_num, strt_num, "All", PS)[2 : (1 + strt_num)];  
  colnames(SNUM) = BBCDname(N, "iter")
  R$SNUM = SNUM; 
  
  R$method = method; 
  R$cov_num = cov_num; 
  R$level_num = level_num; 
  R$n = n; 
  R$iteration = N; 
  R$'Data Type' = "Simulated"; 
  R$DataGeneration = Replace;
  
  DIF0 = RES[5, 1][[1]]; 
  DIF = DIF0[c(1, 1 + order(PS0ordinPS), (1 + strt_num + 1) : (1 + strt_num + sum(level_num))), , drop = FALSE]
  colnames(DIF) = BBCDname(N, "iter"); 
  rownames(DIF) = nameString(cov_num, level_num, strt_num, "All", PS); 
  R$DIF = DIF; 
  
  class(R) = "careval"; 
  return(R);
}

