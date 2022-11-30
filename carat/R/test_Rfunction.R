getData<-function(n,cov_num,level_num,pr,type,beta,mu1,mu2,sigma = 1,method = "HuHuCAR",...){
  FUN = switch(method,"HuHuCAR" = HuHuCAR_getData, "PocSimMIN" = PocSimMIN_getData,
              "StrBCD" = StrBCD_getData, "StrPBR" = StrPBR_getData,
              "DoptBCD" = DoptBCD_getData, "AdjBCD" = AdjBCD_getData)
  data = FUN(n,cov_num,level_num,pr,type,beta,mu1,mu2,sigma,...)
  names = c(paste0("covariate",1:cov_num),"assignment","outcome")
  datafr = data.frame(data, row.names = names, stringsAsFactors = TRUE)
  return(datafr)
}

rand.test<-function(data,Reps = 200,method = c("HuHuCAR","PocSimMIN","StrBCD","StrPBR","DoptBCD","AdjBCD"),conf = 0.95, binwidth = 30,...){
  dataverify = unique(apply(data[, 1, drop = FALSE], 2, class)); 
  if(!(length(dataverify) == 1 && dataverify[1] == "integer"||length(dataverify) == 1 && dataverify[1] == "numeric")){
    stop("Please numeric data!")
  }
  dname = deparse(substitute(data))
  method = match.arg(method)
  FUN = switch(method,"HuHuCAR" = HuHuCAR_RT, "PocSimMIN" = PocSimMIN_RT,
               "StrBCD" = StrBCD_RT, "StrPBR" = StrPBR_RT,
               "DoptBCD" = DoptBCD_RT, "AdjBCD" = AdjBCD_RT)
  result = FUN(data,Reps,conf,...)
  pval = result$pval
  testmethod<-"Randomization test"
  estimate = result$estimate
  names(estimate)<-"difference for treatment effect"
  rval<-list(p.value = pval,estimate = estimate,conf.int = c(result$CIl,result$CIu),
             method = testmethod,data.name = dname,data = result$Randata,binwidth = binwidth)
  class(rval)<-c("plotrandtest","htest")
  return(rval)
}

boot.test<-function(data,B=200,method = c("HuHuCAR","PocSimMIN","StrBCD","StrPBR","DoptBCD","AdjBCD"),conf = 0.95,...){
  dataverify = unique(apply(data[, 1, drop = FALSE], 2, class)); 
  if(!(length(dataverify) == 1 && dataverify[1] == "integer"||length(dataverify) == 1 && dataverify[1] == "numeric")){
    stop("Please numeric data!")
  }
  dname = deparse(substitute(data))
  method = match.arg(method)
  FUN = switch(method,"HuHuCAR" = HuHuCAR_BT, "PocSimMIN" = PocSimMIN_BT,
               "StrBCD" = StrBCD_BT, "StrPBR" = StrPBR_BT,
               "DoptBCD" = DoptBCD_BT, "AdjBCD" = AdjBCD_BT)
  result = FUN(data,B,...)
  testmethod<-"Bootstrap t-test"
  estimate = result[1]
  stderr = result[2]
  tstat = result[3]
  pval = result[4]*2
  cint = c(estimate + stderr*qnorm((1-conf)/2),estimate - stderr*qnorm((1-conf)/2))
  attr(cint,"conf.level")<-conf
  names(tstat)<-"t"
  names(estimate)<-"difference for treatment effect"
  rval<-list(statistic = tstat,p.value = pval,conf.int = cint,
             estimate = estimate,stderr = stderr,
             method = testmethod,data.name = dname)
  class(rval)<-"htest"
  return(rval)
}

corr.test<-function(data,conf = 0.95){
  dataverify = unique(apply(data[, 1, drop = FALSE], 2, class)); 
  if(!(length(dataverify) == 1 && dataverify[1] == "integer"||length(dataverify) == 1 && dataverify[1] == "numeric")){
    stop("Please numeric data!")
  }
  dname = deparse(substitute(data))
  result = CTT(data)
  testmethod<-"Corrected t-test"
  estimate = result[1]
  stderr = result[2]
  tstat = result[3]
  pval = result[4]*2
  cint = c(estimate + stderr*qnorm((1-conf)/2),estimate - stderr*qnorm((1-conf)/2))
  attr(cint,"conf.level")<-conf
  names(tstat)<-"t"
  names(estimate)<-"difference for treatment effect"
  rval<-list(statistic = tstat,p.value = pval,conf.int = cint,
             estimate = estimate,stderr = stderr,
             method = testmethod,data.name = dname)
  class(rval)<-"htest"
  return(rval)
}

evalPower<-function(n,cov_num,level_num,pr,type,beta,di = seq(0,0.5,0.1),sigma = 1,Iternum,sl = 0.05,
                    method = c("HuHuCAR","PocSimMIN","StrBCD","StrPBR","DoptBCD","AdjBCD"),
                    test = c("boot.test","corr.test","rand.test"),plot = TRUE,...){
  a = Sys.time()
  if(!is.logical(plot)){
    print("Please specify whether to plot or not!")
    return(NULL)
  }
  else{
    test = match.arg(test)
    mu2 = rep(0,length(di))
    method = match.arg(method)
    if(test == "rand.test"){
      FUN = switch(method,"HuHuCAR" = HuHuCAR_RT_power, "PocSimMIN" = PocSimMIN_RT_power,
                   "StrBCD" = StrBCD_RT_power, "StrPBR" = StrPBR_RT_power,
                   "DoptBCD" = DoptBCD_RT_power, "AdjBCD" = AdjBCD_RT_power)
    }
    else if(test == "boot.test"){
      FUN = switch(method,"HuHuCAR" = HuHuCAR_BT_power, "PocSimMIN" = PocSimMIN_BT_power,
                   "StrBCD" = StrBCD_BT_power, "StrPBR" = StrPBR_BT_power,
                   "DoptBCD" = DoptBCD_BT_power, "AdjBCD" = AdjBCD_BT_power)
    }
    else if(test == "corr.test"){
      FUN = switch(method,"HuHuCAR" = HuHuCAR_CT_power, "PocSimMIN" = PocSimMIN_CT_power,
                   "StrBCD" = StrBCD_CT_power, "StrPBR" = StrPBR_CT_power,
                   "DoptBCD" = DoptBCD_CT_power, "AdjBCD" = AdjBCD_CT_power)
    }
    else{stop("Please enter a valid test! rand.test, boot.test or corr.test")}
    result = FUN(n,cov_num,level_num,pr,type,beta,di,mu2,sigma,Iternum,sl,...)
    if(plot == TRUE){
      diff = di 
      value = result[1:length(di)]
      se = round(result[(length(di)+1):(2*length(di))], digits = 3)
      tgg=data.frame(diff, value, se, stringsAsFactors = TRUE)
      pic = ggplot2::ggplot(tgg, ggplot2::aes(x=di, y=value)) + ggplot2::geom_line() + geom_point(size=4, shape=20)+
        xlab("Difference in means")+ylab("Power")
      b = Sys.time()
      result = list(Powers = tgg,Plot = pic,Time = paste("Execute time:",round(as.numeric(b-a), digits = 2),units(b-a)))
      return(result)
    }
    else{
      diff = di
      value = result[1:length(di)]
      se = round(result[(length(di)+1):(2*length(di))], digits = 3)
      tgg=data.frame(diff, value, se, stringsAsFactors = TRUE)
      b = Sys.time()
      result = list(Powers = tgg,Time = paste("Execute time:",round(as.numeric(b-a), digits = 2),units(b-a)))
      return(result)
    }
  }
}



compPower<-function(powers,diffs,testname){
  if(!is.vector(testname) || !is.character(testname)){
    stop("Input of testname must be a vector of character!")
  }
  if(!is.vector(diffs) || !is.numeric(diffs)){
    stop("Input of testname must be a vector!")
  }
  if(!is.list(powers)){
    stop("Input of powers must be a list!")
  }else if(length(powers) != length(testname)){
    stop("The length of powers must match that of testname!")
  }else{
    k = length(powers)
    l = length(powers[[1]]$Powers$value)
    for(i in 2:k){
      if(length(powers[[k]]$Powers$value)!=l){
        stop("The length of power vectors must match!")
      }
    }
    if(length(diffs)!=l){
      stop("The length of powers and diffs must match!")
    }
    Lines = NULL
    popp = NULL
    popp_out = NULL
    power_temp = rep('',l)
    for(i in 1:k){
      Lines = c(Lines,rep(testname[i],l))
      for(j in 1:l){
        power_temp[j] = paste(powers[[i]]$Powers$value[j], "(", round(powers[[i]]$Powers$se[j],digits = 3),")",sep = '')
      }
      popp_out = c(popp_out,power_temp)
      popp = c(popp,powers[[i]]$Powers$value)
    }
    diffp = rep(diffs,k)
  }
  letit = "Sample Size"
  if(grepl("corr",testname)[1]||grepl("rand",testname)[1]||grepl("boot",testname)[1]||grepl("Simple",testname)[1]){
    letit = "Tests"
  }
  Lines = factor(Lines,levels = testname)
  tgg = data.frame(Lines,diffp,popp, stringsAsFactors = TRUE)
  pic = ggplot2::ggplot(tgg,ggplot2::aes(x = diffp,y = popp,color = Lines,shape = Lines)) + 
    ggplot2::geom_line() + ggplot2::geom_point(size=4) + ggplot2::xlab("Difference in means") + 
    ggplot2::ylab("Power") + ggplot2::scale_colour_hue(name = letit) + 
    ggplot2::scale_shape_discrete(name = letit) + ggplot2::theme(legend.position="bottom")
  tpp = t(matrix(popp_out,nrow = l))
  rownames(tpp) = testname
  tpp = data.frame(tpp, stringsAsFactors = TRUE)
  colnames(tpp) = diffs
  result = list(powers = tpp,plot = pic)
  return(result)
}


