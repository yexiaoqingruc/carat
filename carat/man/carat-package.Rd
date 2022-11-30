\name{carat-package}
\alias{carat-package}
\alias{carat}
\docType{package}
\title{
  carat-package: Covariate-Adaptive Randomization for Clinical Trials
}
\description{
  Provides functions and a command-line user interface to generate allocation sequences for clinical trials with covariate-adaptive randomization methods. It currently supports six different covariate-adaptive randomization procedures, including stratified randomization, minimization, and a general family of designs proposed by Hu and Hu (2012) <doi:10.1214/12-AOS983>. Three hypothesis testing methods, all valid and robust under covariate-adaptive randomization are also included in the package to facilitate the inference for treatment effects under the included randomization procedures. Additionally, the package provides comprehensive and efficient tools for the performance evaluation and comparison of randomization procedures and tests based on various criteria.
}

\section{Acknowledgement}{
This work was supported by the Fundamental Research Funds for the Central Universities, and the Research Funds of Renmin University of China [grant number 20XNA023].
}
\references{
Atkinson A C. \emph{Optimum biased coin designs for sequential clinical trials with prognostic factors}[J]. Biometrika, 1982, 69(1): 61-67. <doi:10.2307/2335853>

Baldi Antognini A, Zagoraiou M. \emph{The covariate-adaptive biased coin design for balancing clinical trials in the presence of prognostic factors}[J]. Biometrika, 2011, 98(3): 519-535. <doi:10.1093/biomet/asr021>

Hu Y, Hu F. \emph{Asymptotic properties of covariate-adaptive randomization}[J]. The Annals of Statistics, 2012, 40(3): 1794-1815. <doi:10.1214/12-AOS983>

Ma W, Hu F, Zhang L. \emph{Testing hypotheses of covariate-adaptive randomized clinical trials}[J]. Journal of the American Statistical Association, 2015, 110(510): 669-680. <doi:10.1080/01621459.2014.922469>

Ma W, Qin Y, Li Y, et al. \emph{Statistical Inference for Covariate-Adaptive Randomization Procedures}[J]. Journal of the American Statistical Association, 2020, 115(531): 1488-1597. <doi:10.1080/01621459.2019.1635483>

Pocock S J, Simon R. \emph{Sequential treatment assignment with balancing for prognostic factors in the controlled clinical trial}[J]. Biometrics, 1975: 103-115. <doi:10.2307/2529712>

Rosenberger W F, Lachin J M. \emph{Randomization in clinical trials: theory and practice}[M]. John Wiley & Sons, 2015. <doi:10.1002/9781118742112>

Shao J., Yu, X. \emph{Validity of tests under covariate-adaptive biased coin randomization and generalized linear models}[J]. Biometrics, 2013, 69(4), 960-969. <doi:10.1111/biom.12062>

Shao J, Yu X, Zhong B. \emph{A theory for testing hypotheses under covariate-adaptive randomization}[J]. Biometrika, 2010, 97(2): 347-360. <doi:10.1093/biomet/asq014>

Zelen M. \emph{The randomization and stratification of patients to clinical trials}[J]. Journal of chronic diseases, 1974, 27(7): 365-375. <doi:10.1016/0021-9681(74)90015-0>
}
\author{
Fuyi Tu \email{fuyi.tu@ruc.edu.cn};Xiaoqing Ye \email{ye_xiaoq@163.com}; 
Wei Ma \email{mawei@ruc.edu.cn};  Feifang Hu \email{feifang@gwu.edu}.
}
\keyword{ carat-package }