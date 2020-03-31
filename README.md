# Overview

## Covariate-Adaptive Randomization for Clinical Trials supporting OpenMP

Provides functions and command-line user interface to generate allocation sequence by covariate-adaptive randomization for clinical trials. It currently supports six covariate-adaptive randomization procedures. Three hypothesis testing methods that are valid and robust under covariate-adaptive randomization are also available in the package to facilitate the inference for treatment effect under the included randomization procedures. Additionally, the package provides comprehensive and efficient tools to allow one to evaluate and compare the performance of randomization procedures and tests based on various criteria.

## Note

Due to starting with 4.0 R no longer supports **OpenMP**, which results from the limitation of system compiler on MacOS. It is removed from the current version of **carat** on CRAN. However, empirical studies indicate that **OpenMP**, which was used in our early version, is a more efficient parallel computing tool than **doSNOW**. Thus, we also provide a parallel version of **carat** using **OpenMP** via GitHub repository. 

# Installation 

## Install the current release from CRAN:

The package without **OpenMP** is available for Windows, Linux and Mac OS X from CRAN: 

```{R, eval = FALSE}
install.packages("carat")
```
## Install OpenMP-supported version from GitHub: 

The parallel version of **carat** implemented by **OpenMP** named **caratOMP** is available via:

```{R, eval = FALSE}
#install.packages('devtools')
devtools::install_github('yexiaoqingruc/carat/caratOMP')
```

