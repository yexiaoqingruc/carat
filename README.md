
# Covariate-Adaptive Randomization for Clinical Trials

Provides functions and command-line user interface to generate allocation sequence by covariate-adaptive randomization for clinical trials. It currently supports six covariate-adaptive randomization procedures. Three hypothesis testing methods that are valid and robust under covariate-adaptive randomization are also available in the package to facilitate the inference for treatment effect under the included randomization procedures. Additionally, the package provides comprehensive and efficient tools to allow one to evaluate and compare the performance of randomization procedures and tests based on various criteria.

# Installation 

## Install OpenMP-supported version from GitHub: 

The OpenMP-supported version of carat, named caratOMP, can be installed from GitHub:

```{R, eval = FALSE}
#install.packages('devtools')
devtools::install_github('yexiaoqingruc/carat/caratOMP')
```

and the package without **OpenMP**, names carat, can be installed from Github: 

```{R, eval = FALSE}
#install.packages('devtools')
devtools::install_github('yexiaoqingruc/carat/carat')
```

## Install the current release from CRAN:

The package without **OpenMP** is available for Windows, Linux and Mac OS X from CRAN: 

```{R, eval = FALSE}
install.packages("carat")
```


