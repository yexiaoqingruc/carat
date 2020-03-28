# Covariate-Adaptive Randomization for Clinical Trials supporting OpenMP

Provides functions and command-line user interface to generate allocation sequence by covariate-adaptive randomization for clinical trials. It currently supports six covariate-adaptive randomization procedures. Three hypothesis testing methods that are valid and robust under covariate-adaptive randomization are also available in the package to facilitate the inference for treatment effect under the included randomization procedures. Additionally, the package provides comprehensive and efficient tools to allow one to evaluate and compare the performance of randomization procedures and tests based on various criteria.

## Install the current release from CRAN:

The package without OpenMP is available for Windows, Linux and Mac OS X from CRAN: 

```{R, eval = FALSE}
install.packages("carat")
```
## Install OpenMP-supported version from GitHub: 

Note that, starting with 4.0 R no longer support OpenMP because it is platform-dependent. It is not supported by the system compiler on macOS. Empirical studies however indicate that OpenMP makes our algorithms in the package more efficient. Hence, we provide OpenMP-supported version on GitHub and it can be loaded via: 

```{R, eval = FALSE}
devtools::install_github('yexiaoqingruc/caratOMP/caratOMP')
```

