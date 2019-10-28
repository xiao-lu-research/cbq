# cbq: An R Package for Conditional Binary Quantile Models


The package `cbq` provides basic functionalities of conditional binary quantile models using MCMC methods. The estimation is conducted through pre-compiled stan codes.

Caution: The package is still under initial development and there is no guarantee about the functionality of the package.

## Installation

```r
# Make sure that the following packages have been installed in your local R environment
if(!require(rstan)) install.packages("rstan")

# Install cirque from github
if(!require(devtools)) install.packages("devtools")
devtools::install_github("xiao-lu-research/cbq")
```


## Usage

```r

# Load the package
library(cbq)

# Get help
?cbq


```

## References

Lu, X. (forthcoming) Discrete Choice Data with Unobserved Heterogeneity: A Conditional Binary Quantile Model. Political Analysis.
