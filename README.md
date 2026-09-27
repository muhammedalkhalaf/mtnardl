# mtnardl

**Bootstrap Multiple Threshold Nonlinear ARDL**

The `mtnardl` R package implements the Multiple Threshold Nonlinear ARDL (MTNARDL) model following Pal and Mitra (2016). It decomposes regressors into regime-specific partial sums based on quantile or custom partitions and tests for cointegration using PSS bounds testing with optional bootstrap critical values.

## Features

- Quintile, quartile, tercile, and binary partitions
- Automatic lag selection via AIC or BIC
- PSS bounds test with Kripfganz & Schneider (2020) asymptotic CVs
- Bootstrap cointegration test (McNown et al., 2018)
- Long-run coefficients via delta method
- Dynamic multipliers per regime

## Installation

```r
install.packages("mtnardl")
```

## Usage

```r
library(mtnardl)

set.seed(42)
TT <- 80
y  <- cumsum(rnorm(TT))
x  <- cumsum(rnorm(TT))

res <- mtnardl(y = y, x = matrix(x, ncol = 1),
               decompose = 1L, partition = "quintile")
summary(res)
```

## References

- Pal, D., & Mitra, S. K. (2016). *Economic Modelling*, 59, 314-328. https://doi.org/10.1016/j.econmod.2016.08.003
- Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). *Journal of Applied Econometrics*, 16(3), 289-326. https://doi.org/10.1002/jae.616
- McNown, R., Sam, C. Y., & Goh, S. K. (2018). *Applied Economics*, 50(13), 1509-1521. https://doi.org/10.1080/00036846.2017.1366643

## License

GPL-3
