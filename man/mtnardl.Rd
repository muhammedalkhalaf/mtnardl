\name{mtnardl}
\alias{mtnardl}
\title{Bootstrap Multiple Threshold Nonlinear ARDL}
\usage{
mtnardl(
  y,
  x,
  decompose = 1L,
  partition = c("quintile", "quartile", "tercile", "binary"),
  cutpoints = NULL,
  maxlag = 4L,
  ic = c("aic", "bic"),
  bootstrap = FALSE,
  reps = 499L,
  horizon = 20L
)
}
\arguments{
\item{y}{Numeric vector. Dependent variable (time series).}

\item{x}{Numeric matrix or data frame. Regressors. Each column is one
independent variable.}

\item{decompose}{Integer vector or character vector giving the column
indices or names of \code{x} to be decomposed into regime partial sums.}

\item{partition}{Character. Partition method for decomposition:
\code{"quintile"} (default), \code{"quartile"}, \code{"tercile"},
or \code{"binary"} (positive / negative, equivalent to standard NARDL).}

\item{cutpoints}{Numeric vector. Custom quantile probabilities (in
\eqn{(0,1)}) used to partition the decomposed variable. Overrides
\code{partition} when supplied.}

\item{maxlag}{Integer. Maximum lag order for both the dependent variable
AR term and regressors. Default is \code{4}.}

\item{ic}{Character. Information criterion for lag selection: \code{"aic"}
(default) or \code{"bic"}.}

\item{bootstrap}{Logical. If \code{TRUE}, compute bootstrap critical values
and p-values (McNown et al., 2018). Default is \code{FALSE}.}

\item{reps}{Integer. Number of bootstrap replications when
\code{bootstrap = TRUE}. Default is \code{499}.}

\item{horizon}{Integer. Number of periods for dynamic multiplier
computation. Default is \code{20}.}
}
\value{
An object of class \code{"mtnardl"}, which is a list with
  components:
  \item{coefficients}{Named numeric vector. OLS coefficients from the
    selected ARDL specification.}
  \item{se}{Named numeric vector. OLS standard errors.}
  \item{tstat}{t-statistics.}
  \item{pval}{p-values.}
  \item{lr_coefs}{Data frame. Long-run coefficients per regime, computed
    via the delta method.}
  \item{ecm_coef}{Numeric. Speed-of-adjustment coefficient.}
  \item{bounds_test}{List with PSS F-statistic, t-statistic, and
    asymptotic critical values.}
  \item{bootstrap_test}{Data frame of bootstrap critical values and
    p-values (when \code{bootstrap = TRUE}); \code{NULL} otherwise.}
  \item{multipliers}{Data frame of dynamic multipliers per regime.}
  \item{partial_sums}{Data frame of regime partial-sum variables.}
  \item{best_p}{Selected AR lag order.}
  \item{nq}{Number of regimes.}
  \item{partition}{Partition type used.}
  \item{ic}{Information criterion used.}
  \item{ic_val}{Information criterion value of the selected model.}
  \item{fit}{Object from \code{\link[stats]{lm}} for the final model.}
  \item{call}{The matched call.}
}
\description{
Estimates a Multiple Threshold Nonlinear ARDL (MTNARDL) model following
Pal and Mitra (2016). Each variable listed in \code{decompose} is split
into regime-specific partial sums according to quantile cut-points or
user-supplied thresholds. The error-correction form is estimated by OLS,
and cointegration is assessed via the PSS (2001) bounds test.
}
\examples{
\donttest{
set.seed(42)
TT <- 60
y  <- cumsum(rnorm(TT))
x  <- cumsum(rnorm(TT))
res <- mtnardl(y = y, x = matrix(x, ncol = 1),
               decompose = 1L,
               partition = "binary",
               maxlag = 2L, ic = "aic")
print(res)
}
}
\references{
Pal, D., & Mitra, S. K. (2016). Asymmetric oil product pricing in India:
Evidence from a multiple threshold nonlinear ARDL
model. \emph{Economic Modelling}, 59, 314-328.
\doi{10.1016/j.econmod.2016.08.003}

Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing
approaches to the analysis of level relationships.
\emph{Journal of Applied Econometrics}, 16(3), 289-326.
\doi{10.1002/jae.616}

McNown, R., Sam, C. Y., & Goh, S. K. (2018). Bootstrapping the
autoregressive distributed lag test for cointegration.
\emph{Applied Economics}, 50(13), 1509-1521.
\doi{10.1080/00036846.2017.1366643}
}
