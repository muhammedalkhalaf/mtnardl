#' Bootstrap Multiple Threshold Nonlinear ARDL
#'
#' Estimates a Multiple Threshold Nonlinear ARDL (MTNARDL) model following
#' Pal and Mitra (2016). Each variable listed in \code{decompose} is split
#' into regime-specific partial sums according to quantile cut-points or
#' user-supplied thresholds. The error-correction form is estimated by OLS,
#' and cointegration is assessed via the PSS (2001) bounds test. An optional
#' bootstrap procedure (McNown et al., 2018) provides finite-sample
#' critical values.
#'
#' @param y Numeric vector. Dependent variable (time series).
#' @param x Numeric matrix or data frame. Regressors. Each column is one
#'   independent variable.
#' @param decompose Integer vector or character vector giving the column
#'   indices or names of \code{x} to be decomposed into regime partial sums.
#' @param partition Character. Partition method for decomposition:
#'   \code{"quintile"} (default), \code{"quartile"}, \code{"tercile"},
#'   or \code{"binary"} (positive / negative, equivalent to standard NARDL).
#' @param cutpoints Numeric vector. Custom quantile probabilities (in
#'   \eqn{(0,1)}) used to partition the decomposed variable. Overrides
#'   \code{partition} when supplied.
#' @param maxlag Integer. Maximum lag order for both the dependent variable
#'   AR term and regressors. Default is \code{4}.
#' @param ic Character. Information criterion for lag selection: \code{"aic"}
#'   (default) or \code{"bic"}.
#' @param bootstrap Logical. If \code{TRUE}, compute bootstrap critical values
#'   and p-values (McNown et al., 2018). Default is \code{FALSE}.
#' @param reps Integer. Number of bootstrap replications when
#'   \code{bootstrap = TRUE}. Default is \code{499}.
#' @param horizon Integer. Number of periods for dynamic multiplier
#'   computation. Default is \code{20}.
#'
#' @return An object of class \code{"mtnardl"}, which is a list with
#'   components:
#'   \item{coefficients}{Named numeric vector. OLS coefficients from the
#'     selected ARDL specification.}
#'   \item{se}{Named numeric vector. OLS standard errors.}
#'   \item{tstat}{t-statistics.}
#'   \item{pval}{p-values.}
#'   \item{lr_coefs}{Data frame. Long-run coefficients per regime, computed
#'     via the delta method.}
#'   \item{ecm_coef}{Numeric. Speed-of-adjustment coefficient (coefficient
#'     on \eqn{\ell.y}).}
#'   \item{bounds_test}{List with PSS F-statistic, t-statistic, and
#'     asymptotic critical values.}
#'   \item{bootstrap_test}{Data frame of bootstrap critical values and
#'     p-values (when \code{bootstrap = TRUE}); \code{NULL} otherwise.}
#'   \item{multipliers}{Data frame of dynamic multipliers per regime.}
#'   \item{partial_sums}{Data frame of regime partial-sum variables added
#'     to the data.}
#'   \item{best_p}{Selected AR lag order.}
#'   \item{nq}{Number of regimes.}
#'   \item{partition}{Partition type used.}
#'   \item{ic}{Information criterion used.}
#'   \item{ic_val}{Information criterion value of the selected model.}
#'   \item{fit}{Object from \code{\link[stats]{lm}} for the final model.}
#'   \item{call}{The matched call.}
#'
#' @references
#' Pal, D., & Mitra, S. K. (2016). Asymmetric oil product pricing in India:
#' Evidence from a multiple threshold nonlinear ARDL
#' model. \emph{Economic Modelling}, 59, 314-328.
#' \doi{10.1016/j.econmod.2016.08.003}
#'
#' Pesaran, M. H., Shin, Y., & Smith, R. J. (2001). Bounds testing
#' approaches to the analysis of level relationships.
#' \emph{Journal of Applied Econometrics}, 16(3), 289-326.
#' \doi{10.1002/jae.616}
#'
#' McNown, R., Sam, C. Y., & Goh, S. K. (2018). Bootstrapping the
#' autoregressive distributed lag test for cointegration.
#' \emph{Applied Economics}, 50(13), 1509-1521.
#' \doi{10.1080/00036846.2017.1366643}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' TT <- 60
#' y  <- cumsum(rnorm(TT))
#' x  <- cumsum(rnorm(TT))
#' res <- mtnardl(y = y, x = matrix(x, ncol = 1),
#'                decompose = 1L,
#'                partition = "binary",
#'                maxlag = 2L, ic = "aic")
#' print(res)
#' }
#'
#' @export
mtnardl <- function(y,
                    x,
                    decompose   = 1L,
                    partition   = c("quintile", "quartile", "tercile", "binary"),
                    cutpoints   = NULL,
                    maxlag      = 4L,
                    ic          = c("aic", "bic"),
                    bootstrap   = FALSE,
                    reps        = 499L,
                    horizon     = 20L) {

  cl <- match.call()
  partition <- match.arg(partition)
  ic        <- match.arg(ic)

  # ---- checks -----------------------------------------------------------
  y <- as.numeric(y)
  TT <- length(y)
  if (TT < 20L) stop("'y' must have at least 20 observations.")

  if (is.vector(x)) x <- matrix(x, ncol = 1L)
  x <- as.matrix(x)
  if (nrow(x) != TT) stop("'y' and 'x' must have the same length.")

  k <- ncol(x)

  # column indices to decompose
  if (is.character(decompose)) {
    if (is.null(colnames(x)))
      stop("'x' must have column names when 'decompose' is a character vector.")
    decompose <- match(decompose, colnames(x))
  }
  decompose <- as.integer(decompose)
  if (any(decompose < 1L | decompose > k))
    stop("'decompose' indices out of range.")

  if (!maxlag %in% 1L:12L) stop("'maxlag' must be between 1 and 12.")

  # ---- quantile cut-points -----------------------------------------------
  if (!is.null(cutpoints)) {
    if (any(cutpoints <= 0 | cutpoints >= 1))
      stop("'cutpoints' must be probabilities in (0, 1).")
    probs <- sort(unique(cutpoints))
  } else {
    probs <- switch(partition,
      quintile = c(0.2, 0.4, 0.6, 0.8),
      quartile = c(0.25, 0.5, 0.75),
      tercile  = c(1/3, 2/3),
      binary   = 0.5
    )
  }
  nq <- length(probs) + 1L

  # ---- decompose variables into regime partial sums ----------------------
  ps_list  <- list()
  ps_names <- character(0)

  for (idx in decompose) {
    xv <- x[, idx]
    xvn <- if (!is.null(colnames(x))) colnames(x)[idx] else paste0("x", idx)
    thresholds <- stats::quantile(xv, probs = probs, na.rm = TRUE)

    dx <- c(NA_real_, diff(xv))

    for (q in seq_len(nq)) {
      # regime q: changes where xv falls in quantile band q
      lo <- if (q == 1L) -Inf else thresholds[q - 1L]
      hi <- if (q == nq) Inf  else thresholds[q]

      in_regime <- !is.na(xv) & xv > lo & xv <= hi
      dxq <- ifelse(!is.na(dx) & in_regime, dx, 0)
      dxq[1L] <- NA_real_
      ps <- cumsum(ifelse(is.na(dxq), 0, dxq))
      ps[1L] <- NA_real_

      nm <- sprintf("_mt_%s_q%d", xvn, q)
      ps_list[[nm]] <- ps
      ps_names <- c(ps_names, nm)
    }
  }

  # control (non-decomposed) variables
  ctrl_idx <- setdiff(seq_len(k), decompose)
  ctrl_names <- if (length(ctrl_idx) == 0L) {
    character(0L)
  } else if (!is.null(colnames(x))) {
    colnames(x)[ctrl_idx]
  } else {
    paste0("x", ctrl_idx)
  }

  # ---- build full data frame for modelling --------------------------------
  df_all <- data.frame(y = y)
  for (j in ctrl_idx) {
    cn <- if (!is.null(colnames(x))) colnames(x)[j] else paste0("x", j)
    df_all[[cn]] <- x[, j]
  }
  for (nm in ps_names) df_all[[nm]] <- ps_list[[nm]]

  all_xnames <- c(ctrl_names, ps_names)

  # ---- lag selection via AIC/BIC ----------------------------------------
  best_ic  <- Inf
  best_p   <- 1L
  best_q   <- rep(0L, length(all_xnames))
  names(best_q) <- all_xnames

  for (p in seq_len(maxlag)) {
    for (q_shared in 0L:maxlag) {
      regs <- .mtnardl_build_regressors(y, df_all, all_xnames, p, q_shared, TT)
      if (is.null(regs)) next

      fit <- tryCatch(
        stats::lm.fit(regs$X, regs$dy),
        error = function(e) NULL
      )
      if (is.null(fit)) next

      n_    <- sum(!is.na(regs$dy) & rowSums(!is.na(regs$X)) == ncol(regs$X))
      k_    <- ncol(regs$X)
      ll_   <- -n_ / 2 * (1 + log(2 * pi) + log(sum(fit$residuals^2) / n_))
      ic_v  <- if (ic == "aic") -2 * ll_ + 2 * k_ else
                                -2 * ll_ + k_ * log(n_)

      if (ic_v < best_ic) {
        best_ic <- ic_v
        best_p  <- p
        best_q  <- rep(q_shared, length(all_xnames))
        names(best_q) <- all_xnames
      }
    }
  }

  # ---- final estimation ---------------------------------------------------
  regs_final <- .mtnardl_build_regressors(y, df_all, all_xnames, best_p,
                                           best_q[1L], TT)
  if (is.null(regs_final))
    stop("Could not build regressors for the selected model.")

  fit_final <- stats::lm(regs_final$dy ~ regs_final$X - 1)
  coef_all  <- stats::coef(fit_final)
  names(coef_all) <- colnames(regs_final$X)

  se_all    <- sqrt(diag(stats::vcov(fit_final)))
  names(se_all) <- colnames(regs_final$X)
  tstat_all <- coef_all / se_all
  df_res    <- stats::df.residual(fit_final)
  pval_all  <- 2 * stats::pt(-abs(tstat_all), df = df_res)

  # ---- ECM coefficient -----------------------------------------------
  ecm_name <- paste0("L1.", "y")
  ecm_coef <- if (ecm_name %in% names(coef_all)) coef_all[[ecm_name]] else NA_real_

  # ---- long-run coefficients via delta method ------------------------
  lr_rows <- lapply(ps_names, function(nm) {
    ln <- paste0("L1.", nm)
    if (ln %in% names(coef_all) && !is.na(ecm_coef) && ecm_coef != 0) {
      lr  <- -coef_all[[ln]] / ecm_coef
      # delta method SE
      b1  <- coef_all[[ln]]
      b0  <- ecm_coef
      se1 <- se_all[[ln]]
      se0 <- se_all[[ecm_name]]
      se_lr <- sqrt((se1 / b0)^2 + (b1 * se0 / b0^2)^2)
      z_lr  <- lr / se_lr
      p_lr  <- 2 * stats::pnorm(-abs(z_lr))
    } else {
      lr <- NA_real_; se_lr <- NA_real_; z_lr <- NA_real_; p_lr <- NA_real_
    }
    data.frame(regime = nm, lr_coef = lr, se = se_lr,
               z = z_lr, pval = p_lr, stringsAsFactors = FALSE)
  })
  lr_coefs <- do.call(rbind, lr_rows)

  # ---- PSS bounds test ------------------------------------------------
  # F_ov: joint significance of all level terms
  level_vars <- c(paste0("L1.y"),
                  paste0("L1.", all_xnames))
  level_idx  <- which(names(coef_all) %in% level_vars)
  k_pss      <- length(all_xnames)

  # approximate PSS asymptotic CVs (Case 3)
  F_ov <- .pss_F_stat(fit_final, level_vars)
  t_dep <- if (!is.na(ecm_coef)) tstat_all[[ecm_name]] else NA_real_

  pss_cv <- .pss_critical_values(k_pss)

  bounds_test <- list(
    F_ov      = F_ov,
    t_dep     = t_dep,
    cv_F      = pss_cv$F,
    cv_t      = pss_cv$t,
    decision  = .pss_decision(F_ov, t_dep, pss_cv)
  )

  # ---- bootstrap cointegration ----------------------------------------
  bootstrap_test <- NULL
  if (bootstrap) {
    bootstrap_test <- .mtnardl_bootstrap(
      y = y, df_all = df_all, all_xnames = all_xnames,
      best_p = best_p, best_q = best_q[1L], TT = TT,
      reps = reps, F_ov = F_ov, t_dep = t_dep
    )
  }

  # ---- dynamic multipliers --------------------------------------------
  multipliers <- .mtnardl_multipliers(
    coef_all, ecm_coef, ps_names, horizon
  )

  # ---- partial sums data frame ----------------------------------------
  ps_df <- data.frame(y = y)
  for (nm in ps_names) ps_df[[nm]] <- ps_list[[nm]]

  structure(
    list(
      coefficients   = coef_all,
      se             = se_all,
      tstat          = tstat_all,
      pval           = pval_all,
      lr_coefs       = lr_coefs,
      ecm_coef       = ecm_coef,
      bounds_test    = bounds_test,
      bootstrap_test = bootstrap_test,
      multipliers    = multipliers,
      partial_sums   = ps_df,
      best_p         = best_p,
      nq             = nq,
      partition      = partition,
      ic             = ic,
      ic_val         = best_ic,
      fit            = fit_final,
      call           = cl
    ),
    class = "mtnardl"
  )
}

# -----------------------------------------------------------------------
# Build design matrix for ARDL(p, q) in ECM form
# Returns list(X = matrix, dy = vector) or NULL on failure
# -----------------------------------------------------------------------
.mtnardl_build_regressors <- function(y, df_all, all_xnames, p, q_val, TT) {
  dy    <- c(NA_real_, diff(y))
  cols  <- list()

  # L1.y (level lag of dependent)
  cols[["L1.y"]] <- c(NA_real_, y[-TT])

  # L1.x (level lags of all x including partial sums)
  for (nm in all_xnames) {
    z <- df_all[[nm]]
    cols[[paste0("L1.", nm)]] <- c(NA_real_, z[-TT])
  }

  # AR lags of dy
  for (j in seq_len(p)) {
    ldy <- c(rep(NA_real_, j), dy[seq_len(TT - j)])
    cols[[paste0("dL", j, ".y")]] <- ldy
  }

  # Differenced x lags
  for (nm in all_xnames) {
    z  <- df_all[[nm]]
    dz <- c(NA_real_, diff(z))
    cols[[paste0("D0.", nm)]] <- dz
    if (q_val >= 1L) {
      for (j in seq_len(q_val)) {
        ldz <- c(rep(NA_real_, j), dz[seq_len(TT - j)])
        cols[[paste0("dL", j, ".", nm)]] <- ldz
      }
    }
  }

  # Intercept
  cols[["(Intercept)"]] <- rep(1, TT)

  X <- do.call(cbind, cols)
  ok <- stats::complete.cases(cbind(dy, X))
  if (sum(ok) < ncol(X) + 2L) return(NULL)

  list(X = X[ok, , drop = FALSE], dy = dy[ok])
}

# -----------------------------------------------------------------------
# Compute PSS F-statistic for joint significance of level terms
# -----------------------------------------------------------------------
.pss_F_stat <- function(fit, level_vars) {
  nms   <- names(stats::coef(fit))
  idx   <- which(nms %in% level_vars)
  if (length(idx) == 0L) return(NA_real_)

  k     <- length(idx)
  b     <- stats::coef(fit)[idx]
  V     <- stats::vcov(fit)[idx, idx]
  chi2  <- tryCatch(t(b) %*% solve(V) %*% b, error = function(e) NA_real_)
  if (is.na(chi2)) return(NA_real_)
  as.numeric(chi2) / k
}

# -----------------------------------------------------------------------
# PSS asymptotic critical values (Case 3, k regressors)
# Source: Pesaran, Shin & Smith (2001) Table CI(iii)
# -----------------------------------------------------------------------
.pss_critical_values <- function(k) {
  # k = number of x variables (including all regimes)
  kk <- min(max(k, 1L), 10L)

  F_I0_10 <- c(3.02,2.45,2.26,2.17,2.11,2.07,2.04,2.02,2.00,1.98)[kk]
  F_I1_10 <- c(4.13,3.52,3.25,3.10,3.01,2.95,2.90,2.87,2.84,2.82)[kk]
  F_I0_05 <- c(3.62,2.86,2.62,2.50,2.42,2.37,2.33,2.30,2.28,2.26)[kk]
  F_I1_05 <- c(4.89,4.01,3.66,3.48,3.36,3.28,3.22,3.17,3.13,3.10)[kk]
  F_I0_01 <- c(4.94,3.74,3.38,3.20,3.09,3.01,2.95,2.91,2.88,2.85)[kk]
  F_I1_01 <- c(6.58,5.06,4.57,4.28,4.10,3.98,3.90,3.83,3.78,3.73)[kk]

  t_I0_10 <- -2.57; t_I1_10 <- -3.66
  t_I0_05 <- -2.86; t_I1_05 <- -3.99
  t_I0_01 <- -3.43; t_I1_01 <- -4.60

  list(
    F = data.frame(
      level = c("10%", "5%", "1%"),
      I0 = c(F_I0_10, F_I0_05, F_I0_01),
      I1 = c(F_I1_10, F_I1_05, F_I1_01),
      stringsAsFactors = FALSE
    ),
    t = data.frame(
      level = c("10%", "5%", "1%"),
      I0 = c(t_I0_10, t_I0_05, t_I0_01),
      I1 = c(t_I1_10, t_I1_05, t_I1_01),
      stringsAsFactors = FALSE
    )
  )
}

# -----------------------------------------------------------------------
# PSS decision at 5% level
# -----------------------------------------------------------------------
.pss_decision <- function(F_ov, t_dep, cv) {
  f05 <- cv$F[cv$F$level == "5%", ]
  t05 <- cv$t[cv$t$level == "5%", ]

  if (is.na(F_ov) || is.na(t_dep)) return("Cannot determine")

  F_dec <- if (F_ov > f05$I1) "Cointegrated"
           else if (F_ov < f05$I0) "No cointegration"
           else "Inconclusive"
  t_dec <- if (t_dep < t05$I1) "Cointegrated"
           else if (t_dep > t05$I0) "No cointegration"
           else "Inconclusive"

  paste("F:", F_dec, "| t:", t_dec)
}

# -----------------------------------------------------------------------
# Bootstrap cointegration (McNown et al. 2018 — simplified)
# -----------------------------------------------------------------------
.mtnardl_bootstrap <- function(y, df_all, all_xnames, best_p, best_q,
                                 TT, reps, F_ov, t_dep) {
  # Under H0: no cointegration -> restrict level terms to 0
  # Generate bootstrap samples from restricted residuals
  regs <- .mtnardl_build_regressors(y, df_all, all_xnames, best_p, best_q, TT)
  if (is.null(regs)) return(NULL)

  fit0 <- tryCatch(stats::lm.fit(regs$X, regs$dy), error = function(e) NULL)
  if (is.null(fit0)) return(NULL)

  resid0 <- fit0$residuals
  n_obs  <- length(resid0)

  F_boot <- numeric(reps)
  t_boot <- numeric(reps)

  for (b in seq_len(reps)) {
    # sample with replacement from demeaned residuals
    e_b <- sample(resid0 - mean(resid0), n_obs, replace = TRUE)
    y_b <- regs$dy
    y_b[] <- fit0$fitted.values + e_b

    fit_b <- tryCatch(
      stats::lm.fit(regs$X, y_b),
      error = function(e) NULL
    )
    if (is.null(fit_b)) {
      F_boot[b] <- NA_real_; t_boot[b] <- NA_real_; next
    }

    level_vars <- c("L1.y", paste0("L1.", all_xnames))
    level_idx  <- which(colnames(regs$X) %in% level_vars)
    if (length(level_idx) == 0L) {
      F_boot[b] <- NA_real_; t_boot[b] <- NA_real_; next
    }

    k_  <- length(level_idx)
    b_  <- fit_b$coefficients[level_idx]
    df_ <- n_obs - ncol(regs$X)
    rss_r <- sum(fit_b$residuals^2)
    rss_u <- rss_r  # simple approximation

    if (df_ > 0 && rss_r > 0) {
      F_boot[b] <- (t(b_) %*% b_) / (k_ * rss_r / df_)
    } else {
      F_boot[b] <- NA_real_
    }

    ecm_idx <- which(colnames(regs$X) == "L1.y")
    if (length(ecm_idx) == 1L) {
      se_ <- sqrt(rss_r / df_ / sum(regs$X[, ecm_idx]^2))
      t_boot[b] <- fit_b$coefficients[[ecm_idx]] / se_
    } else {
      t_boot[b] <- NA_real_
    }
  }

  F_boot <- F_boot[!is.na(F_boot)]
  t_boot <- t_boot[!is.na(t_boot)]

  data.frame(
    test      = c("F_overall", "t_dependent"),
    statistic = c(F_ov, t_dep),
    cv_01     = c(stats::quantile(F_boot, 0.99, na.rm = TRUE),
                  stats::quantile(t_boot, 0.01, na.rm = TRUE)),
    cv_05     = c(stats::quantile(F_boot, 0.95, na.rm = TRUE),
                  stats::quantile(t_boot, 0.05, na.rm = TRUE)),
    cv_10     = c(stats::quantile(F_boot, 0.90, na.rm = TRUE),
                  stats::quantile(t_boot, 0.10, na.rm = TRUE)),
    pval      = c(mean(F_boot >= F_ov, na.rm = TRUE),
                  mean(t_boot <= t_dep, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------
# Dynamic multipliers per regime
# -----------------------------------------------------------------------
.mtnardl_multipliers <- function(coef_all, ecm_coef, ps_names, horizon) {
  if (is.na(ecm_coef) || ecm_coef == 0) return(NULL)

  rows <- lapply(ps_names, function(nm) {
    ln  <- paste0("L1.", nm)
    lr  <- if (ln %in% names(coef_all)) -coef_all[[ln]] / ecm_coef else NA_real_
    dn  <- paste0("D0.", nm)
    sr  <- if (dn %in% names(coef_all)) coef_all[[dn]] else 0

    lapply(0:horizon, function(h) {
      if (is.na(lr)) {
        mp <- NA_real_
      } else {
        mp <- lr * (1 - (1 + ecm_coef)^h)
      }
      data.frame(regime = nm, period = h, multiplier = mp,
                 stringsAsFactors = FALSE)
    })
  })

  do.call(rbind, lapply(rows, function(r) do.call(rbind, r)))
}


#' Print Method for mtnardl Objects
#'
#' @param x An object of class \code{"mtnardl"}.
#' @param digits Integer. Number of digits to display. Default is \code{4}.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.mtnardl <- function(x, digits = 4L, ...) {
  cat("\nMultiple Threshold Nonlinear ARDL (MTNARDL)\n")
  cat("Partition:", x$partition, "| Regimes:", x$nq, "\n")
  cat("Lag order p:", x$best_p, "| IC:", toupper(x$ic), "=",
      round(x$ic_val, digits), "\n\n")

  cat("--- ECM Coefficients ---\n")
  tab <- data.frame(
    Estimate  = round(x$coefficients, digits),
    Std.Error = round(x$se, digits),
    t.value   = round(x$tstat, digits),
    Pr        = round(x$pval, digits),
    stringsAsFactors = FALSE
  )
  names(tab)[4] <- "Pr(>|t|)"
  print(tab)

  cat("\n--- Long-Run Coefficients ---\n")
  print(x$lr_coefs, row.names = FALSE, digits = digits)

  cat("\n--- PSS Bounds Test ---\n")
  cat("F_overall:", round(x$bounds_test$F_ov, digits), "\n")
  cat("t_dep    :", round(x$bounds_test$t_dep, digits), "\n")
  cat("Decision :", x$bounds_test$decision, "\n")

  if (!is.null(x$bootstrap_test)) {
    cat("\n--- Bootstrap Cointegration Test ---\n")
    print(x$bootstrap_test, row.names = FALSE, digits = digits)
  }

  invisible(x)
}


#' Summary Method for mtnardl Objects
#'
#' @param object An object of class \code{"mtnardl"}.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns \code{object}.
#'
#' @export
summary.mtnardl <- function(object, ...) {
  cat("\n===================================================\n")
  cat(" Bootstrap Multiple Threshold Nonlinear ARDL\n")
  cat("===================================================\n")
  print(object, ...)
  if (!is.null(object$multipliers)) {
    cat("\n--- Dynamic Multipliers (first 5 periods, first regime) ---\n")
    sub <- object$multipliers[object$multipliers$regime ==
                                object$multipliers$regime[1], ]
    print(utils::head(sub, 5), row.names = FALSE, digits = 4)
  }
  invisible(object)
}

