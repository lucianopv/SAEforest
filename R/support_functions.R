# Small functions that support functionality for point and MSE estimates.
#' @importFrom stats aggregate as.formula complete.cases dist formula median na.omit optimize predict quantile rnorm sd


# Calculation of indicators from smearing or MC pseudo-populations ------------------------
calc_indicat <- function(Y, threshold, custom, select.indicator = NULL) {

  want <- function(name) is.null(select.indicator) || name %in% select.indicator

  Ys <- sort(Y)  # sort() drops NA by default (na.last = NA), matching quantile(na.rm = TRUE)
  n <- length(Ys)

  q_at <- function(p) {
    h <- (n - 1) * p + 1
    lo <- floor(h); hi <- ceiling(h)
    Ys[lo] + (h - lo) * (Ys[hi] - Ys[lo])
  }

  mean_est <- mean(Y, na.rm = TRUE)
  Hcr_est <- mean(Y < threshold, na.rm = TRUE)

  quant_preds <- c(Quant10 = NA_real_, Quant25 = NA_real_, Median = NA_real_,
                    Quant75 = NA_real_, Quant90 = NA_real_)
  if (want("Quant10") || want("Quant25") || want("Median") || want("Quant75") || want("Quant90")) {
    quant_preds <- c(Quant10 = q_at(0.1), Quant25 = q_at(0.25), Median = q_at(0.5),
                      Quant75 = q_at(0.75), Quant90 = q_at(0.9))
  }

  Gini_est <- NA_real_
  if (want("Gini")) {
    # ineq::Gini(y, na.rm = TRUE) does y <- as.numeric(na.omit(y)) internally after
    # y[y < 0] <- NA -- reuse the already-sorted Ys directly when there are no
    # negatives to remove (the common case, e.g. after exp() back-transformation),
    # otherwise re-sort the filtered subset.
    Yg <- if (any(Y < 0, na.rm = TRUE)) sort(Y[Y >= 0 & !is.na(Y)]) else Ys
    ng <- length(Yg)
    i <- seq_len(ng)
    Gini_est <- (2 * sum(i * Yg) / (ng * sum(Yg))) - (ng + 1) / ng
  }

  Pgap_est <- NA_real_
  if (want("Pgap")) {
    y2 <- Y
    y2[y2 < 0] <- NA
    Pgap_est <- mean((y2 < threshold) * (threshold - y2) / threshold, na.rm = TRUE)
  }

  Qsr_est <- NA_real_
  if (want("Qsr")) {
    q80 <- q_at(0.8); q20 <- q_at(0.2)
    # Uses the ORIGINAL Y (not the NA-stripped Ys) so NA propagates identically to
    # the original qsr_function(y) -- do not change this to Ys, it would silently
    # produce a different (non-NA) result when Y contains NA.
    Qsr_est <- sum(Y[(Y > q80)]) / sum(Y[(Y < q20)])
  }

  indicators <- cbind(mean_est, t(quant_preds), Gini_est, Hcr_est, Pgap_est, Qsr_est)

  colnames(indicators) <- c(
    "Mean", "Quant10", "Quant25", "Median", "Quant75",
    "Quant90", "Gini", "Hcr", "Pgap", "Qsr"
  )

  if (!is.null(custom)) {
    custom_ind <- unlist(lapply(custom, function(f) f(Y, threshold)))
    indicators <- cbind(indicators, t(custom_ind))
  }

  if (!is.null(select.indicator)) {
    invalid <- setdiff(select.indicator, colnames(indicators))
    if (length(invalid) > 0L) {
      stop(
        "Invalid 'select.indicator' value(s): ", paste(invalid, collapse = ", "),
        ". Valid indicators are: ", paste(colnames(indicators), collapse = ", "), ".",
        call. = FALSE
      )
    }
    indicators <- indicators[, select.indicator]
  }

  return(indicators)
}


# Function providing SAE-specific information ---------------------------------------------
sae_specs <- function(dName, cns, smp) {
  in_dom <- unique(smp[[dName]])
  total_dom <- unique(cns[[dName]])
  OOsamp <- !total_dom %in% in_dom

  return(list(
    N_surv = length(smp[[dName]]),
    N_pop = length(cns[[dName]]),
    D_out = sum(OOsamp),
    D_in = length(in_dom),
    D_total = length(total_dom),
    ni_smp = table(smp[[dName]]),
    ni_pop = table(cns[[dName]]),
    dName = dName
  ))
}


# Draws random survey samples in the MSE procedures ----------------------------------------
sample_select <- function(pop, smp, dName) {
  smpSizes <- table(smp[dName])
  smpSizes <- data.frame(
    smpidD = as.character(names(smpSizes)), n_smp = as.numeric(smpSizes),
    stringsAsFactors = FALSE
  )

  smpSizes <- dplyr::left_join(data.frame(idD = as.character(unique(pop[[dName]]))),
    smpSizes,
    by = c("idD" = "smpidD")
  )

  smpSizes$n_smp[is.na(smpSizes$n_smp)] <- 0

  splitPop <- split(pop, pop[[dName]])

  stratSamp <- function(dfList, ns) {
    do.call(rbind, mapply(dfList, ns, FUN = function(df, n) {
      popInd <- seq_len(nrow(df))
      sel <- base::sample(popInd, n, replace = FALSE)
      df[sel, ]
    }, SIMPLIFY = F))
  }

  samples <- stratSamp(splitPop, smpSizes$n_smp)

  rm(splitPop)
  return(samples)
}


# Computes REB random components in the MSE procedures ------------------------------------
ran_comp <- function(mod, smp_data, Y, dName, ADJsd) {
  scale_to <- function(x, target) {
    s <- sd(x)
    if (!is.finite(s) || s == 0) return(x - mean(x))
    (x / s) * target - mean((x / s) * target)
  }

  forest_res1 <- Y - predict(mod$Forest, smp_data)$predictions
  smp_data$forest_res <- forest_res1

  # Random Effects
  formRF <- formula(paste("forest_res ~", paste0(dName)))
  ran_effs1 <- aggregate(data = smp_data, formRF, FUN = mean)
  colnames(ran_effs1) <- c(dName, "r_bar")

  smp_data <- dplyr::left_join(smp_data, ran_effs1, by = dName)
  smp_data$forest_eij <- smp_data$forest_res - smp_data$r_bar

  # prepare for sampling (guarded against zero/non-finite SD, then centered)
  forest_res <- smp_data$forest_eij
  forest_res <- scale_to(forest_res, ADJsd)

  # prepare for sampling (guarded against zero/non-finite SD, then centered)
  ran_effs <- ran_effs1$r_bar
  ran_effs <- scale_to(ran_effs, mod$RanEffSD)

  return(list(
    forest_res = forest_res,
    ran_effs = ran_effs,
    smp_data = smp_data
  ))
}


# Computes wild random components in the MSE procedures -----------------------------------
ran_comp_wild <- function(mod, smp_data, Y, dName, ADJsd) {
  forest_res <- Y - mod$Forest$predictions - predict(mod$EffectModel, smp_data)
  forest_res <- forest_res - mean(forest_res)

  # Random Effects
  ran_effs <- unique(predict(mod$EffectModel, smp_data))
  ran_effs <- (ran_effs / sd(ran_effs)) * mod$RanEffSD
  ran_effs <- ran_effs - mean(ran_effs)

  return(list(
    forest_res = forest_res,
    ran_effs = ran_effs
  ))
}


# Block-samples errors in the MSE procedures ----------------------------------------------
block_sample <- function(domains, in_samp, smp_data, dName, pop_data, forest_res) {
  block_err <- vector(mode = "list", length = length(domains))

  for (idd in which(in_samp)) {
    block_err[[idd]] <- sample(forest_res[smp_data[dName] == domains[idd]], size = sum(pop_data[dName] == domains[idd]), replace = TRUE)
  }
  if (sum(in_samp) != length(domains)) {
    for (idd in which(!in_samp)) {
      block_err[[idd]] <- sample(forest_res, size = sum(pop_data[dName] == domains[idd]), replace = TRUE)
    }
  }
  return(unlist(block_err))
}


# Samples wild errors in the MSE procedures -----------------------------------------------
wild_errors <- function(x, mod, smp_data, forest_res) {
  fitted_s <- predict(mod$Forest, smp_data)$predictions + predict(mod$EffectModel, smp_data)
  indexer <- vapply(x, function(x) {
    which.min(abs(x - fitted_s))
  },
  FUN.VALUE = integer(1)
  )
  # superpopulation individual errors
  eps <- forest_res[indexer]
  wu <- sample(c(-1, 1), size = length(eps), replace = TRUE)
  eps <- abs(eps) * wu

  return(eps)
}


# Processing the threshold input for various nonlinear functions --------------------------
get_thresh <- function(x, threshold) {
  if (is.numeric(threshold)) {
    thresh <- threshold
  }
  if (is.null(threshold)) {
    thresh <- 0.6 * median(x, na.rm = TRUE)
  }
  if (is.function(threshold)) {
    thresh <- threshold(x)
  }
  return(thresh)
}


# Sort and process final function outputs -------------------------------------------------
sortAlpha <- function(x, dName) {
  orderIndicat <- order(as.character(x[[dName]]))
  return_val <- x[orderIndicat, ]
  rownames(return_val) <- NULL
  return(return_val)
}


# Wrapper to compute the empirical likelihood ---------------------------------------------
elm_wrapper <- function(X_input_elm, mu_input_elm) {
  n <- nrow(X_input_elm)
  p <- ncol(X_input_elm)

  elm_results <- elm(
    x = X_input_elm, mu = mu_input_elm, lam = rep(0, p), maxit = 25, gradtol = 1e-7,
    svdtol = 1e-9, itertrace = 0
  )

  logelr <- elm_results$logelr
  lambda <- elm_results$lambda

  zz <- X_input_elm - t(matrix(mu_input_elm, nrow = p, ncol = n))
  tmp3 <- n * (zz %*% lambda + 1)
  prob <- 1 / tmp3

  return(list(prob = prob, experWeight = elm_results$wts))
}


# Dependecies for elm Wrapper -------------------------------------------------------------
elm <- function(x, mu = rep(0, ncol(x)), lam = t(rep(0, ncol(x))), maxit = 25,
                gradtol = 1e-7, svdtol = 1e-9, itertrace = 0) {
  n <- nrow(x)
  p <- ncol(x)

  if (ncol(x) != length(mu)) {
    stop("Need size(mu) = size(x(i,:)")
  }

  if (length(lam) != length(mu)) {
    stop("Lam must be sized like mu")
  }

  if (gradtol < 1e-16) {
    gradtol <- 1e-16
  }
  if (svdtol < 1e-8) {
    gradtol <- 1e-8
  }

  z <- x - t(matrix(mu, nrow = p, ncol = n))
  newton_wts <- c(1 / 3^c(0:3), rep(0, 12))
  gradient_wts <- .5^c(0:15)
  gradient_wts <- (gradient_wts^2 - newton_wts^2)^0.5

  gradient_wts[12:16] <- gradient_wts[12:16] / (.1^-(1:5))

  nits <- 0
  gsize <- gradtol + 1.0

  while (nits < maxit & gsize > gradtol) {
    arg <- 1 + z %*% lam
    wts1 <- plog(arg, 1 / n, 1)
    wts2 <- (-plog(arg, 1 / n, 2))^0.5

    grad <- -z * matrix(wts1, nrow = nrow(wts1), ncol = p)
    grad <- colSums(grad)
    gsize <- mean(abs(grad))

    hess <- z * matrix(wts2, nrow = nrow(wts2), ncol = p)

    hu <- svd(hess)$u
    hs <- diag(length(svd(hess)$d)) * svd(hess)$d
    hv <- svd(hess)$v
    dhs <- svd(hess)$d

    if (min(dhs) < max(dhs) * svdtol + 1e-128) {
      dhs <- dhs + svdtol * max(dhs) + 1e-128
    }

    dhs <- 1 / dhs
    hs <- diag(dhs)

    nstep <- hv %*% hs %*% t(hu) %*% (wts1 / wts2)
    gstep <- -grad
    if (sum(nstep^2) < sum(gstep^2)) {
      gstep <- gstep * sum(nstep^2)^0.5 / sum(gstep^2)^0.5
    }

    ologelr <- -sum(plog(arg, 1 / n))
    ninner <- 0

    for (i in 1:length(newton_wts)) {
      nlam <- lam + newton_wts[i] * nstep + gradient_wts[i] * gstep
      nlogelr <- -sum(plog(1 + z %*% nlam, 1 / n))

      if (nlogelr < ologelr) {
        lam <- nlam
        ninner <- i
        break
      }
    }

    nits <- nits + 1
    if (ninner == 0) {
      nits <- maxit
    }

    if (itertrace == 1) {
      print(list(lam, nlogelr, gsize, ninner, nits))
    }
  }

  logelr <- nlogelr
  lambda <- lam
  hess <- t(hess) %*% hess
  wts <- wts1

  Ergebnisse <- list(logelr, lambda, grad, hess, wts, nits)
  names(Ergebnisse) <- c("logelr", "lambda", "grad", "hess", "wts", "nits")

  return(Ergebnisse)
}


# Pseudo-logarithm
plog <- function(z, eps, d = 0) {
  zsize <- dim(z)
  out <- c(z)
  low <- out < eps

  if (d == 0) {
    out[!low] <- log(out[!low])
    out[low] <- log(eps) - 1.5 + 2 * out[low] / eps - 0.5 * (out[low] / eps)^2
  } else if (d == 1) {
    out[!low] <- 1 / out[!low]
    out[low] <- 2 / eps - out[low] / eps^2
  } else if (d == 2) {
    out[!low] <- -1 / out[!low]^2
    out[low] <- -1 / eps^2
  } else {
    stop("Unknown option d for pseudologarithm")
  }
  llogz <- matrix(out, zsize[1], zsize[2])
  return(llogz)
}


# Run a bootstrap loop optionally in parallel across `cores` fork workers with
# reproducible L'Ecuyer streams, scoped to this call. `expr` is function(cl) that
# invokes the pbapply call with the given cl. cores == 1 -> serial (cl = NULL), no
# RNG side effects (byte-identical to a bare pbapply(..., cl = NULL) call). cores > 1
# -> L'Ecuyer-CMRG seeded from a value drawn from the current (already-seeded) stream,
# so the parallel bootstrap is reproducible for a fixed (top-level seed, cores); the
# caller's RNGkind and .Random.seed are restored on exit so neither the session nor
# the already-computed point estimates are affected.
with_parallel_rng <- function(cores, expr) {
  if (!(cores > 1L)) {
    return(expr(NULL))
  }
  par_seed <- sample.int(.Machine$integer.max, 1L)
  old_kind <- RNGkind()
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(par_seed)
  expr(cores)
}

# Extra estimator arguments (the MSE driver's `...`) with ranger's num.threads forced
# to 1 when running in parallel, so each fork worker uses one thread (no
# oversubscription). cores == 1 returns the args unchanged. Warns once if the caller
# set a conflicting num.threads.
force_serial_threads <- function(cores, ...) {
  dots <- list(...)
  if (cores > 1L) {
    if (!is.null(dots$num.threads) && !identical(as.integer(dots$num.threads), 1L)) {
      warning("num.threads forced to 1 inside parallel MSE workers (cores > 1) to avoid oversubscription.")
    }
    dots$num.threads <- 1L
  }
  dots
}
