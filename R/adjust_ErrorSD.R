# Auxiliary function for MSEs to perform bias correction ---------------------------------

adjust_ErrorSD <- function(Y, X, smp_data, mod, B = 100, ...) {
  pred_OOB <- matrix(mod$Forest$predictions, ncol = B,
                     nrow = length(mod$Forest$predictions), byrow = FALSE)

  e_ij <- Y - predict(mod$Forest, smp_data)$predictions
  e_ij <- e_ij - mean(e_ij)

  y_star_OOB <- pred_OOB + sample(e_ij, size = length(pred_OOB), replace = TRUE)

  my_estim_f2 <- function(x) {
    ranger::ranger(y = x, x = X, data = smp_data, ...)
  }
  my_f_n2 <- pbapply::pbapply(y_star_OOB, 2, my_estim_f2)

  my_pred_f <- function(x) {
    x$predictions
  }
  pred_OOB_star <- sapply(my_f_n2, my_pred_f)

  Adjustment <- (pred_OOB_star - pred_OOB)^2

  bias <- mean(rowMeans(Adjustment))
  if (mod$ErrorSD^2 < bias) {
    warning("Variance bias correction exceeds ErrorSD^2; clamping adjusted SD to 0.")
  }
  outvar <- sqrt(max(0, mod$ErrorSD^2 - bias))

  return(outvar)
}

# Auxiliary function for MERFranger's per-iteration variance bias correction.
# Returns list(K, m): K = mean over units and inner-bootstrap forests of the squared
# gap between bootstrap-refit OOB predictions and the original OOB predictions
# (Mendez & Lohr, 2011); m = number of inner forests actually fit.
#
# adj_tol = 0 (or NULL): legacy exact path — fixed B forests, vectorised
#   (byte-identical to upstream / pre-adaptive numbers).
# adj_tol > 0: adaptive early-stop — accumulate per-forest gaps g_b (K = mean(g_b))
#   in batches, stop once the relative Monte-Carlo SE of K, sd(g)/sqrt(m)/mean(g),
#   drops below adj_tol (after at least min(B_MIN, B) forests), or the cap B is hit.
#
# No seed argument: this runs inside MERFranger's iteration (and inside the MSE
# bootstrap). It inherits the RNG stream seeded once at the top of SAEforest_model,
# so every refit is reproducible run-to-run yet distinct across bootstrap replicates.
adjust_ErrorSD_ <- function(Y, X, smp_data, rf, B = 100, adj_tol = 0, ...) {
  if (is.null(adj_tol)) adj_tol <- 0

  # centred OOB residuals (paper step 3a) — RNG-free, shared by both paths.
  # ranger returns NA for an observation that happens to be in-bag in EVERY tree,
  # so that observation contributes no OOB residual. Drop those before centring: a
  # single NA would otherwise make mean(e_ij) NA and poison every bootstrap draw
  # built from the pool. No-op when all residuals are finite, so the legacy
  # fixture stays byte-identical.
  e_ij <- Y - rf$predictions
  e_ij <- e_ij[is.finite(e_ij)]
  if (!length(e_ij)) {
    stop("adjust_ErrorSD_: no finite OOB residuals to resample from ",
         "(the forest returned no out-of-bag predictions).")
  }
  e_ij <- e_ij - mean(e_ij)

  # A row whose OOB prediction is itself missing (in-bag in every tree -> 0/0
  # NaN) cannot anchor a synthetic response: one such cell poisons y_star_b and
  # every inner forest dies with ranger's "Missing data in dependent variable."
  # (observed for real: GridSAE merf_coarse task 1443, seed 511100059).
  # Restrict the inner bootstrap to rows with a finite OOB prediction -- the
  # same drop the residual pool above already applies, and the same estimand
  # argument: a row with no OOB prediction carries no information about the
  # OOB gap. No-op (byte-identical, same RNG draws) when all rows are finite.
  ok <- is.finite(rf$predictions)
  pred_ok <- rf$predictions[ok]
  n_ok <- length(pred_ok)
  X_ok <- X[ok, , drop = FALSE]
  smp_ok <- smp_data[ok, , drop = FALSE]

  if (adj_tol <= 0) {
    ## ---- legacy exact path: fixed B, vectorised ----
    pred_OOB <- matrix(pred_ok, ncol = B, nrow = n_ok, byrow = FALSE)
    y_star_OOB <- pred_OOB + sample(e_ij, size = length(pred_OOB), replace = TRUE)
    fits <- apply(y_star_OOB, 2, function(x) {
      ranger::ranger(y = x, x = X_ok, data = smp_ok, ...)
    })
    pred_OOB_star <- sapply(fits, function(x) x$predictions)
    # na.rm at BOTH levels: this matrix is n x B (1.3M cells in the GridSAE
    # consumption fit), and a cell is NA whenever that inner forest left that
    # observation with no OOB prediction. Without na.rm one such cell makes its
    # whole row NA and the row makes K NA, which then killed the caller's
    # `if (K > naive_unadj^2)`. Missingness here is a resampling coincidence,
    # independent of the residual value, so averaging the observed cells is the
    # right estimator rather than a patch. No-op when nothing is missing.
    return(list(K = mean(rowMeans((pred_OOB_star - pred_OOB)^2, na.rm = TRUE),
                         na.rm = TRUE), m = B))
  }

  ## ---- adaptive early-stop path ----
  B_MIN   <- min(20L, B)   # never stop below this many forests
  B_BATCH <- min(10L, B)   # forests fit per convergence check
  g <- numeric(0)          # per-forest mean squared OOB gap; K = mean(g)
  repeat {
    for (k in seq_len(B_BATCH)) {
      if (length(g) >= B) break
      y_star_b <- pred_ok + sample(e_ij, size = n_ok, replace = TRUE)
      fit_b <- ranger::ranger(y = y_star_b, x = X_ok, data = smp_ok, ...)
      g <- c(g, mean((fit_b$predictions - pred_ok)^2, na.rm = TRUE))
    }
    Kbar   <- mean(g, na.rm = TRUE)
    # length(g) >= 2 guard: sd() of a length-1 vector is NA, which would make the
    # stop condition below error on a degenerate B = 1. No-op for realistic B >= 2
    # (first batch yields length(g) = min(B_BATCH, B) >= 2), so it changes no numbers.
    # is.finite(Kbar) guard: same failure as the legacy path -- a non-finite Kbar
    # made this very `if` throw "missing value where TRUE/FALSE needed". Treat it
    # as "not yet converged" and keep fitting rather than dying.
    rel_se <- if (length(g) >= 2 && is.finite(Kbar) && Kbar > 0)
      (stats::sd(g, na.rm = TRUE) / sqrt(length(g))) / Kbar else Inf
    if (length(g) >= B_MIN && rel_se < adj_tol) break
    if (length(g) >= B) break
  }
  list(K = mean(g, na.rm = TRUE), m = length(g))
}
