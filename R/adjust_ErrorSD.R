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
  n <- length(rf$predictions)

  # centred OOB residuals (paper step 3a) — RNG-free, shared by both paths
  e_ij <- Y - rf$predictions
  e_ij <- e_ij - mean(e_ij)

  if (adj_tol <= 0) {
    ## ---- legacy exact path: fixed B, vectorised ----
    pred_OOB <- matrix(rf$predictions, ncol = B, nrow = n, byrow = FALSE)
    y_star_OOB <- pred_OOB + sample(e_ij, size = length(pred_OOB), replace = TRUE)
    fits <- apply(y_star_OOB, 2, function(x) {
      ranger::ranger(y = x, x = X, data = smp_data, ...)
    })
    pred_OOB_star <- sapply(fits, function(x) x$predictions)
    return(list(K = mean(rowMeans((pred_OOB_star - pred_OOB)^2)), m = B))
  }

  ## ---- adaptive early-stop path ----
  B_MIN   <- min(20L, B)   # never stop below this many forests
  B_BATCH <- min(10L, B)   # forests fit per convergence check
  g <- numeric(0)          # per-forest mean squared OOB gap; K = mean(g)
  repeat {
    for (k in seq_len(B_BATCH)) {
      if (length(g) >= B) break
      y_star_b <- rf$predictions + sample(e_ij, size = n, replace = TRUE)
      fit_b <- ranger::ranger(y = y_star_b, x = X, data = smp_data, ...)
      g <- c(g, mean((fit_b$predictions - rf$predictions)^2))
    }
    Kbar   <- mean(g)
    # length(g) >= 2 guard: sd() of a length-1 vector is NA, which would make the
    # stop condition below error on a degenerate B = 1. No-op for realistic B >= 2
    # (first batch yields length(g) = min(B_BATCH, B) >= 2), so it changes no numbers.
    rel_se <- if (length(g) >= 2 && Kbar > 0) (stats::sd(g) / sqrt(length(g))) / Kbar else Inf
    if (length(g) >= B_MIN && rel_se < adj_tol) break
    if (length(g) >= B) break
  }
  list(K = mean(g), m = length(g))
}
