# Fixed-weight variant of MERFranger's per-iteration variance bias correction.
#
# Krennmair & Schmid (2022, step 3) recompute the bootstrap OOB predictions
# "using a random forest", i.e. a new forest per draw, which is what
# adjust_ErrorSD_() does. Mendez & Lohr (2011), whose correction this is, keep
# the forest's weights FIXED in their Procedures 1 and 2 and grow no new trees;
# their Section 3.1 reports that the refitting variant "did not perform better
# ... in fact, the full bootstrap may slightly overcorrect for the bias because
# new trees are grown for each iteration".
#
# This file implements their Procedure 2 (nonparametric bootstrap, equation 11),
# which is the one that matches the resampled-residual scheme already used here:
# the pseudo-errors are resampled from the centred OOB residuals, and the
# bootstrap OOB predictions are W %*% y_star with W held fixed.
#
# Opt-in only, via MERFranger(var.adjust = TRUE, adjust_method = "weights").
# The default path is untouched.

# The forest as a linear smoother (Mendez & Lohr, equation 8).
#
# For a tree k in which observation i is out-of-bag, the tree's prediction for i
# is the in-bag mean of its terminal node, each observation counted with its
# bootstrap multiplicity b_jk:
#
#   w_ijk = b_jk / n(tau_k(x_j))   if x_i falls in the same node as x_j, else 0,
#
# where n(tau) sums those multiplicities, so the row sums to one. Averaging over
# the K_i trees in which i is out-of-bag gives w_bar_ij, and the vector of OOB
# predictions is exactly W %*% y. Requires the forest to have been grown with
# keep.inbag = TRUE.
#
# Returns an n x n matrix with attribute "n_oob" (the K_i). Rows for
# observations that were never out-of-bag are all zero, matching ranger, which
# returns NA for their OOB prediction; callers must drop those rows.
#
# Memory note: dense, hence O(n^2). Fine for the simulation sizes this was
# written for (n ~ 1e3); a sparse representation would be needed before this
# could run at GridSAE's n ~ 1.3e4.
oob_weight_matrix_ <- function(rf, X) {
  if (is.null(rf$inbag.counts)) {
    stop("oob_weight_matrix_: the forest was grown without keep.inbag = TRUE, ",
         "so the bootstrap multiplicities needed for the weights are unavailable.")
  }
  nodes <- stats::predict(rf, X, type = "terminalNodes")$predictions
  inbag <- do.call(cbind, rf$inbag.counts)
  n <- nrow(nodes)
  W <- matrix(0, n, n)
  n_oob <- integer(n)

  for (k in seq_len(ncol(nodes))) {
    b <- inbag[, k]
    leaf <- nodes[, k]
    oob <- which(b == 0)
    if (!length(oob)) next
    n_oob[oob] <- n_oob[oob] + 1L
    leaf_total <- tapply(b, leaf, sum)
    for (lf in unique(leaf[oob])) {
      total <- leaf_total[[as.character(lf)]]
      if (!is.finite(total) || total == 0) next
      members <- which(leaf == lf & b > 0)
      rows <- oob[leaf[oob] == lf]
      W[rows, members] <- W[rows, members] +
        rep(b[members] / total, each = length(rows))
    }
  }

  W <- W / pmax(n_oob, 1)
  attr(W, "n_oob") <- n_oob
  W
}

# K by Mendez & Lohr Procedure 2 with the weights held fixed.
#
# Same contract as adjust_ErrorSD_(): returns list(K, m). Rows without an OOB
# prediction are dropped, exactly as adjust_ErrorSD_() drops them, for the same
# reason -- they carry no information about the OOB gap.
#
# No seed argument, for the same reason as adjust_ErrorSD_(): it inherits the
# stream seeded once by the caller. It consumes less RNG than the refit route
# (n * B draws and no forests), so the two routes' streams differ by design.
adjust_ErrorSD_weights_ <- function(Y, X, rf, B = 100) {
  if (!(is.numeric(B) && length(B) == 1 && !is.na(B) && B >= 1)) {
    stop("adjust_ErrorSD_weights_: B must be a single numeric value >= 1.")
  }
  W <- oob_weight_matrix_(rf, X)
  ok <- attr(W, "n_oob") > 0 & is.finite(rf$predictions)
  if (!any(ok)) {
    stop("adjust_ErrorSD_weights_: no out-of-bag predictions to work from ",
         "(every observation was in-bag in every tree).")
  }

  f_oob <- rf$predictions
  e_ij <- Y - f_oob                      # marginal OOB residuals
  e_ij <- e_ij[is.finite(e_ij)]
  e_ij <- e_ij - mean(e_ij)              # centred, as in Procedure 2 step 1
  if (!length(e_ij)) {
    stop("adjust_ErrorSD_weights_: no finite OOB residuals to resample from.")
  }

  W_ok <- W[ok, , drop = FALSE]
  f_ok <- f_oob[ok]
  f_full <- ifelse(is.finite(f_oob), f_oob, 0)  # anchor for the rows W reads from
  n <- length(f_oob)

  g <- vapply(seq_len(B), function(b) {
    y_star <- f_full + sample(e_ij, n, replace = TRUE)   # step 2
    mean((as.vector(W_ok %*% y_star) - f_ok)^2)          # steps 3 and 4
  }, numeric(1))

  list(K = mean(g), m = as.integer(B))
}
