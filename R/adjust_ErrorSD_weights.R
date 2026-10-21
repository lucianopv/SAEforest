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
# Memory note: dense, hence O(n^2) (about 1.15 GB at n = 12,000). Kept as the
# reference the tests hold oob_predict_boot_() to; adjust_ErrorSD_weights_() no
# longer builds it.
oob_weight_matrix_ <- function(rf, X, num.threads = 1) {
  if (is.null(rf$inbag.counts)) {
    stop("oob_weight_matrix_: the forest was grown without keep.inbag = TRUE, ",
         "so the bootstrap multiplicities needed for the weights are unavailable.")
  }
  # num.threads is pinned: predict.ranger defaults to NULL, which ranger reads as
  # "use every core". Inside a forked worker (the MSE bootstrap, or any parallel
  # driver) that oversubscribes the machine, so the caller's threading choice is
  # honoured instead, defaulting to serial.
  nodes <- stats::predict(rf, X, type = "terminalNodes",
                          num.threads = num.threads)$predictions
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

# W %*% Ystar without forming W.
#
# Same quantity as oob_weight_matrix_(rf, X) %*% Ystar on rows with an
# out-of-bag prediction (other rows 0), computed tree by tree: within tree k the
# in-bag-weighted sum of Ystar by terminal node, divided by the node's in-bag
# count, is the tree's prediction for every observation falling in that node, and
# only the observations out-of-bag in k receive it. Averaging over the K_i trees
# in which i is out-of-bag gives row i of W %*% Ystar. Memory is O(n B).
#
# Ystar is an n x B matrix. Returns an n x B matrix with attribute "n_oob".
# `nodes` lets a caller that has already run predict(type = "terminalNodes")
# pass its result (predict.ranger can touch R's RNG, so the caller may need to
# fix where in its stream that happens).
oob_predict_boot_ <- function(rf, X, Ystar, num.threads = 1, nodes = NULL) {
  if (is.null(rf$inbag.counts)) {
    stop("oob_predict_boot_: the forest was grown without keep.inbag = TRUE, ",
         "so the bootstrap multiplicities needed for the weights are unavailable.")
  }
  Ystar <- as.matrix(Ystar)
  if (is.null(nodes)) {
    nodes <- stats::predict(rf, X, type = "terminalNodes",
                            num.threads = num.threads)$predictions
  }
  inbag <- do.call(cbind, rf$inbag.counts)
  n <- nrow(nodes)
  P <- matrix(0, n, ncol(Ystar))
  n_oob <- integer(n)

  for (k in seq_len(ncol(nodes))) {
    b <- inbag[, k]
    leaf <- nodes[, k]
    oob <- which(b == 0)
    if (!length(oob)) next
    n_oob[oob] <- n_oob[oob] + 1L
    S <- rowsum(b * Ystar, leaf, reorder = FALSE)   # in-bag weighted sums per leaf
    tot <- rowsum(b, leaf, reorder = FALSE)[, 1]    # in-bag counts per leaf
    idx <- match(leaf[oob], as.numeric(rownames(S)))
    t_oob <- tot[idx]
    keep <- is.finite(t_oob) & t_oob != 0           # as the dense code skips them
    if (!any(keep)) next
    rows <- oob[keep]
    P[rows, ] <- P[rows, , drop = FALSE] +
      S[idx[keep], , drop = FALSE] / t_oob[keep]
  }

  P <- P / pmax(n_oob, 1)
  attr(P, "n_oob") <- n_oob
  P
}

# The forest's in-sample prediction for `rows`, from its leaves and no RNG.
#
# In each tree the prediction for an observation is its terminal node's in-bag
# mean of Y, each member counted with its bootstrap multiplicity (ranger's own
# leaf value); averaging over all trees gives predict(rf, X)$predictions without
# calling predict.ranger, which can touch R's RNG. Y must be the outcome the
# forest was grown on and `nodes` its terminal nodes for X.
inbag_predict_ <- function(rf, nodes, Y, rows) {
  inbag <- do.call(cbind, rf$inbag.counts)
  acc <- numeric(length(rows))
  for (k in seq_len(ncol(nodes))) {
    b <- inbag[, k]
    leaf <- nodes[, k]
    S <- rowsum(b * Y, leaf, reorder = FALSE)[, 1]   # in-bag weighted sums per leaf
    tot <- rowsum(b, leaf, reorder = FALSE)[, 1]     # in-bag counts per leaf
    idx <- match(leaf[rows], as.numeric(names(S)))
    acc <- acc + S[idx] / tot[idx]
  }
  unname(acc / ncol(nodes))
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
adjust_ErrorSD_weights_ <- function(Y, X, rf, B = 100, num.threads = 1) {
  if (!(is.numeric(B) && length(B) == 1 && !is.na(B) && B >= 1)) {
    stop("adjust_ErrorSD_weights_: B must be a single numeric value >= 1.")
  }
  if (is.null(rf$inbag.counts)) {
    stop("adjust_ErrorSD_weights_: the forest was grown without keep.inbag = TRUE, ",
         "so the bootstrap multiplicities needed for the weights are unavailable.")
  }
  # The terminal nodes come first, as the dense route's weight matrix did, so the
  # draws below start from the same point of the RNG stream (predict.ranger
  # consumes it).
  nodes <- stats::predict(rf, X, type = "terminalNodes",
                          num.threads = num.threads)$predictions
  # n_oob (the K_i) from the in-bag counts alone: the number of trees with b == 0
  n_oob <- rowSums(do.call(cbind, rf$inbag.counts) == 0)
  ok <- n_oob > 0 & is.finite(rf$predictions)
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

  f_ok <- f_oob[ok]
  # The anchor of Ystar, which W reads for every in-bag row. A row in-bag in every
  # tree has no OOB prediction; it takes the forest's in-sample one, so the anchor
  # moves with Y as the OOB predictions do (anchored at 0 instead, each OOB
  # neighbour's gap picked up about -w * f_j and K depended on where Y sits).
  # Computed without RNG, so the draws below are unchanged.
  f_full <- f_oob
  miss <- which(!is.finite(f_oob))
  if (length(miss)) f_full[miss] <- inbag_predict_(rf, nodes, Y, miss)
  n <- length(f_oob)

  # Step 2, B pseudo-outcome vectors, drawn in the order the dense route drew them
  Ystar <- matrix(0, n, B)
  for (b in seq_len(B)) Ystar[, b] <- f_full + sample(e_ij, n, replace = TRUE)

  # Steps 3 and 4, without forming W (O(n B) memory)
  P <- oob_predict_boot_(rf, X, Ystar, num.threads = num.threads, nodes = nodes)
  g <- colMeans((P[ok, , drop = FALSE] - f_ok)^2)

  list(K = mean(g), m = as.integer(B))
}
