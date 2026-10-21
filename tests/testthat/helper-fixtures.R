# Tiny deterministic subset of eusilcA for fast, reproducible tests.
tiny_saef_data <- function(n_dom = 4, per_dom = 8, seed = 2026) {
  testthat::skip_if_not_installed("emdi")
  set.seed(seed)
  data("eusilcA_smp", package = "emdi", envir = environment())
  data("eusilcA_pop", package = "emdi", envir = environment())
  doms <- sort(unique(as.character(eusilcA_smp$district)))[seq_len(n_dom)]

  smp <- do.call(rbind, lapply(doms, function(d) {
    rows <- which(eusilcA_smp$district == d)
    eusilcA_smp[sample(rows, min(per_dom, length(rows))), ]
  }))
  pop <- eusilcA_pop[eusilcA_pop$district %in% doms, ]
  # cap population per domain to keep smearing cheap
  pop <- do.call(rbind, lapply(doms, function(d) {
    rows <- which(pop$district == d)
    pop[rows[seq_len(min(40, length(rows)))], ]
  }))

  smp <- droplevels(smp); pop <- droplevels(pop)
  list(
    Y = smp$eqIncome,
    X = smp[, c("gender", "eqsize", "cash", "self_empl", "unempl_ben")],
    smp = smp[, c("district", "gender", "eqsize", "cash",
                  "self_empl", "unempl_ben", "eqIncome")],
    pop = pop[, c("district", "gender", "eqsize", "cash",
                  "self_empl", "unempl_ben")],
    dName = "district"
  )
}

# Same as tiny_saef_data(), plus a "subunit" column (2 subunits per district)
# nested within district, added to BOTH pop and smp -- for testing aggregate_to
# in either direction: finer than dName (dName="district", aggregate_to="subunit")
# or coarser than dName (dName="subunit", aggregate_to="district").
tiny_saef_data_with_subunit <- function(n_dom = 4, per_dom = 8, seed = 2026) {
  base <- tiny_saef_data(n_dom = n_dom, per_dom = per_dom, seed = seed)
  base$pop$subunit <- paste0(base$pop$district, "_", rep(1:2, length.out = nrow(base$pop)))
  base$smp$subunit <- paste0(base$smp$district, "_", rep(1:2, length.out = nrow(base$smp)))
  base
}

# Tiny deterministic aggData=TRUE fixture (aggregated-covariate / ELM-calibration path),
# using SAEforest's own bundled data (no emdi dependency). Restricted to in-sample-only
# domains so the out-of-sample augmentation branch of point_meanAGG() is never exercised.
tiny_agg_data <- function(seed = 2026, per_dom = 15) {
  data("eusilcA_smp", "eusilcA_popAgg", "popNsize", package = "SAEforest", envir = environment())

  doms <- c("Graz (Stadt)", "Linz (Stadt)", "Salzburg-Umgebung")
  covars <- c("eqsize", "cash", "self_empl", "unempl_ben")

  set.seed(seed)
  smp <- do.call(rbind, lapply(doms, function(d) {
    rows <- which(as.character(eusilcA_smp$district) == d)
    eusilcA_smp[sample(rows, min(per_dom, length(rows))), ]
  }))
  smp <- droplevels(smp)

  # popAgg and popN must have identical row counts AND identical dName ordering
  # (error_checkFunctions.R checks all.equal(popnsize[[dName]], pop_data[[dName]])).
  popAgg <- eusilcA_popAgg[match(doms, as.character(eusilcA_popAgg$district)), ]
  popAgg <- droplevels(popAgg)
  rownames(popAgg) <- NULL

  popN <- popNsize[match(doms, as.character(popNsize$district)), ]
  rownames(popN) <- NULL

  list(
    Y = smp$eqIncome,
    X = smp[, covars],
    smp = smp[, c("district", covars, "eqIncome")],
    popAgg = popAgg[, c("district", covars)],
    popN = popN,
    dName = "district"
  )
}

# Fixture in which sample_select()'s two orderings diverge.
#
# tiny_saef_data() hides the divergence twice over: its `pop` is built in sorted
# domain order (so first-appearance order == factor-level order), and every
# domain gets the same number of sampled rows (so permuting the sizes is a
# no-op). This fixture breaks both, reproducing the GridSAE geometry -- pop_data
# comes from merge(..., by = "cell_id"), so it is ordered by cell, not by the
# random-effect domain.
tiny_saef_data_unsorted <- function(seed = 2026) {
  d <- tiny_saef_data(n_dom = 4, per_dom = 8, seed = seed)
  doms <- sort(unique(as.character(d$smp[[d$dName]])))
  stopifnot(length(doms) == 4L)

  # Unequal per-domain sample sizes, ascending in sorted-domain order.
  keep <- unlist(Map(function(dom, k) {
    utils::head(which(as.character(d$smp[[d$dName]]) == dom), k)
  }, doms, c(2L, 4L, 6L, 8L)), use.names = FALSE)
  d$smp <- d$smp[keep, , drop = FALSE]
  d$Y   <- d$Y[keep]
  d$X   <- d$X[keep, , drop = FALSE]

  # Population in REVERSE domain order: first-appearance order is now the
  # reverse of split()'s factor-level order, deterministically.
  d$pop <- d$pop[order(as.character(d$pop[[d$dName]]), decreasing = TRUE), ,
                 drop = FALSE]

  rownames(d$smp) <- NULL
  rownames(d$pop) <- NULL
  d
}

# The dense n x n route of adjust_ErrorSD_weights_(), kept verbatim as the
# reference the matrix-free route is held to (same RNG draws, same K).
dense_K_reference <- function(Y, X, rf, B = 100, num.threads = 1) {
  W <- oob_weight_matrix_(rf, X, num.threads = num.threads)
  ok <- attr(W, "n_oob") > 0 & is.finite(rf$predictions)
  f_oob <- rf$predictions
  e_ij <- Y - f_oob
  e_ij <- e_ij[is.finite(e_ij)]
  e_ij <- e_ij - mean(e_ij)
  W_ok <- W[ok, , drop = FALSE]
  f_ok <- f_oob[ok]
  f_full <- ifelse(is.finite(f_oob), f_oob, 0)
  n <- length(f_oob)
  g <- vapply(seq_len(B), function(b) {
    y_star <- f_full + sample(e_ij, n, replace = TRUE)
    mean((as.vector(W_ok %*% y_star) - f_ok)^2)
  }, numeric(1))
  list(K = mean(g), m = as.integer(B))
}

# Synthetic toy for the var.adjust routes: n units in D balanced areas, a
# non-linear signal and an area effect. No real data involved.
make_toy <- function(n = 120, D = 6, seed = 1) {
  set.seed(seed)
  area <- factor(rep(seq_len(D), each = n / D))
  x1 <- rnorm(n); x2 <- rnorm(n)
  v <- rnorm(D, 0, 0.5)[as.integer(area)]
  y <- 2 * x1 + x2^2 + v + rnorm(n, 0, 1)
  list(Y = y, X = data.frame(x1 = x1, x2 = x2),
       data = data.frame(area = area, x1 = x1, x2 = x2, y = y))
}

# Non-degenerate var.adjust fits on make_toy(), one per route the study uses,
# reduced to what an EM-loop or anchor change could move. Held against
# fixtures/varadjust_fits_ref.rds (frozen at 3250638 by
# make_varadjust_fits_reference.R): the degenerate-case guards must leave every
# one of these fits byte-identical.
varadjust_guard_fits <- function() {
  d <- make_toy()
  common <- list(Y = d$Y, X = d$X, random = "(1|area)", data = d$data,
                 var.adjust = TRUE, num.threads = 1, seed = 3)
  cfgs <- list(
    weights  = list(adjust_method = "weights", B_adj = 10, num.trees = 100,
                    min.node.size = 5),
    refit    = list(B_adj = 5, adj_tol = 0, num.trees = 50),
    plugin   = list(K_fixed = c(0.4, 0.3), num.trees = 50),
    crossfit = list(adjust_method = "weights", re_method = "crossfit",
                    re_folds = 3, B_adj = 10, num.trees = 50)
  )
  lapply(cfgs, function(cfg) {
    m <- do.call(MERFranger, c(common, cfg))
    list(K = m$K, sigma_e = m$sigma_e, sigma_u = m$sigma_u, LogLik = m$LogLik,
         IterationsUsed = m$IterationsUsed, ErrorSD = m$ErrorSD,
         RanEffSD = m$RanEffSD, OOBresiduals = m$OOBresiduals,
         predictions = m$Forest$predictions, ranef = unlist(m$RandomEffects),
         naive = m$sd_Mendez_Lohr_naive_unad, ErrorSD_unadj = m$ErrorSD_unadj,
         RanEffSD_unadj = m$RanEffSD_unadj)
  })
}
