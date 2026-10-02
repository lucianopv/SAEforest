# One-off, 2026-10-02: regenerate ONLY the `nonlin` element of
# mse_memory_parity_ref.rds after the smearing-pool centring fix in
# MERFranger.R.
#
# Why the golden moved. The fix centres OOBresiduals on the var.adjust path,
# because `r_ij_new - AllEffects` left the observation-weighted BLUP mean in the
# pool and that is a pure level error in every smeared indicator. The parity
# goldens were frozen for edits that were meant to be numerically INERT, so a
# deliberate numerical change legitimately breaks them -- the guard did its job
# and this records the change rather than loosening the guard.
#
# Why this is safe to accept. On THIS fixture the change is floating-point only:
# max relative difference 8.637e-15 across all cells and indicators, about 39x
# machine epsilon. The reason is that tiny_saef_data() has 4 domains of exactly
# 8 observations (CV = 0.000), and with balanced domains the observation-weighted
# BLUP mean equals the domain-weighted one, which lme4 sets to zero -- so the
# subtraction is a mathematical no-op and only rounding differs. The fixture
# therefore could never have caught this bug; tests/testthat/test-pool-centring.R
# was added to cover it with deliberately size-correlated unbalanced domains.
#
# Why ONLY `nonlin`, and not a full re-run of
# make_mse_memory_parity_reference.R:
#   * `mean`       (test 4) runs var.adjust at its default FALSE, so the fix
#                  cannot reach it; it still passes and must not be disturbed.
#   * `nonlin_agg` (test 5) was ALREADY failing before this fix, for an
#                  unrelated golden row-order reason recorded in
#                  PROJECT_NOTES.md. Regenerating it would make that test pass
#                  and silently bury a real, still-undiagnosed discrepancy. It
#                  is deliberately left frozen so the suite keeps flagging it.
# Expected state after this script: test 3 passes, test 4 passes, test 5 still
# fails. A green test 5 means this script over-reached.

devtools::load_all(quiet = TRUE)
source(testthat::test_path("helper-fixtures.R"))

path <- testthat::test_path("fixtures", "mse_memory_parity_ref.rds")
old <- readRDS(path)
stopifnot(identical(sort(names(old)), c("mean", "nonlin", "nonlin_agg")))

d <- tiny_saef_data()
# Identical call to make_mse_memory_parity_reference.R's r1 -- same fixture,
# same seed, same arguments. Only the package code underneath has changed.
r1 <- SAEforest_model(
  Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
  meanOnly = FALSE, MSE = "nonparametric", B = 3,
  threshold = median(d$Y), var.adjust = TRUE, B_adj = 5, adj_tol = 0.05,
  num.trees = 50, mtry = 3, seed = 5
)$MSE_Estimates

# Refuse to write if the change is larger than floating-point noise: that would
# mean something other than the centring moved, and the golden should not be
# quietly replaced on this script's say-so.
num <- names(r1)[vapply(r1, is.numeric, logical(1))]
mx <- max(vapply(intersect(num, names(old$nonlin)), function(k) {
  max(abs(r1[[k]] - old$nonlin[[k]]) /
        pmax(abs(old$nonlin[[k]]), .Machine$double.eps), na.rm = TRUE)
}, numeric(1)))
cat(sprintf("max relative change vs frozen golden: %.3e\n", mx))
if (mx > 1e-10) {
  stop("change is larger than floating-point noise (", signif(mx, 3),
       "); investigate before regenerating this golden.")
}

new <- old
new$nonlin <- r1                      # mean and nonlin_agg preserved untouched
saveRDS(new, path)
cat("rewrote only $nonlin in ", path, "\n", sep = "")
