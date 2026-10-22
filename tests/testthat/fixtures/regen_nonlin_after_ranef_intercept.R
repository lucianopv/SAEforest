# One-off, 2026-10-22 (Task 13a): regenerate ONLY the `nonlin` element of
# mse_memory_parity_ref.rds after var.adjust predictions stopped carrying the
# effect model's intercept (ranef_part_() in support_functions.R).
#
# Why the golden moved. On the var.adjust path the effect model is
# lmer(r_ij_new ~ (1 | area)); its fixed intercept used to be added to every
# unit prediction (point estimates, the MSE replicates' refits, the wild
# fitted values). It is now left out, as in the authors' code. The parity
# golden is held at tolerance = 0, so any arithmetic change breaks it.
#
# Why this is safe to accept. tiny_saef_data() has 4 domains of exactly 8
# observations, and r_ij_new is centred, so the REML intercept is the grand mean
# of a centred vector: 1.2e-12 on an income scale of 1e4. Subtracting it is a
# no-op up to rounding (max relative change 7.6e-15 when this was written). The
# fixture could not have caught the bug; test-ranef-intercept.R covers it with
# unbalanced areas whose intercept is clearly non-zero.
#
# Why ONLY `nonlin` (as regen_nonlin_after_pool_centring.R): `mean` runs
# var.adjust = FALSE and must not move; `nonlin_agg` already fails for a golden
# row-order reason and is deliberately left frozen so the suite keeps
# flagging it.

devtools::load_all(quiet = TRUE)
source(testthat::test_path("helper-fixtures.R"))

path <- testthat::test_path("fixtures", "mse_memory_parity_ref.rds")
old <- readRDS(path)
stopifnot(identical(sort(names(old)), c("mean", "nonlin", "nonlin_agg")))

d <- tiny_saef_data()
# Identical call to make_mse_memory_parity_reference.R's r1.
r1 <- SAEforest_model(
  Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
  meanOnly = FALSE, MSE = "nonparametric", B = 3,
  threshold = median(d$Y), var.adjust = TRUE, B_adj = 5, adj_tol = 0.05,
  num.trees = 50, mtry = 3, seed = 5
)$MSE_Estimates

# Refuse to write if the change is larger than floating-point noise.
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
