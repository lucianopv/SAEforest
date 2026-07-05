# Golden reference for the var.adjust = TRUE pipeline. Run ONCE, after the
# var.adjust flag has been threaded end-to-end through SAEforest_model, to
# freeze the adjusted point estimates that test-regression-snapshot.R checks
# against. Unlike make_upstream_reference.R, this file is expected to be
# re-run whenever the variance-adjustment methodology intentionally changes.
# adj_tol is left unset below, so this captures output under the adaptive
# default (adj_tol = 0.05); the byte-identical legacy path (adj_tol = 0) is
# pinned separately in adj_nonLin_none_legacy.rds.
devtools::load_all(quiet = TRUE)
source(testthat::test_path("helper-fixtures.R"))
d <- tiny_saef_data()
set.seed(2026)
adj <- SAEforest_model(
  Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
  meanOnly = FALSE, MSE = "none", var.adjust = TRUE, B_adj = 2,
  num.trees = 50, mtry = 3
)
saveRDS(adj$Indicators,
        testthat::test_path("fixtures", "adj_nonLin_none.rds"))
