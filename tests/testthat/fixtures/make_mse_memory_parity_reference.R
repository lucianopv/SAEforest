# Golden reference for the two numerically inert MSE-bootstrap edits: the
# `groups` hoist (Task 3, hoisting the per-domain row split out of the
# `for (i in 1:B)` loop in MSE_SAEforest_nonLin.R / MSE_SAEforest_mean.R) and
# the pseudo-population memory edits (Task 4). Run ONCE, on the code as it
# stood after Task 2's sample_select() fix and BEFORE either edit, to freeze
# MSE_Estimates that both tasks are held to reproduce exactly.
devtools::load_all(quiet = TRUE)
source(testthat::test_path("helper-fixtures.R"))
d <- tiny_saef_data()

r1 <- SAEforest_model(
  Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
  meanOnly = FALSE, MSE = "nonparametric", B = 3,
  threshold = median(d$Y), var.adjust = TRUE, B_adj = 5, adj_tol = 0.05,
  num.trees = 50, mtry = 3, seed = 5
)$MSE_Estimates

r2 <- SAEforest_model(
  Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp, pop_data = d$pop,
  meanOnly = TRUE, MSE = "nonparametric", B = 3,
  num.trees = 50, mtry = 3, seed = 5
)$MSE_Estimates

# r3: the `aggregate_to` branch of `my_agg` in MSE_SAEforest_nonLin.R. This is
# the branch every GridSAE fit actually exercises (dName finer than
# aggregate_to, e.g. dName="subunit"/aggregate_to="district"), but r1/r2 above
# only cover aggregate_to = NULL. Final review 2026-08-24 (finding M2):
# verified against commit 00c553c's version of MSE_SAEforest_nonLin.R (the
# state before the `groups` hoist and the memory edits) via a throwaway
# worktree -- identical at tolerance = 0, so capturing this golden at current
# HEAD is valid.
d3 <- tiny_saef_data_with_subunit()

r3 <- SAEforest_model(
  Y = d3$Y, X = d3$X, dName = "subunit", smp_data = d3$smp, pop_data = d3$pop,
  meanOnly = FALSE, MSE = "nonparametric", B = 3,
  threshold = median(d3$Y), var.adjust = TRUE, B_adj = 5, adj_tol = 0.05,
  num.trees = 50, mtry = 3, seed = 5,
  aggregate_to = "district", select.indicator = c("Mean", "Hcr")
)$MSE_Estimates

saveRDS(list(nonlin = r1, mean = r2, nonlin_agg = r3),
        testthat::test_path("fixtures", "mse_memory_parity_ref.rds"))
