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

saveRDS(list(nonlin = r1, mean = r2),
        testthat::test_path("fixtures", "mse_memory_parity_ref.rds"))
