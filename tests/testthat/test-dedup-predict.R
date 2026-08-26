# Helpers for exact cell-deduplicated population predictions (spec 2026-08-26
# design A). Both helpers must be BYTE-identical to full-row prediction.

test_that("predict_forest_dedup equals row-wise predict on duplicated data", {
  skip_on_cran(); skip_if_not_installed("emdi")
  d <- tiny_saef_data()
  set.seed(11)
  rf <- ranger::ranger(y = d$Y, x = d$X, num.trees = 25)
  pop2 <- d$pop[rep(seq_len(nrow(d$pop)), each = 3), , drop = FALSE]
  pop2$cellkey <- rep(seq_len(nrow(d$pop)), each = 3)
  rownames(pop2) <- NULL
  full <- predict(rf, pop2)$predictions
  expect_identical(predict_forest_dedup(rf, pop2, "cellkey"), full)
  # NULL key = passthrough, exactly the old call
  expect_identical(predict_forest_dedup(rf, pop2, NULL), full)
})

test_that("predict_ranef_dedup equals predict.merMod incl. unseen levels", {
  skip_on_cran(); skip_if_not_installed("emdi")
  d <- tiny_saef_data()
  set.seed(12)
  d$smp$r <- rnorm(nrow(d$smp))
  fit <- lme4::lmer(r ~ (1 | district), data = d$smp, REML = TRUE)
  nd <- d$pop[rep(seq_len(nrow(d$pop)), each = 2), , drop = FALSE]
  nd$district <- as.character(nd$district)
  nd$district[1:2] <- "UNSEEN_LEVEL"   # exercises allow.new.levels = TRUE
  full <- predict(fit, newdata = nd, allow.new.levels = TRUE)
  dd <- predict_ranef_dedup(fit, nd, "district")
  expect_identical(unname(dd), unname(full))
})

# End-to-end: dedup_by must be BYTE-identical through the full pipeline, on
# both the production transformation paths, with var.adjust + MSE on.
dedup_e2e_fixture <- function() {
  d <- tiny_saef_data()
  pop2 <- d$pop[rep(seq_len(nrow(d$pop)), each = 3), , drop = FALSE]
  pop2$cellkey <- rep(seq_len(nrow(d$pop)), each = 3)
  rownames(pop2) <- NULL
  d$pop2 <- pop2
  d
}

test_that("dedup_by is byte-identical end-to-end (none + log transformations)", {
  skip_on_cran(); skip_if_not_installed("emdi")
  d <- dedup_e2e_fixture()
  for (tr in c("none", "log")) {
    run <- function(dd) SAEforest_model(Y = d$Y, X = d$X, dName = d$dName,
      smp_data = d$smp, pop_data = d$pop2, meanOnly = FALSE,
      MSE = "nonparametric", B = 2, num.trees = 50, threshold = median(d$Y),
      var.adjust = TRUE, B_adj = 5, adj_tol = 0.05, transformation = tr,
      seed = 1, dedup_by = dd)
    base <- run(NULL); dedup <- run("cellkey")
    expect_identical(base$Indicators, dedup$Indicators)
    expect_identical(base$MSE_Estimates, dedup$MSE_Estimates)
  }
})

test_that("dedup_by errors loudly when a covariate varies within a level", {
  skip_on_cran(); skip_if_not_installed("emdi")
  d <- dedup_e2e_fixture()
  d$pop2$cash[1] <- d$pop2$cash[1] + 1   # break constancy in level 1 only
  expect_error(
    SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
      pop_data = d$pop2, meanOnly = FALSE, MSE = "none", num.trees = 25,
      threshold = median(d$Y), seed = 1, dedup_by = "cellkey"),
    "vary within")
})

test_that("dedup_by naming a missing column errors loudly", {
  skip_on_cran(); skip_if_not_installed("emdi")
  d <- dedup_e2e_fixture()
  expect_error(
    SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
      pop_data = d$pop2, meanOnly = FALSE, MSE = "none", num.trees = 25,
      threshold = median(d$Y), seed = 1, dedup_by = "no_such_column"),
    "single column of pop_data")
})

test_that("dedup_by is refused on meanOnly and aggData paths", {
  skip_on_cran(); skip_if_not_installed("emdi")
  d <- dedup_e2e_fixture()
  expect_error(
    SAEforest_model(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
      pop_data = d$pop2, meanOnly = TRUE, num.trees = 25, seed = 1,
      dedup_by = "cellkey"),
    "unit-level non-linear")
})
