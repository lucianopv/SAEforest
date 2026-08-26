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
