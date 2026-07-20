# Behaviour of the MSE random-component helpers when the model reports
# RanEffSD == 0 (a degenerate/singular random-effect fit).
#
# This is reachable in practice: a sufficiently flexible forest with
# area-identifying covariates drives the estimated random-effect variance to
# exactly zero, after which lme4 predicts one constant for every unit.
#
# An lm(y ~ 1) fit is used as the EffectModel stand-in wherever a degenerate
# random-effect model is needed -- predict() returns the same value for every
# row, which is exactly what a singular (1 | area) fit does.

mk_smp <- function(D = 6, per = 5, seed = 1) {
  set.seed(seed)
  data.frame(
    area = factor(rep(seq_len(D), each = per)),
    x = stats::rnorm(D * per),
    y = stats::rnorm(D * per, mean = 10, sd = 2)
  )
}

mk_mod <- function(smp, ran_eff_sd, degenerate_effects = TRUE) {
  rf <- ranger::ranger(y ~ x, data = smp, num.trees = 20, num.threads = 1)
  list(
    Forest = rf,
    EffectModel = if (degenerate_effects) stats::lm(y ~ 1, data = smp)
                  else stats::lm(y ~ x, data = smp),
    RanEffSD = ran_eff_sd,
    ErrorSD = 1
  )
}

test_that("ran_comp warns when the random-effect SD is zero", {
  smp <- mk_smp()
  mod <- mk_mod(smp, ran_eff_sd = 0)
  expect_warning(
    ran_comp(mod = mod, smp_data = smp, Y = smp$y, dName = "area", ADJsd = 1),
    "random.effect|between-area|RanEffSD"
  )
})

test_that("ran_comp does not warn when the random-effect SD is positive", {
  smp <- mk_smp()
  mod <- mk_mod(smp, ran_eff_sd = 0.5)
  expect_no_warning(
    ran_comp(mod = mod, smp_data = smp, Y = smp$y, dName = "area", ADJsd = 1)
  )
})

test_that("ran_comp still returns all-zero random effects when the SD is zero", {
  # the arithmetic is correct and must not change -- only the silence does
  smp <- mk_smp()
  mod <- mk_mod(smp, ran_eff_sd = 0)
  out <- suppressWarnings(
    ran_comp(mod = mod, smp_data = smp, Y = smp$y, dName = "area", ADJsd = 1)
  )
  expect_true(all(out$ran_effs == 0))
  expect_false(anyNA(out$ran_effs))
})

test_that("ran_comp_wild returns finite random effects when the SD is zero", {
  # regression: unique() collapses the constant predictions to length 1, and
  # sd() of a length-1 vector is NA, which propagated all the way into
  # wild_errors() and surfaced as an unintelligible vapply error.
  smp <- mk_smp()
  mod <- mk_mod(smp, ran_eff_sd = 0, degenerate_effects = TRUE)
  # it warns for the same reason ran_comp does
  expect_warning(
    ran_comp_wild(mod = mod, smp_data = smp, Y = smp$y, dName = "area", ADJsd = 1),
    "random.effect|between-area|RanEffSD"
  )
  out <- suppressWarnings(
    ran_comp_wild(mod = mod, smp_data = smp, Y = smp$y, dName = "area", ADJsd = 1)
  )
  expect_false(anyNA(out$ran_effs))
  expect_true(all(out$ran_effs == 0))
})

test_that("ran_comp_wild is numerically unchanged when the SD is positive", {
  # guards the fix against altering the non-degenerate path
  smp <- mk_smp()
  mod <- mk_mod(smp, ran_eff_sd = 0.7, degenerate_effects = FALSE)
  out <- ran_comp_wild(mod = mod, smp_data = smp, Y = smp$y,
                       dName = "area", ADJsd = 1)
  # the original formula, verbatim
  ref <- unique(stats::predict(mod$EffectModel, smp))
  ref <- (ref / stats::sd(ref)) * mod$RanEffSD
  ref <- ref - mean(ref)
  expect_equal(out$ran_effs, ref)
})

test_that("MSE = 'wild' completes rather than erroring when the forest absorbs the area effects", {
  skip_on_cran()
  # area-constant covariates let the forest recover area identity, which drives
  # RanEffSD to exactly zero -- the condition that used to abort the wild path.
  set.seed(11)
  D <- 10; per_pop <- 20; per_smp <- 12; n_z <- 10
  area <- rep(seq_len(D), each = per_pop)
  z <- lapply(seq_len(n_z), function(k) stats::runif(D)[area])
  names(z) <- paste0("z", seq_len(n_z))
  pop <- data.frame(area = factor(area), x1 = stats::rnorm(D * per_pop))
  for (nm in names(z)) pop[[nm]] <- z[[nm]]
  pop$y <- 2 * pop$x1 + stats::rnorm(D, sd = 0.6)[area] + stats::rnorm(D * per_pop)
  idx <- unlist(lapply(split(seq_len(nrow(pop)), area),
                       function(r) sample(r, per_smp)), use.names = FALSE)
  smp <- pop[sort(idx), ]
  covs <- c("x1", names(z))

  fit <- suppressWarnings(suppressMessages(
    SAEforest_model(Y = smp$y, X = smp[, covs], dName = "area",
                    smp_data = smp[, c("area", covs, "y")],
                    pop_data = pop[, c("area", covs)],
                    meanOnly = FALSE, MSE = "wild", B = 2,
                    num.trees = 50, num.threads = 1,
                    select.indicator = c("Mean", "Hcr"))
  ))
  expect_equal(fit$MERFmodel$RanEffSD, 0)
  expect_s3_class(fit$MSE_Estimates, "data.frame")
  expect_gt(nrow(fit$MSE_Estimates), 0)
})
