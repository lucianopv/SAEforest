# var.adjust predictions leave out the effect model's intercept (Task 13a).
#
# On the var.adjust path MERFranger fits lmer(r_ij_new ~ (1 | area)), which has
# a fixed intercept, and trains the forest on Target - (predict - fixef): the
# area effects alone. A unit prediction is therefore the forest's plus the
# area's BLUP, as the authors' code computes it
# (_attic/SAEforest_files/pointestim_wrapper.R: predict(EffectModel, ...) -
# fixef(EffectModel)). Keeping the intercept shifts every prediction by a
# replicate-random constant. make_toy_unbalanced() makes that intercept clearly
# non-zero; with balanced areas it is zero up to rounding and nothing here could
# fail. The upstream model (Target ~ -1 + (1 | area)) has no intercept, and its
# outputs are held byte-identical to a snapshot from e261df1.

va_routes <- list(
  weights = list(adjust_method = "weights", B_adj = 10),
  refit   = list(B_adj = 5, adj_tol = 0)
)

fit_va <- function(d, route) {
  do.call(MERFranger, c(list(Y = d$Y, X = d$X, random = "(1|area)",
                             data = d$smp, var.adjust = TRUE, num.trees = 50,
                             num.threads = 1, seed = 3), route))
}

# The area's BLUP for every row of `data`; 0 for an area the fit never saw.
blup_of <- function(m, data, dName = "area") {
  re <- lme4::ranef(m$EffectModel)[[dName]]
  u <- re[as.character(data[[dName]]), 1]
  u[is.na(u)] <- 0
  u
}

# What every point estimate is built from: forest prediction plus area BLUP.
units_of <- function(m, data, dName = "area") {
  predict(m$Forest, data)$predictions + blup_of(m, data, dName)
}

expect_intercept_nonzero <- function(m) {
  b <- lme4::fixef(m$EffectModel)
  expect_length(b, 1L)
  expect_gt(abs(b[[1]]), 1e-3)          # vacuous if the intercept were ~0
}

point_args <- function(d, nonlin = TRUE) {
  a <- list(Y = d$Y, X = d$X, dName = d$dName, smp_data = d$smp,
            pop_data = d$pop, initialRandomEffects = 0, ErrorTolerance = 1e-4,
            MaxIterations = 25, aggregate_to = NULL, num.trees = 50,
            num.threads = 1, seed = 3)
  if (nonlin) a <- c(a, list(threshold = NULL, custom_indicator = NULL))
  a
}

test_that("ranef_part_ subtracts a fixed intercept and passes a model without one through", {
  d <- make_toy_unbalanced()
  nd <- d$pop                                  # includes the unseen area
  d$smp$r <- d$Y
  with_b <- lme4::lmer(r ~ (1 | area), data = d$smp, REML = TRUE)
  expect_gt(abs(lme4::fixef(with_b)[[1]]), 1e-3)
  got <- ranef_part_(with_b, newdata = nd, allow.new.levels = TRUE)
  expect_lt(max(abs(unname(got) - blup_of(list(EffectModel = with_b), nd))), 1e-10)

  # The upstream shape: no fixed effects at all, so not one operation is added.
  no_b <- lme4::lmer(r ~ -1 + (1 | area), data = d$smp, REML = FALSE)
  expect_length(lme4::fixef(no_b), 0L)
  expect_identical(ranef_part_(no_b, newdata = nd, allow.new.levels = TRUE),
                   predict(no_b, newdata = nd, allow.new.levels = TRUE))
  # positional newdata, as predict.MERFranger forwards it (in-sample areas,
  # since allow.new.levels is then FALSE)
  nd_in <- nd[nd$area != "new", , drop = FALSE]
  expect_identical(ranef_part_(no_b, nd_in), predict(no_b, nd_in))
})

test_that("predict.MERFranger on the var.adjust path is the forest plus the area's BLUP", {
  d <- make_toy_unbalanced()
  for (route in names(va_routes)) {
    m <- fit_va(d, va_routes[[route]])
    expect_intercept_nonzero(m)
    got <- unname(predict(m, d$pop))
    expect_lt(max(abs(got - units_of(m, d$pop))), 1e-10)
  }
})

test_that("point_nonLin's unit predictions are the forest plus the area's BLUP (with and without dedup_by)", {
  d <- make_toy_unbalanced()
  seen <- new.env()
  testthat::local_mocked_bindings(smear_indicators = function(preds, ...) {
    seen$preds[[length(seen$preds) + 1L]] <- preds
    c(Mean = mean(preds))
  })
  for (route in names(va_routes)) {
    for (dd in list(NULL, "cell")) {
      seen$preds <- list()
      out <- do.call(point_nonLin, c(point_args(d), va_routes[[route]],
                                     list(var.adjust = TRUE, dedup_by = dd)))
      m <- out[[2]]
      expect_intercept_nonzero(m)
      # one smear_indicators() call per domain, in table() order
      want <- split(units_of(m, d$pop), d$pop$area)[out[[1]]$area]
      expect_lt(max(abs(unlist(seen$preds) - unlist(want))), 1e-10)
    }
  }
})

test_that("point_MC_nonLin and point_mean build on the forest plus the area's BLUP", {
  d <- make_toy_unbalanced()
  dom_mean <- function(m) tapply(units_of(m, d$pop), d$pop$area, mean)
  # Zero bootstrap components: y_star is then exactly the unit predictions.
  testthat::local_mocked_bindings(ran_comp = function(mod, smp_data, ...) {
    list(forest_res = c(0, 0), ran_effs = c(0, 0), smp_data = smp_data)
  })
  for (route in names(va_routes)) {
    set.seed(5)
    mc <- do.call(point_MC_nonLin, c(point_args(d), va_routes[[route]],
                                     list(var.adjust = TRUE, B_point = 2)))
    expect_intercept_nonzero(mc$model)
    want <- dom_mean(mc$model)[as.character(mc[[1]]$area)]
    expect_lt(max(abs(mc[[1]]$Mean - want)), 1e-10)

    # point_mean is never given var.adjust by SAEforest_model (meanOnly
    # ignores it), but it forwards its dots to MERFranger like the others.
    pm <- do.call(point_mean, c(point_args(d, nonlin = FALSE),
                                va_routes[[route]], list(var.adjust = TRUE)))
    expect_intercept_nonzero(pm[[2]])
    want <- dom_mean(pm[[2]])[as.character(pm[[1]]$area)]
    expect_lt(max(abs(pm[[1]]$Mean - want)), 1e-10)
  }
})

test_that("the wild MSE helpers use the area effects without the intercept", {
  d <- make_toy_unbalanced()
  for (route in names(va_routes)) {
    m <- fit_va(d, va_routes[[route]])
    expect_intercept_nonzero(m)
    u <- blup_of(m, d$smp)

    # ran_comp_wild: sample residuals net of the area's BLUP, and the BLUPs
    # themselves. Both are centred, so this holds with or without the
    # intercept; it pins the definition.
    rcw <- ran_comp_wild(mod = m, smp_data = d$smp, Y = d$Y, dName = "area",
                         ADJsd = m$ErrorSD)
    r <- d$Y - m$Forest$predictions - u
    expect_lt(max(abs(rcw$forest_res - (r - mean(r)))), 1e-10)
    expect_lt(max(abs(rcw$ran_effs - scale_to(unique(u), m$RanEffSD))), 1e-10)

    # wild_errors matches each population value to the nearest sample fitted
    # value. A sample unit's own forest-plus-BLUP value must find that unit;
    # with the intercept in the fitted values it finds a neighbour instead.
    x <- units_of(m, d$smp)
    idx <- as.numeric(seq_along(x))
    set.seed(1)
    eps <- wild_errors(x = x, mod = m, smp_data = d$smp, forest_res = idx)
    expect_identical(unname(abs(eps)), idx)
  }
})

test_that("var.adjust = FALSE outputs are unchanged at every effect-model site", {
  skip_on_cran()
  ref <- readRDS(testthat::test_path("fixtures", "upstream_effects_ref.rds"))
  expect_identical(upstream_effect_outputs(), ref)
})
