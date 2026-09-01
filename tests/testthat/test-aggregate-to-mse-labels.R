# Regression tests for the aggregate_to MSE label alignment.
#
# Aggregated MSE values are produced in tapply()/table() level order, but the
# label column used to be read from unique(pop_data[col]) -- first-appearance
# order. pop_data is sorted by dName before aggregation, so the two orders
# agree when aggregate_to is dName itself, or a coarser grouping that the dName
# sort happens to order. They do NOT agree when aggregate_to is finer than
# dName (the cell-level-output direction), and every MSE row was then stapled
# to another domain's value -- with sortAlpha() sorting on the wrong label
# afterwards, so the shipped table looked perfectly ordinary.

test_that("labelsInValueOrder returns labels in tapply/table level order", {
  pop <- data.frame(
    P = rep(c("P1", "P2"), each = 4),
    D = c("Dc", "Dc", "Dd", "Dd", "Da", "Da", "Db", "Db"),
    stringsAsFactors = FALSE)
  # first appearance is Dc Dd Da Db; the values are computed in sorted order
  expect_equal(as.character(unique(pop$D)), c("Dc", "Dd", "Da", "Db"))
  expect_equal(as.character(labelsInValueOrder(pop, "D")$D),
               names(table(pop$D)))
  expect_equal(as.character(labelsInValueOrder(pop, "D")$D),
               c("Da", "Db", "Dc", "Dd"))
})

test_that("labelsInValueOrder is the identity when the orders already agree", {
  # The pre-existing paths: aggregate_to == dName, and a coarser aggregate_to
  # that the dName sort already orders. Both must come back byte-identical, so
  # no previously-correct result changes.
  pop <- data.frame(
    D = rep(c("Da", "Db", "Dc"), each = 3),
    Region = rep(c("R1", "R1", "R2"), each = 3),
    stringsAsFactors = FALSE)
  for (col in c("D", "Region")) {
    plain <- unique(pop[col]); rownames(plain) <- NULL
    expect_identical(labelsInValueOrder(pop, col), plain, info = col)
  }
})

test_that("labelsInValueOrder follows factor level order, not alphabetical", {
  # tapply()/table() order a factor by its LEVELS, so the labels must too --
  # forcing alphabetical order here would break the dName path for a factor
  # domain whose levels are not alphabetical.
  pop <- data.frame(D = factor(c("b", "a", "c"), levels = c("c", "b", "a")),
                    stringsAsFactors = FALSE)
  expect_equal(as.character(labelsInValueOrder(pop, "D")$D), c("c", "b", "a"))
  expect_s3_class(labelsInValueOrder(pop, "D")$D, "factor")
})

test_that("nonparametric MSE attaches each value to its own aggregate_to domain", {
  # End-to-end pairing check. The bootstrap generates errors from ONE global
  # error distribution, so per-domain noise does not propagate -- domain
  # POPULATION SIZE does: a domain's true bootstrap mean has variance
  # sigma_e^2 / N_i, so a tiny domain must carry a far larger MSE than a large
  # one. District "Da" is given a tiny population and the rest large ones.
  #
  # The districts are laid out so that sorting pop_data by dName ("P") leaves
  # them in a NON-sorted first-appearance order -- the configuration that used
  # to permute the table, sending Da's MSE to whichever district appeared first.
  skip_on_cran()
  set.seed(11)
  prov <- paste0("P", 1:4)
  dist <- c("Dg", "Dh", "De", "Df", "Dc", "Dd", "Da", "Db")   # backwards vs prov
  tiny <- "Da"
  N_i  <- ifelse(dist == tiny, 6L, 1500L)
  map  <- data.frame(P = rep(prov, each = 2), D = dist, stringsAsFactors = FALSE)
  pop  <- map[rep(seq_len(nrow(map)), times = N_i), ]; rownames(pop) <- NULL
  expect_false(identical(as.character(unique(pop$D)), sort(unique(pop$D))))

  pop$x1 <- rnorm(nrow(pop)); pop$x2 <- rnorm(nrow(pop))
  u <- setNames(rnorm(length(prov), sd = 0.2), prov)
  pop$y <- 12 + 1.5 * pop$x1 - 0.8 * pop$x2 + u[pop$P] + rnorm(nrow(pop), sd = 3)
  smp <- do.call(rbind, lapply(split(pop, pop$P), function(g) g[sample(nrow(g), 60), ]))
  rownames(smp) <- NULL

  fit <- suppressWarnings(SAEforest_model(
    Y = smp$y, X = smp[, c("x1", "x2")], dName = "P",
    smp_data = smp, pop_data = pop, aggregate_to = "D",
    meanOnly = FALSE, MSE = "nonparametric", B = 8, B_adj = 5, num.trees = 30,
    var.adjust = TRUE, threshold = 12, select.indicator = c("Mean", "Hcr"),
    seed = 42, na.rm = TRUE))

  mse <- fit$MSE_Estimates
  expect_setequal(as.character(mse$D), sort(unique(pop$D)))
  # the tiny domain must own the largest MSE, by a wide margin -- not a coin flip
  expect_equal(as.character(mse$D[which.max(mse$Mean)]), tiny)
  expect_gt(max(mse$Mean), 10 * max(mse$Mean[as.character(mse$D) != tiny]))
})
