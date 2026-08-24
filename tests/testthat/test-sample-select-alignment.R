# Regression suite for the sample_select() domain/size misalignment.
# See GridSAE docs/superpowers/specs/2026-08-24-sample-select-fix-and-hpc-efficiency-design.md

drawn_counts <- function(x, dName, doms) {
  as.integer(table(factor(as.character(x[[dName]]), levels = doms)))
}

test_that("the unsorted fixture really does diverge (guards the test itself)", {
  d <- tiny_saef_data_unsorted()
  first_appearance <- as.character(unique(d$pop[[d$dName]]))
  level_order <- names(split(seq_len(nrow(d$pop)), d$pop[[d$dName]]))
  expect_false(identical(first_appearance, level_order))
  # and the sample sizes must not all be equal, or a permutation is invisible
  expect_gt(length(unique(as.integer(table(d$smp[[d$dName]])))), 1L)
})

test_that("sample_select draws each domain's OWN sample size when pop is unsorted", {
  d <- tiny_saef_data_unsorted()
  doms <- sort(unique(as.character(d$smp[[d$dName]])))

  set.seed(11)
  got <- SAEforest:::sample_select(d$pop, smp = d$smp, dName = d$dName)

  expect_equal(drawn_counts(got, d$dName, doms),
               drawn_counts(d$smp, d$dName, doms))
})

test_that("a population domain with no sampled units draws nothing", {
  d <- tiny_saef_data_unsorted()
  doms <- sort(unique(as.character(d$smp[[d$dName]])))
  drop_dom <- doms[[1]]
  d$smp <- d$smp[as.character(d$smp[[d$dName]]) != drop_dom, , drop = FALSE]

  set.seed(12)
  got <- SAEforest:::sample_select(d$pop, smp = d$smp, dName = d$dName)

  expect_equal(sum(as.character(got[[d$dName]]) == drop_dom), 0L)
  expect_equal(nrow(got), nrow(d$smp))
})
