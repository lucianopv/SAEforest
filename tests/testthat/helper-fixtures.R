# Tiny deterministic subset of eusilcA for fast, reproducible tests.
tiny_saef_data <- function(n_dom = 4, per_dom = 8, seed = 2026) {
  set.seed(seed)
  data("eusilcA_smp", package = "emdi", envir = environment())
  data("eusilcA_pop", package = "emdi", envir = environment())
  doms <- sort(unique(as.character(eusilcA_smp$district)))[seq_len(n_dom)]

  smp <- do.call(rbind, lapply(doms, function(d) {
    rows <- which(eusilcA_smp$district == d)
    eusilcA_smp[sample(rows, min(per_dom, length(rows))), ]
  }))
  pop <- eusilcA_pop[eusilcA_pop$district %in% doms, ]
  # cap population per domain to keep smearing cheap
  pop <- do.call(rbind, lapply(doms, function(d) {
    rows <- which(pop$district == d)
    pop[rows[seq_len(min(40, length(rows)))], ]
  }))

  smp <- droplevels(smp); pop <- droplevels(pop)
  list(
    Y = smp$eqIncome,
    X = smp[, c("gender", "eqsize", "cash", "self_empl", "unempl_ben")],
    smp = smp[, c("district", "gender", "eqsize", "cash",
                  "self_empl", "unempl_ben", "eqIncome")],
    pop = pop[, c("district", "gender", "eqsize", "cash",
                  "self_empl", "unempl_ben")],
    dName = "district"
  )
}
