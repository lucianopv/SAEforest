# MSE bootstrap for nonlinear indicators under unit-level covariate data ------------------

MSE_SAEforest_nonLin <- function(Y,
                                 X,
                                 dName,
                                 threshold,
                                 smp_data,
                                 mod,
                                 ADJsd,
                                 pop_data,
                                 B = 100,
                                 B_point,
                                 initialRandomEffects = 0,
                                 ErrorTolerance = 0.0001,
                                 MaxIterations = 25,
                                 custom_indicator,
                                 wild,
                                 MC,
                                 aggregate_to,
                                 var.adjust = FALSE,
                                 B_adj = 100,
                                 adj_tol = 0,
                                 adj_tol_mse = NULL,
                                 mse_tol = 0,
                                 transformation = c("none", "log"),
                                 select.indicator = NULL,
                                 cores = 1,
                                 worker.threads = 1L,
                                 dedup_by = NULL,
                                 ...) {

  transformation <- match.arg(transformation, c("none", "log"))
  if (is.null(adj_tol_mse)) adj_tol_mse <- adj_tol

  rand_struc <- paste0(paste0("(1|", dName), ")")
  domains <- t(unique(pop_data[dName]))
  in_samp <- domains %in% t(unique(smp_data[dName]))
  N_i <- as.numeric(table(pop_data[[dName]]))

  pred_vals <- predict_forest_dedup(mod$Forest, pop_data, dedup_by)
  pred_mat <- matrix(pred_vals, nrow = length(pred_vals), ncol = B)

  # Prepare data for sampling
  if (wild == TRUE) {
    ran_obj <- ran_comp_wild(Y = Y, smp_data = smp_data, mod = mod, ADJsd = ADJsd, dName = dName)
    ran_effs <- ran_obj$ran_effs
    forest_res <- ran_obj$forest_res

    sample_e <- function(x) {
      wild_errors(x = x, mod = mod, smp_data = smp_data, forest_res = forest_res)
    }
  }

  if (wild == FALSE) {
    ran_obj <- ran_comp(Y = Y, smp_data = smp_data, mod = mod, ADJsd = ADJsd, dName = dName)
    ran_effs <- ran_obj$ran_effs
    forest_res <- ran_obj$forest_res
    smp_data <- ran_obj$smp_data

    sample_e <- function(x) {
      sample(forest_res, size = sum(N_i), replace = TRUE)
    }
  }

  sample_ui <- function(x) {
    rep(sample(ran_effs,
      size = length(N_i),
      replace = TRUE
    ), N_i)
  }

  u_i <- apply(pred_mat, 2, sample_ui)

  smp_data$forest_res <- NULL

  # combine to y_star
  mu_ij <- pred_mat + u_i
  # (the matrix(NA, ...) that used to sit here allocated a full N x B matrix --
  # ~2.74 GB at 6.85M rows x B=50 -- and was overwritten on the next line)
  e_ij <- apply(mu_ij, 2, sample_e)

  y_star <- mu_ij + e_ij

  if (transformation == "log") {
    y_star <- exp(y_star); y_star[!is.finite(y_star)] <- NA
  }

  # get tau_star
  y_star_L <- split(y_star, col(y_star))
  thresh_L <- lapply(y_star_L, function(x) get_thresh(x, threshold = threshold))

  # The threshold is a scalar per replicate. It used to be broadcast into a
  # second N-length column via Map(cbind, ...) and recovered with unique(x[, 2]),
  # which doubled the largest object in this function (an extra ~5.5 GB at
  # 6.85M rows x B=50) to carry one number. Passed as an argument instead.
  if (is.null(aggregate_to)) {
    my_agg <- function(y, th) {
      tapply(y, pop_data[[dName]], calc_indicat, threshold = th,
             custom = custom_indicator, select.indicator = select.indicator)
    }
  } else {
    my_agg <- function(y, th) {
      tapply(y, pop_data[[aggregate_to]], calc_indicat, threshold = th,
             custom = custom_indicator, select.indicator = select.indicator)
    }
  }

  tau_star <- Map(my_agg, y_star_L, thresh_L)

  if(is.null(aggregate_to)){
    comb <- function(x) {
    matrix(unlist(x), nrow = length(N_i), byrow = TRUE)
    }
  } else{
    N_i_agg <- as.numeric(table(pop_data[[aggregate_to]]))
    comb <- function(x) {
      matrix(unlist(x), nrow = length(N_i_agg), byrow = TRUE)
    }
  }

  tau_star <- sapply(tau_star, comb, simplify = FALSE)

  # get bootstrap samples
  boots_sample <- vector(mode = "list", length = B)

  # The domain column is invariant across replicates -- the loop below rewrites
  # only pop_data$y_star -- so the per-domain row index split is computed once
  # here instead of B times inside sample_select().
  grp <- split(seq_len(nrow(pop_data)), pop_data[[dName]])

  for (i in 1:B) {
    pop_data$y_star <- y_star[, i]
    boots_sample[[i]] <- sample_select(pop_data, smp = smp_data, dName = dName,
                                       groups = grp)
  }

  # Everything above is dead once boots_sample and tau_star exist, but all of it
  # stays reachable from this frame when with_parallel_rng() forks below. R's GC
  # marks every object in each forked worker, which writes to every page and
  # defeats copy-on-write -- so these ~20 GB (at 6.85M rows x B=50) are
  # materialised once per core. Freeing them here is what makes an 8-core task
  # fit in ~80 GB instead of ~350 GB.
  pop_data$y_star <- NULL
  rm(pred_mat, u_i, mu_ij, e_ij, y_star, y_star_L, thresh_L, grp)
  gc()

  # uses sample to estimate tau_b
  if (MC == TRUE) {
    dots <- force_serial_threads(cores, ..., worker.threads = worker.threads)
    my_estim_f <- function(x) {
      if (cores > 1L) pin_blas_threads()
      do.call(point_MC_nonLin, c(list(
        Y = x$y_star, X = x[, colnames(X)], dName = dName, threshold = threshold,
        smp_data = x, pop_data = pop_data, initialRandomEffects = initialRandomEffects,
        ErrorTolerance = ErrorTolerance, B_point = B_point, MaxIterations = MaxIterations,
        custom_indicator = custom_indicator, aggregate_to = aggregate_to,
        var.adjust = var.adjust, B_adj = B_adj, adj_tol = adj_tol_mse,
        select.indicator = select.indicator,
        dedup_by = dedup_by), dots))[[1]][, -1]
    }
  }

  if (MC == FALSE) {
    dots <- force_serial_threads(cores, ..., worker.threads = worker.threads)
    my_estim_f <- function(x) {
      if (cores > 1L) pin_blas_threads()
      do.call(point_nonLin, c(list(
        Y = x$y_star, X = x[, colnames(X)], dName = dName, threshold = threshold,
        smp_data = x, pop_data = pop_data, initialRandomEffects = initialRandomEffects,
        ErrorTolerance = ErrorTolerance, MaxIterations = MaxIterations,
        custom_indicator = custom_indicator, aggregate_to = aggregate_to,
        var.adjust = var.adjust, B_adj = B_adj, adj_tol = adj_tol_mse,
        transformation = transformation, select.indicator = select.indicator,
        dedup_by = dedup_by), dots))[[1]][, -1]
    }
  }

  mean_square <- function(x, y) {
    (x - y)^2
  }

  if (is.null(mse_tol)) mse_tol <- 0

  if (mse_tol <= 0) {
    ## ---- exact path: fixed B (byte-identical to previous behaviour) ----
    tau_b <- with_parallel_rng(cores, function(cl) {
      pbapply::pbsapply(boots_sample, my_estim_f, cl = cl, simplify = FALSE)
    })
    Mean_square_B <- mapply(mean_square, tau_b, tau_star, SIMPLIFY = FALSE)
    MSE_estimates <- Reduce("+", Mean_square_B) / length(Mean_square_B)
  } else {
    ## ---- adaptive early-stop: batched refits, stop on median relative MC-SE ----
    B_MIN   <- min(20L, B)
    B_BATCH <- min(max(as.integer(cores), 10L), B)
    res_mse <- with_parallel_rng(cores, function(cl) {
      sq <- vector("list", B)
      m <- 0L
      repeat {
        idx <- (m + 1L):min(m + B_BATCH, B)
        batch <- pbapply::pbsapply(boots_sample[idx], my_estim_f, cl = cl, simplify = FALSE)
        for (j in seq_along(idx)) sq[[idx[j]]] <- mean_square(batch[[j]], tau_star[[idx[j]]])
        m <- max(idx)
        if (m >= B_MIN && m >= 2L) {
          # sq[[i]] may be a data.frame (multi-indicator) or a matrix
          # (single-indicator, via select.indicator). simplify2array() on a
          # list of data.frames uses ncol() as "length" and silently drops the
          # domain dimension, so coerce to plain matrices first -- this always
          # yields a domains x indicators x m array regardless of shape.
          arr <- simplify2array(lapply(sq[seq_len(m)], as.matrix))
          mse <- apply(arr, c(1, 2), mean)
          sdv <- apply(arr, c(1, 2), stats::sd)
          rel <- (sdv / sqrt(m)) / mse
          rel <- rel[is.finite(rel) & mse > 0]
          if (length(rel) > 0 && stats::median(rel) < mse_tol) break
        }
        if (m >= B) break
      }
      list(sq = sq[seq_len(m)], m = m)
    })
    Mean_square_B <- res_mse$sq
    MSE_estimates <- Reduce("+", Mean_square_B) / length(Mean_square_B)
  }

  mse_col_names <- colnames(MSE_estimates)
  if (is.null(mse_col_names)) mse_col_names <- select.indicator

  if(is.null(aggregate_to)){
    MSE_estimates_out <- data.frame(unique(pop_data[dName]), MSE_estimates)
    colnames(MSE_estimates_out) <- c(dName, mse_col_names)
  } else{
    MSE_estimates_out <- data.frame(unique(pop_data[aggregate_to]), MSE_estimates)
    colnames(MSE_estimates_out) <- c(aggregate_to, mse_col_names)
  }

  rownames(MSE_estimates_out) <- NULL

  # __________________________
  return(MSE_estimates_out)
}
