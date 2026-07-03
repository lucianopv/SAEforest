# Auxiliary function for MSEs to perform bias correction ---------------------------------

adjust_ErrorSD <- function(Y, X, smp_data, mod, B = 100, ...) {
  pred_OOB <- matrix(mod$Forest$predictions, ncol = B,
                     nrow = length(mod$Forest$predictions), byrow = FALSE)

  e_ij <- Y - predict(mod$Forest, smp_data)$predictions
  e_ij <- e_ij - mean(e_ij)

  y_star_OOB <- pred_OOB + sample(e_ij, size = length(pred_OOB), replace = TRUE)

  my_estim_f2 <- function(x) {
    ranger::ranger(y = x, x = X, data = smp_data, ...)
  }
  my_f_n2 <- pbapply::pbapply(y_star_OOB, 2, my_estim_f2)

  my_pred_f <- function(x) {
    x$predictions
  }
  pred_OOB_star <- sapply(my_f_n2, my_pred_f)

  Adjustment <- (pred_OOB_star - pred_OOB)^2

  bias <- mean(rowMeans(Adjustment))
  if (mod$ErrorSD^2 < bias) {
    warning("Variance bias correction exceeds ErrorSD^2; clamping adjusted SD to 0.")
  }
  outvar <- sqrt(max(0, mod$ErrorSD^2 - bias))

  return(outvar)
}

# Auxiliary function for MERFranger's per-iteration variance bias correction.
# Returns the raw bias term K = mean over units of the mean squared gap between
# bootstrap-refit OOB predictions and the original OOB predictions (Mendez & Lohr, 2011).
adjust_ErrorSD_ <- function(Y, X, smp_data, rf, B = 100, ...) {
  # No seed argument: this runs inside MERFranger's iteration (and inside the MSE
  # bootstrap). It inherits the RNG stream seeded once at the top of SAEforest_model,
  # so every refit is reproducible run-to-run yet distinct across bootstrap replicates.
  pred_OOB <- matrix(rf$predictions, ncol = B,
                     nrow = length(rf$predictions), byrow = FALSE)

  e_ij <- Y - rf$predictions
  e_ij <- e_ij - mean(e_ij)

  y_star_OOB <- pred_OOB + sample(e_ij, size = length(pred_OOB), replace = TRUE)

  my_estim_f2 <- function(x) {
    ranger::ranger(y = x, x = X, data = smp_data, ...)
  }
  my_f_n2 <- apply(y_star_OOB, 2, my_estim_f2)
  pred_OOB_star <- sapply(my_f_n2, function(x) x$predictions)

  mean(rowMeans((pred_OOB_star - pred_OOB)^2))
}
