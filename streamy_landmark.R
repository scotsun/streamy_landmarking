#' Streamy Landmarking
#' 
#' @author Scott Sun
#' 
B <- 100
BATCH_SIZE <- 50

ipl_star_formula <- ~ x * (poly(LM/7, degree = 1, raw = TRUE) + I(exp(LM/7) - 1))
p <-  ncol(model.matrix(ipl_star_formula, data = data.frame(x = NA, LM = NA))[,-1])
batch_ref <- numeric(p)
betas <- numeric(p)
cum_nrow_superdata <- 0
for (b in seq_len(B)) {
  batch_data <- gen_df(BATCH_SIZE) # or read from file
  batch_superdata <- gen_superdata(
    dat = batch_data,
    w = 3,
    s = seq(0, 7, 0.2),
    timevar = "t_c", eventvar = "delta", fixed_x = "x"
  )
  # adaptive mean-centering
  batch_superX <- model.matrix(ipl_star_formula, data = batch_superdata)[,-1]
  p <- ncol(batch_superX)
  for (k in seq_len(p)) {
    if (all(batch_superX[,k] %in% c(0,1))) {
      next
    } else {
      batch_ref[k] <- (batch_ref[k] * cum_nrow_superdata + sum(batch_superX[,k])) / (cum_nrow_superdata + nrow(batch_superX))
    }
  }
  cum_nrow_superdata <- cum_nrow_superdata + nrow(batch_superX)
  batch_superX <- sweep(batch_superX, 1, batch_ref, "-")
  # calculate gradient
  
}

