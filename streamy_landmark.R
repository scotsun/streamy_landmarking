#' Streamy Landmarking
#' 
#' computation efficiency is bad; solve memory issue
#' @author Scott Sun
library(progress)
source("./landmark_prediction.R")

B <- 500
BATCH_SIZE <- 100

ipl_star_formula <- ~ x * (poly(LM/7, degree = 1, raw = TRUE) + I(exp(LM/7) - 1))
p <-  ncol(model.matrix(ipl_star_formula, data = data.frame(x = NA, LM = NA))[,-1])
batch_ref <- numeric(p)

betas <- prev_betas <- numeric(p)

cum_nrow_superdata <- 0

BETAS <- matrix(NA, nrow = B, ncol = p)
colnames(BETAS) <- colnames(model.matrix(ipl_star_formula, data = data.frame(x = 1, LM = NA))[,-1])
ALPHAS <- numeric(B)
pb <- progress_bar$new(
  format = "  downloading [:bar] :percent eta: :eta",
  total = B, clear = FALSE, width= 60)
for (b in seq_len(B)) {
  pb$tick()
  batch_data <- gen_df(BATCH_SIZE) # or read from file
  batch_superdata <- gen_superdata(
    dat = batch_data,
    w = 3,
    s = seq(0, 7, 0.2),
    timevar = "t_c", eventvar = "delta", fixed_x = "x"
  )
  # adaptive mean-centering (useful when calculate streaming h0's)
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
  
  # calculate gradient
  grad <- prev_grad <- numeric(p)
  for (i in seq_len(nrow(batch_superX))) {
    Xbar_numl <- batch_superX * (batch_superdata$t_c >= batch_superdata$t_c[i]) * exp(drop(batch_superX %*% betas))
    Xbar_denl <- (batch_superdata$t_c >= batch_superdata$t_c[i]) * exp(drop(batch_superX %*% betas))
    grad <- grad - batch_superdata$delta[i] * (batch_superX[i,] - apply(Xbar_numl, 2, sum) / sum(Xbar_denl))
  }
  # update
  step_size <- 2e-3
  betas <- betas - step_size * grad
  BETAS[b,] <- betas; ALPHAS[b] <- step_size
}

BETAS <- data.frame(BETAS) 
BETASlong <- BETAS %>% 
  mutate(iter = 1:B) %>% 
  pivot_longer(cols = colnames(BETAS)[1:5])

ggplot(data = BETASlong,
       aes(x = iter, y = value, color = name)) +
  geom_line()

