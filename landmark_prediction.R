#' Landmark dynamic prediction for survival
#' @author Scott Sun
#' @param
#' @param
#' @return
#' @examples
library(tidyverse)
library(survival)
library(dynpred)


gen_df <- function(n) {
  u <- runif(n, 0, 1)
  x <- sample(c(-1,1), size = n, prob = c(0.5, 0.5), replace = TRUE)
  # this makes the effect time-dependent we cannot factor out a exp(LP) where LP is independent of t
  t <- 5 * (-log(u))^(1/exp(0.3*x))
  delta <- t <= 10
  t_c <- pmin(t, 10)
  
  df <- cbind(delta, t_c, x) %>% as.data.frame()
  df$x <- as.factor(x)
  return(df)
}

#' generate a super dataset
#' 
#' @examples 
#' gen_superdata(dat = df, w = 3, s = seq(0, 7, 1), timevar = "t_c", eventvar = "delta", fixed_x = "x")
gen_superdata <- function(dat, w, s, timevar, eventvar, fixed_x) {
  output <- lapply(s, function(t_lm)
    cutLM(
      dat,
      outcome = list(time = timevar, status = eventvar),
      LM = t_lm,
      horizon = t_lm + w,
      covs = list(fixed = fixed_x)
    ))
  return(bind_rows(output))
}


#' landmark CS prediction
#' 
#' @examples
#' newdata <- data.frame(x = factor(-1, levels = c(-1, 1)), LM = seq(0, 7))
#' landmarking_predict(ipl_star, newdata, 3, seq(0, 7))
landmarking_predict <- function(ipl_star_model, newdata, w, news) {
  L <- length(news)
  if (is.null(ipl_star_model[["bhazard"]])) {
    return(numeric(L))
  }
  # H0
  bhazard <- ipl_star_model[["bhazard"]]
  h0 <- bhazard[,"h0"]
  t <- bhazard[,"t"]
  H0_cond <- numeric(L)
  for (l in 1:L) {
    H0_cond[l] <- cumulative_baseline_hazard(news[l] + w, t, h0)
  }
  # prepare newdata
  newformula <- update(ipl_star_model$formula, NULL ~ .)
  newdata <- model.matrix(newformula, data = newdata)[,-1]
  newdata <- sweep(newdata, 2, ipl_star_model$means)
  # get LP and convert to CS
  lp <- newdata %*% ipl_star_model$coefficients
  cs <- exp(-exp(lp) * H0_cond)
  return(cs)
}


#' calculate baseline hazard rate
#' 
#' 
#' @examples baseline_hazard(ipl_dat, ipl_star, "t_c", "delta")
baseline_hazard <- function(superdata, ipl_star_model, timevar, eventvar) {
  eventtime_tbl <- data.frame(table(superdata[superdata[[eventvar]] == 1, timevar]))
  t <- as.numeric(levels(eventtime_tbl[, 1]))[eventtime_tbl[, 1]]
  d <- eventtime_tbl[, 2]
  h0 <- numeric(length(t))
  x_centered <- sweep(ipl_star_model$x, 2, ipl_star_model$means, "-")
  for (i in 1:length(t)) {
    h0[i] <- d[i] / sum(exp(x_centered %*% ipl_star_model$coefficients)[superdata[[timevar]] >= t[i]])
  }
  return(cbind(h0, t))
}

#' H0
#' 
#' @param
#' @param
#' @return
#' @examples
cumulative_baseline_hazard <- function(t, observed_event_time, hazard) {
  H0 <- sum(hazard[which(observed_event_time <= t)])
  return(H0)
}

#'utility plot wrapper
#'CS(3|s) in true values & standard cox predictions
simulation_cond_surv_plot <- function(df) {
  mod <- coxph(Surv(t_c, delta) ~ x, data = df, model = TRUE)
  
  t <- seq(0, 7, 0.05)
  CS_w_est1 <- summary(
    survfit(mod, newdata = data.frame(x = as.factor(1)), se.fit = FALSE), 
    times = t + 3
  )[["surv"]] / summary(
    survfit(mod, newdata = data.frame(x = as.factor(1)), se.fit = FALSE), 
    times = seq(0, 7, 0.05)
  )[["surv"]]
  CS_w_true1 <- exp(- ((t + 3)/5)^(exp(0.3)) + (t/5)^(exp(0.3)))
  
  CS_w_est2 <- summary(
    survfit(mod, newdata = data.frame(x = as.factor(-1)), se.fit = FALSE), 
    times = t + 3
  )[["surv"]] / summary(
    survfit(mod, newdata = data.frame(x = as.factor(-1)), se.fit = FALSE), 
    times = seq(0, 7, 0.05)
  )[["surv"]]
  CS_w_true2 <- exp(- ((t + 3)/5)^(exp(-0.3)) + (t/5)^(exp(-0.3)))
  
  plot(t, CS_w_est1, 
       col = "darkolivegreen3", 
       type = "l", lty = "dotted", lwd = 3,
       ylab = "CS(3|t)",
       xlab = "t (time)", ylim = c(0,1))
  lines(t, CS_w_true1, col = "chartreuse4", lwd = 3)
  lines(t, CS_w_est2, col = "coral", lty = "dotted", lwd = 3)
  lines(t, CS_w_true2, col = "coral3", lwd = 3)
  legend("topleft",
         legend = c("Case (X = -1)", "Control (X = 1)", "Predicted Case", "Predicted Control"),
         col = c("chartreuse4", "coral3", "darkolivegreen3", "coral"),
         lty = c("solid", "solid", "dotted", "dotted"),
         lwd = 3)
}


