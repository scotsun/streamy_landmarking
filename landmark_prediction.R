#' Landmark dynamic prediction for survival
#' @author Scott Sun
#' @param
#' @param
#' @return
#' @examples

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

# basis function:
# scaled according to s and g(0) = 0 (e.x. log(s + 1) )

ipl_dat <- gen_superdata(dat = df, w = 3, s = seq(0, 7, 1), timevar = "t_c", eventvar = "delta", fixed_x = "x")
ipl_star <- coxph(Surv(t_c, delta) ~ x*(poly(LM/7, degree = 1, raw = TRUE) + I(exp(LM/7) - 1)),
                  data = ipl_dat, model = TRUE, method = "breslow",
                  x = TRUE)
ipl_star$bhazard <- baseline_hazard(ipl_dat, ipl_star, "t_c", "delta")

simulation_cond_surv_plot()
points(seq(0, 7, 0.5),
       landmarking_predict(
         ipl_star,
         data.frame(x = factor(1, levels = c(-1, 1)), LM = seq(0, 7, 0.5)),
         3,
         seq(0, 7, 0.5)
       ),
       pch = 19, col = "darkolivegreen3")
points(seq(0, 7, 0.5),
       landmarking_predict(
         ipl_star,
         data.frame(x = factor(-1, levels = c(-1, 1)), LM = seq(0, 7, 0.5)),
         3,
         seq(0, 7, 0.5)
       ),
       pch = 19, col = "coral")

