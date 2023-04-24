library(tidyverse)
library(survival)
library(dynpred)


set.seed(561)
n <- 5000
u <- runif(n, 0, 1)
x <- sample(c(-1,1), size = n, prob = c(0.5, 0.5), replace = TRUE)
# this makes the effect time-dependent we cannot factor out a exp(LP) where LP is independent of t
t <- 5 * (-log(u))^(1/exp(0.3*x))
delta <- t <= 10
t_c <- pmin(t, 10)

df <- cbind(delta, t_c, x) %>% as.data.frame()
df$x <- as.factor(x)

rm(u, x, t, delta, t_c)

mod <- coxph(Surv(t_c, delta) ~ x, data = df, model = TRUE)


# use "survfit", then "summary" to estimate probability
plot(survfit(mod, newdata = data.frame(x = as.factor(1)), se.fit = FALSE), 
     col = "darkolivegreen3", lty = "dotted", lwd = 3,
     ylab = "S(t)",
     xlab = "time",
     main = "Cox-PH model gives incorrect survival predictions with time-varying effects")
lines(sort(t_c), 
      exp(-(sort(t_c)/5)^(4/3)), col = "chartreuse4",
      lwd = 3)
lines(sort(t_c),
      summary(survfit(mod, newdata = data.frame(x = as.factor(-1)), se.fit = FALSE), times = sort(t_c))[["surv"]], 
      col = "coral", lty = "dotted", lwd = 3)
lines(sort(t_c), exp(-(sort(t_c)/5)^(3/4)), col = "coral3", lwd = 3)
legend(7, 1,
       legend = c("Case", "Control", "Predicted Case", "Predicted Control"),
       col = c("chartreuse4", "coral3", "darkolivegreen3", "coral"),
       lty = c("solid", "solid", "dotted", "dotted"),
       lwd = 3)

# conditional survival, time-window w = 3
# the simple/regular cox model does not pick up the dynamic variation at all

simulation_cond_surv_plot <- function() {
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
       ylab = "CS(3|s)",
       xlab = "s (landmarking time)", ylim = c(0,1))
  lines(t, CS_w_true1, col = "chartreuse4", lwd = 3)
  lines(t, CS_w_est2, col = "coral", lty = "dotted", lwd = 3)
  lines(t, CS_w_true2, col = "coral3", lwd = 3)
  legend(4.6, 1,
         legend = c("Case", "Control", "Predicted Case", "Predicted Control"),
         col = c("chartreuse4", "coral3", "darkolivegreen3", "coral"),
         lty = c("solid", "solid", "dotted", "dotted"),
         lwd = 3)
}






