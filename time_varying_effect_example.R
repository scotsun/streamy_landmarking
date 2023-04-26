library(tidyverse)
library(survival)
library(dynpred)
source("./landmark_prediction.R")

set.seed(561)
df <- lapply(1:10, function(i) gen_df(500)) %>% 
  bind_rows()
# the density of the landmark grid will not significantly affect the prediction
ipl_dat <- gen_superdata(dat = df, w = 3, s = seq(0, 7, 1), timevar = "t_c", eventvar = "delta", fixed_x = "x")
# basis function:
# scaled according to s and g(0) = 0 (e.x. log(s + 1), exp(s) - 1)
ipl_star <- coxph(Surv(t_c, delta) ~ x*(poly(LM/7, degree = 1, raw = TRUE) + I(exp(LM/7) - 1)),
                  data = ipl_dat, model = TRUE, method = "breslow",
                  x = TRUE)
ipl_star$bhazard <- baseline_hazard(ipl_dat, ipl_star, "t_c", "delta")
# conditional survival, time-window w = 3
# the simple/regular cox model does not pick up the dynamic variation at all
simulation_cond_surv_plot(df)
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

# use "survfit", then "summary" to estimate probability
mod <- coxph(Surv(t_c, delta) ~ x, data = df, model = TRUE)
plot(survfit(mod, newdata = data.frame(x = as.factor(1)), se.fit = FALSE), 
     col = "darkolivegreen3", lty = "dotted", lwd = 3,
     ylab = "S(t)",
     xlab = "t (time)")
lines(sort(df$t_c), 
      exp(-(sort(df$t_c)/5)^(4/3)), col = "chartreuse4",
      lwd = 3)
lines(sort(df$t_c),
      summary(survfit(mod, newdata = data.frame(x = as.factor(-1)), se.fit = FALSE), times = sort(df$t_c))[["surv"]], 
      col = "coral", lty = "dotted", lwd = 3)
lines(sort(df$t_c), exp(-(sort(df$t_c)/5)^(3/4)), col = "coral3", lwd = 3)
legend("topright",
       legend = c("Case (X = -1)", "Control (X = 1)", "Predicted Case", "Predicted Control"),
       col = c("chartreuse4", "coral3", "darkolivegreen3", "coral"),
       lty = c("solid", "solid", "dotted", "dotted"),
       lwd = 3)






