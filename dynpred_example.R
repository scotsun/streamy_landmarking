library(tidyverse)
library(dynpred)

test0 <- data.frame(
  id = c(1, 1, 1, 2, 2, 2),
  survyrs = c(2.3, 2.3, 2.3, 2.7, 2.7, 2.7),
  survstat = c(1, 1, 1, 0, 0, 0),
  age = c(76, 76, 76, 68, 68, 68),
  gender = c(1, 1, 1, 2, 2, 2),
  bp = c(80, 84, 88, 92, 90, 89),
  bptime = c(1, 2, 2.2, 0, 1, 2),
  sbp = c(80, 84, 88, 92, 90, 89),
  sbptime = c(0.5, 1.5, 2, 0, 0.5, 2.2)
)
cutLM(test0, outcome=list(time="survyrs", status="survstat"),
      LM=1, horizon=2.5, covs=list(fixed=c("age","gender"),varying="bp"),
      format="long", id="id", rtime="bptime")

cutLM(test0, outcome=list(time="survyrs", status="survstat"),
      LM=1, horizon=2.5, covs=list(fixed=c("age","gender"),varying=c("bp", "sbp")),
      format="long", id="id", rtime=c("bptime"), right=FALSE)

test1 <-
  data.frame(
    id = 1:4,
    survyrs = c(7.6, 8.4, 5.3, 2.6),
    survstat = c(0, 1, 1, 0),
    age = c(48, 52, 76, 18),
    gender = c(1, 2, 2, 1),
    recyrs = c(7.6, 5.2, 0.8, 2.6),
    recstat = c(0, 1, 1, 0)
  )
cutLM(test1, outcome=list(time="survyrs", status="survstat"),
      LM=3, horizon=8, covs=list(fixed=c("id","age","gender"),varying="recyrs"))

test2 <-
  data.frame(
    id = c(1, 2, 2, 3, 3, 4),
    survyrs = c(7.6, 8.4, 8.4, 5.3, 5.3, 2.6),
    survstat = c(0, 1, 1, 1, 1, 0),
    age = c(48, 52, 52, 76, 76, 18),
    gender = c(1, 2, 2, 2, 2, 1),
    rec = c(0, 0, 1, 0, 1, 0),
    rectime = c(0, 0, 5.2, 0, 0.8, 0)
  )
cutLM(test2, outcome=list(time="survyrs", status="survstat"),
      LM=1, horizon=2, covs=list(fixed=c("age","gender"),varying="rec"),
      format="long", id="id", rtime="rectime")

## supermodel landmark ----------------------
# a single step
## --- --- --- --- --- --- --- --- --- --- --- --- --- 
library(survival)

# landmark sL in [0, 6], t_hor = 10

superdata1 <- lapply(seq(0, 6, 0.2), function(s)
  cutLM(
    df,
    outcome = list(time = "t_c", status = "delta"),
    LM = s,
    horizon = 10,
    covs = list(fixed = c("x"))
  )) 
superdata1 <- bind_rows(superdata1)
  
lm_mod1 <- coxph(Surv(t_c, delta) ~ x*poly(LM/7, degree = 5, raw = TRUE), data = superdata1, model = TRUE)
lm_mod1



