# It is necessary to have the simulated dataset saved in the directory
setwd("C:/Users/albaf/OneDrive/Escritorio/Prácticas IEO/Simulación/GitHub_Simulation")

# Then you could load the data
load("C:/Users/albaf/OneDrive/Escritorio/Prácticas IEO/Simulación/GitHub_Simulation/dataframes.RData")

# Transform response variable
### Abundance
dataset1$abut <- dataset1$abu + 0.1
### CPUE
dataset1$CPUE <- dataset1$abu / dataset1$eff
dataset1$CPUEt <- dataset1$abut / dataset1$eff

# MODELS WITHOUT ZEROS

## Model 1: Biomass + eff + bath (Gamma)

for.mod1 <-
  abut ~ eff + bath # formula (response variable and lineal predictor)

set.seed(12345)
mod1 <-
  inla(
    for.mod1,
    family = "gamma",
    data = dataset1,
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )

## Model 2: Biomass + eff + bath (ln Normal)

for.mod2 <-
  abut ~ eff + bath # formula (response variable and lineal predictor)

set.seed(12345)
mod2 <-
  inla(
    for.mod1,
    family = "lognormal",
    data = dataset1,
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )

## Model 3: CPUE + bath (Gamma)

for.mod3 <-
  CPUEt ~ bath # formula (response variable and lineal predictor)

set.seed(12345)
mod3 <-
  inla(
    for.mod3,
    family = "gamma",
    data = dataset1,
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )

## Model 4: CPUE + bath (ln Normal)

for.mod4 <-
  CPUEt ~ bath # formula (response variable and lineal predictor)

set.seed(12345)
mod4 <-
  inla(
    for.mod3,
    family = "lognormal",
    data = dataset1,
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
  )


c(mod1$dic$dic, mod2$dic$dic, mod3$dic$dic, mod4$dic$dic)

# Models with saptial term

## Simulated coordinates
coords <- as.matrix(dataset1[, 4:5])
## Specify the boundary parameters
bound = inla.nonconvex.hull(
  as.matrix(dataset1[, 4:5]),
  convex = -0.09,
  eps = 0.05,
  resolution = 60
)
## Create the mesh
mesh <-
  inla.mesh.2d(
    loc = coords,
    boundary = bound,
    max.edge = c(0.01, 0.3),
    offset = c(0.09, 0.09),
    cutoff = 0.035,
    min.angle = 0.05
  )

## Define the weihting factors a_ik
A.random <-
  inla.spde.make.A(mesh, loc = matrix(c(dataset1$x, dataset1$y), ncol = 2))

## Function in R to define the SPDE approach
spde <- inla.spde2.matern(mesh)

## Define spatial field
w.index <- inla.spde.make.index(name = "w", n.spde = spde$n.spde)

## Data structure
Xm <- model.matrix(~ -1 + eff + bath, data = dataset1)
Xm1 <- model.matrix( ~ -1 + bath, data = dataset1)
## Covariates
X <- data.frame(eff = Xm[, 1],
                bath = Xm[, 2])
X1 <- data.frame(bath = Xm1)
N <- nrow(dataset1)

## Combine the data
StackFit <- inla.stack(
  tag = "Fit",
  # stack name
  data = list(y = dataset1$abut),
  # response variable
  A = list(1, 1, A.random),
  # the 1 correspond to an argument in the effects
  effects = list(
    Intercept = rep(1, N),
    X = X,
    w = w.index
  )
)

StackFit1 <- inla.stack(
  tag = "Fit",
  # stack name
  data = list(y1 = dataset1$CPUEt),
  # response variable
  A = list(1, 1, A.random),
  # the 1 correspond to an argument in the effects
  effects = list(
    Intercept = rep(1, N),
    X = X1,
    w = w.index
  )
)


## Model 5: Biomass ~ eff + bath + u (Gamma)

for.mod5 <- y ~ -1 + Intercept + eff + bath + f(w, model = spde)

set.seed(12345)
mod5 <-
  inla(
    for.mod5,
    family = "gamma",
    data = inla.stack.data(StackFit),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
    control.predictor = list(A = inla.stack.A(StackFit))
  )

## Model 6: Biomass ~ eff + bath + u (ln Normal)

for.mod6 <- y ~ -1 + Intercept + eff + bath + f(w, model = spde)

set.seed(12345)
mod6 <-
  inla(
    for.mod6,
    family = "lognormal",
    data = inla.stack.data(StackFit),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
    control.predictor = list(A = inla.stack.A(StackFit))
  )

## Model 7: CPUE ~ bath + u (Gamma)

for.mod7 <- y1 ~ -1 + Intercept + bath + f(w, model = spde)

set.seed(12345)
mod7 <-
  inla(
    for.mod7,
    family = "gamma",
    data = inla.stack.data(StackFit1),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
    control.predictor = list(A = inla.stack.A(StackFit1))
  )

## Model 8: CPUE ~ bath + u (Gamma)

for.mod8 <- y1 ~ -1 + Intercept + bath + f(w, model = spde)

set.seed(12345)
mod8 <-
  inla(
    for.mod8,
    family = "lognormal",
    data = inla.stack.data(StackFit1),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
    control.predictor = list(A = inla.stack.A(StackFit1))
  )

c(mod5$dic$dic, mod6$dic$dic, mod7$dic$dic, mod8$dic$dic)

## Model 9: Biomass ~ eff + f(bath, model = "rw1") + u (Gamma)

for.mod9 <-
  y ~ -1 + Intercept + eff + f(inla.group(bath, n = 30), model = "rw1") + f(w, model = spde)

set.seed(12345)
mod9 <-
  inla(
    for.mod9,
    family = "gamma",
    data = inla.stack.data(StackFit),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
    control.predictor = list(A = inla.stack.A(StackFit))
  )

## Model 10: CPUE ~ f(bath, model = "rw1")  + u (Gamma)

for.mod10 <-
  y1 ~ -1 + Intercept + f(inla.group(bath, n = 30), model = "rw1") + f(w, model = spde)

set.seed(12345)
mod10 <-
  inla(
    for.mod10,
    family = "gamma",
    data = inla.stack.data(StackFit1),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
    control.predictor = list(A = inla.stack.A(StackFit1))
  )

## Model 11: Biomass ~ eff + f(bath, model = "rw2") + u (Gamma)

for.mod11 <-
  y ~ -1 + Intercept + eff + f(inla.group(bath, n = 30), model = "rw2") + f(w, model = spde)

set.seed(12345)
mod11 <-
  inla(
    for.mod11,
    family = "gamma",
    data = inla.stack.data(StackFit),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
    control.predictor = list(A = inla.stack.A(StackFit))
  )

## Model 12: CPUE ~ f(bath, model = "rw2") + u (Gamma)

for.mod12 <-
  y1 ~ -1 + Intercept + f(inla.group(bath, n = 30), model = "rw2") + f(w, model = spde)

set.seed(12345)
mod12 <-
  inla(
    for.mod12,
    family = "gamma",
    data = inla.stack.data(StackFit1),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
    control.predictor = list(A = inla.stack.A(StackFit1))
  )

c(mod9$dic$dic, mod10$dic$dic, mod11$dic$dic, mod12$dic$dic)

## Model 13: Biomass + offset(eff) + f(bath, model = "rw1") + u (Gamma)

for.mod13 <-
  y ~ -1 + Intercept + offset(log(eff)) + f(inla.group(bath, n = 30), model = "rw1") + f(w, model = spde)

set.seed(12345)
mod13 <-
  inla(
    for.mod13,
    family = "gamma",
    data = inla.stack.data(StackFit),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
    control.predictor = list(A = inla.stack.A(StackFit))
  )

## Model 14: Biomass ~ eff + f(bath, model = "rw1") + u (ln Normal)

for.mod14 <-
  y ~ -1 + Intercept + eff + f(inla.group(bath, n = 30), model = "rw1") + f(w, model = spde)

set.seed(12345)
mod14 <-
  inla(
    for.mod14,
    family = "lognormal",
    data = inla.stack.data(StackFit),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
    control.predictor = list(A = inla.stack.A(StackFit))
  )

## Model 15: CPUE ~ f(bath, model = "rw1")  + u (ln Normal)

for.mod15 <-
  y1 ~ -1 + Intercept + f(inla.group(bath, n = 30), model = "rw1") + f(w, model = spde)

set.seed(12345)
mod15 <-
  inla(
    for.mod15,
    family = "lognormal",
    data = inla.stack.data(StackFit1),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
    control.predictor = list(A = inla.stack.A(StackFit1))
  )

## Model 16: Biomass ~ eff + f(bath, model = "rw2") + u (ln Normal)

for.mod16 <-
  y ~ -1 + Intercept + eff + f(inla.group(bath, n = 30), model = "rw2") + f(w, model = spde)

set.seed(12345)
mod16 <-
  inla(
    for.mod16,
    family = "lognormal",
    data = inla.stack.data(StackFit),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
    control.predictor = list(A = inla.stack.A(StackFit))
  )

## Model 17: CPUE ~ f(bath, model = "rw2") + u (ln Normal)

for.mod17 <-
  y1 ~ -1 + Intercept + f(inla.group(bath, n = 30), model = "rw2") + f(w, model = spde)

set.seed(12345)
mod17 <-
  inla(
    for.mod17,
    family = "lognormal",
    data = inla.stack.data(StackFit1),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
    control.predictor = list(A = inla.stack.A(StackFit1))
  )

## DIC and LCPO

DIC <-
  c(
    mod1$dic$dic,
    mod2$dic$dic,
    mod3$dic$dic,
    mod4$dic$dic,
    mod5$dic$dic,
    mod6$dic$dic,
    mod7$dic$dic,
    mod8$dic$dic,
    mod9$dic$dic,
    mod10$dic$dic,
    mod11$dic$dic,
    mod12$dic$dic,
    mod13$dic$dic,
    mod14$dic$dic,
    mod15$dic$dic,
    mod16$dic$dic,
    mod17$dic$dic
  )

LCPO <-
  c(
    -mean(log(mod1$cpo$cpo)),
    -mean(log(mod2$cpo$cpo)),
    -mean(log(mod3$cpo$cpo)),
    -mean(log(mod4$cpo$cpo)),
    -mean(log(mod5$cpo$cpo)),
    -mean(log(mod6$cpo$cpo)),
    -mean(log(mod7$cpo$cpo)),
    -mean(log(mod8$cpo$cpo)),
    -mean(log(mod9$cpo$cpo)),
    -mean(log(mod10$cpo$cpo)),
    -mean(log(mod11$cpo$cpo)),
    -mean(log(mod12$cpo$cpo)),
    -mean(log(mod13$cpo$cpo)),
    -mean(log(mod14$cpo$cpo)),
    -mean(log(mod15$cpo$cpo)),
    -mean(log(mod16$cpo$cpo)),
    -mean(log(mod17$cpo$cpo))
  )

CPO <-
  c(
    sum(mod1$cpo$cpo),
    sum(mod2$cpo$cpo),
    sum(mod3$cpo$cpo),
    sum(mod4$cpo$cpo),
    sum(mod5$cpo$cpo),
    sum(mod6$cpo$cpo),
    sum(mod7$cpo$cpo),
    sum(mod8$cpo$cpo),
    sum(mod9$cpo$cpo),
    sum(mod10$cpo$cpo),
    sum(mod11$cpo$cpo),
    sum(mod12$cpo$cpo),
    sum(mod13$cpo$cpo),
    sum(mod14$cpo$cpo),
    sum(mod15$cpo$cpo),
    sum(mod16$cpo$cpo),
    sum(mod17$cpo$cpo)
  )


mods <- cbind(DIC, LCPO)

rownames(mods) <- c(
  "mod1",
  "mod2",
  "mod3",
  "mod4",
  "mod5",
  "mod6",
  "mod7",
  "mod8",
  "mod9",
  "mod10",
  "mod11",
  "mod12",
  "mod13",
  "mod14",
  "mod15",
  "mod16",
  "mod17"
)



# MODELS WITH ZEROS: HURDLE MODELS
## Hurdle models fit two process: (1) binary process, which means that we need to add a variable with 0 and 1 (presence and absence)
## and (2) continuous process. 

dataset1$CPUE_01 <- ifelse(dataset1$CPUE > 0, 1, 0)
dataset1$abu_01 <- ifelse(dataset1$abu > 0, 1, 0)

## MESH
### coordinates
coords <- as.matrix(dataset1[, 4:5])
### Boundary
bound = inla.nonconvex.hull(
  as.matrix(dataset1[, 4:5]),
  convex = -0.09,
  eps = 0.05,
  resolution = 60
)

mesh <-
  inla.mesh.2d(
    loc = coords,
    boundary = bound,
    max.edge = c(0.01, 0.3),
    offset = c(0.09, 0.09),
    cutoff = 0.035,
    min.angle = 0.05
  )

## SPDE
spde <- inla.spde2.matern(mesh)

## make.A.matrix
est.temp <-
  inla.spde.make.A(mesh, loc = matrix(c(dataset1$x, dataset1$y), ncol = 2))


## INDEX for both process bernouilli/continuous (w)
mesh.index.bin <- inla.spde.make.index("i.bin", n.spde = spde$n.spde)
mesh.index.con <- inla.spde.make.index("i.con", n.spde = spde$n.spde)

## Stack 1

X <- data.frame(bath.bin = dataset1$bath) # in this case we only hace one covariate
N <- nrow(dataset1)

est.bin <- inla.stack(
  tag = "est.bin", # name of the stack
  # stack name
  data = list(y = cbind(dataset1$CPUE_01, NA)), # response variable (We need to add NA because after we will combine both process (Stack 1 (est.bin) and 2 (est.con)))
  A = list(1, 1, est.temp), # 1 for the intercept and the covariates and the est.temp is the a (u = a x w)
  effects = list(
    Intercept = rep(1, N),
    X = X,
    w = mesh.index.bin # specifies the w related with the process
  )
)


X1 <- data.frame(bath.con = dataset1$bath)

est.con <- inla.stack(
  tag = "est.con",
  # stack name
  data = list(y = cbind(
    NA, ifelse(dataset1$CPUE > 0, dataset1$CPUE, NA)
  )),
  # response variable: in this case we only add the positve part
  A = list(1, 1, est.temp),
  # the 1 correspond to an argument in the effects
  effects = list(
    Intercept = rep(1, N),
    X = X1,
    w = mesh.index.con
  )
)


est.stack <- inla.stack(est.bin, est.con)

## MODH 1: CPUE ~ f(bath) (the effect is share for both process "copy") + u_1Bernouilli + u_2Gamma
modh1 <-
  y ~ -1 + Intercept + f(inla.group(bath.bin), model = "rw1") + f(inla.group(bath.con), copy =
                                                                    "inla.group(bath.bin)", fixed = F) +
  f(i.bin, model = spde) +
  f(i.con, model = spde) 

set.seed(12345)
modh1 <- inla(
  modh1,
  family = c('binomial', "gamma"),
  data = inla.stack.data(est.stack),
  control.compute = list(
    dic = TRUE,
    cpo = TRUE,
    waic = T,
    config = FALSE
  ),
  control.predictor = list(A = inla.stack.A(est.stack), compute =
                             TRUE)
) 

modh1dic <- modh1$dic$dic
modh1lcpo <- -mean(log(modh1$cpo$cpo), na.rm = T)

## MODH 2: f(bath) (the effect is share for both process "copy") + u (the spatial effect is share)
modh2 <-
  y ~ -1 + Intercept + f(inla.group(bath.bin), model = "rw1") + f(inla.group(bath.con), copy =
                                                                    "inla.group(bath.bin)", fixed = F) +
  f(i.bin, model = spde) +
  f(i.con, copy = "i.bin", fixed = F)

set.seed(12345)
modh2 <- inla(
  modh2,
  family = c('binomial', "gamma"),
  data = inla.stack.data(est.stack),
  control.compute = list(
    dic = TRUE,
    cpo = TRUE,
    waic = T,
    config = FALSE
  ),
  control.predictor = list(A = inla.stack.A(est.stack), compute =
                             TRUE)
) 


modh2dic <- modh2$dic$dic
modh2lcpo <- -mean(log(modh2$cpo$cpo), na.rm = T)

## MODH 3: f(bath) (the effect is share for both process "copy") + u (the spatial effect is share)
modh3 <-
  y ~ -1 + Intercept + f(inla.group(bath.bin), model = "rw1") + f(inla.group(bath.con), copy =
                                                                    "inla.group(bath.bin)", fixed = F) +
  f(i.bin, model = spde) +
  f(i.con, copy = "i.bin", fixed = F)

set.seed(12345)
modh3 <- inla(
  modh3,
  family = c('binomial', "lognormal"),
  data = inla.stack.data(est.stack),
  control.compute = list(
    dic = TRUE,
    cpo = TRUE,
    waic = T,
    config = FALSE
  ),
  control.predictor = list(A = inla.stack.A(est.stack), compute =
                             TRUE)
) 


modh3dic <- modh3$dic$dic
modh3lcpo <- -mean(log(modh3$cpo$cpo), na.rm = T)

## MODH 4: f(bath) (the effect is share for both process "copy") + u_1Bernouilli + u_2lognormal
modh4 <-
  y ~ -1 + Intercept + f(inla.group(bath.bin), model = "rw1") + f(inla.group(bath.con), copy =
                                                                    "inla.group(bath.bin)", fixed = F) +
  f(i.bin, model = spde) +
  f(i.con, model = spde)

set.seed(12345)
modh4 <- inla(
  modh4,
  family = c('binomial', "lognormal"),
  data = inla.stack.data(est.stack),
  control.compute = list(
    dic = TRUE,
    cpo = TRUE,
    waic = T,
    config = FALSE
  ),
  control.predictor = list(A = inla.stack.A(est.stack), compute =
                             TRUE)
) # the fit is better with lognormal family in terms of DIC and LCPO but the validation is better with gamma distribution


modh4dic <- modh4$dic$dic
modh4lcpo <- -mean(log(modh4$cpo$cpo), na.rm = T)

## MODH5: response variable (Catches/Biomass)

## Stack2

X <- data.frame(bath.bin = dataset1$bath)
N <- nrow(dataset1)

est.bin1 <- inla.stack(
  tag = "est.bin1",
  # stack name
  data = list(y1 = cbind(dataset1$abu_01, NA)),
  # response variable
  A = list(1, 1, est.temp),
  # the 1 correspond to an argument in the effects
  effects = list(
    Intercept1 = rep(1, N),
    X = X,
    w = mesh.index.bin
  )
)


X1 <- data.frame(bath.con = dataset1$bath)

est.con1 <- inla.stack(
  tag = "est.con1",
  # stack name
  data = list(y1 = cbind(
    NA, ifelse(dataset1$abu > 0, dataset1$abu, NA)
  )),
  # response variable
  A = list(1, 1, est.temp),
  # the 1 correspond to an argument in the effects
  effects = list(
    Intercept = rep(1, N),
    X = X1,
    w = mesh.index.con
  )
)

est.stack1 <- inla.stack(est.bin1, est.con1)

modh5 <-
  y1 ~ -1 + Intercept + f(inla.group(bath.bin), model = "rw1") + f(inla.group(bath.con), copy =
                                                                     "inla.group(bath.bin)", fixed = F) +
  f(i.bin, model = spde) +
  f(i.con, model = spde)

set.seed(12345)
modh5 <- inla(
  modh5,
  family = c('binomial', "gamma"),
  data = inla.stack.data(est.stack1),
  control.compute = list(
    dic = TRUE,
    cpo = TRUE,
    waic = T,
    config = FALSE
  ),
  control.predictor = list(A = inla.stack.A(est.stack1), compute =
                             TRUE)
) # the fit is better with lognormal family in terms of DIC and LCPO but the validation is better with gamma distribution



modh5dic <- modh5$dic$dic
modh5lcpo <- -mean(log(modh5$cpo$cpo), na.rm = T)

## MODH6

## Stack 3

X2 <- data.frame(bath = dataset1$bath)
N <- nrow(dataset1)

est.bin2 <- inla.stack(
  tag = "est.bin2",
  # stack name
  data = list(y2 = cbind(dataset1$CPUE_01, NA)),
  # response variable
  A = list(1, 1, est.temp),
  # the 1 correspond to an argument in the effects
  effects = list(
    Intercept1 = rep(1, N),
    X = X2,
    w = mesh.index.bin
  )
)

est.con2 <- inla.stack(
  tag = "est.con2",
  # stack name
  data = list(y2 = cbind(
    NA, ifelse(dataset1$CPUE > 0, dataset1$CPUE, NA)
  )),
  # response variable
  A = list(1, 1, est.temp),
  # the 1 correspond to an argument in the effects
  effects = list(
    Intercept = rep(1, N),
    X = X2,
    w = mesh.index.con
  )
)

est.stack2 <- inla.stack(est.bin2, est.con2)

modh6 <- y2 ~ -1 + Intercept + I(bath ^ 2) +
  f(i.bin, model = spde) +
  f(i.con, model = spde)

set.seed(12345)
modh6 <- inla(
  modh6,
  family = c('binomial', "gamma"),
  data = inla.stack.data(est.stack2),
  control.compute = list(
    dic = TRUE,
    cpo = TRUE,
    waic = T,
    config = FALSE
  ),
  control.predictor = list(A = inla.stack.A(est.stack2), compute =
                             TRUE)
) # the fit is better with lognormal family in terms of DIC and LCPO but the validation is better with gamma distribution


modh6dic <- modh6$dic$dic
modh6lcpo <- -mean(log(modh6$cpo$cpo), na.rm = T)

## DIC and LCPO

modshdic <-
  c(modh1dic, modh2dic, modh3dic, modh4dic, modh5dic, modh6dic)
modshlcpo <-
  c(modh1lcpo,
    modh2lcpo,
    modh3lcpo,
    modh4lcpo,
    modh5lcpo,
    modh6lcpo)

modsh <- cbind(modshdic, modshlcpo)
rownames(modsh) <-
  c("modh1", "modh2", "modh3", "modh4", "modh5", "modh6")
modsh
