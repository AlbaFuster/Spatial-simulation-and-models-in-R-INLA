# ----------------- PREFERENTIAL SAMPLING

setwd("C:/Users/albaf/OneDrive/Escritorio/Prácticas IEO/Simulación/GitHub_Simulation")
source(file = "funcion_dual_mesh.R") # function book.mesh.dual
load("C:/Users/albaf/OneDrive/Escritorio/Prácticas IEO/Simulación/GitHub_Simulation/dataframes.RData")

# ----------------- PACKAGES

library(INLA)
library(inlabru)
library(ggplot2)
library(spatstat)
library(RandomFields)
library(rgeos)
library(sp)
library(maptools)
library(rgdal)
library(rgeos)
library(raster)
library(rasterVis)

# ----------------- TRANSFORM RESPONSE VARIABLE

dataset2$abut <- dataset2$abu + 0.1
dataset2$CPUE <- dataset2$abu / dataset2$eff
dataset2$CPUEt <- dataset2$abut / dataset2$eff

# ----------------- MESH

coords <- as.matrix(dataset2[, 4:5])
bound = inla.nonconvex.hull(
  as.matrix(dataset2[, 4:5]),
  convex = -0.085,
  eps = 0.01,
  resolution = 60
)
mesh <-
  inla.mesh.2d(
    loc = coords,
    boundary = bound,
    max.edge = c(0.050, 0.25),
    offset = c(0.05, 0.3),
    cutoff = 0.01,
    min.angle = 0.05
  )
plot(mesh)
points(dataset2[, 4:5])

# ----------------- SPDE

spde <- inla.spde2.matern(mesh)

# ----------------- DUAL MESH

nv <- mesh$n
n <- length(dataset2$CPUEt)
loc.d <- cbind(c(0,1,1,0,0), c(0,0,1,1,0))

dmesh <- book.mesh.dual(mesh)
domain.polys <- Polygons(list(Polygon(loc.d)), '0')
domainSP <- SpatialPolygons(list(domain.polys))

w<- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i, ], domainSP))
    return(gArea(gIntersection(dmesh[i, ], domainSP)))
  else return(0)
})

plot(mesh$loc, asp = 1, col = (w == 0) + 1, pch = 19,main =" ") 
plot(dmesh, add=TRUE)
lines(loc.d, col = 3)
points(dataset2[,4:5],col = "yellow")

# ----------------- Vector representing the points (1) and the mesh nodes (0)

y.pp <- rep(0:1, c(nv,n))

# ----------------- Vector exposure (E)

e.pp <- c(w,rep(0,n))

# ----------------- Projector matrix
# ----------------- For the integration points this is just a diagonal matrix because these locations are just the mesh vertices

imat <- Diagonal(nv, rep(1,nv))

# ----------------- For the observed points, another projection matrix is defined

lmat <- inla.spde.make.A(mesh, loc=coords)

# ----------------- The entire projection matrix

A.pp <- rbind(imat, lmat)

# ----------------- MODELS 

# MOD1: CPUE ~ b0.p + b0.c + u1 + alphau1 + f(bath, rw1)

stk1 <- inla.stack(
  data = list(y = cbind(dataset2$CPUEt, NA), e = rep(0, n)), 
  A = list(lmat, 1,1), 
  effects = list(i = 1:nv, b0.y = rep(1, n), bath = dataset2$bath), 
  tag = 'resp2')

stk.pp1 <- inla.stack(data = list(y = cbind(NA, y.pp), e = e.pp), 
                      A = list(A.pp, 1),
                      effects = list(j = 1:nv, b0.pp = rep(1, nv + n)),
                      tag = 'pp2')

j.stk1 <- inla.stack(stk1, stk.pp1)

jform1 <- y ~ 0 + b0.pp + b0.y + f(i, model = spde) +
  f(j, copy = 'i', fixed = FALSE) + f(inla.group(bath), model = "rw1")

set.seed(12345)
j.res1 <- inla(jform1, family = c('lognormal', 'poisson'), 
               data = inla.stack.data(j.stk1),
               E = inla.stack.data(j.stk1)$e,
               control.predictor = list(A = inla.stack.A(j.stk1)), control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE), num.threads = 2, verbose = TRUE)
summary(j.res1)

# MOD2: CPUE ~ b0.p + b0.c + u1 + alphau1 + f(bath, rw2)

jform2 <- y ~ 0 + b0.pp + b0.y + f(i, model = spde) +
  f(j, copy = 'i', fixed = FALSE) + f(inla.group(bath), model = "rw2")

set.seed(12345)
j.res2 <- inla(jform5, family = c('lognormal', 'poisson'), 
               data = inla.stack.data(j.stk1),
               E = inla.stack.data(j.stk1)$e,
               control.predictor = list(A = inla.stack.A(j.stk1)), control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE), num.threads = 2, verbose = TRUE)
summary(j.res2)

# MOD3: Catches ~ b0.p + b0.c + u1 + alphau1 + f(bath, rw1) + f(eff, rw1) 

Xm  <- model.matrix(~ -1 + bath + eff, data = dataset2) 
X <- data.frame(bath = Xm[,1],
                eff = Xm[,2]) 

stk2 <- inla.stack(
  data = list(y = cbind(dataset2$abut, NA), e = rep(0, n)), 
  A = list(lmat, 1,1),
  effects = list(i = 1:nv, b0.y = rep(1, n), X = X),
  tag = 'resp2')

stk.pp2 <- inla.stack(data = list(y = cbind(NA, y.pp), e = e.pp), 
                      A = list(A.pp, 1),
                      effects = list(j = 1:nv, b0.pp = rep(1, nv + n)),
                      tag = 'pp2')

j.stk2 <- inla.stack(stk2, stk.pp2)

jform2 <- y ~ 0 + b0.pp + b0.y + f(i, model = spde) + f(j, copy = "i", fixed = FALSE) + f(inla.group(eff), model = "rw1") + f(inla.group(bath), model = "rw1") 

set.seed(12345)
j.res2 <- inla(jform2, family = c('lognormal', 'poisson'), 
               data = inla.stack.data(j.stk2),
               E = inla.stack.data(j.stk2)$e,
               control.predictor = list(A = inla.stack.A(j.stk2)), control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE), num.threads = 2, verbose = TRUE)


summary(j.res2)

# MOD3: Catches ~ b0.p + b0.c + u1 + alphau1 + f(bath, rw2) + f(eff, rw2)

jform3 <- y ~ 0 + b0.pp + b0.y + f(i, model = spde) + f(j, copy = "i", fixed = FALSE) + f(inla.group(eff, n = 35), model = "rw2") + f(inla.group(bath), model = "rw2") 

set.seed(12345)
j.res3 <- inla(jform3, family = c('lognormal', 'poisson'), 
               data = inla.stack.data(j.stk2),
               E = inla.stack.data(j.stk2)$e,
               control.predictor = list(A = inla.stack.A(j.stk2)), control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE), num.threads = 2, verbose = TRUE)


summary(j.res3)

