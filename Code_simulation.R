# GEOREFERENCED DATA: RANDOM AND PREFERENTIAL SAMPLING

## Simulate Log-Gaussian Cox Process

### Code of random sampling

#### Function parameters for data simulation
beta0.r <- 6.3
win.r <- owin(c(0, 1), c(0, 1))
nu.r <- 1
variance.r <- 0.0005
scale.r <- 0.1

#### Apply the function
set.seed(12345)
lg.s.r <-
  rLGCP("exp", #### homogeneous LGCP with exponential covariance function
        beta0.r,
        var = variance.r,
        scale = scale.r,
        win = win.r) #### returns a list with the coordinates of the points
#### and the value of lambda

Lam.r <- attr(lg.s.r, 'Lambda')
rf.s.r <- log(Lam.r$v)

xy.r <- cbind(lg.s.r$x, lg.s.r$y)[, 2:1]
n.r <- length(xy.r[, 1])
x.r <- xy.r[, 1]
y.r <- xy.r[, 2]

data.r <- data.frame(x = x.r, y = y.r)
coords.r <- rbind(data.r$x, data.r$y)

### Code of preferential sampling

#### Function parameters for data simulation
beta0.p <- 4  
win.p <- owin(c(0, 1), c(0, 1))
nu.p <- 1
var.p <- 7.5
scale.p <- 0.15

#### Apply the function
set.seed(12345)
lg.s.p <- rLGCP(
  'matern', #### inhomogeneous LGCP with Matern covariance function
  beta0.p,
  var = var.p,
  scale = scale.p,
  nu = nu.p,
  win = win.p
) #### returns a list with the coordinates of the points
#### and the value of lambda

Lam.p <- attr(lg.s.p, 'Lambda')
rf.s.p <- log(Lam.p$v)

xy.p <- cbind(lg.s.p$x, lg.s.p$y)[, 2:1]
x.p <- xy.p[, 1]
y.p <- xy.p[, 2]

data.p <- data.frame(x = x.p, y = y.p)
coords.p <- rbind(data.p$x, data.p$y)

## Simulation of the study variables 

### Response variable

#### Abundance
Lam.bat <- as.vector(rf.s.p) ### transform log lambda in vector
set.seed(12345)
abu <- vector()
for (i in 1:length(Lam.bat)) {
  if (Lam.bat[i] > 0 & Lam.bat[i] <= 4.469) {
    abu[i] <- rpois(1, 0.2) ### A poisson is used to simulate zeros.
  } else if (Lam.bat[i] > 4.469 & Lam.bat[i] <= 6.423) {
    abu[i] <- rgamma(1, 8, 0.05)
  } else {
    abu[i] <- rnorm(1, 350, 40)
  }
}
abu.image <- matrix(abu, ncol = 128, nrow = 128)

### Explanatory variables

#### Bathymetry
set.seed(12345)
bath <- matrix(ncol = 128, nrow = 128)
for (i in 1:nrow(Lam.p$v)) {
  for (j in 1:ncol(Lam.p$v)) {
    a <- abu.image[i, j]
    if (a > 100 &a <=  250) {
      bath[i, j] <- rnorm(1, 250, 50)
      ### Assign values according to coordinates
    } else if (a > 250 &a <= 350) {
      bath[i,j] <- rnorm(1,400, 50)
    } else if (a >= 350) {
      bath[i,j] <- rnorm(1,500, 50)
    } else if (Lam.p$yrow[i] < 0.3) {
      bath[i, j] <- rgamma(1, 1, 0.04)
    } else {
      bath[i, j] <- rnorm(1, 850.05, 100)
    }
  }
}
bath.image <- matrix(bath, ncol = 128, nrow = 128)

### Extraction of variable values related to point patterns


#### Extract for preferential sampling
##### Variable response: abundance
r.abu.p <- raster(abu.image) ##### transform values of abudance in raster format
r.abu.p <- flip(t(r.abu.p), direction = 'y') 
sp.p <- SpatialPoints(data.p) ##### transform coords in a SpatialPoints
val1 <- extract(r.abu.p,sp.p) ##### Extract values from the raster
##### that coincide with the SpatialPoints.

##### Variable explanatory: bathymetry
r.bath.p <- raster(bath.image)
r.bath.p <- flip(t(r.bath.p), direction='y')
val2 <- extract(r.bath.p,sp.p)

#### Extract for random sampling
r.abu.r <- raster(abu.image)
r.abu.r <- flip(t(r.abu.r), direction = 'y')
sp.r <- SpatialPoints(data.r)
val3 <- extract(r.abu.r,sp.r)

r.bath.r <- raster(bath.image)
r.bath.r <- flip(t(r.bath.r), direction='y')
val4 <- extract(r.bath.r,sp.r)

##### Include the variables in a dataframe 
##### Preferential sampling
data.variables.p <- data.frame(bath.p = val2, abu.p = val1)
##### Random sampling
data.variables.r <- data.frame(bath.r = val4, abu.r = val3)

#### Effort

##### Effort for random sampling
set.seed(12345)
effort.r <- rnorm(535, 30, 0.5)

##### Effort for preferential sampling
set.seed(12345)
effort.p <- vector()
for (i in 1:length(data.variables.p$abu.p)) {
  if (data.variables.p$abu.p[i] >= 0 &
      data.variables.p$abu.p[i] <= 100) {
    effort.p[i] <- rnorm(1, 25, 3)
  } else if (data.variables.p$abu.p[i] > 100 &
             data.variables.p$abu.p[i] <= 200.4) {
    effort.p[i] <- rnorm(1, 35, 3.5)
  } else {
    effort.p[i] <- rnorm(1, 45, 2.5)
  }
}

data.variables.r$effort.r <- effort.r
data.variables.p$effort.p <- effort.p

# AREAL DATA 

## Generate a 10 x 10 matrix 
mat <- matrix(1:100, nrow = 10, byrow = TRUE)
## Transform the matrix into a raster 
raster <- raster(mat)
map <- rasterToPolygons(raster)

## Variable response: Abundance 

### Sum of the abundance value for each of the squares
areal_data <- seq(0, 1, 0.1) ### square sides
sum_abu.p <- vector()
n <- 1

for (i in (length(areal_data)):2) {
  for (j in 1:(length(areal_data) - 1)) {
    ### Conditional to know where you are on the map 
    cond <- which(
      data.p$x >= areal_data[j] &
        data.p$x < areal_data[j + 1] &
        data.p$y < areal_data[i] &
        data.p$y >= areal_data[i - 1] 
    )
    cap.grid <- data.variables.p$abu.p[cond]
    sum_abu.p[n] <- sum(cap.grid) ### Sum of abundace 
    n <- n + 1
  }
}

map$abu.p <- sum_abu.p
map$abu.p[is.na(map$abu.p)] <- 0 ### Change NA to 0
map$layer <- NULL 

## Explanatory variables: Effort and bathymetry

### Effort
sum_effort.p <- vector()
n <- 1
for (i in (length(areal_data)):2) {
  for (j in 1:(length(areal_data) - 1)) {
    cond <- which(
      data.p$x >= areal_data[j] &
        data.p$x < areal_data[j + 1] &
        data.p$y < areal_data[i] &
        data.p$y >= areal_data[i - 1]
    )
    cap.grid <- data.variables.p$effort.p[cond]
    sum_effort.p[n] <- sum(cap.grid)
    n <- n + 1
  }
}

map$effort.p <- sum_effort.p
map$effort.p[is.na(map$effort.p)] <- 0

### Bathymetry
mean_bat.p <- vector()
n <- 1
for (i in (length(areal_data)):2) {
  for (j in 1:(length(areal_data) - 1)) {
    cond <- which(
      data.p$x >= areal_data[j] &
        data.p$x < areal_data[j + 1] &
        data.p$y < areal_data[i] &
        data.p$y >= areal_data[i - 1]
    )
    cap.grid <- data.variables.p$bath.p[cond]
    mean_bat.p[n] <- mean(cap.grid)
    n <- n + 1
  }
}

map$bat.p <- mean_bat.p
map$bat.p[is.na(map$bat.p)] <- 0

## Resize the square

mat <- matrix(1:200, nrow = 20, byrow = TRUE)
raster <- raster(mat)
map <- rasterToPolygons(raster)

### Abundance
areal_data <- seq(0,1,0.1) 
sum_abu.p <- vector()
n <- 1 

for(i in (length(areal_data)):2) {
  for(j in 1:(length(areal_data)-1)) {
    cond <- which(
      data.p$x >= areal_data[j] & 
        data.p$x < areal_data[j + 1] &
        data.p$y < areal_data[i] &
        data.p$y >= areal_data[i - 1]
    )
    cap.grid <- data.variables.p$abu.p[cond]
    sum_abu.p[n] <- sum(cap.grid)
    n <- n + 1
  }
}

map$abu.p <- sum_abu.p
map$abu.p[is.na(map$abu.p)] <- 0
map$layer <- NULL

### Effort
areal_data <- seq(0,1,0.1) # square sides
sum_effort.p <- vector()
n <- 1 
for(i in (length(areal_data)):2) {
  for(j in 1:(length(areal_data)-1)) {
    cond <- which(
      data.p$x >= areal_data[j] & 
        data.p$x < areal_data[j + 1] &
        data.p$y < areal_data[i] &
        data.p$y >= areal_data[i - 1]
    )
    cap.grid <- data.variables.p$effort.p[cond]
    sum_effort.p[n] <- sum(cap.grid)
    n <- n + 1
  }
}

map$effort.p <- sum_effort.p
map$effort.p[is.na(map$effort.p)] <- 0

### Bathymetry
areal_data <- seq(0,1,0.1) 
mean_bat.p <- vector()
n <- 1 
for(i in (length(areal_data)):2) {
  for(j in 1:(length(areal_data)-1)) {
    cond <- which(
      data.p$x >= areal_data[j] & 
        data.p$x < areal_data[j + 1] &
        data.p$y < areal_data[i] &
        data.p$y >= areal_data[i - 1]
    )
    cap.grid <- data.variables.p$bath.p[cond]
    mean_bat.p[n] <- mean(cap.grid)
    n <- n + 1
  }
}

map$bat.p <- mean_bat.p
map$bat.p[is.na(map$bat.p)] <- 0

# Prepare data 

## Random sampling 

dataset1 <- data.frame(abu = data.variables.r$abu.r, bath = data.variables.r$bath.r, eff = data.variables.r$effort.r, x = data.r$x, y = data.r$y)

## Preferential sampling

dataset2 <- data.frame(abu = data.variables.p$abu.p, bath = data.variables.p$bath.p, eff = data.variables.p$effort.p, x = data.p$x, y = data.p$y)
