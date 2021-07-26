# Scenario 3: Areal_data

## This scenario depends on the others scenarios 

mat <- matrix(1:100, nrow = 10, byrow = TRUE)
raster <- raster(mat)
map <- rasterToPolygons(raster)

plot(map)
title("Map areal data")
points(xy)

areal_data <- seq(0,1,0.1) # square sides
mean_catches <- vector()
n <- 1 

for(i in (length(areal_data)):2) {
  for(j in 1:(length(areal_data)-1)) {
    cond <- which(
      data$x >= areal_data[j] & 
        data$x < areal_data[j + 1] &
        data$y < areal_data[i] &
        data$y >= areal_data[i - 1]
    )
    cap.grid <- data$catch[cond]
    mean_catches[n] <- mean(cap.grid)
    n <- n + 1
  }
}

map$catches <- mean_catches
map$catches[is.na(map$catches)] <- 0
map$layer <- NULL


## Effort

areal_data <- seq(0,1,0.1) # square sides
mean_effort <- vector()
n <- 1 
for(i in (length(areal_data)):2) {
  for(j in 1:(length(areal_data)-1)) {
    cond <- which(
      data$x >= areal_data[j] & 
        data$x < areal_data[j + 1] &
        data$y < areal_data[i] &
        data$y >= areal_data[i - 1]
    )
    cap.grid <- data$effort[cond]
    mean_effort[n] <- mean(cap.grid)
    n <- n + 1
  }
}

map$effort <- mean_effort
map$effort[is.na(map$effort)] <- 0

## Bathymetry

areal_data <- seq(0,1,0.1) # square sides
mean_bat <- vector()
n <- 1 
for(i in (length(areal_data)):2) {
  for(j in 1:(length(areal_data)-1)) {
    cond <- which(
      data$x >= areal_data[j] & 
        data$x < areal_data[j + 1] &
        data$y < areal_data[i] &
        data$y >= areal_data[i - 1]
    )
    cap.grid <- data$bat[cond]
    mean_bat[n] <- mean(cap.grid)
    n <- n + 1
  }
}

map$bat <- mean_bat
map$bat[is.na(map$bat)] <- 0

## Results

spplot(map[,1], main="Catches")
spplot(map[,2], main="Effort")
spplot(map[,3], main = "Bathymetry")
