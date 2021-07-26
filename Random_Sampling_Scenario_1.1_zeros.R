# Scenario 1.1: RANDOM SAMPLING (georeferencing data)

## parameters for the function rLGCP

beta0 <- 5.5 # determinate the number of points 
win <- owin(c(0,1),c(0,1)) # window 0,1 and 0,1
nu <- 1 # nugget
variance <- 0.1 # the variance
scale <- 0.1 # the scale

## Simulate the Random Point Pattern

set.seed(12345)

pattern <- rLGCP("exp", beta0, var = variance, scale = scale, win=win) # Use the function rLGCP to generate a random point pattern. Homogenuous LGCP with exponential covariance function. 

Lam <- attr(pattern, 'Lambda') # log(lambda(s)) =beta_0 + S(s)
rf.s <- log(Lam$v)

xy <- cbind(pattern$x, pattern$y)[,2:1] # studied coords
image.plot(rf.s, col=mako(20)) # Random Gaussian Field
points(xy, type="p", pch=19) # Marks
n <- length(xy[,1]) # number of marks
x <- xy[,1] 
y <- xy[,2]

data <- data.frame(x = x, y =y)
coords <- rbind(data$x, data$y)

## Simulate the response variable (catches)

set.seed(12345)
a <- rpois(n = 1, lambda = 2.1) # simulate 0
b <- rgamma(n = a, shape = 5, rate = 0.01) # simulate positive values
c <- sum(b)
dat1 <- data.frame(a,c)


for (count in 1:251){
  a <- rpois(n = 1, lambda = 2.1)
  b <- rgamma(n = a, shape = 5, rate = 0.01) 
  c <- sum(b)
  dat1[nrow(dat1) + 1,] <- list(a, c)
}

## Calculate distance between two coordinates

library(geosphere)
d <- vector()
for(i in 1:252) { 
  d[i] <- distm(c(0,0), coords[,i], fun = distHaversine)
}

data$d <- d # distance between coordinates 
data <- data[order(d),]  
catches.sort <- sort(dat1$c) # order max to min the response variable
data$catch <- catches.sort # 

## Plotting of the marks with the value of the associated catches

ggplot(data, aes(x=x, y=y, size = catch)) +
  geom_point(alpha=0.1) + theme_classic()

## Simulate the effort and bathymetry

set.seed(12345)
effort <- rnorm(252, 30, 0.5)
data$effort <- effort

set.seed(12345)
bat1 <- rnorm(84, 50, 17)
bat2 <- rnorm(84, 150, 20)
bat3 <- rnorm(84,350, 60)

data.bat1 <- data.frame(bat = bat3)
new.data.bat <- data.frame(bat = bat2)
new.data.bat1 <- data.frame(bat=bat1)

data.bat <- rbind(data.bat1, new.data.bat, new.data.bat1)

data$bat <- data.bat$bat

## Relation catches and bathymetry
ggplot(data, aes(catch, bat)) + geom_point() + geom_smooth(method = "loess")

ggplot(data, aes(catch, effort)) + geom_point() + geom_smooth(method = "loess")




