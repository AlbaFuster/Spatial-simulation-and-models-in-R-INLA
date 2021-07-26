# Scenario 2: Preferential_Sampling

## Parameters

beta0 <- 5.5
win <- owin(c(0,1),c(0,1))
nu <- 1

## Simulate the Random Point Pattern

set.seed(12345)
lg.s <- rLGCP('matern', beta0, var = 2,
              scale = 0.15, nu = nu, win = win)

Lam <- attr(lg.s, 'Lambda')
rf.s <- log(Lam$v)
summary(as.vector(rf.s))
mean(rf.s)

xy <- cbind(lg.s$x, lg.s$y)[, 2:1]
x <- xy[,1] 
y <- xy[,2]

data <- data.frame(x = x, y =y)
coords <- rbind(data$x, data$y)

image.plot(rf.s, col=mako(20)) # GRF
points(xy, type="p", pch=19) # Simulated marks

## Simulate response variable catches

set.seed(12345)
a <- rpois(n = 1, lambda = 3.5)
b <- rgamma(n = a, shape = 4, rate = 0.01) 
c <- sum(b)
dat1 <- data.frame(a,c)


for (count in 1:487){
  a <- rpois(n = 1, lambda = 3.5)
  b <- rgamma(n = a, shape = 4, rate = 0.01) 
  c <- sum(b)
  dat1[nrow(dat1) + 1,] <- list(a, c)
}

## Function to found nearest coordinate in the grid to my point

closest <- function(coord, x) {
  dt <- data.table(x, val = x)
  setattr(dt, "sorted", "x")
  setkey(dt, "x")
  point <- dt[J(coord), roll="nearest"][,2]
  as.numeric(point)
} 

## Function to found the closest Lam$v value to the point

find.grid <- function(coords, grid) {
  coord.x <- closest(coords[1], grid$xcol)
  coord.y <- closest(coords[2], grid$yrow)
  
  coord.x <- match(coord.x, grid$xcol)  
  coord.y <- match(coord.y, grid$yrow) 
  
  grid$v[coord.x,coord.y]
} 

val <- vector()
for(i in 1:488) {
  val[i] <- find.grid(coords[,i], Lam)
}

data <- data.frame(x = data$x, y = data$y, val = val)
data <- data[order(val),] 

catch.sort <- sort(dat1$c)

data$catch <- catch.sort


ggplot(data, aes(x=x, y=y, size = catch)) +
  geom_point(alpha=0.5) + theme_classic() # result 

## Simulate effort 

set.seed(12345)
effort <- rnorm(488, 30, 5)
effort <- sort(effort)

data$effort <- effort 

ggplot(data, aes(catch, effort)) + geom_point() + geom_smooth(method = "loess")

## Simulate bathymethry

set.seed(12345)
bat1 <- rnorm(163, 50, 17)
bat2 <- rnorm(162, 150, 20)
bat3 <- rnorm(163,350, 60)

data.bat1 <- data.frame(bat = bat3)
new.data.bat <- data.frame(bat = bat2)
new.data.bat1 <- data.frame(bat=bat1)

data.bat <- rbind(data.bat1, new.data.bat, new.data.bat1)

data$bat <- data.bat$bat

ggplot(data, aes(catch, bat)) + geom_point() + geom_smooth(method = "loess")

data <- data[,-3]