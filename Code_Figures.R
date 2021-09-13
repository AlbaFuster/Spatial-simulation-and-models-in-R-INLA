# PLOTS

## Figure 1

image.plot(rf.s.p, col=mako(20), main = "Figure 1. LGCP of random and preferential sampling") 
points(xy.p, type="p", pch=19, col = "yellow") 
points(xy.r, type="p", pch=19, col = 2)

## Figure 2

par(mfrow=c(2,2))
image.plot(abu.image, col=mako(20))
title(main = "Figure 2.1. Abundance")
image.plot(bath.image, col=mako(20))
title(main = "Figure 2.2. Bathymetry")
image.plot(abu.image, col=mako(20))
title(main = "Figure 2.3. Abundance and point patterns")
points(xy.p, type="p", pch=19, col = "yellow") 
points(xy.r, type="p", pch=19, col = "red")
image.plot(bath.image, col=mako(20))
title(main = "Figure 2.4. Bathymetry and point patterns")
points(xy.p, type="p", pch=19, col = "yellow") 
points(xy.r, type="p", pch=19, col = "red")

## Figure 3

p1.1 <- ggplot(data.variables.r, aes(x = bath.r)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "red") +
  geom_density() + labs(x = "Bathymetry", y = "Frequency") + ggtitle("Figure 3.1. Bathymetry (random sampling)")

p1.2 <- ggplot(data.variables.r, aes(x = abu.r)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "red") +
  geom_density() + labs(x = "Abundance", y = "Frequency")  + ggtitle("Figure 3.2. Abundance (random sampling)")

p1.3 <- ggplot(data.variables.p, aes(x = bath.p)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "yellow") +
  geom_density() +  labs(x = "Bathymetry", y = "Frequency") + ggtitle("Figure 3.3. Bathymetry (preferential sampling)")

p1.4 <- ggplot(data.variables.p, aes(x = abu.p)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "yellow") +
  geom_density() +  labs(x = "Abundance", y = "Frequency")  + ggtitle("Figure 3.4. Abundance (preferential sampling)")


grid.arrange(p1.2,p1.1,p1.4,p1.3, ncol=2, nrow=2)

## Figure 4

p1 <-
  ggplot(data.variables.p, aes(bath.p, abu.p)) + geom_point() + geom_smooth(method = "loess") + labs(x = "bathymetry", y = "abundance") + ggtitle("Figure 4.1. Bath vs abu (preferential sampling)")
p2 <-
  ggplot(data.variables.r, aes(bath.r, abu.r)) + geom_point() + geom_smooth(method = "loess") + labs(x = "bathymetry", y = "abundance") + ggtitle("Figure 4.2. Bath vs abu (random sampling)")
p3 <-
  ggplot(data.variables.p, aes(effort.p, abu.p)) + geom_point() + geom_smooth(method = "loess") + labs(x = "effort", y = "abundance") + ggtitle("Figure 4.3. Effort vs abu (preferential sampling)")
p4 <-
  ggplot(data.variables.r, aes(effort.r, abu.r)) + geom_point() + geom_smooth(method = "loess") + labs(x = "effort", y = "abundance") + ggtitle("Figure 4.4. Effort vs abu (random sampling)")

grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

## Figure 5

data.variables.p$coords.x.p <- data.p$x
data.variables.p$coords.y.p <- data.p$y
data.variables.r$coords.x.r <- data.r$x
data.variables.r$coords.y.r <- data.r$y

p1 <- ggplot(data.variables.p, aes(x=coords.x.p, y=coords.y.p, size = abu.p)) +
  geom_point(alpha=0.1) + theme_classic() + labs(x = "coords x", y = "coords y") + ggtitle("Figure 5.1. Distribution of abundance (preferential sampling)")

p2 <- ggplot(data.variables.r, aes(x=coords.x.r, y=coords.y.r, size = abu.r)) +
  geom_point(alpha=0.1) + theme_classic() + labs(x = "coords x", y = "coords y") + ggtitle("Figure 5.2. Distribution of abundance (random sampling)")

grid.arrange(p1, p2, ncol=2)

## Figure 6

plot(map)
title("Figure 6. Map for areal data")
points(xy.p, pch=19)

## Figure 7

p1.1 <- spplot(map[,1], main="Abundance")
p2.2 <- spplot(map[,2], main="Effort")
p3.3 <- spplot(map[,3], main = "Bathymetry")

grid.arrange(arrangeGrob(p2.2,p3.3, ncol=1, nrow=2),
             arrangeGrob(p1.1, ncol=1, nrow=1), heights=c(20,1), widths=c(1,1), top = textGrob("Figure 7. Response and explanatory variables (Areal data)",gp=gpar(fontsize=15)))

## Figure 8

p1.1 <- spplot(map[,1], main="Abundance")
p2.2 <- spplot(map[,2], main="Effort")
p3.3 <- spplot(map[,3], main = "Bathymetry")

grid.arrange(arrangeGrob(p2.2,p3.3, ncol=1, nrow=2),
             arrangeGrob(p1.1, ncol=1, nrow=1), heights=c(20,1), widths=c(1,1), top = textGrob("Figure 8. Resize the square",gp=gpar(fontsize=15)))