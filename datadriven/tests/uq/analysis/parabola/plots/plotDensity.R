plot2dNormal <- function() {
  par(mar=c(4.5, 5, 3, 10))

  ## get data
  n <- 20
  x <- seq(0, 1, length.out = n)
  data <- matrix(, nrow = n * n, ncol = 3)
  k <- 1
  for (i in seq(along = x)) {
    for (j in seq(along = x)) {
      data[k, 1] <- x[i]
      data[k, 2] <- x[j]
      data[k, 3] <- dnorm(x[i], 0.2, 0.1) * dnorm(x[j], 0.2, 0.1)
      k <- k + 1
    }
  }

  ## plot function
  plotFunction(data, 
               xlab = "$\\theta_1$",
               ylab = "$\\theta_2$")
  
  image.plot(legend.only=TRUE,
             zlim=c(min(data[, 3]), max(data[, 3])),
             legend.width = 0.8,
             legend.shrink = 0.8,
             nlevel = 256,
             legend.mar = 10.5,
             axis.args = list(cex.axis = 2),
             legend.args = list(text='$p(\\theta_1, \\theta_2)$',
               side=4, font=2, line=6, cex=2))
}

plot2dNormalReduced <- function() {
  par(mar=c(0,0,0,0))

  ## get data
  n <- 20
  x <- seq(-3, 3, length.out = n)
  data <- matrix(, nrow = n * n, ncol = 3)
  k <- 1
  for (i in seq(along = x)) {
    for (j in seq(along = x)) {
      data[k, 1] <- (x[i] + 3) / 6
      data[k, 2] <- (x[j] + 3) / 6
      data[k, 3] <- dnorm(x[i], sd = 1.5) * dnorm(x[j], sd = 1.5)
      k <- k + 1
    }
  }
}
