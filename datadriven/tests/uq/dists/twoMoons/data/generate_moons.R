library("clusterSim")
## source("~/Rlibs/devplot.R")

## generate two moons data set
show("start generating moons...")
moons <- shapes.two.moon(numObjects=20000,
                         shape1a=-1.0,
                         shape2b=1.0,
                         shape1rFrom=1.4,
                         shape1rTo=2.0,
                         shape2rFrom=1.4,
                         shape2rTo=2.0,
                         outputCsv="",
                         outputCsv2="",
                         outputColNames=TRUE,
                         outputRowNames=TRUE)$data

show(ncol(moons))
## normalize them to [-1, 1]
moons[, 1] <- (moons[, 1] + 2) / 3
moons[, 2] <- (moons[, 2] + 3) / 5

## write result to file
write.table(moons, paste("twomoons.csv", sep = ""),
            append=FALSE,
            row.names = FALSE,
            col.names = FALSE,
            fileEncoding="UTF-8")

## ## plot results
## moons <- read.table("moons.csv")
## devplot(paste("moons.jpg", sep = ""), function() {
##     ## ------------------------------------------------------------------------
##     main.ranges <- paste("[[", round(min(moon[, 1]), digits=2), ", ",
##                          round(max(moon[, 1]), digits=2), "]",
##                          " [", round(min(moon[, 2]), digits=2), ", ",
##                          round(max(moon[, 2]), digits=2), "]]", sep = "")
##     main.samples <- paste(", n = ", nrow(moon), sep = "")
##     main <- paste("two moons, ", main.ranges, ", ", main.samples, sep = "")
##     ## ------------------------------------------------------------------------
##     plot(moons[, 1], moons[, 2], main=main)
## }, dev = "jpg")
