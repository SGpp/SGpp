library("clusterSim")
source("~/Rlibs/devplot.R")

## generate two moons data set
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
                         outputRowNames=TRUE)

## rotate the result by 90 degrees
clusters <- list("train" =
                     list("key" = 1,
                          "R" = matrix(c(0, 1, 1, 0), ncol=2),
                          "b" = c(0, 0)
                          ),
                 "test" =
                     list("key" = 2,
                          "R" = matrix(c(0, 1, 1, 0), ncol=2),
                          "b" = c(1.0, 1.0)
                          )
                 )

## write the two clusters as separate files
for (name in c("train")) {
    key <- clusters[[name]]$key
    R <- clusters[[name]]$R
    b <- clusters[[name]]$b

    ## just select one of the two moons
    ixs <- which(moons$cluster == key)
    moon <- R %*% t(moons$data[ixs, ]) + b

    ## normalize them to [-1, 1]
    moon[1, ] <- moon[1, ] / 2

    ## write result to file
    write.table(moon, paste("moon.csv", sep = ""),
                append=FALSE,
                row.names = FALSE,
                col.names = FALSE,
                fileEncoding="UTF-8")
}

## plot results
moon <- read.table("moon.csv")
devplot(paste("moon.jpg", sep = ""),
        function() {
            ## ------------------------------------------------------------------------
            main.name <- paste("cluster = ", key, sep = "")
            main.ranges <- paste("[[", round(min(moon[1, ]), digits=2), ", ",round(max(moon[1, ]), digits=2), "]",
                                 " [", round(min(moon[2, ]), digits=2), ", ", round(max(moon[2, ]), digits=2), "]]", sep = "")
            main.samples <- paste(", n = ", ncol(moon), sep = "")
            main <- paste(main.name, ", ", main.ranges, ", ", main.samples, sep = "")
            ## ------------------------------------------------------------------------
            plot(t(moon), main=main)
        }, dev = "jpg")
