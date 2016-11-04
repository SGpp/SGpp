source('~/.Rprofile')
source('readData.R')
source('plotTube.R')
source('plotSurpluses.R')
source('plotErrorDecay.R')
source('plotDensity.R')

basepath <- "plots"

## devplot(paste(basepath, 'dnorm.pdf', sep="/"),
##         function() {
##           plot2dNormal()
##         }, dev = 'latex', width = 5.03, height = 5)

## devplot(paste(basepath, 'expectation-value-data-6.pdf', sep="/"),
##         function() {
##           par(mar=c(5.1, 5, 4.1, 16), xpd=FALSE)
##           plotExpectationValue(stats$mc$moments,
##                                col='black', pch=0, ylim=c(0, 0.2))
##           ## plotExpectationValue(stats$fgb0deg2$sg_l1$moments, add=TRUE, col='orange', pch=1)
##           plotExpectationValue(stats$sgb0deg2$sg_l4$moments, add=TRUE, col='purple', pch=2)
##           plotExpectationValue(stats$fgb0deg2$sg_l3$moments, add=TRUE, col='red', pch=3)
##           ## plotExpectationValue(stats$sgb0deg2$sg_l4$moments, add=TRUE, col='darkgreen', pch=4)
##           plotExpectationValue(stats$csr2b0deg2$i1$moments, add=TRUE, col='blue', pch=5)
##           ## plotExpectationValue(stats$csr2b0deg2$moments, add=TRUE, col='black', pch=6)
##           ## plotExpectationValue(stats$csrSGDE2b0deg2i4$moments, add=TRUE, col='purple', pch=3)

##           par(xpd=TRUE)
##           legend('topright',
##                  c('Monte Carlo (n=1066)',
##                    'regular SG level 1 (n=1)',
##                    'regular SG level 2 (n=7)',
##                    'regular SG level 3 (n=31)',
##                    'regular SG level 4 (n=111)'),
##                    ## 'adaptive SG (n=17)',
##                    ## 'adaptive SG (n=25)',
##                    ## 'adaptive SG (n=53)',
##                    ## 'adaptive SG (n=74)',
##                    ## 'aSG data (n=499)'),
##                    ## 'aSG SGDE (n=74506)'),
##                  pch = c(0, 1, 2, 3, 4, 5, 6),
##                  col = c('black', 'orange', 'purple', 'red',
##                    'darkgreen', 'blue', 'black'),
##                  cex = 1.3,
##                  box.lwd = 1, box.col = "black", bg = "white",
##                  inset=c(-0.73,0))
##         }, dev='latex', width = 7.03, height = 5)


## devplot(paste(basepath, 'variance-value-0.pdf', sep=""),
##         function() {
##           par(mar=c(5.1, 5, 4.1, 15), xpd=FALSE)
##           plotVariance(stats$mc$moments, col='black', pch=0, ylim=c(0, 0.012))
##           plotVariance(stats$sgb0deg2$sg_l4$moments, add=TRUE, col='blue', pch=5)
##           plotVariance(stats$fgb0deg2$sg_l3$moments, add=TRUE, col='orange', pch=1)

##           par(xpd=TRUE)
##           legend('topright',
##                  c('Monte Carlo (n=1066)'),
##                  ## 'adaptive SG (n=74)'),
##                  ## 'regular SG (n=111)'),
##                  pch = c(0, 5, 1),
##                  col = c('black', 'blue', 'orange'),
##                  cex = 1.3,
##                  box.lwd = 1, box.col = "black", bg = "white",
##                  inset=c(-0.65,0))

##         }, dev='latex', width = 7.03, height = 5)


## for (i in seq(along = stats$mc$iterative)) {
##   devplot(paste(basepath, '/MC-expectation-i', i, '.png', sep=''),
##           function() {
##             data <- stats$mc$iterative[[i]]
##             n <- data[, 'grid_size']
##             mu <- data[, 'expectation_value']
##             alpha <- 0.01
##             z <- qnorm(1 - alpha / 2)
##             sigma <- z * sqrt(data[, 'variance'] / n)
##             ylim <- c(min(mu - sigma), max(mu + sigma))

##             plot(n, mu,
##                  type='b', ylim = ylim,
##                  xlab = 'samples',
##                  ylab = '$\\E(u)$',
##                  main = paste('t = ', data[1, 'iteration'] / (24 * 60^2), sep=''))
##             lines(n, mu - sigma, col='red')
##             lines(n, mu + sigma, col='red')
##           }, dev='png', width = 5.03, height = 5)
## }

data <- list('fgb0deg1' = stats$uniform$fgb0deg1$squared.stats,
             'sgb0deg1' = stats$uniform$sgb0deg1$squared.stats,
             'sccb0deg1' = stats$uniform$sccb0deg1$squared.stats)
             ## 'rss2b0deg1' = stats$rss2b0deg1$squared.stats,
             ## 'rev2b0deg1' = stats$rev2b0deg1$squared.stats,
             ## ## 'ava1b0deg1' = stats$ava1b0deg1$squared.stats,
             ## ## 'ava2b0deg1' = stats$ava2b0deg1$squared.stats,
             ## ## 'avc2b0deg1' = stats$avc2b0deg1$squared.stats,
             ## 'ava2b0deg1\\_2' = stats$ava2b0deg1_2$squared.stats[c(seq(1, 14), 18, 25, 37), ])

devplot(paste(basepath, 'l2ErrorDecay_squared.pdf', sep="/"),
        function() {
          par(mar=c(5.1, 6.6, 4.1, 13), xpd=FALSE)
          plotErrorDecay(data, dtype = 'testL2Error',
                         log='xy',
                         cols=stats$cols,
                         pchs=stats$pchs,
                         xaxis.base = 10)
          par(xpd=TRUE)
          title('$L_2$ error' )
          legend("topright",
                 names(data),
                 col = stats$cols,
                 pch = stats$pchs,
                 box.lwd = 1, box.col = "black", bg = "white",
                 cex = 1.3, lty = 1,
                 inset=c(-0.42,0))
          box()
        }, dev="latex", output.format="eps", width = 8, height = 5)


data <- list('fgb0deg1' = stats$uniform$fgb0deg1$moments,
             'sgb0deg1' = stats$uniform$sgb0deg1$moments,
             'sccb0deg1' = stats$uniform$sccb0deg1$moments)
             ##'rss2b0deg1' = stats$rss2b0deg1$moments,
             ## 'ava2b0deg1' = stats$ava2b0deg1$moments)
             ## 'ava1b0deg1' = stats$ava1b0deg1$moments,
             ## 'ava2b0deg1' = stats$ava2b0deg1$moments,
             ## 'avc2b0deg1' = stats$avc2b0deg1$moments,
             ## 'ava2b0deg1_2' = stats$ava2b0deg1_2$moments[c(seq(1, 14), 18, 25, 37), ])


devplot(paste(basepath, 'meanErrorDecay.pdf', sep="/"),
        function() {
          par(mar=c(5.1, 6.6, 4.1, 13), xpd=FALSE)
          plotErrorDecay(data, dtype = 'meanRelativeDifference',
                         cols=stats$cols,
                         pchs=stats$pchs,
                         log='xy',
                         xaxis.base = 10)
          par(xpd=TRUE)
          title('$(\\mathbb{E}[g] - \\mathbb{E}[f_N]) / \\mathbb{E}[g]$')
          legend("topright",
                 names(data),
                 col = stats$cols,
                 pch = stats$pchs,
                 box.lwd = 1, box.col = "black", bg = "white",
                 cex = 1.3, lty = 1,
                 inset=c(-0.42,0))
          box()
        }, dev="latex", output.format="eps", width = 8, height = 5)


devplot(paste(basepath, 'varianceErrorDecay.pdf', sep="/"),
        function() {
          par(mar=c(5.1, 6.6, 4.1, 13), xpd=FALSE)
          plotErrorDecay(data, dtype = 'varRelativeDifference',
                         cols=stats$cols,
                         pchs=stats$pchs,
                         log='xy',
                         xaxis.base = 10)
          par(xpd=TRUE)
          title('$(\\mathbb{V}[g] - \\mathbb{V}[f_N]) / \\mathbb{V}[g]$')
          legend("topright",
                 names(data),
                 col = stats$cols,
                 pch = stats$pchs,
                 box.lwd = 1, box.col = "black", bg = "white",
                 cex = 1.3, lty = 1,
                 inset=c(-0.42,0))
          box()
        }, dev="latex", output.format="eps", width = 8, height = 5)

data <- list('fgb0deg1' = stats$uniform$fgb0deg1$simple.stats,
             'sgb0deg1' = stats$uniform$sgb0deg1$simple.stats,
             'sccb0deg1' = stats$uniform$sccb0deg1$simple.stats)

devplot(paste(basepath, 'l2normDecay.pdf', sep="/"),
        function() {
          par(mar=c(5.1, 6.6, 4.1, 13), xpd=FALSE)
          plotErrorDecay(data, dtype = 'testL2Error',
                         cols=stats$cols,
                         pchs=stats$pchs,
                         log='xy',
                         xaxis.base = 10)
          par(xpd=TRUE)
          title('$||f - f^{\\text{SG}}||_{L_2}$')
          legend("topright",
                 names(data),
                 col = stats$cols,
                 pch = stats$pchs,
                 box.lwd = 1, box.col = "black", bg = "white",
                 cex = 1.3, lty = 1,
                 inset=c(-0.42,0))
          box()
        }, dev="latex", output.format="eps", width = 8, height = 5)



## plotSurpluses(stats$uniform$sgb0deg1)
