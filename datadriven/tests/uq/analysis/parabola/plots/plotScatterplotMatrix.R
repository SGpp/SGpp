plotScatterPlotMatrix <- function(data, grid, alpha, names, t) {
  show('Start Testing Pairs Plot...')

  scatterPlotMatrix(data, grid, alpha,
                    labels = labels,
                    lower=insertHeatmap,
                    upper=insertGrid)

}

##   dim <- ncol(data)
  
##   if (dim > 2) {
##   } else {
##     # plot grid
##     devplot(paste("UQ", cap, "-", ref, "-", level, "-sg.pdf", sep=""), function() {
##       par(mar=c(4.5, 4, 3, 8))
##       plotGrid(grid, alpha,
##                xlab = "$\\phi$",
##                ylab = "$K_a$")

##       image.plot(legend.only=TRUE,
##                  zlim=c(min(alpha), max(alpha)),
##                  legend.width = 0.8,
##                  legend.shrink = 0.8,
##                  nlevel = 256,
##                  legend.mar = 7.5,
##                  axis.args = list(cex.axis = 2),
##                  legend.args = list(text='surpluses',
##                    side=4, font=2, line=4, cex=2))
##     }, dev='latex', width = 5.75, height = 5)

##     # plot function
##     devplot(paste("UQ", cap, "-", ref, "-", level, "-sg-function.pdf", sep=""), function() {
##       par(mar=c(4.5, 4, 3, 8))
##       plotFunction(data, alpha,
##                    xlab = "$\\Delta x$",
##                    ylab = "$\\alpha$")

##       image.plot(legend.only=TRUE,
##                  zlim=c(min(data[, 3]), max(data[, 3])),
##                  legend.width = 0.8,
##                  legend.shrink = 0.8,
##                  nlevel = 256,
##                  legend.mar = 7.5,
##                  axis.args = list(cex.axis = 2),
##                  legend.args = list(text='damage',
##                    side=4, font=2, line=4, cex=2))
##     }, dev='latex', width = 5.75, height = 5)
##   }
## }
