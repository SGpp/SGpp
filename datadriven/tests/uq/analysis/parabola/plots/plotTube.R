plotTube <- function(moments) {
  show('plotting tube')

  time <- moments[, 1, drop=FALSE] / (24 * 60 * 60)
  E <- moments[, 2]
  sigma <- moments[, 3, drop=FALSE]

  plot(time, E, type = 'l',
       ylim = c(min(E-sigma), max(E+sigma)),
       xaxs = 'i', yaxs = 'i',
       xlab = 'Time [days]',
       ylab = 'CO$_2$ Leakage (\\%)',
       cex.axis = 1.7,
       cex.lab = 1.7)
  
  polygon(c(time, rev(time)), c(E + sigma, rev(E)),
          col = '#6F8BB550', border = NA)
  polygon(c(time, rev(time)), c(E - sigma, rev(E)),
          col = '#6F8BB550', border = NA)
  lines(time, E - sigma, lwd = 2, col = '#00488490')
  lines(time, E + sigma, lwd = 2, col = '#00488490')
  lines(time, E, lwd = 2, col = 'black')
  legend("bottomright", legend=c(
                          '$\\mathbb{E}(f_\\mathcal{I})$',
                          '$\\mathbb{E}(f_\\mathcal{I}) \\pm \\mathbb{V}(f_\\mathcal{I})$'),
         lty = 1, col = c('black', '#00488490'), lwd=4,
         box.lwd = 1, box.col = "black", bg = "white",
         cex = 1.5)
  box()
}

plotExpectationValue <- function(moments, add=FALSE, ...) {
  show('plotting expectation value')
  
  time <- moments[, 1, drop=FALSE] / (24 * 60 * 60)
  time.sorted <- sort(time, index.return=TRUE)
  time <- time.sorted$x
  E <- moments[time.sorted$ix, 'mean']

  if (!add) {
    plot(time, E,
         xaxs = 'i', yaxs = 'i',
         xlab = 'Time [days]',
         ylab = '$\\mathbb{E}(u)$ of CO$_2$ Leakage (\\%)',
         cex.axis = 1.7,
         cex.lab = 1.7, ...)
  } else {
    points(time, E, ...)
  }
  
  lines(time, E, lwd = 2, ...)
  ## legend("bottomright", legend=c('$\\mathbb{E}(f_\\mathcal{I})$'),
  ##        lty = 1, col = c('black', '#00488490'), lwd=4,
  ##        box.lwd = 1, box.col = "black", bg = "white",
  ##        cex = 1.5)
  box()  
}


plotVariance <- function(moments, add=FALSE, ...) {
  show('plotting var')

  time <- moments[, 1, drop=FALSE] / (24 * 60 * 60)
  time.sorted <- sort(time, index.return=TRUE)
  time <- time.sorted$x
  sigma <- moments[time.sorted$ix, 'var', drop=FALSE]  

  if (!add) {
    plot(time, sigma,
         xaxs = 'i', yaxs = 'i',
         xlab = 'Time [days]',
         ylab = '$\\mathbb{V}(u)$ of CO$_2$ Leakage (\\%)',
         cex.axis = 1.7,
         cex.lab = 1.7, ...)
  } else {
    points(time, sigma, ...)
  }

  lines(time, sigma, lwd = 2, ...)
  box()  
}
