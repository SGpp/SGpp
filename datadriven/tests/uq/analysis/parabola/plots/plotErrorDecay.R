##----------------------------------------------------------
## Error decay
##----------------------------------------------------------

plotErrorDecay <- function(stats, dtype = 'meanRelativeDifference',
                           xlab='', ylab='',
                           cols = NULL, pchs = NULL, xlim = NULL,
                           ylim = NULL, log = 'y',
                           xaxis.base = 2, yaxis.base = 10,
                           ...) {
  options(scipen=99)
  gridsizes <- unlist(lapply(stats, function(x) {
    return (x[, 'grid_size', drop = FALSE])
  }))

  errors <- unlist(lapply(stats, function(x) {
    return (x[, dtype, drop = FALSE])
  }))
  yaxis.pot <- c(trunc(log(min(errors), base=yaxis.base)) - 1,
                ceiling(log(max(errors), base=yaxis.base)))
  xaxis.pot <- c(0, ceiling(log(max(gridsizes), base=xaxis.base)))

  ## par(mar=c(5, 6, 4, 2))
  ## par(mar=c(5, 6, 2, 2))

  if (log == 'xy') {
    xaxt <- 'n'
    if (is.null(xlim))
      xlim = c(xaxis.base^xaxis.pot[1], xaxis.base^xaxis.pot[2])
  } else if (log == 'y') {
    xaxt <- 's'
    if (is.null(xlim))
      xlim <- c(min(gridsizes), max(gridsizes))
  }

  if (is.null(ylim)) {
    ylim <- c(yaxis.base^yaxis.pot[1], yaxis.base^yaxis.pot[2])
  }
  
  plot(NULL, NULL, type = "n",
       yaxt = "n",
       xaxt = xaxt,
       xlab = "",
       ylab = "",
       log = log,
       ylim = ylim,
       xlim = xlim,
       cex.lab = 1.7, , cex.axis = 1.7,
       ...)

  mtext(xlab, side = 1, lin = 3, cex = 1.7)
  mtext(ylab, side = 2, lin = 5, cex = 1.7)

  ## yaxis
  s <- seq(yaxis.pot[1], yaxis.pot[2], by = 1)
  axis(2,
       at = sapply(s, function(x) {yaxis.base^x}),
       label = sapply(s, function(x) {
         paste('$', yaxis.base, '^{', x, '}$', sep = '') }),
       las = 1, cex.axis = 1.7)

  # draw grid
  for (si in s) {
    for (i in seq(1, yaxis.base, 2)) {
      abline(h=i * yaxis.base^si, lty='dashed', lwd = .3, col='lightgrey')
    }
  }

  ## xaxis
  if (log == 'xy') {
    s <- seq(xaxis.pot[1], xaxis.pot[2], by = 1)
    axis(1,
         at = sapply(s, function(x) {xaxis.base^x}),
         label = sapply(s, function(x) {
           paste('$', xaxis.base, '^{', x, '}$', sep = '') }),
         las = 1, cex.axis = 1.7)

    ## draw grid
    for (si in s) {
      for (i in seq(1, xaxis.base)) {
        abline(v=i * xaxis.base^si, lty='dashed', lwd = .3, col='lightgrey')
      }
    }
  } else {
    ## draw grid
    for (si in seq(xlim[1], xlim[2], 100)) {
      abline(v=si, lty='dashed', lwd = .3, col='lightgrey')
    }
  }

  if (is.null(cols) || length(cols) < length(stats)) {
    cols <- seq(along = stats)
  }

  if (is.null(pchs) || length(pchs) < length(stats)) {
    pchs <- seq(along = stats)
  }

  i <- 1
  for (items in stats) {
    lines(items[, 'grid_size'],
          items[, dtype],
          col = cols[i], type='l', lwd = 2)
    if (nrow(items) < 100) {
      points(items[, 'grid_size'],
             items[, dtype],
             pch = pchs[i], col=cols[i])
    } else {
      ix <- seq(1, nrow(items), 5)
      points(items[ix, 'grid_size'],
             items[ix, dtype],
             pch = pchs[i], col=cols[i])
    }
    i <- i + 1
  }
}
