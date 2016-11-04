##----------------------------------------------------------
## Surpluses
##----------------------------------------------------------

collectFGInfos <- function(stats, time, levels) {
  s <- 0
  ans <- c()
  for (i in levels) {
    key <- paste('sg_l', i, sep='')
    data <- stats[[key]]$infos
    data <- data[which(data$time == time), ]
    ## set iteration
    data[, 'iteration'] <- i - 1
    ## search for unseen points
    plevels <- data[, c('l_0', 'l_1', 'l_2')]

    ixs <- c()
    ## collect grid points on the border
    for (j in seq(1, nrow(plevels))) {
      level <- plevels[j, ]
      if (any(level == i) || (i == 1 && any(level == 0))) {
        ixs <- c(ixs, j)
      }
    }

    ans <- rbind(ans, data[ixs, ])
    s <- s + length(ixs)
    show(paste(nrow(data), "==", s, "(", i, ")"))
  }
  return (ans)
}

collectSGInfos <- function(stats, time, levels) {
  s <- 0
  m <- 1
  ans <- c()
  for (i in levels) {
    key <- paste('sg_l', i, sep='')
    data <- stats[[key]]$infos
    data <- data[which(data[, 'time'] == time), ]
    
    ## set iteration
    data[, 'iteration'] <- i - 1
    ## search for unseen points
    plevels <- data[, c('l_0', 'l_1', 'l_2')]
    levelsum <- apply(plevels, 1, sum)
    ## collect inner grid points
    ixs <- which(levelsum > m)
    
    ## collect grid points on the border
    for (j in seq(1, nrow(plevels))) {
      level <- plevels[j, ]
        
      if (any(level == 0)) {
        if (max(level == i) || (i == 1 && any(level == 0))) {        
          ixs <- c(ixs, j)
        }
      }
    }

    ans <- rbind(ans, data[ixs, ])
    m <- max(levelsum)
    s <- s + length(ixs)
    show(paste(length(levelsum), "==", s, "(", m, ")"))
  }
  return (ans)
}


plotSurpluses <- function(stats, time=864000, levels=c(3)) {
  data <- collectSGInfos(stats, time, levels)

  d <- data.frame()
  
  x <- c()
  ym <- c()
  yu <- c()
  yl <- c()

  for (iteration in unique(data[, 2])) {
    ix <- which(data[, 2] == iteration)
    d <- rbind(d, cbind(rep(iteration + 1, length(ix)), data[ix, ncol(data)]))

    x <- c(x, iteration + 1)
    ym <- c(ym, mean(data[ix, ncol(data)]))
    yu <- c(yu, max(data[ix, ncol(data)]))
    yl <- c(yl, min(data[ix, ncol(data)]))
  }
  names(d) <- c('iteration', 'surplus')
  
  boxplot(surplus~iteration, data = d,
          ylim = c(-0.5, 1),#min(data[, 8]), max(data[, 8])),
          xlab = 'iteration',
          ylab = 'hierarchical coefficients $v_{\\mathbf{l}, \\mathbf{i}}$',
          col = 'grey80',
          cex.axis = 1.7, cex.lab = 1.7)
  ## lines(x, yu)
  ## lines(x, yl)  
}

