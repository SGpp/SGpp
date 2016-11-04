source('~/.Rprofile')

##----------------------------------------------------------
## Prepare the data sets
##----------------------------------------------------------
basepath <- '../results'

stats <- list()
for (param.setting in list.dirs(basepath, recursive = FALSE,
                                full.names = FALSE)) {
    ## reference values
    param.basepath <- paste(basepath, "/", param.setting, sep="")
    stats[[param.setting]] <- list()
    stats[[param.setting]]$ref <- readARFFData(paste(param.basepath, '/',
                                    param.setting, '.ref.moments.arff', sep = ''))

    ## adaptive sparse grid
    available <- c()

    ## 'rss2b0deg1',
    ##                'rev2b0deg1',
    ##                'ava1b0deg1',
    ##                'ava2b0deg1',
    ##                'ava2b0deg1_2',
    ##                'avc2b0deg1')

    for (key in available) {
        path <- paste(param.basepath, '/', key, '/', sep='')
        moments <- readARFFData(paste(path, 'sg', '.moments.arff', sep=''))
        simpleStats <- readARFFData(paste(path, 'sg.SIMPLE.stats.arff', sep=''))
        squaredStats <- readARFFData(paste(path, 'sg.SQUARED.stats.arff', sep=''))
        for (iteration in unique(moments[, 'iteration'])) {
            ix <- which(moments[, 'iteration'] == iteration)
            itkey <- paste('i', iteration, sep = '')
            stats[[param.setting]][[key]][[itkey]] <- list(
                'moments' = moments[ix, , drop=FALSE],
                'simple.stats' = simpleStats[ix, , drop=FALSE],
                'squared.stats' = squaredStats[ix, , drop=FALSE]
            )
            ## add absolute relative error of the estimated
            ## values if there is a reference value given
            i <- which('ref' == names(stats[[param.setting]]))
            if (length(i) > 0) {
                ## some reference values found
                itkeyMoments <- stats[[param.setting]][[key]][[itkey]]$moments
                if (stats[[param.setting]]$ref[, 'mean'] != 0) {
                    err <- abs((stats[[param.setting]]$ref[, 'mean'] - itkeyMoments[, 'mean']) / stats[[param.setting]]$ref[, 'mean'])
                    itkeyMoments <- cbind(itkeyMoments, err)
                    colnames(itkeyMoments)[ncol(itkeyMoments)] <- 'meanRelativeDifference'
                }
                if (stats[[param.setting]]$ref[, 'var'] != 0) {
                    err <- abs((stats[[param.setting]]$ref[, 'var'] - itkeyMoments[, 'var']) / stats[[param.setting]]$ref[, 'var'])
                    itkeyMoments <- cbind(itkeyMoments, err)
                    colnames(itkeyMoments)[ncol(itkeyMoments)] <- 'varRelativeDifference'
                }
                stats[[param.setting]][[key]][[itkey]]$moments <- itkeyMoments
            }
        }
        ## gather moments for the regular parts
        for (items in stats[[param.setting]][[key]]) {
            for (name in names(items)) {
                stats[[param.setting]][[key]][[name]] <- rbind(stats[[param.setting]][[key]][[name]], items[[name]])
            }
        }
    }

    ## regular sparse grid
    available <- list('sgb0deg1' = list('levels' = seq(1, 9)),
                      'fgb0deg1' = list('levels' = seq(1, 6)),
                      'sccb0deg1' = list('levels' = seq(1, 9)))

    for (key in names(available)) {
        ## gather moments
        radix = strsplit(key, "b")[[1]][1]
        for (level in available[[key]]$levels) {
            path <- paste(param.basepath, '/', key, '/', radix, '_l', level, '/', radix, sep='')
            levelkey <- paste('l', level, sep='')
            stats[[param.setting]][[key]][[levelkey]] <- list(
                'moments' = readARFFData(paste(path, '.moments.arff', sep='')),
                'simple.stats' = readARFFData(paste(path, '.SIMPLE.stats.arff', sep='')),
                'squared.stats' = readARFFData(paste(path, '.SQUARED.stats.arff', sep=''))
            )
            ## add absolute relative error of the estimated
            ## values if there is a reference value given
            i <- which('ref' == names(stats[[param.setting]]))
            if (length(i) > 0) {
                ## some reference values found
                moments <- stats[[param.setting]][[key]][[levelkey]]$moments
                if (stats[[param.setting]]$ref[, 'mean'] != 0) {
                    err <- abs((stats[[param.setting]]$ref[, 'mean'] - moments[, 'mean']) / stats[[param.setting]]$ref[, 'mean'])
                    moments <- cbind(moments, err)
                    colnames(moments)[ncol(moments)] <- 'meanRelativeDifference'
                }
                if (stats[[param.setting]]$ref[, 'var'] != 0) {
                    err <- abs((stats[[param.setting]]$ref[, 'var'] - moments[, 'var']) / stats[[param.setting]]$ref[, 'var'])
                    moments <- cbind(moments, err)
                    colnames(moments)[ncol(moments)] <- 'varRelativeDifference'
                }
                stats[[param.setting]][[key]][[levelkey]]$moments <- moments
            }
        }

                                        # gather moments for the regular parts
        stats[[param.setting]]$regular[[key]] <- list()
        for (items in stats[[param.setting]][[key]]) {
            for (name in names(items)) {
                stats[[param.setting]][[key]][[name]] <- rbind(stats[[param.setting]][[key]][[name]], items[[name]])
            }
        }

        ## L2 error surpluses is the difference between levels
        stats[[param.setting]][[key]]$simple.stats[, 'L2ErrorSurpluses'] <-
            c(stats[[param.setting]][[key]]$simple.stats[1, 'L2ErrorSurpluses'],
              diff(stats[[param.setting]][[key]]$simple.stats[, 'L2ErrorSurpluses']))
        stats[[param.setting]][[key]]$squared.stats[, 'L2ErrorSurpluses'] <-
            c(stats[[param.setting]][[key]]$squared.stats[1, 'L2ErrorSurpluses'],
              diff(stats[[param.setting]][[key]]$squared.stats[, 'L2ErrorSurpluses']))
    }

}
# stats$mc <- list('moments' = readARFFData(paste(basepath, '/mc/mc.moments.arff', sep='')))

## ## read mc values per time steps and sample size
## stats$mc$iterative <- list()
## files <- list.files(paste(basepath, '/mc', sep = ''), pattern='_', full.names = TRUE)

## i <- 1
## for (file in files) {
##   stats$mc$iterative[[i]] <- readARFFData(file)
##   i <- i + 1
## }

stats[['cols']] <- c(1, 2, 3, 4, 'yellow3', 6, 'purple3', 'brown1', 'orangered4')
stats[['pchs']] <- seq(4, 4 + length(stats))
stats[['ixs']] <- c(1, 2)
stats[['names']] <- c('regular sg',    #1
                      'adaptive sg',   #2
                      'monte carlo')    #3)
