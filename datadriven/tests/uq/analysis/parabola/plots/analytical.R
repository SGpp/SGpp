library(lattice)

f <- function(x, y) (4 * x * (1 - x)) * (4 * y * (1 - y))
f <- function(x, y) dnorm(x, 0.2, 0.1) * dnorm(y, 0.2, 0.1)

x <- seq(0, 1, len = 40)
y <- seq(0, 1, len = 40)
g <- expand.grid(x = x, y = y)
g$z <- f(g$x, g$y)

devplot('analytical_pdf.pdf', function() {
  cex <- 1.6
  lattice.options(axis.padding = list(numeric = .1))
  print(wireframe(z ~ x * y, data = g,
                  scales = list(
                    arrows = FALSE,
                    cex = cex),
                  drape = TRUE, colorkey = FALSE,
                  col.regions = tim.colors(100),
                  screen = list(z = 30, x = -70),
                  par.settings = list(axis.line = list(col = "transparent"),
                    layout.heights=list(top.padding=-10, bottom.padding=-10),
                    layout.widths=list(right.padding=-5, left.padding=0)),
                  xlab = list(label='$\\theta_1$', cex=cex),
                  ylab = list(label='$\\theta_2$', cex=cex),
                  zlab = list(label='$u(\\theta_1, \\theta_2)$', cex=cex)))
}, dev='latex')
