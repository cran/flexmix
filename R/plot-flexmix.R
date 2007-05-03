#
#  Copyright (C) 2004-2005 Bettina Gruen
#  $Id: plot.R 2669 2006-07-10 13:57:41Z gruen $
#

determine_y <- function(h, root) {
  y <- h$counts
  if (root) y <- sqrt(y)
  return(y)
}

panel.rootogram <- 
function (x, breaks, equal.widths = TRUE, nint = max(round(log2(length(x)) + 1), 3), alpha = plot.polygon$alpha,
          col = plot.polygon$col, border = plot.polygon$border, 
          lty = plot.polygon$lty, lwd = plot.polygon$lwd, subscripts, groups, mark, root = TRUE, markcol, ...) 
{
    x <- as.numeric(x)
    plot.polygon <- trellis.par.get("plot.polygon")
    require(grid)
    grid.lines(x = c(0.05, 0.95), y = unit(c(0, 0), "native"), 
        gp = gpar(col = border, lty = lty, lwd = lwd, alpha = alpha),
               default.units = "npc")
    if (length(x) > 0) {
        if (is.null(breaks)) {
            breaks <- if (equal.widths) 
                do.breaks(range(x, finite = TRUE), nint)
            else quantile(x, 0:nint/nint, na.rm = TRUE)
        }
        h <- lattice:::hist.constructor(x, breaks = breaks, plot = FALSE, ...)
        y <- determine_y(h, root)
        if (!is.null(mark)) {
          h1 <- lattice:::hist.constructor(x[groups[subscripts] == mark], breaks = h$breaks, plot = FALSE, ...)
          y1 <- determine_y(h1, root)
        }
        nb <- length(breaks)
        if (length(y) != nb - 1) 
            warning("problem with hist computations")
        if (nb > 1) {
            panel.rect(x = breaks[-nb], y = 0, height = y, width = diff(breaks), 
                col = col, alpha = alpha, border = border, lty = lty, 
                lwd = lwd, just = c("left", "bottom"))
            if (!is.null(mark)) panel.rect(x = breaks[-nb], y = 0, height = y1, width = diff(breaks),
                                           col = markcol, alpha = alpha, border = border, lty = lty, 
                                           lwd = lwd, just = c("left", "bottom"))
        }
    }
}

prepanel.rootogram <- 
function (x, breaks = NULL, equal.widths = TRUE, nint = max(round(log2(length(x)) + 1), 3), root = TRUE, ...) 
{
  if (length(x) < 1) 
    list(xlim = NA, ylim = NA, dx = NA, dy = NA)
  else {
    if (is.factor(x)) {
      isFactor <- TRUE
      xlimits <- levels(x)
    }
    else isFactor <- FALSE
    if (!is.numeric(x)) 
      x <- as.numeric(x)
    if (is.null(breaks)) {
      breaks <- if (equal.widths) 
        do.breaks(range(x, finite = TRUE), nint)
      else quantile(x, 0:nint/nint, na.rm = TRUE)
    }
    h <- lattice:::hist.constructor(x, breaks = breaks, plot = FALSE, ...)
    y <- determine_y(h, root)
    xlim <- range(x, finite = TRUE)
    list(xlim = if (isFactor) xlimits else range(x, breaks, 
           finite = TRUE), ylim = range(0, y, finite = TRUE), 
         dx = 1, dy = 1)
  }
}


setMethod("plot", signature(x="flexmix", y="missing"),
function(x, y, mark=NULL, markcol=NULL, col=NULL, 
         eps=1e-4, root=TRUE, ylim=TRUE, main=NULL, xlab = "", ylab = "",
         as.table = TRUE, breaks = NULL, ...){

    k <- length(x@prior)

    if(is.null(markcol)) markcol <- flexmix:::FullColors[5]
    if(is.null(col)) col <- flexmix:::LightColors[4]

    if(is.null(main)){
        main <- ifelse(root,
                      "Rootogram of posterior probabilities",
                      "Histogram of posterior probabilities")
        main <- paste(main, ">", eps)
    }
    if (is.null(x@weights))
      z <- data.frame(posterior = as.vector(x@posterior$scaled),
                      component = rep(paste("Comp.", 1:x@k), each = nrow(x@posterior$scaled)),
                      cluster = rep(as.vector(x@cluster), k))
    else
      z <- data.frame(posterior = rep(as.vector(x@posterior$scaled/x@weights), rep(x@weights, k)),
                      component = rep(paste("Comp.", 1:x@k), each = sum(x@weights)),
                      cluster = rep(rep(as.vector(x@cluster), x@weights), k))

    require(lattice)
    panel <- function(x, subscripts, groups, breaks, ...)
      panel.rootogram(x, breaks = breaks, root = root, mark = mark, col = col, markcol = markcol,
                      subscripts = subscripts, groups = groups, ...)
    prepanel <- function(x, breaks, ...) prepanel.rootogram(x, breaks = breaks, root = root, ...)
    z <- subset(z, posterior > eps)
    if (is.logical(ylim)) {
      scales <- if (ylim) list() else list(y = list(relation = "free"))
      hh <- histogram(~ posterior | component, data = z, main = main,  ylab = ylab, xlab = xlab, groups = cluster, 
                      panel = panel, prepanel = prepanel, scales = scales, as.table = as.table, breaks = breaks, ...)
    }
    else hh <- histogram(~ posterior | component, data = z, main = main, ylab = ylab, xlab = xlab, groups = cluster, 
                         ylim = ylim, panel = panel, prepanel = prepanel, as.table = as.table, breaks = breaks, ...)
    hh
  })
