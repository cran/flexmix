prepanel.default.coef <- function (x, y, subscripts, groups=NULL, horizontal = TRUE, nlevels, origin = NULL, 
                                   ...) 
{
  if (any(!is.na(x) & !is.na(y))) {
    if (horizontal) {
      if (!is.factor(y)) {
        if (missing(nlevels)) 
          nlevels <- length(unique(y))
        y <- factor(y, levels = 1:nlevels)
      }
      if (!is.null(groups)) {
        if (!is.numeric(x)) stop("x must be numeric")
              x <- rep(x, each = 2) + rep(groups[subscripts], each = 2) *c(-1,1)
      }
      list(xlim = if (is.numeric(x)) range(x, origin, finite = TRUE) else levels(x), 
           ylim = levels(y), yat = sort(unique(as.numeric(y))), 
           dx = 1, dy = 1)
    }
    else {
      if (!is.factor(x)) {
        if (missing(nlevels)) 
          nlevels <- length(unique(x))
        x <- factor(x, levels = 1:nlevels)
      }
      if (!is.null(groups)) {
        if (!is.numeric(y)) stop("y must be numeric")
        y <- rep(as.numeric(y), each = 2) + rep(groups[subscripts], each = 2) *c(-1,1)
      }
      list(xlim = levels(x), xat = sort(unique(as.numeric(x))), 
           ylim = if (is.numeric(y)) range(y, origin, finite = TRUE) else levels(y), 
           dx = 1, dy = 1)
    }
  }
  else list(xlim = c(NA, NA), ylim = c(NA, NA), dx = 1, dy = 1)
}

panel.coef <- function(x, y, subscripts, groups, significant = NULL, horizontal = TRUE, 
                       lwd = 2, col, col.line = c("black", "grey"), ...)
{
  col.sig <- rep(col.line[1], length(x))
  if (!is.null(significant)) {
    if (missing(col)) col <-   c("grey", "white")
    col.fill <- rep(col[1], length(x))
    col.sig[!significant[subscripts]] <- col.line[2]
    col.fill[!significant[subscripts]] <- col[2]
  }
  else if (missing(col)) col.fill <- "grey" else col.fill <- col
  panel.barchart(x, y, border = col.sig, col = col.fill, horizontal = horizontal, ...)
  if (!missing(groups)) {
    if (horizontal) {
      z <- x + rep(c(-1,1), each = length(x)) * matrix(rep(groups[subscripts], 2), ncol = 2)
      for (i in 1:length(x)) {
        panel.xyplot(z[i,], rep(y[i], 2), type = "l", col = col.sig[i], lwd = lwd)
      }
    }
    else {
      z <- y + rep(c(-1,1), each = length(y)) * matrix(rep(groups[subscripts], 2), ncol = 2)
      for (i in 1:length(y)) {
        panel.xyplot(rep(x[i], 2), z[i,], type = "l", col = col.sig[i], lwd = lwd)
      }
    }
  }
}

setGeneric("getCoefs", function(x, ...) standardGeneric("getCoefs"))

setMethod("getCoefs", signature(x = "list"),
function(x, alpha = 0.05, components, ...) {
  Comp <- list()
  for (n in components) {
    Comp[[n]] <- getCoefs(x[[n]], alpha)
    Comp[[n]]$Component <- n
  }
  only <- names(which(table(unlist(sapply(Comp, names))) == length(components)))
  do.call("rbind", lapply(Comp, "[", TRUE, only))
})

setMethod("getCoefs", signature(x = "FLXRMRglm"),
function(x, alpha = 0.05, ...) {
  cm <- x@summary
  data.frame(Value = cm[,1],
             SD = cm[,2] * qt(1-alpha/2, df=x@fitted$df.residual),
             Variable = rownames(cm),
             Significance = cm[,4] <= alpha)
})

setMethod("getCoefs", signature(x = "FLXRM"),
function(x, alpha = 0.05, ...)
{
  cm <- x@summary
  data.frame(Value = cm,
             Variable = names(cm))
})


setMethod("getCoefs", signature(x = "FLXRP"),
function(x, alpha = 0.05, components, ...)
{
  cm <- x@summary
  data.frame(Value = cm[components],
             Variable = "pi",
             Component = components)
})

setMethod("getCoefs", signature(x = "FLXRPmultinom"),
function(x, alpha = 0.05, components, ...) {
  cm <- x@summary
  cm <- cm[names(cm) %in% components]
  Comp <- lapply(names(cm), function(n) 
                 data.frame(Value = cm[[n]][,1],
                            SD = cm[[n]][,2] * qt(1-alpha/2, df=x@fitted$df.residual),
                            Variable = rownames(cm[[n]]),
                            Component = n,
                            Significance = cm[[n]][,4] <= alpha))
  do.call("rbind", Comp)
})

setMethod("plot", signature(x="FLXR", y="missing"),
function(x, y, bycluster=TRUE, alpha=0.05, components, labels=NULL,
         scale = bycluster, xlab = NULL, ylab = NULL, ci = TRUE,
         scales = list(), as.table = TRUE, horizontal = TRUE, ...)
{
    if (missing(components)) components <- 1:x@k
    plot.data <- getCoefs(summary(x)@refit, alpha, components)
    if (!is.null(labels)) plot.data$Variable <- factor(plot.data$Variable, labels = labels)
    if (!"SD" %in% colnames(plot.data)) ci <- FALSE
    if (scale) {
      ord <- order(plot.data$Variable)
      m <- with(plot.data, tapply(abs(Value), Variable, max))
      if (ci) plot.data[,c("Value", "SD")] <- plot.data[,c("Value", "SD")]/m[plot.data$Variable]
      else plot.data[,"Value"] <- plot.data[,"Value"]/m[plot.data$Variable]
      if (horizontal) scales$x$draw <- FALSE else scales$y$draw <- FALSE
    }
    plot.data$Component <- with(plot.data, factor(Component, sort(unique(Component)), labels = paste("Comp.", sort(unique(Component)))))
    if (bycluster) {
      formula <- if (horizontal) Variable ~ Value | Component else Value ~ Variable | Component
      plot.data$Variable <- with(plot.data, factor(Variable, levels = rev(Variable)))
    }
    else {
      formula <- if (horizontal) Component ~ Value | Variable else Value ~ Component | Variable
      plot.data$Component <- with(plot.data, factor(Component, levels = rev(levels(Component))))
    }
    require(lattice)
    groups <- if (ci) plot.data$SD else NULL
    significant <- if (ci) plot.data$Significance else NULL
    xyplot(formula, data = plot.data, xlab = xlab, ylab = ylab, origin = 0, horizontal = horizontal,
           scales = scales, as.table = as.table, significant = significant,
           groups = groups, prepanel = function(...) prepanel.default.coef(...),
           panel = function(x, y, subscripts, groups, ...)
           panel.coef(x, y, subscripts, groups, ...), ...)
})

