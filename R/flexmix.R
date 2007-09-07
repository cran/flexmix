#
#  Copyright (C) 2004-2005 Friedrich Leisch and Bettina Gruen
#  $Id: flexmix.R 3691 2007-09-07 08:30:31Z gruen $
#


## The following two methods only fill in and rearrange the model argument
setMethod("flexmix",
          signature(formula = "formula", model="missing"),
function(formula, data=list(), k=NULL, cluster=NULL,
         model=NULL, concomitant=NULL, control=NULL, weights=NULL)
{
    mycall = match.call()
    z <- flexmix(formula=formula, data=data, k=k, cluster=cluster,
                 model=list(FLXMRglm()), concomitant=concomitant,
                 control=control, weights = weights)
    z@call <- mycall
    z
})

setMethod("flexmix",
          signature(formula = "formula", model="FLXM"),
function(formula, data=list(), k=NULL, cluster=NULL, 
         model=NULL, concomitant=NULL, control=NULL, weights=NULL)
{
  mycall = match.call()
  z <- flexmix(formula=formula, data=data, k=k, cluster=cluster, 
               model=list(model), concomitant=concomitant,
               control=control, weights=weights)
  z@call <- mycall
  z
})


## This is the real thing
setMethod("flexmix",
          signature(formula = "formula", model="list"),
function(formula, data=list(), k=NULL, cluster=NULL,
         model=NULL, concomitant=NULL, control=NULL, weights=NULL)
{
    mycall = match.call()
    control = as(control, "FLXcontrol")
    if (!is(concomitant, "FLXP")) concomitant <- FLXPconstant()
    
    groups <- .FLXgetGrouping(formula, data)
    model <- lapply(model, FLXcheckComponent, k, cluster)
    k <- unique(unlist(sapply(model, FLXgetK, k)))
    if (length(k) > 1) stop("number of clusters not specified correctly")
    
    model <- lapply(model, FLXgetModelmatrix, data, groups$formula)
    
    groups$groupfirst <-
        if (length(groups$group)) groupFirst(groups$group)
        else rep(TRUE, FLXgetObs(model[[1]]))
    
    concomitant <- FLXgetModelmatrix(concomitant, data = data,
                                     groups = groups)
    
    if (is(weights, "formula")) {
      weights <- model.frame(weights, data = data, na.action = NULL)[,1]
    }
    ## check if the weights are integer
    ## if non-integer weights are wanted modifications e.g.
    ## for classify != weighted and
    ## plot,flexmix,missing-method are needed
    if (!is.null(weights) & !identical(weights, as.integer(weights)))
      stop("only integer weights allowed")
    
    postunscaled <- initPosteriors(k, cluster, FLXgetObs(model[[1]]))
    
    z <- FLXfit(model=model, concomitant=concomitant, control=control,
                postunscaled=postunscaled, groups=groups, weights = weights)
    
    z@formula = formula
    z@call = mycall
    z@k0 = as.integer(k)
    z
})

###**********************************************************

setMethod("FLXgetK", signature(model = "FLXM"), function(model, k, ...) k)
setMethod("FLXgetObs", signature(model = "FLXM"), function(model) nrow(model@x))
setMethod("FLXcheckComponent", signature(model = "FLXM"), function(model, ...) model)
setMethod("FLXremoveComponent", signature(model = "FLXM"), function(model, ...) model)

setMethod("FLXmstep", signature(model = "FLXM"), function(model, weights, ...) {
  apply(weights, 2, function(w) model@fit(model@x, model@y, w))
})

setMethod("FLXdeterminePostunscaled", signature(model = "FLXM"), function(model, components, ...) {
  sapply(components, function(x) x@logLik(model@x, model@y))
})

###**********************************************************
setMethod("FLXfit", signature(model="list"),
function(model, concomitant, control, postunscaled=NULL, groups, weights)
{
  ### initialize
  k <- ncol(postunscaled)
  N <- nrow(postunscaled)
  control <- allweighted(model, control, weights)
  if(control@verbose>0)
    cat("Classification:", control@classify, "\n")
  if (control@classify=="random") iter.rm <- 0
  group <- groups$group
  groupfirst <- groups$groupfirst
  if(length(group)>0) postunscaled <- groupPosteriors(postunscaled, group)
  postscaled <- postunscaled/rowSums(postunscaled)
  
  llh <- -Inf
  if (control@classify=="random") llh.max <- -Inf
  converged <- FALSE   
  ### EM
  for(iter in 1:control@iter.max) {
      ### M-Step
      postscaled = .FLXgetOK(postscaled, control, weights)
      prior <- if (is.null(weights))
        ungroupPriors(concomitant@fit(concomitant@x, postscaled[groupfirst,,drop=FALSE]),
                      group, groupfirst)
      else ungroupPriors(concomitant@fit(concomitant@x, (postscaled/weights)[groupfirst & weights > 0,,drop=FALSE], weights[groupfirst & weights > 0]),
                         group, groupfirst)
      # Check min.prior
      nok <- which(colSums(prior)/sum(prior) < control@minprior)
      if(length(nok)) {
        if(control@verbose>0)
          cat("*** Removing",length(nok), "component(s) ***\n")
        prior <- prior[,-nok,drop=FALSE]
        prior <- prior/sum(prior)
        postscaled <- postscaled[,-nok,drop=FALSE]
        postscaled[rowSums(postscaled) == 0,] <- prior
        postscaled <- postscaled/rowSums(postscaled)
        if (!is.null(weights)) postscaled <- postscaled * weights
        k <- ncol(prior)
        if (k == 0) stop("all components removed")
        if (control@classify=="random") {
          llh.max <- -Inf
          iter.rm <- iter
        }
        model <- lapply(model, FLXremoveComponent, nok)
      }
      components <- lapply(model, FLXmstep, postscaled)
      postunscaled <- matrix(0, nrow = N, ncol = k)
      for (n in 1:length(model)) 
        postunscaled <- postunscaled + FLXdeterminePostunscaled(model[[n]], components[[n]])
      if(length(group)>0)
        postunscaled <- groupPosteriors(postunscaled, group)
      ### E-Step     
      postunscaled <- exp(postunscaled)
      for(m in 1:k)
        postunscaled[,m] <- prior[,m] * postunscaled[,m]
      ##<FIXME>: wenn eine beobachtung in allen Komonenten extrem
      ## eine postunscaled-werte hat, ist exp(-postunscaled)
      ## numerisch Null, und damit postscaled NaN
      ## log(rowSums(postunscaled)) ist -Inf
      ##</FIXME>
      postscaled <- postunscaled/rowSums(postunscaled)
      if (any(is.nan(postscaled))) {
        index <- which(as.logical(rowSums(is.nan(postscaled))))
        postscaled[index,] <- rep(prior, each=length(index))
        postunscaled[index,] <- .Machine$double.xmin
      }
      ### check convergence
      llh.old <- llh
      llh <- if (is.null(weights)) sum(log(rowSums(postunscaled[groupfirst,,drop=FALSE])))
             else sum(log(rowSums(postunscaled[groupfirst,,drop=FALSE]))*weights[groupfirst])
      if(is.na(llh) | is.infinite(llh))
        stop(paste(formatC(iter, width=4),
                   "Log-likelihood:", llh))
      if (abs(llh-llh.old)/(abs(llh)+0.1) < control@tolerance){
        if(control@verbose>0){
          printIter(iter, llh)
          cat("converged\n")
        }
        converged <- TRUE
        break
      }
      if (control@classify=="random") {
        if (llh.max < llh) {
          components.max <- components
          prior.max <- prior
          postscaled.max <- postscaled
          postunscaled.max <- postunscaled
          llh.max <- llh
        }
      }
      if(control@verbose && (iter%%control@verbose==0))
        printIter(iter, llh)
    }

  ### Construct return object
  if (control@classify=="random") {
    components <- components.max
    prior <- prior.max
    postscaled <- postscaled.max
    postunscaled <- postunscaled.max
    llh <- llh.max
    iter <- control@iter.max - iter.rm
  }

  components <- lapply(1:k, function(i) lapply(components, function(x) x[[i]]))
  names(components) <- paste("Comp", 1:k, sep=".") 
  cluster <- max.col(postscaled)
  size <-  if (is.null(weights)) tabulate(cluster, nbins=k) else tabulate(rep(cluster, weights), nbins=k)
  names(size) <- 1:k
  concomitant <- fillConcomitant(concomitant, postscaled[groupfirst,,drop=FALSE], weights)
  df <- concomitant@df(concomitant@x, k) + sum(sapply(components, sapply, slot, "df"))
  control@nrep <- 1
  
  retval <- new("flexmix", model=model, prior=colMeans(prior),
                posterior=list(scaled=postscaled,
                  unscaled=postunscaled),
                weights = weights,
                iter=iter, cluster=cluster, size = size,
                logLik=llh, components=components,
                concomitant=concomitant,
                control=control, df=df, group=group, k=as(k, "integer"),
                converged=converged)
  retval
})

###**********************************************************
.FLXgetOK = function(p, control, weights){
    n = ncol(p)
    N = 1:n
    if (is.null(weights)) {
      if (control@classify == "weighted")
        return(p)
      else {
        z = matrix(FALSE, nrow = nrow(p), ncol = n)
        if(control@classify %in% c("CEM", "hard")) 
            m = max.col(p)
        else if(control@classify %in% c("SEM", "random")) 
            m = apply(p, 1, function(x) sample(N, size = 1, prob = x))
        else stop("Unknown classification method")
        z[cbind(1:nrow(p), m)] = TRUE
      }
    }else {
      if(control@classify=="weighted")
        z <- p * weights
      else{
        z = matrix(FALSE,  nrow=nrow(p), ncol=n)       
        if(control@classify %in% c("CEM", "hard")) {
          m = max.col(p)
          z[cbind(1:nrow(p), m)] = TRUE
          z <- z * weights
        }
        else if(control@classify %in% c("SEM", "random")) 
          z = t(sapply(1:nrow(p), function(i) table(factor(sample(N, size=weights[i], prob=p[i,], replace=TRUE), N))))
        else stop("Unknown classification method")
      }
    }
    z
}    

###**********************************************************

.FLXgetGrouping <- function(formula, data) {
  group = factor(integer(0))
  lf = length(formula)
  formula1 = formula
  if(length(formula[[lf]])>1 &&
     deparse(formula[[lf]][[1]]) == "|"){
    group = as.factor(eval(formula[[lf]][[3]], data))
    formula1[[lf]] = formula[[lf]][[2]]
  }
  return(list(group=group, formula=formula1))
}

setMethod("FLXgetModelmatrix", signature(model="FLXM"),
function(model, data, formula, lhs=TRUE, ...)
{
  if(is.null(model@formula))
    model@formula = formula
  
  ## model@fullformula = update.formula(formula, model@formula)
  ## <FIXME>: ist das der richtige weg, wenn ein punkt in beiden
  ## formeln ist?
  model@fullformula = update(terms(formula, data=data), model@formula)
  ## </FIXME>
  
  if (lhs) {
    mf <- model.frame(model@fullformula, data=data, na.action = NULL)
    model@terms <- attr(mf, "terms")
    model@y <- as.matrix(model.response(mf))
    model@y <- model@preproc.y(model@y)
  }
  else {
    mt1 <- terms(model@fullformula, data=data)
    mf <- model.frame(delete.response(mt1), data=data, na.action = NULL)
    model@terms<- attr(mf, "terms")
    ## <FIXME>: warum war das da???
    ## attr(mt, "intercept") <- attr(mt1, "intercept")
    ## </FIXME>
  }
  X <- model.matrix(model@terms, data=mf)
  model@contrasts <- attr(X, "contrasts")
  model@x <- X
  model@x <- model@preproc.x(model@x)
  model@xlevels <- .getXlevels(model@terms, mf)
  model
})

## groupfirst: for grouped observation we need to be able to use
## the posterior of each group, but for computational simplicity
## post(un)scaled has N rows (with mutiple identical rows for each
## group). postscaled[groupfirst,] extracts posteriors of each
## group ordered as the appear in the data set.
groupFirst <- function(x) !duplicated(x)

## if we have a group variable, set the posterior to the product
## of all density values for that group (=sum in logarithm)
groupPosteriors <- function(x, group)
{    
    for(g in levels(group)){
        gok <- group==g
        if(any(gok)){
            x[gok,] <- matrix(colSums(x[gok,,drop=FALSE]),
                              nrow=sum(gok), ncol=ncol(x), byrow=TRUE)
        }
    }
    x
}

ungroupPriors <- function(x, group, groupfirst) {
  if (!length(group)) group <- 1:length(groupfirst)
  if (nrow(x) > 1) {
    x <- x[order(as.integer(group[groupfirst])),,drop=FALSE]
    x <- x[as.integer(group),,drop=FALSE]
  }
  x
}

allweighted <- function(model, control, weights) {
  allweighted = all(sapply(model, function(x) x@weighted))
  if(allweighted){
    if(control@classify=="auto")
      control@classify="weighted"
  }
  else{
    if(control@classify %in% c("auto", "weighted"))
      control@classify="hard"
    if(!is.null(weights))
      stop("it is not possible to specify weights for models without weighted ML estimation")
  }
  control
}

initPosteriors <- function(k, cluster, N) {
  if(is(cluster, "matrix")){
    postunscaled <- cluster
    if (!is.null(k)) if (k != ncol(postunscaled)) stop("specified k does not match the number of columns of cluster")
  }
  else{
    if(is.null(cluster)){
      if(is.null(k))
        stop("either k or cluster must be specified")
      else
        cluster <- sample(1:k, size=N, replace=TRUE)
    }
    else{
      cluster <- as(cluster, "integer")
      if (!is.null(k)) if (k != max(cluster)) stop("specified k does not match the values in cluster")
      k <- max(cluster)
    }
    postunscaled <- matrix(0.1, nrow=N, ncol=k)
    for(K in 1:k){
      postunscaled[cluster==K, K] <- 0.9
    }
  }
  postunscaled
}

###**********************************************************

setMethod("predict", signature(object="FLXdist"),
function(object, newdata=list(), aggregate=FALSE, ...){
    if (missing(newdata)) return(fitted(object, aggregate=aggregate, drop=FALSE))
    x = list()
    for(m in 1:length(object@model)) {
      comp <- lapply(object@components, "[[", m)
      x[[m]] <- predict(object@model[[m]], newdata, comp, ...)
    }
    if (aggregate) {
      z <- lapply(x, function(z) matrix(rowSums(matrix(sapply(1:object@k, function(K) z[[K]] * object@prior[K]), ncol = object@k)),
                                        nrow = nrow(z[[1]])))
    }
    else {
      z <- list()
      for (k in 1:object@k) {
        z[[k]] <- do.call("cbind", lapply(x, "[[", k))
      }
      names(z) <- paste("Comp", 1:object@k, sep=".")
    }
    z
})

###**********************************************************

setMethod("parameters", signature(object="FLXdist"),
function(object, component=NULL, model=NULL, simplify=TRUE, drop=TRUE)
{
    if (is.null(component)) component <- 1:object@k
    if (is.null(model)) model <- 1:length(object@model)
    
    if (simplify) {
      parameters <- sapply(model, function(m) sapply(object@components[component], function(x) unlist(x[[m]]@parameters), simplify=TRUE),
                           simplify = FALSE)
    }
    else {
      parameters <- sapply(model, function(m) sapply(object@components[component], function(x) x[[m]]@parameters, simplify=FALSE),
                           simplify = FALSE)
    }
    if (drop) {
      if (length(component) == 1 && !simplify) parameters <- lapply(parameters, "[[", 1)
      if (length(model) == 1) parameters <- parameters[[1]]
    }
    return(parameters)
})

setMethod("prior", signature(object="FLXdist"),
function(object)
{
  object@prior
})

setMethod("posterior", signature(object="flexmix", newdata="missing"),
function(object)
{
    object@posterior$scaled
})

setMethod("posterior", signature(object="FLXdist", newdata="listOrdata.frame"),
          function(object, newdata, unscaled=FALSE,...) {
            comp <- lapply(object@components, "[[", 1)
            postunscaled <- posterior(object@model[[1]], newdata, comp, ...)
            if (length(object@model) > 1) {
              for (m in 2:length(object@model)) {
                comp <- lapply(object@components, "[[", m)
                postunscaled <- postunscaled + posterior(object@model[[m]], newdata, comp, 
                                                         ...)
              }
            }
            groups <- .FLXgetGrouping(object@formula, newdata)
            if(length(groups$group)>0)
              postunscaled <- groupPosteriors(postunscaled, groups$group)
            for(m in 1:object@k)
              postunscaled[,m] <- object@prior[m] * exp(postunscaled[,m])
            if (unscaled) return(postunscaled)
            else return(postunscaled/rowSums(postunscaled))
})            

setMethod("posterior", signature(object="FLXM", newdata="listOrdata.frame"),
function(object, newdata, components, ...) {
            mf <- model.frame(terms(object@fullformula, data=newdata),
                              data = newdata, na.action = NULL)
            x <- model.matrix(attr(mf, "terms"), data = mf)
            y <- as.matrix(model.response(mf))
            matrix(sapply(components, function(z) z@logLik(x, y, ...)),
                   nrow = nrow(y))
})
    
setMethod("cluster", signature(object="flexmix", newdata="missing"),
function(object)
{
    object@cluster
})
    
setMethod("cluster", signature(object="flexmix", newdata="ANY"),
function(object, newdata, ...)
{
    max.col(posterior(object, newdata, ...))
})

###**********************************************************

setMethod("summary", "flexmix",
function(object, eps=1e-4, ...){    
    z <- new("summary.flexmix",
             call = object@call,
             AIC = AIC(object),
             BIC = BIC(object),
             logLik = logLik(object))

    TAB <- data.frame(prior=object@prior,
                      size=object@size)
    rownames(TAB) <- paste("Comp.", 1:nrow(TAB), sep="")
    TAB[["post>0"]] <- colSums(object@posterior$scaled > eps)
    TAB[["ratio"]] <- TAB[["size"]]/TAB[["post>0"]]
    
    z@comptab = TAB
    z
    
})

###**********************************************************

