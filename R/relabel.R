#
#  Copyright (C) 2008 Friedrich Leisch and Bettina Gruen
#  $Id: relabel.R 4142 2008-10-01 13:19:19Z leisch $
#
###*********************************************************

setMethod("relabel", signature(object="flexmix", by="character"),
function(object, by, data=NULL, ...)
{
    MYCALL <- match.call()
    p <- parameters(object)
    n <- grep(by, rownames(p))
    if(length(n)==0)
        stop("No parameters match argument", sQuote("by"))
    if(length(n)!=1)
        stop("More than one parameter matches argument", sQuote("by"))
        
    o <- order(p[n,])
    p <- posterior(object)[,o]

    args <- object@call

    if(is.null(data)){
        mfc <- call("model.frame", formula = args$formula, data=args$data)
        data <- eval(mfc, parent.frame())
    }
    args[[1]] <- NULL
    for(n in c("k", "nrep", "control", "verbose"))
        if(n %in% names(args)) args[[n]] <- NULL

    control <- object@control
    control@iter.max=1
    control@verbose=0
    
    args <- as.list(args)
    args$data <- data
    args$control <- control
    args$cluster <- p

    z <- do.call("flexmix", args)
    z@call <- MYCALL
    z@converged <- TRUE
    z
})
